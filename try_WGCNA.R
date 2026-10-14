library(dplyr)
library(tidyr)
library(tibble)
library(WGCNA)
library(mixOmics)

hormones <- c("PROG", "X17BES", "FSH", "LH", "PRL")
cov_names <- c("storage_months", "Age", "BMI", "batch", "from")
phases <- c("F", "O", "EL", "LL")

# TRUE: retain cycle-associated variation
# FALSE: remove average phase effects as well
retain_phase <- TRUE

# Assumes all other columns in d_wide are the shared proteins
proteins <- setdiff(names(d_wide), c("SampleID", "ID", "phase"))

stopifnot(
  !anyDuplicated(d_wide$SampleID),
  !anyDuplicated(covariates$SampleID),
  !anyDuplicated(metadata$SampleID)
)

# Match all three tables by SampleID
samples <- Reduce(
  intersect,
  list(d_wide$SampleID, covariates$SampleID, metadata$SampleID)
)

meta <- covariates[
  match(samples, covariates$SampleID),
  c("SampleID", "ID", "phase", cov_names)
]

P <- as.matrix(
  d_wide[match(samples, d_wide$SampleID), proteins, drop = FALSE]
)

H <- as.matrix(
  metadata[match(samples, metadata$SampleID), hormones, drop = FALSE]
)

# Require complete metadata and hormone measurements
keep <- complete.cases(meta) & rowSums(!is.finite(H)) == 0

meta <- meta[keep, , drop = FALSE]
P <- P[keep, , drop = FALSE]
H <- H[keep, , drop = FALSE]

# Within-woman analysis requires at least two observations
keep <- ave(
  rep(1, nrow(meta)), meta$ID, FUN = length
) >= 2

meta <- meta[keep, , drop = FALSE]
P <- P[keep, , drop = FALSE]
H <- H[keep, , drop = FALSE]

meta$ID <- factor(meta$ID)
meta$phase <- factor(meta$phase, levels = phases)
meta$batch <- factor(meta$batch)
meta$from <- factor(meta$from)

rownames(P) <- rownames(H) <- meta$SampleID

# Transparent complete-protein analysis: no imputation
keep_protein <- colSums(!is.finite(P)) == 0
message("Removing ", sum(!keep_protein), " proteins with missing/nonfinite values.")
P <- P[, keep_protein, drop = FALSE]

# Subtract each woman's mean, separately for every variable
within_center <- function(M, id) {
  M <- as.matrix(M)
  
  for (person in unique(as.character(id))) {
    idx <- which(as.character(id) == person)
    M[idx, ] <- sweep(
      M[idx, , drop = FALSE],
      2,
      colMeans(M[idx, , drop = FALSE]),
      "-"
    )
  }
  
  M
}

Pw <- within_center(P, meta$ID)
Hw <- within_center(H, meta$ID)

# Within-woman phase design
Dphase <- model.matrix(~ phase, data = meta)[, -1, drop = FALSE]
Dphase <- within_center(Dphase, meta$ID)

# Exclude covariates with only one observed value
active_covs <- cov_names[
  vapply(meta[cov_names], function(x) length(unique(x)) > 1, logical(1))
]

Dcov <- model.matrix(reformulate(active_covs), data = meta)
Dcov <- within_center(Dcov, meta$ID)

# Remove intercept and columns constant within women
Dcov <- Dcov[, apply(abs(Dcov), 2, max) > 1e-10, drop = FALSE]

D <- cbind(Dphase, Dcov)

if (qr(D)$rank < ncol(D)) {
  stop("Phase and within-woman covariates are confounded/redundant; inspect D.")
}

adjust_within <- function(M) {
  fit <- lm.fit(x = D, y = M)
  
  if (!retain_phase) {
    out <- fit$residuals
  } else if (ncol(Dcov) > 0) {
    nuisance_idx <- ncol(Dphase) + seq_len(ncol(Dcov))
    out <- M - Dcov %*% fit$coefficients[nuisance_idx, , drop = FALSE]
  } else {
    out <- M
  }
  
  dimnames(out) <- dimnames(M)
  out
}

X <- adjust_within(Pw)
Y <- adjust_within(Hw)

# Remove proteins with no remaining variation
X <- X[, apply(X, 2, sd) > 1e-8, drop = FALSE]

if (any(apply(Y, 2, sd) <= 1e-8)) {
  stop("At least one hormone has no remaining variation.")
}

cat(
  "Samples:", nrow(X),
  "\nWomen:", n_distinct(meta$ID),
  "\nProteins:", ncol(X), "\n"
)






# Inspect candidate soft-thresholding powers
powers <- c(1:10, 12, 14, 16, 18, 20)

sft <- WGCNA::pickSoftThreshold(
  X,
  powerVector = powers,
  networkType = "signed",
  verbose = 0
)

print(sft$fitIndices)

par(mfrow = c(1, 2))

plot(
  sft$fitIndices$Power,
  sft$fitIndices$SFT.R.sq,
  type = "b",
  xlab = "Power",
  ylab = "Scale-free topology fit"
)

plot(
  sft$fitIndices$Power,
  sft$fitIndices$mean.k.,
  type = "b",
  xlab = "Power",
  ylab = "Mean connectivity"
)

par(mfrow = c(1, 1))

# Exploratory starting value: inspect diagnostics and sensitivity
soft_power <- 6

net <- WGCNA::blockwiseModules(
  X,
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  corType = "pearson",
  minModuleSize = 20,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  maxBlockSize = ncol(X) + 1,
  randomSeed = 123,
  verbose = 0
)

# Protein-to-module assignments
module_membership <- tibble(
  protein = colnames(X),
  module = unname(net$colors)
)

module_membership %>% count(module, sort = TRUE)

# Grey = proteins not assigned to a module
module_scores <- as.data.frame(net$MEs)
module_scores <- module_scores[
  , names(module_scores) != "MEgrey", drop = FALSE
]

rownames(module_scores) <- meta$SampleID

# Visualise dendrogram and module assignments
WGCNA::plotDendroAndColors(
  net$dendrograms[[1]],
  net$colors[net$blockGenes[[1]]],
  groupLabels = "Module",
  dendroLabels = FALSE,
  hang = 0.03
)

# Descriptive module–hormone correlations
module_hormone_cor <- cor(as.matrix(module_scores), Y)

pheatmap::pheatmap(
  module_hormone_cor,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks = seq(-1, 1, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = nrow(module_hormone_cor) > 1,
  cluster_cols = FALSE,
  main = "Within-woman module–hormone correlations"
)

# Mean module profiles across the four phases
module_phase_profiles <- bind_cols(
  meta[, c("SampleID", "ID", "phase")],
  module_scores
) %>%
  pivot_longer(
    cols = starts_with("ME"),
    names_to = "module",
    values_to = "score"
  ) %>%
  group_by(module, phase) %>%
  summarise(mean_score = mean(score), .groups = "drop")

ggplot2::ggplot(
  module_phase_profiles,
  ggplot2::aes(phase, mean_score, group = module)
) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ module, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Phase", y = "Mean module score")




ncomp <- 2
keepX <- rep(min(20, ncol(X)), ncomp)
keepY <- rep(ncol(Y), ncomp)

set.seed(123)

spls_fit <- mixOmics::spls(
  X = X,
  Y = Y,
  ncomp = ncomp,
  mode = "canonical",
  keepX = keepX,
  keepY = keepY,
  scale = TRUE
)

# Selected proteins and their component loadings
selected_proteins <- bind_rows(
  lapply(seq_len(ncomp), function(k) {
    loadings <- spls_fit$loadings$X[, k]
    
    tibble(
      component = k,
      protein = rownames(spls_fit$loadings$X),
      loading = as.numeric(loadings)
    ) %>%
      filter(abs(loading) > 1e-10) %>%
      arrange(desc(abs(loading)))
  })
)

selected_proteins

# Hormone loadings
hormone_loadings <- bind_rows(
  lapply(seq_len(ncomp), function(k) {
    tibble(
      component = k,
      hormone = rownames(spls_fit$loadings$Y),
      loading = as.numeric(spls_fit$loadings$Y[, k])
    )
  })
)

hormone_loadings

# Sample component scores
spls_scores <- bind_cols(
  meta[, c("SampleID", "ID", "phase")],
  as.data.frame(spls_fit$variates$X) %>%
    rename_with(~ paste0("protein_", .x)),
  as.data.frame(spls_fit$variates$Y) %>%
    rename_with(~ paste0("hormone_", .x))
)

mixOmics::plotIndiv(
  spls_fit,
  comp = c(1, 2),
  group = meta$phase,
  ind.names = FALSE,
  legend = TRUE
)

mixOmics::plotVar(
  spls_fit,
  comp = c(1, 2),
  var.names = c(FALSE, TRUE)
)