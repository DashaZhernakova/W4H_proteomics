# Mediation analysis
library(mediation)
all_combined <- full_join(d_wide, pheno, by = c("SampleID","ID", "phase"))
all_combined$SampleID.x <- NULL
all_combined$SampleID.y <- NULL

all_combined_1 <- all_combined[all_combined$phase == 'F',]

phase = 'F'
exposure = 'PROG'
mediator = 'TFPI'
outcome = 'LDL'
res <- run_simple_mediation(all_combined[all_combined$phase == phase,], exposure, mediator, outcome)
summary(res)

res <- run_simple_mediation(all_combined[all_combined$phase == "EL",], "PROG", "TFPI", "LDL")
summary(res)

trios <- data.frame(matrix(c("HDL", "FAP", "COL",
                    "LDL", "MEGF9", "COL",
                    "COL", "LDLR", "LDL",
                    "COL", "GAS6", "LDL",
                    "COL", "TFPI", "LDL",
                    "ALT", "GSTA1", "COL",
                    "COL", "GAS6", "AST",
                    "COL", "ADH4", "ALT",
                    "TRI", "AGRP", "COL",
                    "TRI", "LDLR", "LDL",
                    "PROG", "EPHA1", "AST",
                    "LH", "SPON2", "LDL",
                    "PROG", "TFPI", "LDL",
                    "PROG", "CHRDL2", "LDL",
                    "17BES", "CHRDL2", "LDL",
                    "HOMA_IR", "IGFBP1", "PROG",
                    "CRLF1", "PROG", "CLSTN2",
                    "HYAL1", "HOMA_IR", "IGFBP1"
                    ), byrow = T, ncol = 3))
colnames(trios) <- c("exposure", "mediator", "outcome")

res_med <- data.frame(matrix(nrow = nrow(trios) * length(all_phases), ncol = 7))
colnames(res_med) <- c("exposure", "mediator", "outcome", "phase", "ACME_p", "ADE_p", "prop_mediated")
cnt <- 1
for (i in 1:nrow(trios)){
  exposure <- trios[i,"exposure"]
  mediator <- trios[i, "mediator"]
  outcome <- trios[i, "outcome"]
  for (phase in all_phases){
    res <- summary(run_simple_mediation(all_combined[all_combined$phase == phase,], exposure, mediator, outcome))
    acme_p <- res$d0.p 
    ade_p <- res$z0.p 
    prop_med <- res$n0
    res_med[cnt,] <- c(trios[i,], phase, acme_p, ade_p, prop_med)
    cnt <- cnt + 1
  }
}

res_med2 <- res_med %>%
  mutate(across(c(ACME_p, ADE_p, prop_mediated), as.numeric))

run_simple_mediation <- function(all_combined_1, exposure, mediator, outcome){
  subs <- na.omit(all_combined_1[,c(exposure, mediator, outcome)])
  colnames(subs) <- c("exp", "med", "out")
  model.0 <- lm(out ~ exp, subs)
  #summary(model.0)
  
  model.M <- lm(med ~ exp, subs)
  #summary(model.M)
  
  model.Y <- lm(out ~ exp + med, subs)
  #summary(model.Y)
  
  med_res <- mediate(model.M, model.Y, treat='exp', mediator='med', sims=1000)
  summary(med_res)
  return(med_res)
}

#### regmed ####
library(regmed)

lmm_res <- read.delim("../results/pheno_prot/prot_vs_pheno_withTP_lmm_adj_covar.txt", sep = "\t", as.is = T, check.names = T)
d1 <- lmm_res[lmm_res$BH_pval < 0.1,c("prot", "pheno", "estimate", "pval")]

exposure = 'PROG'
outcome = 'INS'

hormones <- c("PROG","LH","FSH","X17BES","PRL")
lipids <- colnames(pheno)[!colnames(pheno) %in% c("ID", "SampleID", "TP", hormones)]

all_combined_1 <- na.omit(all_combined_1)
all_combined_1$TP <- NULL
all_combined_1 <-all_combined_1 %>%
  remove_rownames %>%
  column_to_rownames(var="ID")
exp_data <- as.matrix(all_combined_1[, hormones])
out_data <- as.matrix(all_combined_1[, lipids])

prots <- unique(d1[, "prot"])
length(prots)

med_data <- as.matrix(all_combined_1[,prots])

row.names(exp_data) <- row.names(out_data) <- row.names(med_data) <- all_combined_1$ID

lambda.grid <- seq(from = 0.4,  to = 0.01, by = -0.1)
#fit.grid <- regmed.grid(exp_data, med_data, out_data, lambda.grid, frac.lasso = 0.8)
#plot.regmed.grid(fit.grid)
#fit.best <- regmed.grid.bestfit(fit.grid)
#summary(fit.best)

fit.grid <- mvregmed.grid(exp_data[,'PROG'], med_data, out_data, lambda.grid)
plot.mvregmed.grid(fit.grid)

fit.best <- mvregmed.grid.bestfit(fit.grid)
summary(fit.best)

edges.any <- regmed.edges(fit.best, type = "any")
edges.any$edges[edges.any$edges == "exp_data"] <-  exposure
edges.any$edges[edges.any$edges == "out_data"] <-  outcome

plot.regmed.edges(edges.any)



