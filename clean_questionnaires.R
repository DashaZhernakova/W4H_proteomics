library(dplyr)
library(purrr)

setwd("/Users/Dasha/work/Sardinia/W4H/phenotypes/batch12/")

d <- read.delim("cleaned_questionnaire_141125.csv", sep = ",", as.is = T, check.names = F, row.names = 1)

d[d == 'NA'] <- NA
d[d == ''] <- NA
d[d == ' '] <- NA

d_per_tp <- d[d$Visit_number != 0,]
d_per_tp <- d_per_tp[,colSums(! is.na(d_per_tp)) > 0]

d <- d[d$Visit_number == 0,]
d <- d[,colSums(! is.na(d)) > 0]

d$t2d_first_second <- 3
d$t2d_first_second[d$t2d_familare == 1 | d$t2d_famil_secondo == 1] <- 1
d$t2d_first_second[d$t2d_familare == 2 & d$t2d_famil_secondo == 2] <- 2
d$t2d_familare <- NULL
d$t2d_famil_secondo <- NULL

d$cholesterol_first_second <- 3
d$cholesterol_first_second[d$colester_famil == 1 | d$colester_famil_secondo == 1] <- 1
d$cholesterol_first_second[d$colester_famil == 2 & d$colester_famil_secondo == 2] <- 2
d$colester_famil <- NULL
d$colester_famil_secondo <- NULL


col_summary <- data.frame(matrix(nrow = ncol(d), ncol = 5))
row.names(col_summary) <- colnames(d)
colnames(col_summary) <- c("count_not_na", "num_levels", "num_1", "num_2", "num_3")
cnt = 1

for(col in colnames(d)) {
  num_not_na <- sum(!is.na(d[[col]]))
  num_levels = length(unique(na.omit(d[[col]])))
  num_1 <- NA
  num_2 <- NA
  num_3 <- NA
  if (num_levels == 2){
    num_1 <- table(d[[col]], useNA = "no")[1]
    num_2 <- table(d[[col]], useNA = "no")[2]
  } else if (num_levels == 3){
    num_1 <- table(d[[col]], useNA = "no")[1]
    num_2 <- table(d[[col]], useNA = "no")[2]
    num_3 <- table(d[[col]], useNA = "no")[3]
  }
  
  col_summary[col,] <- c(num_not_na, num_levels, num_1, num_2, num_3 )
}
View(col_summary)

col_summary$select <- NA

col_summary[col_summary$num_levels == 1, "select"] <- F
col_summary[col_summary$count_not_na < 30, "select"] <- F

col_summary[col_summary$num_levels == 2 & (col_summary$num_1 < 10 | col_summary$num_2 < 10), "select"] <- F
col_summary[col_summary$num_levels == 2 & (col_summary$num_1 >= 10 & col_summary$num_2 >= 10), "select"] <- T

col_summary[col_summary$num_levels == 3 & (col_summary$num_1 < 10 | col_summary$num_2 < 10 | col_summary$num_3 < 10), "select"] <- F
col_summary[col_summary$num_levels == 3 & (col_summary$num_1 >= 10 & col_summary$num_2 >= 10 & col_summary$num_3 >= 10), "select"] <- T

col_summary[grepl("freq",row.names(col_summary)) , "select"] <- "diet"
col_summary[grepl("consumo",row.names(col_summary)) , "select"] <- "diet"
col_summary[grepl("data",row.names(col_summary)) , "select"] <- F
col_summary[grepl("names$",row.names(col_summary)) , "select"] <- F

col_summary['Code', ]$select <- T
col_summary['ID', ]$select <- T
col_summary['Visit_number', ]$select <- T

col_summary['Swab_morning_or_not', ]$select <- F
col_summary['Swab_after_feces', ]$select <- F
col_summary['Bristol_stool_scale', ]$select <- F

col_summary['tipo_patologia', ]$select <- F
col_summary['uso_farmaci', ]$select <- F
col_summary['partners_3_mesi', ]$select <- F

col_summary['tipologia_integratore', ]$select <- F
col_summary['tipo_intolleranza', ]$select <- F
col_summary['altre_note', ]$select <- F

col_summary['fonte', ]$select <- F
col_summary['datavisitaodierna', ]$select <- F
col_summary['giorno_mestruazione', ]$select <- F

col_summary['farmaci_names', ]$select <- F
col_summary['note', ]$select <- F
col_summary['app_check', ]$select <- F
col_summary['app_nouse', ]$select <- F

col_summary['numero_figli', ]$select <- F

phenos_to_select <- row.names(col_summary[col_summary$select == T | is.na(col_summary$select),])
phenos_to_select

d_flt <- d[,phenos_to_select]
write.table(d_flt, file = "questionnaire_201025_selected_visit_0.txt", quote = F, sep = "\t", row.names = FALSE)


d_per_tp$Swab_morning_or_not <- NULL
d_per_tp$Swab_after_feces <- NULL
d_per_tp$numero_gravidanze <- NULL
d_per_tp$datavisitaodierna <- NULL
d_per_tp$giorno_mestruazione <- NULL

d_per_tp$from <- NULL
d_per_tp$app_check <- NULL
d_per_tp$app_nouse <- NULL
d_per_tp$note <- NULL
d_per_tp$raccolta_quando <- NULL

d_per_tp$farmaci_names <- NULL
d_per_tp$integratori_names <- NULL
d_per_tp$ovuli <- NULL
d_per_tp$malattie_names <- NULL
d_per_tp$covidvaccino_data <- NULL
d_per_tp$dataraccolta <- NULL
d_per_tp$partners <- NULL

write.table(d_per_tp, file = "questionnaire_201025_selected_longitudinal.txt", quote = F, sep = "\t", row.names = FALSE)

