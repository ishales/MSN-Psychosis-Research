#Data arrangement
library(dplyr)
library(tidyr)
library(tibble)

#Reading and writing tables
library(readr)
library(openxlsx)

#create data table
install.packages("data.table")
library("data.table")
library("ggplot2")

library(broom)
library(tools)

#String operations
library(stringr)

#Functional programming
library(magrittr)
library(purrr)
library(rlist)

# permutation of linear models
library(lmPerm)
library(tictoc)

#Plotting
library(ggplot2)
library(Cairo)
#library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)
library(ggseg)
#library(ggsegExtra)
library(dplyr)

#Functional programming
library(magrittr)
library(broom)
library(HGNChelper)

#Statistical Analyses
library(effsize)
library(pls)

# permutation of linear models
library(lmPerm)
library(tictoc)

#String operations
library(stringr)
library(tools)

#Overlap analysis
#library(GeneOverlap)
library(forcats)
library(plyr)
library(tidyverse)
library(reshape2)
library(purrrlyr)
library(rlist)
library(psych)

library(janitor)

#install neuro combat function
source("https://neuroconductor.org/neurocLite.R")
neuro_install('neuroCombat', release = "stable", release_repo = "github")
library(neuroCombat)


# set directory
Cdirectory <- "/Users/isabellehales/Desktop/capstonefinal/"
setwd(Cdirectory)
# where to save outputs
outdir <- paste(Cdirectory,sep="")


# load in 2 databases
NAPLS2_demo <- read.csv("ClientInfo_13Jan2016_BL_Sub_ID_fixed.csv")
NAPLS2_cannabis <- read.csv("Cannabis_13Jan2016.csv")

# join 2 databases
NAPLS2_all <- left_join(NAPLS2_demo, NAPLS2_cannabis, by = join_by(SiteNumber, SubjectNumber, SubjectType, VisitNumber, VisitLabel, DataCollectedDate, DataQuality))


# MSN function to work with asymmetric hemispheres
MSN_freesurfer <- function(sub, parc, sub_dir){
  lh <- read.table(paste0(sub_dir,sub,"/stats/lh.",parc,".stats"),comment.char = "#")[,c(3,4,5,7,8,9,10)]
  print(lh)
  rh <- read.table(paste0(sub_dir,sub,"/stats/rh.",parc,".stats"),comment.char = "#")[,c(3,4,5,7,8,9,10)]
  print(rh)
  lhnames <- read.table(paste0(sub_dir,sub,"/stats/lh.",parc,".stats"),comment.char = "#")[,1]
  lhnames <- lapply(lhnames, function (x) paste("lh_", x, sep=""))
  rhnames <- read.table(paste0(sub_dir,sub,"/stats/rh.",parc,".stats"),comment.char = "#")[,1]
  rhnames <- lapply(rhnames, function (x) paste("rh_", x, sep=""))
  rownames(lh) <- lhnames
  rownames(rh) <- rhnames
  lhrh <- rbind(lh,rh)
  msn <- cor(t(scale(lhrh)))
  msn
} # 3=area; 4=volume; 5=thickness; 7=mean curvature; 8=gaussian curvature; 9=folding index; 10=curved index

# get subject lists
final.N.all <- NAPLS2_all[!is.na(NAPLS2_all$cannabis_ever_smoked),]
subjects.all <- final.N.all[,"Sub_ID"]

# add empty column for labeling
final.N.all$Subgroup = NA
# replace subgroup title with appropriate measure
final.N.all$Subgroup[final.N.all$Diagnosis == "HC" & final.N.all$cannabis_ever_smoked == "1"] <- "HC.ca"
final.N.all$Subgroup[final.N.all$Diagnosis == "HC" & final.N.all$cannabis_ever_smoked == "0"] <- "HC.noCa"
final.N.all$Subgroup[final.N.all$Diagnosis == "CHRc" & final.N.all$cannabis_ever_smoked == "0"] <- "CHR.noCa"
final.N.all$Subgroup[final.N.all$Diagnosis == "CHRnc" & final.N.all$cannabis_ever_smoked == "0"] <- "CHR.noCa"
final.N.all$Subgroup[final.N.all$Diagnosis == "CHRc" & final.N.all$cannabis_ever_smoked == "1"] <- "CHR.ca"
final.N.all$Subgroup[final.N.all$Diagnosis == "CHRnc" & final.N.all$cannabis_ever_smoked == "1"] <- "CHR.ca"

# extract subjects of each subgroup
sub.HC.can <- final.N.all[final.N.all[,73] == "HC.ca","Sub_ID"]
sub.HC.nocan <- final.N.all[final.N.all[,73] == "HC.noCa","Sub_ID"]
sub.CHR.can <- final.N.all[final.N.all[,73] == "CHR.ca","Sub_ID"]
sub.CHR.nocan <- final.N.all[final.N.all[,73] == "CHR.noCa","Sub_ID"]

# generate MSNs for updated parcellation 500.aparc
mdir = "/Users/isabellehales/Desktop/capstonefinal/napls2D/"
cp_MSNs <- sapply(subjects.all, function(x) list(MSN_freesurfer(x,"500.aparc",mdir)))


# annotate MSNs
#get rid of unknown rows and columns
for (i in 1:length(cp_MSNs)) {
  temp = cp_MSNs[i]
  #View(temp)
  id = temp[[1]]
  #View(id)
  #erase unknown column 1 (left hem)
  id = id[,-1]
  #erase unknown row 1 (left hem)
  id = id[-1,]
  #erase unknown row 153 (right hem)
  id = id[-153,]
  #erase unknown column 153 (right hem)
  id = id[,-153]
  cp_MSNs[[i]] = id
}

#msn means
cp_msn_means_allsub <- sapply(cp_MSNs, rowMeans, simplify = TRUE)
cp_msn_means_allsub_df <- data.frame(cp_msn_means_allsub)
cp_msn_means_allsub.t <- t(cp_msn_means_allsub)

# neurocombat correction
site <- final.N.all[, colnames(final.N.all)[c(4)]]
site.df <- data.frame(site)
site.t <-as.data.frame(t(site.df))
site.tm <- as.matrix(site.t)
all_msn_means_combat <- neuroCombat(cp_msn_means_allsub_df, site.tm) 

# final corrected MSNs
all_correct_msns <- t(all_msn_means_combat[["dat.combat"]])
all_correct_msns <- as.data.frame(all_correct_msns)
all_correct_msns <- rownames_to_column(all_correct_msns, var = "Sub_ID")
all_correct_msns$Sub_ID<-gsub("X","",as.character(all_correct_msns$Sub_ID))

#subgroup msns (uncorrected with neuro combat)
HC.can.msns <- cp_MSNs[names(cp_MSNs) %in% sub.HC.can]
HC.nocan.msns <- cp_MSNs[names(cp_MSNs) %in% sub.HC.nocan]
CHR.can.msns <- cp_MSNs[names(cp_MSNs) %in% sub.CHR.can]
CHR.nocan.msns <- cp_MSNs[names(cp_MSNs) %in% sub.CHR.nocan]

#subgroup means
HC.c.mean <- all_correct_msns[all_correct_msns$Sub_ID %in% sub.HC.can, ]
HC.nc.mean <- all_correct_msns[all_correct_msns$Sub_ID %in% sub.HC.nocan, ]
CHR.c.mean <- all_correct_msns[all_correct_msns$Sub_ID %in% sub.CHR.can, ]
CHR.nc.mean <- all_correct_msns[all_correct_msns$Sub_ID %in% sub.CHR.nocan, ]



# anova prep
# Convert to data frame if needed
library(tibble)
ctf <- as.data.frame(all_correct_msns)
ctf$Sub_ID<-gsub("X","",as.character(ctf$Sub_ID))

pf <- as.data.frame(final.N.all)
msn_means_pheno_F <- left_join(pf, ctf, by = "Sub_ID")

# save output of msns and pheno data
write.xlsx(all_correct_msns, "Corrected MSNs (for Site) for Cannabis Study - 4 Groups.xlsx")
write.xlsx(msn_means_pheno_F, "Demographic data with Correct MSNs for Cannabis Study - 4 Groups.xlsx")

## anova function
msnROI.list <- rownames(cp_msn_means_allsub)

# filtering difference files

mean_pheno_MSN_HC.can <- filter(msn_means_pheno_F, Subgroup == "HC.ca")
mean_pheno_MSN_HC.nocan <- filter(msn_means_pheno_F, Subgroup == "HC.noCa")
mean_pheno_MSN_CHR.can <- filter(msn_means_pheno_F, Subgroup == "CHR.ca")
mean_pheno_MSN_CHR.nocan <- filter(msn_means_pheno_F, Subgroup == "CHR.noCa")

mean.val.MSN.HC.can <- mean_pheno_MSN_HC.can[, c(74:381)]
mean.val.MSN.HC.nocan <- mean_pheno_MSN_HC.nocan[, c(74:381)]
mean.val.MSN.CHR.can <- mean_pheno_MSN_CHR.can[, c(74:381)]
mean.val.MSN.CHR.nocan <- mean_pheno_MSN_CHR.nocan[, c(74:381)]

HC.can.mean <- apply(mean.val.MSN.HC.can, 2, mean)
HC.nocan.mean <- apply(mean.val.MSN.HC.nocan, 2, mean)
CHR.can.mean <- apply(mean.val.MSN.CHR.can, 2, mean)
CHR.nocan.mean <- apply(mean.val.MSN.CHR.nocan, 2, mean)

HC.nocan.can.DIFF <- HC.nocan.mean - HC.can.mean
CHR.nocan.can.DIFF <- CHR.nocan.mean - CHR.can.mean
HC.nocan.CHR.nocan.DIFF <- HC.nocan.mean - CHR.nocan.mean
HC.nocan.CHR.can.DIFF <- HC.nocan.mean - CHR.can.mean

HC.can.CHR.nocan.DIFF <- HC.can.mean - CHR.nocan.mean
HC.can.CHR.can.DIFF <- HC.can.mean - CHR.can.mean

write.xlsx(list(HC.nocan.can.DIFF), "Diff.Mean.Msn_HC.nocan.vs.HC.can.xlsx", rownames = FALSE)
write.xlsx(list(CHR.nocan.can.DIFF ), "Diff.Mean.Msn_CHR.nocan.vs.CHR.can.xlsx")
write.xlsx(list(HC.can.CHR.can.DIFF), "Diff.Mean.Msn_HC.can.vs.CHR.can.xlsx")
write.xlsx(list(HC.nocan.CHR.nocan.DIFF), "Diff.Mean.Msn_HC.nocan.vs.CHR.nocan.xlsx")
write.xlsx(list(HC.can.CHR.nocan.DIFF), "Diff.Mean.Msn_HC.can.vs.CHR.nocan.xlsx")
write.xlsx(list(HC.nocan.CHR.can.DIFF), "Diff.Mean.Msn_HC.nocan.vs.CHR.can.xlsx")



# original ANOVA function
msnANOVA.wCovar <- function(ROI.vector, data_table) {
  ROI.vector_string <- as.name(ROI.vector)
  print(data_table$ROI.vector)
  print(ROI.vector)
  ROI.vector_tmp = data_table[[ROI.vector]]
  print(head(ROI.vector_tmp))
  model_df <- cbind(ROI.vector_tmp, dplyr::select(data_table, Subgroup, demo_age_ym, demo_sex, SiteNumber))
  print(head(model_df))
  msnANOVA.wCovar <- lm(ROI.vector_tmp ~ factor(Subgroup) + demo_age_ym + factor(demo_sex) + factor(SiteNumber), model_df) %>% anova %>% 
    tidy %>%
    filter(!is.na(p.value))
  #round p-values
  msnANOVA.wCovar$p.value
  #print(model_lm)
  #return(model_lm)
}

## NO FDR corrected parametric ANOVAs
###### HC vs CHRc (can) w/ FDR correction
anova.HC.CHR.Can.wCovar <- map(msnROI.list, msnANOVA.wCovar, msn_means_pheno_F) %>%
  reduce(rbind) %>% as_tibble %>%
  #mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.HC.CHR.Can.wCovar) <- c("Subgroup", "demo_age_ym", "demo_sex", "SiteNumber")
anova.HC.CHR.Can.wCovar$ROI <- msnROI.list
write.xlsx(anova.HC.CHR.Can.wCovar, "Anova_4_Groups_NO_FDR.xlsx")

## FDR corrected parametric ANOVAs
###### HC vs CHRc (can) w/ FDR correction
anova.HC.CHR.Can.wCovar.FDR <- map(msnROI.list, msnANOVA.wCovar, msn_means_pheno_F) %>%
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.HC.CHR.Can.wCovar.FDR) <- c("Subgroup", "demo_age_ym", "demo_sex", "SiteNumber")
anova.HC.CHR.Can.wCovar.FDR$ROI <- msnROI.list
write.xlsx(anova.HC.CHR.Can.wCovar.FDR, "Anova_4_Groups.xlsx")


## posthoc analysis

#### posthoc pairwise comparisons - diagnosis (FDR correction)
sig_ROI_D_fdr <- filter(anova.HC.CHR.Can.wCovar.FDR, Subgroup <= .056)
sig_ROI_names_D_fdr <- sig_ROI_D_fdr$ROI

posthoc_lm <- function(sig_ROI_names, data_table) {
  module_tmp = data_table[[sig_ROI_names]]
  
  model_df <- cbind(module_tmp, dplyr::select(data_table, Subgroup))
  print(model_df)
  
  model_aov <- aov(module_tmp ~ as.factor(Subgroup), model_df)  %>% TukeyHSD
  print(model_aov)
  return(model_aov)
}

posthoc.lm.output.diag.FDR <- map(sig_ROI_names_D_fdr, posthoc_lm, msn_means_pheno_F) 
names(posthoc.lm.output.diag.FDR) <- sig_ROI_names_D_fdr 
capture.output(posthoc.lm.output.diag.FDR, file="posthoc.lm.output.fdr_diag")

file_content <- readLines("posthoc.lm.output.fdr_diag")
writeLines(file_content, file.path(Cdirectory, "posthoc_fdr"))



#### posthoc pairwise comparisons - age (FDR correction)
sig_ROI_A_FDR <- filter(anova.HC.CHR.Can.wCovar.FDR, demo_age_ym <= .056)
sig_ROI_names_A_FDR <- sig_ROI_A_FDR$ROI

ph_age <- function(sig_ROI_names, data_table) {
  module_tmp = data_table[[sig_ROI_names]]
  
  model_df <- cbind(module_tmp, dplyr::select(data_table, demo_age_ym))
  print(model_df)
  
  model_aov <- aov(module_tmp ~ as.factor(demo_age_ym), model_df)  %>% TukeyHSD
  print(model_aov)
  return(model_aov)
}

lm_posthoc_age.FDR.output <- map(sig_ROI_names_A_FDR, ph_age, msn_means_pheno_F) 
names(lm_posthoc_age.FDR.output) <- sig_ROI_names_A_FDR
capture.output(lm_posthoc_age.FDR.output, file="posthoc.lm.output.age.FDR")

file_content <- readLines("posthoc.lm.output.age.FDR")
writeLines(file_content, file.path(Cdirectory, "posthoc_age_FDR_corr"))


# re-assessing data with different cannabis threshold values
level_of_use <- final.N.all$cannabis_lifetime_use
table(level_of_use, useNA = "ifany")

hc.thres.nocan <- final.N.all$Sub_ID[final.N.all$Diagnosis == "HC" & final.N.all$cannabis_lifetime_use == 0 | final.N.all$Diagnosis == "HC" & is.na(final.N.all$cannabis_lifetime_use)]
hc.thres.can3 <- final.N.all$Sub_ID[final.N.all$Diagnosis == "HC" & final.N.all$cannabis_lifetime_use > 3 & !is.na(final.N.all$cannabis_lifetime_use)]
chr.thres.nocan <- final.N.all$Sub_ID[final.N.all$Diagnosis == "CHRc" & final.N.all$cannabis_lifetime_use == 0 | final.N.all$Diagnosis == "CHRc" & is.na(final.N.all$cannabis_lifetime_use) | final.N.all$Diagnosis == "CHRnc" & final.N.all$cannabis_lifetime_use == 0 | final.N.all$Diagnosis == "CHRnc" & is.na(final.N.all$cannabis_lifetime_use)]
chr.thres.can3 <- final.N.all$Sub_ID[final.N.all$Diagnosis == "CHRc" & final.N.all$cannabis_lifetime_use > 3 & !is.na(final.N.all$cannabis_lifetime_use) | final.N.all$Diagnosis == "CHRnc" & final.N.all$cannabis_lifetime_use > 3 & !is.na(final.N.all$cannabis_lifetime_use)]
notused3 <- final.N.all$Sub_ID[!(final.N.all$Sub_ID %in% hc.thres.nocan) & !(final.N.all$Sub_ID %in% hc.thres.can3) & !(final.N.all$Sub_ID %in% chr.thres.nocan) & !(final.N.all$Sub_ID %in% chr.thres.can3)]

t3_data <- msn_means_pheno_F[!(msn_means_pheno_F$Sub_ID %in% notused3), ]
## FDR corrected parametric ANOVAs
###### HC vs CHRc (can) w/ FDR correction
anova.thres.3.wCovar.FDR <- map(msnROI.list, msnANOVA.wCovar, t3_data) %>%
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.thres.3.wCovar.FDR) <- c("Subgroup", "demo_age_ym", "demo_sex", "SiteNumber")
anova.thres.3.wCovar.FDR$ROI <- msnROI.list
write.xlsx(anova.thres.3.wCovar.FDR, "Anova_Threshold_3.xlsx")


## threshold 10
hc.thres.can10 <- final.N.all$Sub_ID[final.N.all$Diagnosis == "HC" & final.N.all$cannabis_lifetime_use > 10 & !is.na(final.N.all$cannabis_lifetime_use)]
chr.thres.can10 <- final.N.all$Sub_ID[final.N.all$Diagnosis == "CHRc" & final.N.all$cannabis_lifetime_use > 10 & !is.na(final.N.all$cannabis_lifetime_use) | final.N.all$Diagnosis == "CHRnc" & final.N.all$cannabis_lifetime_use > 3 & !is.na(final.N.all$cannabis_lifetime_use)]
notused10 <- final.N.all$Sub_ID[!(final.N.all$Sub_ID %in% hc.thres.nocan) & !(final.N.all$Sub_ID %in% hc.thres.can10) & !(final.N.all$Sub_ID %in% chr.thres.nocan) & !(final.N.all$Sub_ID %in% chr.thres.can10)]

t10_data <- msn_means_pheno_F[!(msn_means_pheno_F$Sub_ID %in% notused10), ]
## FDR corrected parametric ANOVAs
###### HC vs CHRc (can) w/ FDR correction
anova.thres.10.wCovar.FDR <- map(msnROI.list, msnANOVA.wCovar, t10_data) %>%
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.thres.10.wCovar.FDR) <- c("Subgroup", "demo_age_ym", "demo_sex", "SiteNumber")
anova.thres.10.wCovar.FDR$ROI <- msnROI.list
write.xlsx(anova.thres.10.wCovar.FDR, "Anova_Threshold_10.xlsx")

## threshold 50
hc.thres.can50 <- final.N.all$Sub_ID[final.N.all$Diagnosis == "HC" & final.N.all$cannabis_lifetime_use > 50 & !is.na(final.N.all$cannabis_lifetime_use)]
chr.thres.can50 <- final.N.all$Sub_ID[final.N.all$Diagnosis == "CHRc" & final.N.all$cannabis_lifetime_use > 50 & !is.na(final.N.all$cannabis_lifetime_use) | final.N.all$Diagnosis == "CHRnc" & final.N.all$cannabis_lifetime_use > 3 & !is.na(final.N.all$cannabis_lifetime_use)]
notused50 <- final.N.all$Sub_ID[!(final.N.all$Sub_ID %in% hc.thres.nocan) & !(final.N.all$Sub_ID %in% hc.thres.can50) & !(final.N.all$Sub_ID %in% chr.thres.nocan) & !(final.N.all$Sub_ID %in% chr.thres.can50)]

t50_data <- msn_means_pheno_F[!(msn_means_pheno_F$Sub_ID %in% notused50), ]
## FDR corrected parametric ANOVAs
###### HC vs CHRc (can) w/ FDR correction
anova.thres.50.wCovar.FDR <- map(msnROI.list, msnANOVA.wCovar, t50_data) %>%
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.thres.50.wCovar.FDR) <- c("Subgroup", "demo_age_ym", "demo_sex", "SiteNumber")
anova.thres.50.wCovar.FDR$ROI <- msnROI.list
write.xlsx(anova.thres.50.wCovar.FDR, "Anova_Threshold_50.xlsx")

## threshold 150
hc.thres.can150 <- final.N.all$Sub_ID[final.N.all$Diagnosis == "HC" & final.N.all$cannabis_lifetime_use > 150 & !is.na(final.N.all$cannabis_lifetime_use)]
chr.thres.can150 <- final.N.all$Sub_ID[final.N.all$Diagnosis == "CHRc" & final.N.all$cannabis_lifetime_use > 150 & !is.na(final.N.all$cannabis_lifetime_use) | final.N.all$Diagnosis == "CHRnc" & final.N.all$cannabis_lifetime_use > 3 & !is.na(final.N.all$cannabis_lifetime_use)]
notused150 <- final.N.all$Sub_ID[!(final.N.all$Sub_ID %in% hc.thres.nocan) & !(final.N.all$Sub_ID %in% hc.thres.can150) & !(final.N.all$Sub_ID %in% chr.thres.nocan) & !(final.N.all$Sub_ID %in% chr.thres.can150)]

t150_data <- msn_means_pheno_F[!(msn_means_pheno_F$Sub_ID %in% notused150), ]
## FDR corrected parametric ANOVAs
###### HC vs CHRc (can) w/ FDR correction
anova.thres.150.wCovar.FDR <- map(msnROI.list, msnANOVA.wCovar, t150_data) %>%
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.thres.150.wCovar.FDR) <- c("Subgroup", "demo_age_ym", "demo_sex", "SiteNumber")
anova.thres.150.wCovar.FDR$ROI <- msnROI.list
write.xlsx(anova.thres.150.wCovar.FDR, "Anova_Threshold_150.xlsx")




##### threshold post hoc analysis

#### THRESHOLD LEVEL 3 - diagnosis (FDR correction)
sr3 <- filter(anova.thres.3.wCovar.FDR, Subgroup <= .056)
sr3_names <- sr3$ROI

posthoc.three.diag.FDR <- map(sr3_names, posthoc_lm, t3_data) 
names(posthoc.three.diag.FDR) <- sr3_names 
capture.output(posthoc.three.diag.FDR, file="posthoc.three.diag")

file_content <- readLines("posthoc.three.diag")
writeLines(file_content, file.path(Cdirectory, "posthoc.three.diag"))


#### THRESHOLD LEVEL 10 - diagnosis (FDR correction)
sr10 <- filter(anova.thres.10.wCovar.FDR, Subgroup <= .056)
sr10_names <- sr10$ROI

posthoc.ten.diag.FDR <- map(sr10_names, posthoc_lm, t10_data) 
names(posthoc.ten.diag.FDR) <- sr10_names 
capture.output(posthoc.ten.diag.FDR, file="posthoc.ten.diag")

file_content <- readLines("posthoc.ten.diag")
writeLines(file_content, file.path(Cdirectory, "posthoc.ten.diag"))


#### THRESHOLD LEVEL 50 - diagnosis (FDR correction)
sr50 <- filter(anova.thres.50.wCovar.FDR, Subgroup <= .056)
sr50_names <- sr50$ROI

posthoc.fifty.diag.FDR <- map(sr50_names, posthoc_lm, t50_data) 
names(posthoc.fifty.diag.FDR) <- sr50_names 
capture.output(posthoc.fifty.diag.FDR, file="posthoc.fifty.diag")

file_content <- readLines("posthoc.fifty.diag")
writeLines(file_content, file.path(Cdirectory, "posthoc.fifty.diag"))


#### THRESHOLD LEVEL 150 - diagnosis (FDR correction)
sr150 <- filter(anova.thres.150.wCovar.FDR, Subgroup <= .056)
sr150_names <- sr150$ROI

posthoc.onefifty.diag.FDR <- map(sr150_names, posthoc_lm, t150_data) 
names(posthoc.onefifty.diag.FDR) <- sr150_names 
capture.output(posthoc.onefifty.diag.FDR, file="posthoc.onefifty.diag")

file_content <- readLines("posthoc.onefifty.diag")
writeLines(file_content, file.path(Cdirectory, "posthoc.onefifty.diag"))




#potential sex effect
# add empty column for labeling
secondFinal <- final.N.all[, c(1:72)]
secondFinal$gender.subgroup = NA
# replace subgroup title with appropriate measure
secondFinal$gender.subgroup[secondFinal$Diagnosis == "HC" & secondFinal$cannabis_ever_smoked == "1" & secondFinal$demo_sex == "1"] <- "HC.ca.m"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "HC" & secondFinal$cannabis_ever_smoked == "0" & secondFinal$demo_sex == "1"] <- "HC.noCa.m"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "CHRc" & secondFinal$cannabis_ever_smoked == "0" & secondFinal$demo_sex == "1"] <- "CHR.noCa.m"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "CHRnc" & secondFinal$cannabis_ever_smoked == "0" & secondFinal$demo_sex == "1"] <- "CHR.noCa.m"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "CHRc" & secondFinal$cannabis_ever_smoked == "1" & secondFinal$demo_sex == "1"] <- "CHR.ca.m"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "CHRnc" & secondFinal$cannabis_ever_smoked == "1" & secondFinal$demo_sex == "1"] <- "CHR.ca.m"

secondFinal$gender.subgroup[secondFinal$Diagnosis == "HC" & secondFinal$cannabis_ever_smoked == "1" & secondFinal$demo_sex == "2"] <- "HC.ca.f"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "HC" & secondFinal$cannabis_ever_smoked == "0" & secondFinal$demo_sex == "2"] <- "HC.noCa.f"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "CHRc" & secondFinal$cannabis_ever_smoked == "0" & secondFinal$demo_sex == "2"] <- "CHR.noCa.f"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "CHRnc" & secondFinal$cannabis_ever_smoked == "0" & secondFinal$demo_sex == "2"] <- "CHR.noCa.f"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "CHRc" & secondFinal$cannabis_ever_smoked == "1" & secondFinal$demo_sex == "2"] <- "CHR.ca.f"
secondFinal$gender.subgroup[secondFinal$Diagnosis == "CHRnc" & secondFinal$cannabis_ever_smoked == "1" & secondFinal$demo_sex == "2"] <- "CHR.ca.f"

# anova prep
pf2 <- as.data.frame(secondFinal)
msn_means_sex_study <- left_join(pf2, ctf, by = "Sub_ID")


msnANOVA.wCovar.sex <- function(ROI.vector, data_table) {
  ROI.vector_string <- as.name(ROI.vector)
  print(data_table$ROI.vector)
  print(ROI.vector)
  ROI.vector_tmp = data_table[[ROI.vector]]
  print(head(ROI.vector_tmp))
  model_df <- cbind(ROI.vector_tmp, dplyr::select(data_table, gender.subgroup, demo_age_ym, demo_sex, SiteNumber))
  print(head(model_df))
  msnANOVA.wCovar <- lm(ROI.vector_tmp ~ factor(gender.subgroup) + demo_age_ym + factor(demo_sex) + factor(SiteNumber), model_df) %>% anova %>% 
    tidy %>%
    filter(!is.na(p.value))
  #round p-values
  msnANOVA.wCovar$p.value
  #print(model_lm)
  #return(model_lm)
}

## FDR corrected parametric ANOVAs
###### HC vs CHRc (can) w/ FDR correction
anova.sex.wCovar.FDR <- map(msnROI.list, msnANOVA.wCovar.sex, msn_means_sex_study) %>%
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.sex.wCovar.FDR) <- c("gender.subgroup", "demo_age_ym", "demo_sex", "SiteNumber")
anova.sex.wCovar.FDR$ROI <- msnROI.list
write.xlsx(anova.sex.wCovar.FDR, "Anova_Sex_Effect.xlsx")



#### posthoc pairwise comparisons - age (FDR correction)
sig_ROI_S_FDR <- filter(anova.sex.wCovar.FDR, gender.subgroup <= .056)
sig_ROI_names_S_FDR <- sig_ROI_S_FDR$ROI

ph_sex <- function(sig_ROI_names, data_table) {
  module_tmp = data_table[[sig_ROI_names]]
  
  model_df <- cbind(module_tmp, dplyr::select(data_table, gender.subgroup))
  print(model_df)
  
  model_aov <- aov(module_tmp ~ as.factor(gender.subgroup), model_df)  %>% TukeyHSD
  print(model_aov)
  return(model_aov)
}

lm_posthoc_sex.FDR.output <- map(sig_ROI_names_S_FDR, ph_sex, msn_means_sex_study) 
names(lm_posthoc_sex.FDR.output) <- sig_ROI_names_S_FDR
capture.output(lm_posthoc_sex.FDR.output, file="posthoc.lm.output.sex.FDR")

file_content <- readLines("posthoc.lm.output.sex.FDR")
writeLines(file_content, file.path(Cdirectory, "posthoc_sex"))















## brain maps with free surfer 

# brain plots
install.packages('fsbrain') # uncomment and run if you haven't installed fsbrain
library(fsbrain)
library(ggplot2)
library(tidyr)
library(dplyr)
install.packages("freesurfer")
library(freesurfer)

# install.packages("devtools")
devtools::install_github("muschellij2/freesurfer")

######### 
path2freesurfer='/Applications/freesurfer/' # set path to freesurfer
options(freesurfer.path=path2freesurfer)

subjects_dir = '/Applications/freesurfer/7.1.1/subjects' # set path to subjects directory (e.g., if you're using fsaverageSubP this should point to where that folder lives)
set_fs_subj_dir(subjects_dir)
#options(freesurfer.subjects.dir = subjects_dir)

default_roi_path = '/Users/isabellehales/Desktop/capstonefinal/308_regions_names.txt' # set path to "308_regions_names.txt" file which should contain the 308 ROI names in the default order (i.e., starts with lh banks)
default_rois <- read.table(default_roi_path, header=F) # Load 500 aparc region names

## Read in summary file (a file containing a vector of 308 values). If you don't have a .csv file to read, and just have a loose vector of MSN values lying around, you can simply pass that vector as the "data" argument granted it is in the expected order (i.e., how it is outlined in the "region_names" variable above; see below for more information)
ROI_diff.hc.nocan.vs.can <- read.xlsx("/Users/isabellehales/Desktop/capstonefinal/Diff.Mean.Msn_HC.nocan.vs.HC.can.xlsx", colNames = FALSE)

## Below is a potentially overcomplicated function that plots 500 aparc values. The simplest way to run it is by ordering a vector of MSN values in the default order (i.e., the order outlined by the region_names file) and passing it as the "data" argument (see line 108 for an example). You should also specify the region names by passing the list of ROI names (from region_names file) as the region_names argument. By default, the ROI values will be projected onto the fsaverageSubP, the folder for which should be within the subjects_dir specified in line 12, but you can change that with the subject_id argument if desired. If you would like to live life more dangerously, you can see the other examples below for how to plot an unordered list of values as long as the values have ROI names attached to them somehow.

visualize.500.roi <- function(data, region_names = NULL, valcol=NULL, subject_id = 'fsaverageSubP', ordered = T, roicol = NULL) {
  if (is.data.frame(data)) {
    if (ordered == T) { # if ordered argument is true (default), we assume that the dataframe is ordered according to the default order of regions (i.e., denoted by region_names variable read above)
      if (is.null(region_names)) {
        print("Please make sure to include the region_names argument which contains the list of region names to superimpose onto the MSN values")
      }
      
      ## Assign region names to vector values. 
      data$ROI = region_names
      
      ## Create two vectors - one containing rh values with corresponding roi name, the other containing lh values and roi names
      rh_indices <- which(grepl('rh_', data$ROI)) # find rh rois
      rh_roivals = data[,valcol][rh_indices] # pull rh MSNs
      names(rh_roivals) = sub('rh_','',data$ROI[rh_indices]) # attach rh ROI names to rh MSNs
      
      # repeat for lh
      lh_indices <- which(grepl('lh_', data$ROI))
      lh_roivals = data[,valcol][lh_indices]
      names(lh_roivals) = sub('lh_','',data$ROI[lh_indices])
      ###########################################
    } else { ## if ordered argument is false, we assume that a column is included in the dataframe that denotes the ROI names that correspond to their values
      if (is.null(roicol)) {
        print("Please include roicol argument as a string denoting the name of the column that is home to the ROI labels")
      } else {
        ## Pull lh values and assign roi names
        lh_indices = which(grepl('lh_', data[,roicol]))
        lh_roivals = data[,valcol][lh_indices]
        names(lh_roivals) = sub('lh_','',data[,roicol][lh_indices])
        
        # Pull rh values and assign roi names
        rh_indices = which(grepl('rh_', data[,roicol]))
        rh_roivals = data[,valcol][rh_indices]
        names(rh_roivals) = sub('rh_','',data[,roicol][rh_indices])
      }
    }
    
  } else { ## If data is not a dataframe, we assume it is a list/vector!
    
    if (ordered == T) { # if ordered argument is T, we assume the provided values are in the default order and assign them names in the order specified by the region_names file
      lh_roivals = data[1:152] # select first 152 values in vector
      names(lh_roivals) = sub('lh_','',region_names[1:152]) ## assign lh names to the first 152 values in the vector
      
      rh_roivals = data[153:length(data)] ## select the last 156 values
      names(rh_roivals) = sub('rh_','',region_names[153:length(region_names)])  ## asign rh names to the last 156 values in the vector
    } else { 
      if (is.null(names(data))) {
        print('You specified that your data was not in the default order; however, your vector does not have names attached so I cannot tell which value corresponds to which ROI. Please provide the names')
      } else { 
        lh_indices = which(grepl('lh_', names(data)))
        lh_roivals = data[lh_indices]
        names(lh_roivals) = sub('lh_', '', names(data)[lh_indices])
        
        rh_indices = which(grepl('rh_', names(data)))
        rh_roivals = data[rh_indices]
        names(rh_roivals) = sub('rh_', '', names(data)[rh_indices])
      }
    }
  }
  
  
  ###########################  Time for some plotting
  ###############
  atlas = '500.aparc' # specify the atlas (i.e., this will pull annot files in format <hemi>.500.aparc.annot)
  
  # initialize color vector
  colFn_diverging = function(n) { grDevices::hcl.colors(n, palette = "Blue-Red 3"); }
  #colFn_diverging = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, name="RdBu")); # this will reverse the color scheme (so red will correspond to negative values)
  
  makecmap_options = list('colFn'=colFn_diverging, 'symm' = T)
  #makecmap_options = list('colFn'=viridis::viridis); # this will use some nice greens/blues but maybe not as interpretable/intuitive. 
  ##########################################################
  ### Adjust plot size so that you can see color bar
  ## Default plot size
  # rgloption = list("windowRect"=c(50,50,800,800));
  
  # Custom plot size to include color bars
  rgloptions=list("windowRect"=c(100, 100, 800, 800)); 
  ##########################################################
  
  # Plot function
  vis.region.values.on.subject(subjects_dir = fs_subj_dir(), subject_id = subject_id,atlas = atlas, lh_region_value_list = lh_roivals, rh_region_value_list = rh_roivals, surface = 'inflated', makecmap_options = makecmap_options, draw_colorbar = 'horizontal', rgloptions = rgloptions);
}

## Example of a call with a vector of values as the main argument
visualize.500.roi(data = HC.nocan.can.DIFF, region_names = default_rois$V1)

# ## Example of a call with a dataframe as the main argument
# visualize.500.roi(data = roi_file, region_names = region_names$V1, valcol='MSN_Value')

## Example of a call with an unordered vector of values (data) with names specified
# visualize.500.roi(data = data, ordered = F)
# 
# Example of a call with an unordered dataframe with a column specifying the ROIs
# visualize.500.roi(data = data, ordered = F, valcol = 'value',roicol = 'roi')




## data exploration
final.N.all

# age plot
ggplot(final.N.all, aes(x = demo_age_ym, fill = Subgroup)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  labs(title = "Age Distribution by Diagnosis", x = "Age", y = "Frequency") +
  theme_minimal()

# race plot
race_table <- table(final.N.all$demo_racial)
print(race_table)

race_diagnosis_table <- xtabs(~ demo_racial + Subgroup, data = final.N.all)
race_mapping <- c("First Nations", "East Asian", "Southeast Asian", "South Asian", "Black", "Central / South American", "Middle Eastern", "White", "Pacific Islander", "Interracial")
rownames(race_diagnosis_table) <- race_mapping
transposed_race <- t(race_diagnosis_table)
print(t(race_diagnosis_table))
print(transposed_race)

## cannabis use
ggplot(final.N.all, aes(x = cannabis_lifetime_use, fill = Subgroup)) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.5) +
  labs(title = "Usage by Diagnosis", x = "Cannabis Use", y = "Frequency") +
  theme_minimal()

hist(final.N.all$cannabis_lifetime_use, breaks=30, main = "Cannabis Lifetime Use", xlab = "Lifetime Count of Usage")
use_table <- table(final.N.all$cannabis_lifetime_use)
print(use_table)

hist(final.N.all$cannabis_age_first_try, breaks = 18, main = "Age of First Cannabis Use", xlab = "Age")















## calculate effect size
## find group means 
row.names(HC.c.mean) <- HC.c.mean$Sub_ID
HC.c.mean.ES <- HC.c.mean[,-1]
HC_groupmeans_c <- t(sapply(1:308, function(x) mean(HC.c.mean.ES[,x])) %>% as.matrix())
colnames(HC_groupmeans_c) <- msnROI.list

row.names(HC.nc.mean) <- HC.nc.mean$Sub_ID
HC.nc.mean.ES <- HC.nc.mean[,-1]
HC_groupmeans_NC <- t(sapply(1:308, function(x) mean(HC.nc.mean.ES[,x])) %>% as.matrix())
colnames(HC_groupmeans_NC) <- msnROI.list

row.names(CHRc.c.mean) <- CHRc.c.mean$Sub_ID
CHRc.c.mean.ES <- CHRc.c.mean[,-1]
CHRc_groupmeans_c <- t(sapply(1:308, function(x) mean(CHRc.c.mean.ES[,x])) %>% as.matrix())
colnames(CHRc_groupmeans_c) <- msnROI.list

row.names(CHRc.nc.mean) <- CHRc.nc.mean$Sub_ID
CHRc.nc.mean.ES <- CHRc.nc.mean[,-1]
CHRc_groupmeans_NC <- t(sapply(1:308, function(x) mean(CHRc.nc.mean.ES[,x])) %>% as.matrix())
colnames(CHRc_groupmeans_NC) <- msnROI.list

# standard deviation
all_msns_SD <-all_correct_msns[,-1]
row.names(all_msns_SD) <- all_correct_msns$Sub_ID
msn_groupSD <- t(sapply(1:308, function(x) SD(all_msns_SD[,x])) %>% as.matrix())

# group mean differences
msn_HC_can_nocan_GD <- HC_groupmeans_c - HC_groupmeans_NC
msn_HC_can_nocan_GD_abs<- abs(HC_groupmeans_c - HC_groupmeans_NC)

## effect size
msn_cohend_HC.can_HC.nocan <- msn_HC_can_nocan_GD/msn_groupSD
msn_cohend_HC.can_HC.nocan_abs <- msn_HC_can_nocan_GD_abs/msn_groupSD
write.xlsx(msn_cohend_HC.can_HC.nocan, "HC.can vs HC.nocan effect size")




#Calculate hedges g and cohens d using the effsize package
all_msns_only_num<-all_correct_msns[,-1]
row.names(all_msns_only_num) <- all_correct_msns$Sub_ID

subgroup_id <- data.frame(Sub_ID = c(final.N.all$Sub_ID), subgroup_title = c(final.N.all$Subgroup))
effect_data <- left_join(all_correct_msns, subgroup_id)

install.packages("effsize")
library(effsize)

# List of diagnosis groups
diagnosis_groups <- unique(effect_data$subgroup_title)

# Use combn to generate all pairs of diagnosis groups
group_pairs <- combn(diagnosis_groups, 2, simplify = TRUE)

# Calculate Cohen's d for each pair
calculate_cohen_d <- function(data, group1, group2) {
  print(group1)
  print(group2)
  x_group1 <- data[data$subgroup_title == group1, 2:309]
  x_group2 <- data[data$subgroup_title == group2, 2:309]
  # Calculate mean for each subgroup across selected columns
  mean_group1 <- colMeans(x_group1, na.rm = TRUE)
  mean_group2 <- colMeans(x_group2, na.rm = TRUE)
  print(mean_group1)
  print(mean_group2)
  # Calculate effect size
  cohen.d(mean_group1, mean_group2)
}

effect_sizes <- apply(group_pairs, 2, function(pair) {
  calculate_cohen_d(effect_data, pair[1], pair[2])
})

#NAPLS2_msn_means_combat_500aparc_allsub_pheno$dx_chr <- as.factor(NAPLS2_msn_means_combat_500aparc_allsub_pheno$dx_chr)
#NAPLS2_msn_means_combat_500aparc_allsub_pheno$dx_chrc <- as.factor(NAPLS2_msn_means_combat_500aparc_allsub_pheno$dx_chrc)
#NAPLS2_msn_means_combat_500aparc_allsub_pheno$dx_chrnc <- as.factor(NAPLS2_msn_means_combat_500aparc_allsub_pheno$dx_chrnc)

HC.c_HC.nc_hedgesg <- sapply(all_correct_msns[2:309], function(x) cohen.d(x ~ [all_correct_msns$Sub_ID %in% sub.HC.can], data= list(all_correct_msns), hedges.correction=TRUE))

HC_CHR_hedgesg_df <- as.data.frame(HC_CHR_hedgesg)
HC_CHR_hedgesg_df.t <- as.data.frame(t(HC_CHR_hedgesg_df))
HC_CHR_hedgesg_df.t <- tibble::rownames_to_column(HC_CHR_hedgesg_df.t, "ROI")
write.xlsx(HC_CHR_hedgesg_df.t, "HC_CHR_hedgesg_df.t_112223.xlsx")

HC_CHR_cohensd <- sapply(NAPLS2_msn_means_combat_500aparc_allsub_pheno[64:371], function(x) cohen.d(x ~ NAPLS2_msn_means_combat_500aparc_allsub_pheno$dx_chr, data= list(NAPLS2_msn_means_combat_500aparc_allsub_pheno)))

HC_CHR_cohensd_df <- as.data.frame(HC_CHR_cohensd)
HC_CHR_cohensd_df.t <- as.data.frame(t(HC_CHR_cohensd_df))
HC_CHR_cohensd_df.t <- tibble::rownames_to_column(HC_CHR_cohensd_df.t, "ROI")
write.xlsx(HC_CHR_cohensd_df.t, "HC_CHR_cohensd_df.t_112223.xlsx")





#final.N.all$Subgroup[final.N.all$Diagnosis == "CHRc" & final.N.all$cannabis_ever_smoked == "1"] <- "CHRc.ca"
#final.N.all$Subgroup[final.N.all$Diagnosis == "CHRc" & final.N.all$cannabis_ever_smoked == "0"] <- "CHRc.noCa"
#final.N.all$Subgroup[final.N.all$Diagnosis == "CHRnc" & final.N.all$cannabis_ever_smoked == "1"] <- "CHRnc.ca"
#final.N.all$Subgroup[final.N.all$Diagnosis == "CHRnc" & final.N.all$cannabis_ever_smoked == "0"] <- "CHRnc.noCa"

#sub.CHRnc.can <- final.N.all[final.N.all[,73] == "CHRnc.ca","Sub_ID"]
#sub.CHRnc.nocan <- final.N.all[final.N.all[,73] == "CHRnc.noCa","Sub_ID"]

#CHRnc.can.msns <- cp_MSNs[names(cp_MSNs) %in% sub.CHRnc.can]
#CHRnc.nocan.msns <- cp_MSNs[names(cp_MSNs) %in% sub.CHRnc.nocan]

#CHRnc.c.mean <- all_correct_msns[all_correct_msns$Sub_ID %in% sub.CHRnc.can, ]
#CHRnc.nc.mean <- all_correct_msns[all_correct_msns$Sub_ID %in% sub.CHRnc.nocan, ]

#msn_means_HC_CHRnc_can <- filter(msn_means_pheno_F, Subgroup == "HC.ca" | Subgroup == "CHRnc.ca") # HC vs CHR-non_con (can)
#msn_means_HC_CHRnc_nocan <- filter(msn_means_pheno_F, Subgroup == "HC.noCa" | Subgroup == "CHRnc.noCa") # HC vs CHR-non-con (no can)
#msn_means_CHRc_CHRnc_can <- filter(msn_means_pheno_F, Subgroup == "CHRc.ca" | Subgroup == "CHRnc.ca") # CHR-c vs CHR-nc (can)
#msn_means_CHRc_CHRnc_nocan <- filter(msn_means_pheno_F, Subgroup == "CHRc.noCa" | Subgroup == "CHRnc.noCa") # CHR-c vs CHR-nc (no can)

#### posthoc pairwise comparisons - diagnosis (no correction)
sig_ROI_D_noFDR <- filter(anova.HC.CHR.Can.wCovar, Subgroup <= .056)
sig_ROI_names_D_noFDR <- sig_ROI_D_noFDR$ROI

posthoc.lm.output.diag <- map(sig_ROI_names_D_noFDR, posthoc_lm, msn_means_pheno_F) 
names(posthoc.lm.output.diag) <- sig_ROI_names_D_noFDR 
capture.output(posthoc.lm.output.diag, file="posthoc.lm.output.noFDR_diag")

file_content <- readLines("posthoc.lm.output.noFDR_diag")
writeLines(file_content, file.path(Cdirectory, "posthoc_noFDR"))


#### posthoc pairwise comparisons - age (no correction)
sig_ROI_A_noFDR <- filter(anova.HC.CHR.Can.wCovar, demo_age_ym <= .056)
sig_ROI_names_A_noFDR <- sig_ROI_A_noFDR$ROI

lm_posthoc_age.output <- map(sig_ROI_names_A_noFDR, ph_age, msn_means_pheno_F) 
names(lm_posthoc_age.output) <- sig_ROI_names_A_noFDR
capture.output(lm_posthoc_age.output, file="posthoc.lm.output.age")

file_content <- readLines("posthoc.lm.output.age")
writeLines(file_content, file.path(Cdirectory, "posthoc_age"))

# post hoc sub sections
sig_thres <- 0.05

# Filter rows 
ph.age.fdr <- list() # diagnosis 
for (i in seq_along(lm_posthoc_age.FDR.output)) {
  # Access the p-value from the current element
  cp_value <- lm_posthoc_age.FDR.output[[i]]$p_adj
  print(cp_value)
  
  # Check if the p-value is less than the significance level
  if (cp_value < sig_thres) {
    # If significant, add the current element to the list of significant results
    ph.age.fdr[[paste0("result_", i)]] <- lm_posthoc_age.FDR.output[[i]]
  }
}
# View the significant results
print(ph.age.fdr)



# group mean msn data frames
msn_means_HC_CHR_can <- filter(msn_means_pheno_F, Subgroup == "HC.ca" | Subgroup == "CHR.ca") # HC vs CHR (can)
msn_means_HC_CHR_nocan <- filter(msn_means_pheno_F, Subgroup == "HC.noCa" | Subgroup == "CHR.noCa") # HC vs CHR(no can)
msn_means_HC_can_nocan <- filter(msn_means_pheno_F, Subgroup == "HC.ca" | Subgroup == "HC.noCa") # HC can vs no can
msn_means_CHR_can_nocan <- filter(msn_means_pheno_F, Subgroup == "CHR.ca" | Subgroup == "CHR.noCa") # CHR can vs no can
msn_means_HC_nocan_CHR_can <- filter(msn_means_pheno_F, Subgroup == "HC.noCa" | Subgroup == "CHR.ca") # HC no can vs CHR can
msn_means_HC_can_CHR_nocan <- filter(msn_means_pheno_F, Subgroup == "HC.ca" | Subgroup == "CHR.noCa") # HC can vs CHR no can

mv.HC.CHR.can <- msn_means_HC_CHR_can[, c(2, 74:381)]             # HC can vs CHR can
mv.HC.CHR.nocan <- msn_means_HC_CHR_nocan[, c(2, 74:381)]            # HC no can vs CHR no can
mv.HC.can.nocan <- msn_means_HC_can_nocan[, c(2, 74:381)]            # HC can vs HC no can
mv.CHR.can.nocan <- msn_means_CHR_can_nocan[, c(2, 74:381)]            # CHR can vs CHR no can
mv.HC.nocan.CHR.can <- msn_means_HC_nocan_CHR_can[, c(2, 74:381)]      # HC no can vs CHR can
mv.HC.can.CHR.nocan <- msn_means_HC_can_CHR_nocan[, c(2, 74:381)]      # HC no can vs CHR can

write.xlsx(mv.HC.CHR.can, "mean MSN HC can vs CHR can.xlsx")
write.xlsx(mv.HC.CHR.nocan, "mean MSN HC no can vs CHR no can.xlsx")
write.xlsx(mv.HC.can.nocan, "mean MSN HC can vs HC no can.xlsx")
write.xlsx(mv.CHR.can.nocan, "mean MSN CHR can vs CHR no can.xlsx")
write.xlsx(mv.HC.nocan.CHR.can, "mean MSN HC no can vs CHR can.xlsx")
write.xlsx(mv.HC.can.CHR.nocan, "mean MSN HC can vs CHR no can.xlsx")





