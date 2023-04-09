#set working directory
setwd("E:/CI_Project")

#load data
load("321.RData")


#Install Packages
install.packages("pacman")
pacman::p_load("BiocManager","glmnet", "caret","ellipsis", "vctrs", "survival", "survminer", "readr", "grpreg","testthat", "pkgload", "devtools", "ROCit", "car", "ggpubr", "gplots", "testthat", "pkgload")

suppressPackageStartupMessages(library("survival","survminer", "grpreg"))
suppressPackageStartupMessages(library(glmnet))
remotes::install_github("cran/plotmo", force = TRUE)
suppressPackageStartupMessages(library(plotmo))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(caret, compareGroups))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(ggkm))
suppressPackageStartupMessages(library(ROCit))

install.packages("varhandle")
suppressPackageStartupMessages(library("varhandle"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library("car"))
suppressPackageStartupMessages(library(gridExtra, grid, "ggpubr"))
suppressPackageStartupMessages(library(lattice, forcats, ellipsis))
suppressPackageStartupMessages(library(vctrs))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))

#import raw data
raw_data1 <- read.csv("raw_data.csv")
raw_data1 <- as.data.frame(raw_data1)

#remove missing clinical data
raw_data_cleaned <- raw_data1[!(raw_data1$ID %in% c("148369","118808","71653","84411","135595",
                                                    "110175","91388","160717","64469","132993","93848","132466","49488","60410","152780","143990","155016","155519")), ]

#import clinical data and clean
clinical_data <- read_excel("clinical_data_all.xlsx")
patient.id <- merged_df_sd$ID
clinical_data <- clinical_data[clinical_data$ID %in% patient.id, ]
clinical_data <- clinical_data[complete.cases(clinical_data$OS.event), ]

os.data <- clinical_data %>% 
  select(ID,Overall.survival..days.,OS.event)

#remove missing mask data
empty_var1 <- raw_data_cleaned$diagnostics_Versions_PyRadiomics == "" | is.na(raw_data_cleaned$diagnostics_Versions_PyRadiomics)
empty_var1_entries <- raw_data_cleaned[empty_var1, ]
empty_ID <- empty_var1_entries$ID
raw_data_cleaned <- raw_data_cleaned[!(raw_data_cleaned$ID %in% empty_ID), ]


#Filter by tissue type 
mask_fat <- grepl("fat", raw_data_cleaned$Mask, ignore.case = TRUE)
fat_filtered <- subset(raw_data_cleaned, mask_fat)

mask_st <- grepl("st", raw_data_cleaned$Mask, ignore.case = TRUE)
st_filtered <- subset(raw_data_cleaned, mask_st)

mask_iv <- grepl("iv", raw_data_cleaned$Mask, ignore.case = TRUE)
iv_filtered <- subset(raw_data_cleaned, mask_iv)

mask_necrosis <- grepl("necrosis", raw_data_cleaned$Mask, ignore.case = TRUE)
necrosis_filtered <- subset(raw_data_cleaned, mask_necrosis)

#remove unwanted columns
fat_filtered <- select(fat_filtered, -(2:33))
fat_filtered <- select(fat_filtered, -(4:5))

st_filtered <- select(st_filtered, -(2:33))
st_filtered <- select(st_filtered, -(4:5))

iv_filtered <- select(iv_filtered, -(2:33))
iv_filtered <- select(iv_filtered, -(4:5))

necrosis_filtered <- select(necrosis_filtered, -(2:33))
necrosis_filtered <- select(necrosis_filtered, -(4:5))

#merge: x = fat, y = necrosis, xx = st, yy = iv

merged_df <- fat_filtered %>%
  inner_join(necrosis_filtered, by = "ID") %>%
  inner_join(st_filtered, by = "ID") %>%
  inner_join(iv_filtered, by = "ID")

merged_df <- clinical_data %>% 
  inner_join(merged_df, by = "ID")
merged_df <- merged_df %>% 
  filter(merged_df$Histology <=3) %>% 
  select(1:3, 5, 9, 15:ncol(merged_df))

#create test and training data frames


merged_df$Stage <- as.numeric(merged_df$Stage)
set.seed(38)
train.index <- createDataPartition(merged_df$Age | merged_df$Overall.survival..days. | merged_df$Stage, p = 4/6, list = FALSE) #age, stage, survival 
train_df <- merged_df[ train.index,]
test_df <- merged_df[-train.index,]

#create covariates
covariates <- c(colnames(merged_df))
covariates <- covariates[6:ncol(merged_df)]

covariates


#standardisation
train_sd <- train_df
for (i in 6:ncol(train_df)){
  if (sd((train_df[[covariates[i-6+1]]]), na.rm = TRUE) != 0){
    train_sd[[covariates[i-6+1]]]<- scale(train_df[[covariates[i-6+1]]])
  } else{
    train_sd[[covariates[i-6+1]]] <- 0
  }
}

test_sd <- test_df
for (i in 6:ncol(test_df)){
  if (sd((test_df[[covariates[i-6+1]]]), na.rm = TRUE) != 0){
    test_sd[[covariates[i-6+1]]]<- scale(test_df[[covariates[i-6+1]]])
  } else{
    test_sd[[covariates[i-6+1]]] <- 0
  }
}

#univ_data_raw <- merge(os.data,merged_df_sd, by = "ID", all=TRUE)

#univariate cox regression
one_level_cols <- which(sapply(train_sd, function(x) length(unique(x))) == 1)
removed_one_levels <- train_sd[,-one_level_cols]
univ_covariates <- colnames(removed_one_levels)
univ_formulas <- sapply(univ_covariates,
                        function(x) as.formula(paste('Surv(Overall.survival..days., OS.event)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = removed_one_levels )})

count_vars <- function(x) {
  length(coef(x))
}

models_to_include <- which(sapply(univ_models, count_vars) == count_vars(univ_models[[1]]))
univ_results <- lapply(univ_models[models_to_include], function(x) {
  x <- summary(x)
  p.value<-signif(x$wald["pvalue"], digits=2)
  wald.test<-signif(x$wald["test"], digits=2)
  beta<-signif(x$coef[1], digits=2)
  HR <-signif(x$coef[2], digits=2)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  res<-c(beta, HR, wald.test, p.value)
  names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", "p.value")
  return(res)
})

res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
coxresult <- as.data.frame(res)
coxresult$FDR<-p.adjust(coxresult$p.value,method="fdr")
coxresult<-coxresult[order(coxresult$FDR, decreasing=F),]

filtered_coxresult <- coxresult %>% 
  select(p.value)

features <- rownames(filtered_coxresult)
filtered_coxresult <- cbind(features, filtered_coxresult)
rownames(filtered_coxresult) <- NULL
filtered_coxresult$p.value <- as.numeric(filtered_coxresult$p.value)
class(filtered_coxresult$p.value)
regular_pval <- format(filtered_coxresult$p.value, scientific = FALSE)
filtered_coxresult$p.value <- regular_pval

#remove os data
filtered_coxresult <- filtered_coxresult[!(filtered_coxresult$features %in% c("Overall.survival..days.", "OS.event","Age","Stage")), ]

#set p value threshold and filter out cox model results
p_thres <- 0.001

coxresult_significant <- filtered_coxresult %>% 
  filter(p.value<p_thres)

#write.table(coxresult, file = "all_coxresults.csv", sep = ",", row.names = TRUE)



#Lasso

train_features <- coxresult_significant$features

train_matching <- names(train_sd)[names(train_sd) %in% train_features]
train_subset <- subset(train_sd, select = train_matching)
train_raw <- train_sd %>% 
  select(ID,Overall.survival..days., OS.event,diagnostics_Mask.interpolated_VolumeNum.x)

test_subset <- subset(test_sd, select = train_matching)

train_subset <- na.omit(train_subset)
removed_rows <- attr(train_subset, "na.action")
removed_rows
train_raw <- train_raw[-c(17), ]
#train_raw <- train_raw[-c(16,67,119,142 ), ]
train_subset <- train_subset[-c(17 ), ]


nFolds <- 50
foldid <- sample(rep(seq(nFolds), length.out = nrow(train_subset)))


forlasso <- train_subset


fit <- glmnet(x=as.matrix(forlasso), y= Surv(train_raw$"Overall.survival..days.", train_raw$"OS.event"), family="cox", alpha=1, nfolds = nFolds,foldid = foldid)
plot_glmnet(fit,xvar='lambda',label=7)


cvfit <- cv.glmnet(x=as.matrix(forlasso), y= Surv(train_raw$"Overall.survival..days.", train_raw$"OS.event"), family="cox",alpha=1,  nfolds = nFolds, foldid = foldid)
plot(cvfit)

#Get cross validated R squared for goodness of fit
#y = Surv(train_raw$"Overall.survival..days.", train_raw$"OS.event")
#rsq = 1 - cvfit$cvm/var(y)
#plot(cvfit$lambda,rsq)
#cvfit$cvm
#var(y)
fit <- glmnet(x=as.matrix(forlasso), y= Surv(train_raw$"Overall.survival..days.", train_raw$"OS.event"), family="cox",alpha=1,nfolds = nFolds, lambda=cvfit$lambda.min, foldid = foldid)
fit$beta[,1]



#Predicting RPV

#training
my_prediction_model_train <- predict(fit, newx = as.matrix(train_subset[,train_features]), cvfit$lambda.min)
#testing
my_prediction_model_test <- predict(fit, newx = as.matrix(test_subset[,train_features]), cvfit$lambda.min)


# K-mean clustering

predict.kmeans <- function(object, newdata){
  centers <- object$centers
  n_centers <- nrow(centers)
  dist_mat <- as.matrix(dist(rbind(centers, newdata)))
  dist_mat <- dist_mat[-seq(n_centers), seq(n_centers)]
  max.col(-dist_mat)
}



set.seed(5)
no_partition <- 5
k_mean_lasso <- cbind(my_prediction_model_train,forlasso)
colnames(k_mean_lasso)[1] <- "RPV"
km.res <- kmeans(k_mean_lasso, no_partition, nstart = 1)

k_mean_lasso$grouping <- km.res$cluster
k_mean_lasso$grouping

class(k_mean_lasso$grouping)

k_mean_lasso$grouping <- as.numeric(k_mean_lasso$grouping)

k_mean_lasso <- cbind(k_mean_lasso,train_raw$OS.event)
k_mean_lasso <- cbind(k_mean_lasso,train_raw$Overall.survival..days.)

colnames(k_mean_lasso)[83] <- "OS.event"
colnames(k_mean_lasso)[84] <- "Overall.survival..days."






for (i in 1:nrow(k_mean_lasso)){
  if (k_mean_lasso$grouping[i] == 3) {
    k_mean_lasso$group_final[i] <- 1
  }
  else{
    k_mean_lasso$group_final[i] <- 2
  }
}

km_fit <- survfit(Surv(Overall.survival..days.,OS.event) ~ group_final, data= k_mean_lasso)



theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.border = element_blank(),
               panel.background = element_blank()
) 

a <- ggsurvplot(
  km_fit,     # survfit object with calculated statistics.
  data = k_mean_lasso,               # data used to fit survival curves. 
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,1100),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 100,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = F, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations
  # in legend of risk table
  legend.labs=c("Low Risk", "High Risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme_survminer(font.main = 12),
  title = "Kaplan-Meier Plots of Patient Groups (Training)\n Stratified Based on RPV",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  #surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent",
  font.main = c(14, "bold"),
  font.x = c(12), 
  font.y = c(12), font.tickslab = c(12),
  font.legend = c(12), risk.table.fontsize = 4
)

a


#TEST
set.seed(69)
k_mean_test <- cbind(my_prediction_model_test,test_subset)
no_partition <-  5
colnames(k_mean_test)[1] <- "RPV"
k_mean_test <- na.omit(k_mean_test)
k_test_removed_rows <- attr(k_mean_test, "na.action")
k_test_removed_rows
#is.na(k_mean_test) 
#k_mean_test <- na.omit(k_mean_test)
km.test <- kmeans(k_mean_test, no_partition, nstart = 1)
k_mean_test$grouping <- km.test$cluster
k_mean_test$grouping

k_mean_test$grouping <- as.numeric(k_mean_test$grouping)
test_raw <- test_sd %>% 
  select(ID,Overall.survival..days., OS.event,diagnostics_Mask.interpolated_VolumeNum.x)
#test_raw <- test_raw[-c(2,23,53), ]
#test_raw <- test_raw[-c(2,18,46),]
k_mean_test <- cbind(k_mean_test,test_raw$OS.event)
k_mean_test <- cbind(k_mean_test,test_raw$Overall.survival..days.)
colnames(k_mean_test)[83] <- "OS.event"
colnames(k_mean_test)[84] <- "Overall.survival..days."


for (i in 1:nrow(k_mean_test)){
  if (k_mean_test$grouping[i] == 3) {
    k_mean_test$group_final[i] <- 1
  }
  else{
    k_mean_test$group_final[i] <- 2
  }
}

k_mean_test$group_final
km_fit_test <- survfit(Surv(Overall.survival..days.,OS.event) ~ group_final, data= k_mean_test)



theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.border = element_blank(),
               panel.background = element_blank()
) 

b <- ggsurvplot(
  km_fit_test,     # survfit object with calculated statistics.
  data = k_mean_test,               # data used to fit survival curves. 
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,1100),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 100,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = F, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations
  # in legend of risk table
  legend.labs=c("Low Risk","High Risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme_survminer(font.main = 12),
  title = "Kaplan-Meier Plots of Patient Groups (Testing) \n Stratified Based on RPV",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  #surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent",
  font.main = c(14, "bold"),
  font.x = c(12), 
  font.y = c(12), font.tickslab = c(12),
  font.legend = c(12), risk.table.fontsize = 4
)

b

mean(merged_df$Age)
sd(merged_df$Age)


cate_dis <- train_raw$OS.event
cate_dis

for (i in 1:length(cate_dis)) {
  if (cate_dis[i] == "1"){
    cate_dis[i] <- '+'
  } else {
    cate_dis[i] <- '-'
  }
}



score_dis <- unname(my_prediction_model_train)



ROCit_obj <- rocit(score=score_dis, class=cate_dis, negref = "+", method ="bin")


AUC_obj <- ciAUC(ROCit_obj, level = 0.95)
p <- plot(ROCit_obj)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj$AUC, 2), ", 95% CI [", round(AUC_obj$lower, 2), ",", round(AUC_obj$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("Prediction of 3-year survival (Training Set)")


p



cate_test <- test_raw$OS.event
cate_test


for (i in 1:length(cate_test)) {
  if (cate_test[i] == "1"){
    cate_test[i] <- '+'
  } else {
    cate_test[i] <- '-'
  }
}


score_test <- unname(my_prediction_model_test)



ROCit_obj_test <- rocit(score=score_test, class=cate_test, negref = "+", method ="bin")


AUC_obj_test <- ciAUC(ROCit_obj_test, level = 0.95)
q <- plot(ROCit_obj_test)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj_test$AUC, 2), ", 95% CI [", round(AUC_obj_test$lower, 2), ",", round(AUC_obj_test$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("Prediction of 3-year survival (Testing Set)")





clinical_filtered <- merge(merged_df,clinical_data,by="ID")
clinical_filtered <- clinical_filtered[, (c(1,4174:ncol(clinical_filtered)))]

sum(clinical_filtered$Performance == 4, na.rm = TRUE)
class(clinical_filtered$Ethnicity)
sum(clinical_filtered$Ethnicity == 1)
sd(clinical_filtered$Age.y)
sum(grepl("", clinical_filtered$Stage.y))
