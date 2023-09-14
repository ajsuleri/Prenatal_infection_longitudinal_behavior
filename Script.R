#################################################################################
##########################Longitudinal behavior project##########################
#################################################################################

#Project: Prenatal infections and the development of internalizing and externalizing problems in offspring 
#Author: Anna Suleri

rm(list = ls()) #clears the environment

set.seed(2022) #set seed for reproducibility 

##########Libraries##########----
libraries <- c('foreign', 'haven', 'dplyr', 'tidyverse', 'mice', 'lme4', 'lattice', 'ggplot2', 'ggeffects', 'splines', 'reshape2', 'rmcorr', 'gridExtra', 'grid', 'ggpubr', 'tidyr', 'broom.mixed', 'xlsx', 'corrplot', 'writexl')

invisible(lapply(libraries, require, character.only = T)) #loading relevant libraries 

##########Setting working directory##########----
setwd('set_path_to_data_folder')
wd <- getwd()

##########Dataset creation##########----
###Loading datasets (exposure and demographics)
df1 <- read.spss('Covariates_MRI_analyses.sav', to.data.frame = T)
ddf1 <- dplyr::select(df1, IDC, EDUCM_3groups,EDUCP_3groups, INCOME)
ddf1$INCOME <- as.factor(ifelse(ddf1$INCOME == 'less than 450' | ddf1$INCOME == '450-600 euro' | ddf1$INCOME == '600-700 euro' | ddf1$INCOME == '700-800 euro' | ddf1$INCOME == '800-900 euro' | ddf1$INCOME == '900-1200 euro' | ddf1$INCOME == '1200-1400 euro' | ddf1$INCOME == '1400-1600 euro' | ddf1$INCOME == '1600-1800 euro' | ddf1$INCOME == '1800-2000 euro' | ddf1$INCOME == '2000-2200 euro', '<2000', ifelse(is.na(ddf1$INCOME), NA, '>2000'))) #recode to 2 categories for above and below 2000 EUR per month 
df2 <- read.spss('Second_hits_DF.sav', to.data.frame = T)
ddf2 <- dplyr::select(df2, IDM, IDC, MOTHER, AGE_M_v2,WISC13_FSIQ,APM_IQ,ETHNMv2,mdrink_updated,f1100101,SMOKE_ALL,GENDER,GSI,INTAKEPERIOD,sumscore_inf_tot,sumscore_inf_tri1,sumscore_inf_tri2,sumscore_inf_tri3, sum_ext_14, cbcl_sum_14,sum_int_14,AGECHILD_GR1093)
ddf2$age_child_cbcl14 <- ddf2$AGECHILD_GR1093 #select vars of interest
ddf2 <- dplyr::select(ddf2, -c('AGECHILD_GR1093')) #drop vars we don't need 

###Loading datasets (behavior data) and renaming repeated behavior variables
df3 <- read.spss("20200111_GR1029_CBCL A1_1.5_incl_Tscores.sav", to.data.frame = T) #cbcl 1.5
ddf3 <- dplyr::select(df3, IDC,AGE18M, cbcl_sum,sum_int,sum_ext)
ddf3$age_child_cbcl1.5 <- ddf3$AGE18M
ddf3$cbcl_sum_1.5 <- ddf3$cbcl_sum
ddf3$sum_int_1.5 <- ddf3$sum_int
ddf3$sum_ext_1.5 <- ddf3$sum_ext
ddf3 <- dplyr::select(ddf3, -c('AGE18M', 'cbcl_sum', 'sum_int', 'sum_ext'))
df4 <- read.spss("CBCL_3_incl_Tscores__GR1065E2_GR1066A1_20201111.sav", to.data.frame = T) #cbcl 3
ddf4 <- dplyr::select(df4, IDC,age_GR1065,cbcl_sum_36m,sum_int_36m,sum_ext_36m)
ddf4$age_child_cbcl3 <- ddf4$age_GR1065
ddf4$cbcl_sum_3 <- ddf4$cbcl_sum_36m
ddf4$sum_int_3 <- ddf4$sum_int_36m
ddf4$sum_ext_3 <- ddf4$sum_ext_36m
ddf4 <- dplyr::select(ddf4, -c('age_GR1065', 'cbcl_sum_36m', 'sum_int_36m', 'sum_ext_36m'))
df5 <- read.spss("CHILDCBCL_6_incl_Tscores_20201111.sav", to.data.frame = T) #cbcl 6
ddf5 <- dplyr::select(df5, IDC, agechild_GR1075,cbcl_sum_5, sum_int_5,sum_ext_5)
ddf5$age_child_cbcl6 <- ddf5$agechild_GR1075
ddf5$cbcl_sum_6<- ddf5$cbcl_sum_5
ddf5$sum_int_6<- ddf5$sum_int_5
ddf5$sum_ext_6<- ddf5$sum_ext_5
ddf5 <- dplyr::select(ddf5, -c('agechild_GR1075', 'cbcl_sum_5', 'sum_int_5', 'sum_ext_5'))
df6 <- read.spss("CHILDCBCL9_incl_Tscores_20201111.sav", to.data.frame = T) #cbcl 9
ddf6 <- dplyr::select(df6, IDC, AgeChild_CBCL9m,cbcl_sum_9m,sum_int_9m,sum_ext_9m)
ddf6$age_child_cbcl9 <- ddf6$AgeChild_CBCL9m
ddf6$cbcl_sum_9<- ddf6$cbcl_sum_9m
ddf6$sum_int_9 <- ddf6$sum_int_9m
ddf6$sum_ext_9 <- ddf6$sum_ext_9m
ddf6 <- dplyr::select(ddf6, -c('AgeChild_CBCL9m', 'cbcl_sum_9m', 'sum_int_9m', 'sum_ext_9m'))

###Loading in PRS and residualizing raw prs
prs <- read.delim("prs_data_file.txt")

# Using the residual score of PRS as covariates, instead of raw PRS, so that you don't need to adjust for PCA in your main prediction models

# 1) For PRS neurodev3.Pt_1
# regress PRS on covariates to get residuals as an adjusted outcome
fit1 <- lm(prs$neurodev3.Pt_1 ~ C1+ C2+ C3+ C4 + C5, data = prs)  

# save and create a new column for the residual score of PRS neurodev3.Pt_1 in your data 
prs$PRS_neurodev_resid <- fit1$residuals

# 2) For PRS moodpsych3.Pt_1 
fit2 <- lm(prs$moodpsych3.Pt_1 ~ C1+ C2+ C3+ C4 + C5, data = prs)  
prs$PRS_moodpsych_resid <- fit2$residuals

# 3) For PRS compuls3.Pt_1
fit3 <- lm(prs$compuls3.Pt_1 ~ C1+ C2+ C3+ C4 + C5, data = prs)   
prs$PRS_compuls_resid <- fit3$residuals

str(prs)

prs_final <- dplyr::select(prs, c('IDC', 'PRS_neurodev_resid', 'PRS_moodpsych_resid', 'PRS_compuls_resid'))

###Loading additional covariates of interest after journal revisions 
secondhits_df <- read.spss('Second_hits_DF.sav', to.data.frame = T)
prenatal_medication <- read.spss('MEDICATIONSELFREPORTPREGNANCY_30112017.sav', to.data.frame = T)
preg_inflam_df1 <- read.spss('MATERNALCOMPLICATIONS_22112016.sav', to.data.frame = T)
preg_inflam_df2 <- read.spss('PARTUS_17072015.sav', to.data.frame = T)
inflam_somatic_df <- read.spss('GR1001-D1-37_16072021.sav', to.data.frame = T) 

# construct scores and select variables of interest (birth complications, pregnancy inflammatory conditions and inflammatory medical conditions including prenatal medication use)
birth_df <- dplyr::select(secondhits_df, 'IDC' ,'WEIGHT', 'GESTBIR', 'UMPI2', 'UtRI2', 'PLAWGHT', 'PLGF_g1','sumscore_childhood_infection')
birth_compl_df <- dplyr::select(preg_inflam_df2, c("IDM", "APGAR5", "PH_NAV"))

prenatal_medication <- dplyr::select(prenatal_medication, c('IDM', 'SSRITOT', 'TRIPTOT', 'PSYTOT', 'TCATOT', 'NSAIDTOT', 'ABIOTOT', 'PMOLTOT', 'CORTTOT', 'MUCOTOT' ,'COUGHTOT')) #selecting medication use of interest 
prenatal_medication$psymed_usage <- ifelse(prenatal_medication$SSRITOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$SSRITOT == 'Early pregnancy only' | prenatal_medication$TRIPTOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$TRIPTOT == 'Early pregnancy only' | prenatal_medication$PSYTOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$PSYTOT == 'Early pregnancy only' | prenatal_medication$TCATOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$TCATOT == 'Early pregnancy only', 1, ifelse(prenatal_medication$SSRITOT == 'no use' | prenatal_medication$TRIPTOT == 'no use' | prenatal_medication$PSYTOT == 'no use' | prenatal_medication$TCATOT == 'no use', 0, NA)) #sum score psy med usage
prenatal_medication$inflam_med_usage <- ifelse(prenatal_medication$NSAIDTOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$NSAIDTOT == 'Early pregnancy only' | prenatal_medication$ABIOTOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$ABIOTOT == 'Early pregnancy only' | prenatal_medication$PMOLTOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$PMOLTOT == 'Early pregnancy only' | prenatal_medication$CORTTOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$MUCOTOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$MUCOTOT == 'Early pregnancy only' | prenatal_medication$COUGHTOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$COUGHTOT == 'Early pregnancy only', 1, ifelse(prenatal_medication$NSAIDTOT == 'no use' | prenatal_medication$ABIOTOT == 'no use' | prenatal_medication$PMOLTOT == 'no use' | prenatal_medication$MUCOTOT == 'no use' | prenatal_medication$COUGHTOT == 'no use', 0, NA)) #sum score anti inflamm or infect med usage
prenatal_medication$cortico_med_usage <- ifelse(prenatal_medication$CORTTOT == 'During pregnancy-first-second-third trimester' | prenatal_medication$CORTTOT == 'Early pregnancy only', 1, ifelse(prenatal_medication$CORTTOT == 'no use', 0, NA)) #sum score cortico med usage

prenatal_med_df <- dplyr::select(prenatal_medication, c('IDM', 'psymed_usage', 'inflam_med_usage', 'cortico_med_usage')) #final df for prenatal medications 

preg_inflam_df1$HELLP <- ifelse(preg_inflam_df1$PE_total == 'HELLP', 1, ifelse(preg_inflam_df1$PE_total == "NA's", NA, 0))
preg_inflam_df1$gest_diab <- ifelse(preg_inflam_df1$DIAB_GRA == 'Yes', 1, ifelse(preg_inflam_df1$DIAB_GRA == 'No', 0, NA))
preg_inflam_df1$preg_hypt <- ifelse(preg_inflam_df1$PIH_v1 == 'Yes', 1, ifelse(preg_inflam_df1$PIH_v1 == 'No', 0, NA))
preg_inflam_df1$preclampsia <- ifelse(preg_inflam_df1$PE == 'Yes', 1, ifelse(preg_inflam_df1$PE == 'No', 0, NA))

preg_inflam_df2$promm <- ifelse(preg_inflam_df2$PROMM == 'Yes', 1, ifelse(preg_inflam_df2$PROMM == 'No', 0, NA))
preg_inflam_df2$sectio <- ifelse(preg_inflam_df2$UITDRIJF == 'prim. sectio' | preg_inflam_df2$UITDRIJF == 'sec. sectio' | preg_inflam_df2$UITDRIJF == 'sectio (onbekend prim./sec.)', 1, ifelse(preg_inflam_df2$UITDRIJF == "NA's", NA, 0))

preg_inflam <- merge(preg_inflam_df1, preg_inflam_df2, by = c('IDM'), all = T)
preg_inflam$preg_inflam_cond_score <- preg_inflam$HELLP + preg_inflam$gest_diab + preg_inflam$preg_hypt + preg_inflam$preclampsia 

preg_inflam$preg_inflam_cond_score_bin <- ifelse(preg_inflam$preg_inflam_cond_score == 1 | preg_inflam$preg_inflam_cond_score > 1, 1, 0) #making  binary score because of low power in sum score

preg_inflam_df <- dplyr::select(preg_inflam, c('IDM', 'preg_inflam_cond_score_bin', 'sectio', 'promm')) #final data df 

inflam_somatic_df$intestinal <- ifelse(inflam_somatic_df$d2100101 == 'Yes', 1, ifelse(inflam_somatic_df$d2100101 == "NA's", NA, 0))
inflam_somatic_df$sle <- ifelse(inflam_somatic_df$d2300101 == 'Yes', 1, ifelse(inflam_somatic_df$d2300101 == "NA's", NA, 0))
inflam_somatic_df$arthritis <- ifelse(inflam_somatic_df$d2400101 == 'Yes', 1, ifelse(inflam_somatic_df$d2400101 == "NA's", NA, 0))
inflam_somatic_df$ms <- ifelse(inflam_somatic_df$d2500101 == 'Yes', 1, ifelse(inflam_somatic_df$d2500101 == "NA's", NA, 0))
inflam_somatic_df$thyroid <- ifelse(inflam_somatic_df$d2600101 == 'Yes', 1, ifelse(inflam_somatic_df$d2600101 == "NA's", NA, 0))
inflam_somatic_df$diabetes <- ifelse(inflam_somatic_df$d1100101 == 'Yes', 1, ifelse(inflam_somatic_df$d1100101 == "NA's", NA, 0))
inflam_somatic_df$inflam_conditions_score <- inflam_somatic_df$intestinal + inflam_somatic_df$sle + inflam_somatic_df$arthritis + inflam_somatic_df$ms + inflam_somatic_df$thyroid + inflam_somatic_df$diabetes
inflam_somatic_df$inflam_conditions_score_bin <- ifelse(inflam_somatic_df$inflam_conditions_score == 1 | inflam_somatic_df$inflam_conditions_score > 1, 1, 0) #making new binary score because of low power in sum score

inflam_somatic_df_final <- dplyr::select(inflam_somatic_df, c('IDM', 'inflam_conditions_score_bin'))

####Merging datasets into one dataframe
df_merge1 <- merge(ddf1, ddf2, by = 'IDC', all.x = TRUE)
df_merge2 <- merge(df_merge1, ddf3, by = 'IDC', all.x = TRUE)
df_merge3 <- merge(df_merge2, ddf4, by = 'IDC', all.x = TRUE)
df_merge4 <- merge(df_merge3, ddf5, by = 'IDC', all.x = TRUE)
df_merge5 <- merge(df_merge4, ddf6, by = 'IDC', all.x = TRUE)
df_merge6 <- merge(df_merge5, prs_final, by = 'IDC', all.x = TRUE)

df_merge7 <- merge(df_merge6, birth_df, by = 'IDC', all.x =TRUE)
df_merge8 <- merge(df_merge7, birth_compl_df, by = 'IDM', all.x = TRUE)
df_merge9 <- merge(df_merge8, prenatal_med_df, by = 'IDM', all.x = TRUE)
df_merge10 <- merge(df_merge9, preg_inflam_df, by = 'IDM', all.x = TRUE)
df_merge11 <- merge(df_merge10, inflam_somatic_df_final, by = 'IDM', all.x = TRUE)

##########Final dataset##########----
###Inclusion and exclusion criteria----
inclu1 <- subset(df_merge11, INTAKEPERIOD == '<18 weeks') #including only mothers who enrolled in tri1, n=7145
inclu2 <- subset(inclu1, complete.cases(sumscore_inf_tot)) #no missing information in prenatal infection, n=4259
inclu3 <- subset(inclu2, complete.cases(cbcl_sum_1.5) | complete.cases(cbcl_sum_3) | complete.cases(cbcl_sum_6) | complete.cases(cbcl_sum_9) | complete.cases(cbcl_sum_14) | complete.cases(sum_int_3) | complete.cases(sum_int_14) & complete.cases(sum_ext_14) | complete.cases(sum_int_6) | complete.cases(sum_ext_3)) #behavior outcome at one of the timepoints available

#'Excluding 1 twin/sibling based on most available information and if equal then at random 
exclu1 <- inclu3[sample(nrow(inclu3)),] #making a random order in the df
exclu1$na_count <- apply(exclu1, 1, function(x) sum(is.na(x))) #add new column to df with na count per row/IDC
exclu2 <- exclu1[order(exclu1$na_count),] #order dataset from least to most missing
exclu3 <- exclu2[!duplicated(exclu2$MOTHER, fromLast = T),] #delete 1 twin/sibling, deleting when a duplicate is found (so second one) which also has the most missings

df_final <- dplyr::select(exclu3, -c('INTAKEPERIOD', 'na_count')) #n=3598 is final sample size

###Baseline table----
setwd('set_path_to_results_folder')

baselinevars <- c("age_child_cbcl1.5","age_child_cbcl3","age_child_cbcl6","age_child_cbcl9","age_child_cbcl14","AGE_M_v2", "GENDER", "ETHNMv2","EDUCM_3groups", "EDUCP_3groups","INCOME","f1100101","SMOKE_ALL","mdrink_updated", "GSI", "sumscore_inf_tot","sumscore_inf_tri1","sumscore_inf_tri2","sumscore_inf_tri3","cbcl_sum_1.5","cbcl_sum_3","cbcl_sum_6","cbcl_sum_9","cbcl_sum_14",'WEIGHT', 'GESTBIR', 'UMPI2', 'UtRI2', 'PLAWGHT', 'PLGF_g1', 'sumscore_childhood_infection', 'APGAR5', 'PH_NAV', 'psymed_usage', 'inflam_med_usage', 'cortico_med_usage', 'preg_inflam_cond_score_bin', 'sectio', 'promm', 'inflam_conditions_score_bin')
df_final$age_child_cbcl1.5 <- df_final$age_child_cbcl1.5/12 #transform scale from months to years 
df_final$age_child_cbcl3 <- df_final$age_child_cbcl3/12 #transform scale from months to years
df_final$age_child_cbcl6 <- df_final$age_child_cbcl6/12 #transform scale from months to years

df_final <- df_final %>% mutate_at(c('psymed_usage', 'inflam_med_usage', 'cortico_med_usage', 'preg_inflam_cond_score_bin', 'sectio', 'promm', 'inflam_conditions_score_bin'), as.factor)

baseline_table_function <- function(baselinevars, df){ #feed baselinevars vector and dataframe 
  
  #create empty dataframe to store results 
  results <- data.frame(Variable = character(0), Output = character(0))
  
  #for each baseline variable as specified in baselinevars decide whether continous or categorical and give output accordingly 
  for(i in baselinevars){ 
    
    #x = i vars that are columns in the dataframe df
    x <- df[, i] 
    
    #show column name as heading per output 
    message(i) 
    
    #function for continuous variables
    summary_continuous <- function(x){
      standev <- sd(x, na.rm = T)
      meanvar <- mean(x, na.rm = T)
      print(paste0(round(meanvar, 1), '(', round(standev, 1), ')'))
    }
    
    #function for categorical variables 
    summary_categorical <- function(x){
      tab1 <- prop.table(table(x, useNA = 'always'))
      tab2 <- table(x, useNA = "always")
      print(paste(round(tab1 * 100, 1), '%', names(tab1), collapse = ','))
      print(paste(tab2, names(tab2)))
    }
    
    #if else to apply correct function for vars type 
    if (class(x) == 'numeric') {
      output <- summary_continuous(x)
    } 
    else 
    {
      output <- summary_categorical(x) 
    }
    results <- rbind(results, data.frame(Variable = i, Output = output))
  }
  write_xlsx(results, 'baseline_table_results.xlsx') #load writexl library at start of script
}

baseline_table_function(baselinevars, df_final)

###Figure: collection of all data----
#figure prenatal vars
variables_prenatal <- c("Maternal education", "Paternal education","Household income","Maternal age at enrollment","Maternal national background","Maternal alcohol consumption","Maternal substance usage","Maternal tobacco consumption","Child sex","Maternal psychopathology","Prenatal infection total sum score","Prenatal infection trimester 1 sum score","Prenatal infection trimester 2 sum score","Prenatal infection trimester 3 sum score", "Polygenic risk scores for psychiatric conditions (3x)", "Placental vascular resistance index", "Placental growth factor", "Premature ruptured membranes", "Maternal psychotropic medication use", "Maternal anti infection or inflammation medication use", "Maternal corticosteroid use", "Pregnancy inflammatory conditions", "Inflammatory medical conditions")
time1 <- c(40, 40, 40, 0, 40, 40, 40, 40, 40, 40, 12, 12, 26, 40, 0, 18, 18, 40, 18, 18, 18, 40, 18)
time2 <- c(rep(NA, 10), 26, rep(NA, 12))
time3<- c(rep(NA, 10), 40, rep(NA, 12))

IDvar <- seq(1:length(variables_prenatal))

prenatal_df <- data.frame(variables_prenatal, time1,
                          time2, time3, IDvar)
prenatal_df_long <- reshape(prenatal_df, 
                            idvar = 'IDvar', 
                            varying = c("time1", 
                                        "time2", 
                                        "time3"),
                            v.names = c("time"),
                            direction ='long')
prenatal_df_long <- prenatal_df_long[order(prenatal_df_long$IDvar),]

prenatal_df_long[, "Measurement"] <- rep(1:3, length.out = nrow(prenatal_df_long)) %>%
  factor()

temp1 <- ggplot(prenatal_df_long, aes(variables_prenatal, time)) + 
  geom_point(aes(group = IDvar,
                 color = Measurement, size = 1.5),
             na.rm = TRUE) + 
  coord_flip() + 
  labs(y= 'Time (weeks gestation)', 
       x = "", 
       title = 'Prenatal period') + theme(axis.text.y=element_text(size=12))

temp1

#figure postnatal vars
variables_postnatal <- c("CBCL total behavioral problems", "CBCL internalizing problems", "CBCL externalizing problems", "Child's age","Maternal IQ","Child IQ", "Birth weight", "Gestational age at birth", "Placental weight at birth", "5-minute Apgar score", "pH Umbilical cord", "Caesarian delivery", "Childhood infections")
time1 <- c(1.5, 1.5, 1.5, 1.5, 6, 13, 40, 40, 40, 40, 40, 40, 0.2)
time2 <- c(3, 3, 3, 3, NA, NA,NA,NA,NA,NA,NA,NA,0.5)
time3 <- c(6, 6, 6, 6, NA, NA, NA, NA, NA, NA, NA, NA, 1)
time4 <- c(10, 10, 10, 10, NA, NA, NA, NA, NA, NA, NA, NA,2)
time5 <- c(14, 14, 14, 14, NA, NA, NA, NA, NA, NA, NA, NA,3)
time6 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 4)
time7 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 5)
time8 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9)

id <- seq(1:length(variables_postnatal))

postnatal_df <- data.frame(id, variables_postnatal, time1,
                          time2, time3, time4, time5, time6, time7, time8)

postnatal_df_long <- reshape(postnatal_df, 
                            idvar = 'id', 
                            varying = c("time1", 
                                        "time2", 
                                        "time3",
                                        "time4",
                                        "time5",
                                        'time6',
                                        'time7',
                                        'time8'),
                            v.names = c("time"),
                            direction ='long')

postnatal_df_long <- postnatal_df_long[order(postnatal_df_long$id),]

postnatal_df_long[, "Measurement"] <- rep(1:5, length.out = nrow(postnatal_df_long)) %>%
  factor()

temp2 <- ggplot(postnatal_df_long, aes(variables_postnatal, time)) + 
  geom_point(aes(group = id,
                 color = Measurement, size = 1.5),
             na.rm = TRUE) + 
  coord_flip() + 
  labs(y= 'Time (child age in years)', 
       x = "", 
       title = 'Postnatal period') + theme(axis.text.y=element_text(size=12))


#combining plots
fig <- ggarrange(temp1, temp2, 
                     labels =c('A', 'B'),
                     ncol = 2, nrow =1)

fig 

###Non-response analysis----
setwd('set_path_to_results_folder') 

nr_vars <- c("INCOME","ETHNMv2","EDUCM_3groups","mdrink_updated","f1100101", "SMOKE_ALL", "AGE_M_v2")

sink('Results_nonresponse_analysis.txt') #function to save output as txt file in wd 

for (i in nr_vars){
  #define df with excluded sample (so who are in genR but not in study sample)
  df_NR1 <- df_merge5[!(df_merge5$IDC %in% df_final$IDC),]
  #for every loop call name i 
  cat(i)
  #if cat vars then chi square test else t test 
  if (class(df_final[,i]) == 'factor') {
    print(chisq.test(cbind(table(df_NR1[,i]), table(df_final[,i]))))
  } 
  else 
  {
    print(t.test(df_NR1[,i], df_final[,i], paired = F))
  }
}

sink()

setwd('set_path_to_data_folder') 

###Transformations----
#'Before transforming, calculating what 1 SD means for each behavior scale
outcomevars <- c("sum_ext_14","sum_int_14","cbcl_sum_14","cbcl_sum_1.5","sum_int_1.5","sum_ext_1.5","cbcl_sum_3","sum_int_3","sum_ext_3","cbcl_sum_6","sum_int_6","sum_ext_6","cbcl_sum_9","sum_ext_9", 'sum_int_9')

sink('results_sd_transforming.txt')

for(i in outcomevars){
  df_final[, i]
  cat(i)
  print(sd(df_final[,i], na.rm = T))
}

sink()

#' Trans 1: sqrt outcome 
sqrt_vars <- c("sum_ext_14","sum_int_14","cbcl_sum_14","cbcl_sum_1.5","sum_int_1.5","sum_ext_1.5","cbcl_sum_3","sum_int_3","sum_ext_3","cbcl_sum_6","sum_int_6","sum_ext_6","cbcl_sum_9","sum_ext_9", 'sum_int_9')

for (i in sqrt_vars){ #function to sqrt multiple variables 
  k <- df_final[i]
  colname <- paste0(colnames(k), "_sqrt")
  df_final[, colname] <- sqrt(k)
}

#' Trans 2: z-scores outcome 
outcomes <- c("sum_ext_1.5_sqrt","sum_int_1.5_sqrt","cbcl_sum_1.5_sqrt","sum_ext_3_sqrt","sum_int_3_sqrt","cbcl_sum_3_sqrt","sum_ext_6_sqrt","sum_int_6_sqrt","cbcl_sum_6_sqrt","sum_ext_9_sqrt","sum_int_9_sqrt","cbcl_sum_9_sqrt","sum_ext_14_sqrt","sum_int_14_sqrt","cbcl_sum_14_sqrt") #vector with all outcomes 

for (x in outcomes){ #function to standardize all continuous exposure and outcome variables at once
  t <- df_final[x]
  colname <- paste0(colnames(t), "_standardized")
  df_final[,colname] <- as.numeric(scale(t))
}

#'Make correlation plot for all the outcomes 
corr_vars <- df_final %>% rename("T1. CBCL TOT" = cbcl_sum_1.5, "T2. CBCL TOT" = cbcl_sum_3,"T3. CBCL TOT" = cbcl_sum_6,"T4. CBCL TOT" = cbcl_sum_9,"T5. CBCL TOT" = cbcl_sum_14, "T1. CBCL INT" =sum_int_1.5,  "T2. CBCL INT" =sum_int_3, "T3. CBCL INT" =sum_int_6, "T4. CBCL INT" =sum_int_9, "T5. CBCL INT" =sum_int_14,"T1. CBCL EXT" = sum_ext_1.5,"T2. CBCL EXT" = sum_ext_3,"T3. CBCL EXT" = sum_ext_6,"T4. CBCL EXT" = sum_ext_9,"T5. CBCL EXT" = sum_ext_14, "Maternal education" = EDUCM_3groups, "Paternal education" = EDUCP_3groups, "Household income" = INCOME, "Maternal age"=AGE_M_v2, "Child IQ"= WISC13_FSIQ, "Maternal IQ" = APM_IQ, "Maternal ethnicity" = ETHNMv2, "Maternal alcohol use" = mdrink_updated, 'Maternal drug use' = f1100101, 'Maternal tobacco use' = SMOKE_ALL, "Child sex" = GENDER, 'Maternal psychopathology' = GSI, 'Prenatal infection total' = sumscore_inf_tot, 'Prenatal infection tri1' = sumscore_inf_tri1, 'Prenatal infection tri2' = sumscore_inf_tri2, 'Prenatal infection tri3' = sumscore_inf_tri3, 'Age T1' = age_child_cbcl1.5, 'Age T2' = age_child_cbcl3, 'Age T3' = age_child_cbcl6, 'Age T4' = age_child_cbcl9, 'Age T5' = age_child_cbcl14, 'PRS Neurodevelopmental' = PRS_neurodev_resid, 'PRS Compulsive' = PRS_compuls_resid, 'PRS Mood-Psychotic' = PRS_moodpsych_resid, 'Birth weight' = WEIGHT, 'Gestational age at birth' = GESTBIR, 'A. umbiliclis pulsatility index' = UMPI2, 'A. uterine resistance index' = UtRI2, 'Placental weight at birth' = PLAWGHT, 'Placental growth factor' = PLGF_g1, '5-minute Apgar score' = APGAR5, 'pH umbilical cord' = PH_NAV, 'Maternal psychotropic medication use' = psymed_usage, 'Maternal anti infection or inflammation medication use' = inflam_med_usage, 'Maternal corticosteroid use' = cortico_med_usage, 'Pregnancy inflammatory conditions score' = preg_inflam_cond_score_bin, 'Caesarian delivery' = sectio, 'Premature ruptured membranes' = promm, 'Inflammatory medical conditions score' = inflam_conditions_score_bin)

corr_vars2 <- dplyr::select(corr_vars, -c(1:2, 6, 59:88))

#recode cat vars to make them numeric for corr plot 
corr_vars2$"Child sex" <- ifelse(corr_vars2$"Child sex" == 'boy', 0, 1) %>% as.numeric
corr_vars2$"Household income" <- ifelse(corr_vars2$"Household income" == '<2000', 0, 1) %>% as.numeric
corr_vars2$"Maternal education" <- ifelse(corr_vars2$'Maternal education' == 'low', 0, 1) %>% as.numeric
corr_vars2$"Paternal education" <- ifelse(corr_vars2$'Paternal education' == 'low', 0, 1) %>% as.numeric
corr_vars2$"Maternal ethnicity" <- ifelse(corr_vars2$'Maternal ethnicity' == 'Not Dutch', 0, 1) %>% as.numeric
corr_vars2$"Maternal alcohol use" <- ifelse(corr_vars2$`Maternal alcohol use` == 'mother never drank in pregnancy', 0, 1) %>% as.numeric
corr_vars2$"Maternal drug use" <- ifelse(corr_vars2$`Maternal drug use` == 'No', 0, 1) %>% as.numeric
corr_vars2$"Maternal tobacco use" <- ifelse(corr_vars2$`Maternal tobacco use` == 'never smoked during pregnancy', 0, 1) %>% as.numeric

correlation <- cor(corr_vars2, use="pairwise.complete.obs")

#write.xlsx(correlation, 'corr_table.xlsx')

corrplot::corrplot(correlation, method = 'color', order = 'FPC', type = 'lower',diag = F, tl.col = 'black', tl.cex = 0.5)

####Wide to long format----
#'This step is for the df to be in the right format for the linear mixed model 
#'Recoding name age and behavior longitudinal to better names to identify different time points
df_final_long <- dplyr::select(df_final, -c(20:22, 25:27, 29:31,33:35,37:39, 59:73)) #drops cols we dont need anymore 
df_final_long['age_child.t1'] <- df_final['age_child_cbcl1.5'] 
df_final_long['age_child.t2'] <- df_final['age_child_cbcl3']
df_final_long['age_child.t3'] <- df_final['age_child_cbcl6']
df_final_long['age_child.t4'] <- df_final['age_child_cbcl9']
df_final_long['age_child.t5'] <- df_final['age_child_cbcl14']
df_final_long['cbcl_sum_sqrt_standardized.t1'] <- df_final['cbcl_sum_1.5_sqrt_standardized']
df_final_long['cbcl_sum_sqrt_standardized.t2'] <- df_final['cbcl_sum_3_sqrt_standardized']
df_final_long['cbcl_sum_sqrt_standardized.t3'] <- df_final['cbcl_sum_6_sqrt_standardized']
df_final_long['cbcl_sum_sqrt_standardized.t4'] <- df_final['cbcl_sum_9_sqrt_standardized']
df_final_long['cbcl_sum_sqrt_standardized.t5'] <- df_final['cbcl_sum_14_sqrt_standardized']
df_final_long['cbcl_int_sqrt_standardized.t1'] <- df_final['sum_int_1.5_sqrt_standardized']
df_final_long['cbcl_int_sqrt_standardized.t2'] <- df_final['sum_int_3_sqrt_standardized']
df_final_long['cbcl_int_sqrt_standardized.t3'] <- df_final['sum_int_6_sqrt_standardized']
df_final_long['cbcl_int_sqrt_standardized.t4'] <- df_final['sum_int_9_sqrt_standardized']
df_final_long['cbcl_int_sqrt_standardized.t5'] <- df_final['sum_int_14_sqrt_standardized']
df_final_long['cbcl_ext_sqrt_standardized.t1'] <- df_final['sum_ext_1.5_sqrt_standardized']
df_final_long['cbcl_ext_sqrt_standardized.t2'] <- df_final['sum_ext_3_sqrt_standardized']
df_final_long['cbcl_ext_sqrt_standardized.t3'] <- df_final['sum_ext_6_sqrt_standardized']
df_final_long['cbcl_ext_sqrt_standardized.t4'] <- df_final['sum_ext_9_sqrt_standardized']
df_final_long['cbcl_ext_sqrt_standardized.t5'] <- df_final['sum_ext_14_sqrt_standardized']

#'Reshaping  df to long format for the timevarying covariates age and outcome, so instead of having multiple columns for y1-y5 (outcome) and t1-t5 (age), we want 2 columns (outcome and time/age) in which the different outcome values per time point are stacked on each other
reshape_df <- reshape(df_final_long, idvar = "IDC", varying = list(c("age_child.t1", "age_child.t2", "age_child.t3", "age_child.t4", "age_child.t5"),c("cbcl_sum_sqrt_standardized.t1", "cbcl_sum_sqrt_standardized.t2", "cbcl_sum_sqrt_standardized.t3", "cbcl_sum_sqrt_standardized.t4", "cbcl_sum_sqrt_standardized.t5"),  c("cbcl_int_sqrt_standardized.t1", "cbcl_int_sqrt_standardized.t2", "cbcl_int_sqrt_standardized.t3", "cbcl_int_sqrt_standardized.t4", "cbcl_int_sqrt_standardized.t5"), c("cbcl_ext_sqrt_standardized.t1", "cbcl_ext_sqrt_standardized.t2", "cbcl_ext_sqrt_standardized.t3", "cbcl_ext_sqrt_standardized.t4", "cbcl_ext_sqrt_standardized.t5")), v.names = c("time", "cbcl_total_score", "cbcl_int_score", "cbcl_ext_score"), direction = 'long',) #this df is ordered by time point 
reshape_df_ordered <- reshape_df[order(reshape_df$IDC),] #order DF by IDC 

write.csv(reshape_df_ordered, 'reshape_df_ordered.csv')

####Multiple imputation----
#' Impute missing covariates, first have an overall glance at missing covariates
missvalues <- cbind("# NA" = sort(colSums(is.na(reshape_df_ordered))),
                     "% NA" = round(sort(colMeans(is.na(reshape_df_ordered))) * 100, 2))
missvalues

#Running setup imputation run 
imp0 <- mice(reshape_df_ordered, maxit = 0, defaultMethod = c("norm", "logreg", "polyreg", "polr"))

imp0$loggedEvents

#Imputation method matrix
meth <- imp0$method

#'Variables that should not be imputed are set to ""
meth[c(1:2,6, 25:27)] <- "" 

#Predictor matrix
pred <- imp0$predictorMatrix

#'Variables that should not be used as predictors set to 0 
pred[, c(1:2,6, 25:27)] <- 0 

#Visit sequence
visSeq <- imp0$visitSequence

#Performing the imputation
#'30 iterations (generally recommended) with 30 datasets for reaching convergence
imp.test <- mice(reshape_df_ordered, method = meth, predictorMatrix = pred, visitSequence = visSeq, maxit = 30, m = 30, printFlag = TRUE, seed = 2022) #imputing long df

imp.test$loggedEvents
plot(imp.test) #check convergence 

#\ (save imputed df once and from then on open saved workspace)
saveRDS(imp.test, 'imputed_df_after_revisions.rds')

save.image("imputed_df.RData")

load("imputed_df.RData")

##########Analyses##########----
##Setting WD for results analyses
setwd('set_path_to_results_folder')

##STEP 1: model with age and behavior to see how this relationship is (linear or non-linear)
#'Spaghetti plots to visualize individual trajectories for each outcome over time (non imp data)
xyplot(cbcl_int_score ~ time, data = reshape_df_ordered, group = IDC, type = 'b')
xyplot(cbcl_ext_score ~ time, data = reshape_df_ordered, group = IDC, type = 'b')
xyplot(cbcl_total_score ~ time | GENDER, data = reshape_df_ordered, group = IDC, type = 'b') #visualizing sex-specific spaghetti plot 

#'Determine boundary knots for spline of time 
quantile(complete.cases(reshape_df_ordered$time), probs = c(0, 0.05, 0.95, 1)) #B is 0 at 5% and 1 at 95% --> this is to be used for the spline for time = ns(time, df=3, B = c(0, 1))

#"Fitting marginal model on non imp df to get rho
marginal_model <- gls(cbcl_total_score ~ sumscore_inf_tot + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101, data=reshape_df_ordered, correlation=corCompSymm(form = ~time | IDC), na.action = na.exclude) 
summary(marginal_model) #rho = 0.47 [this is the average correlation estimate between repeated measurements]

##STEP 2: Prenatal infection & child behavior over time (mixed model)----
#Extend random effects; test if adding random slope improves model fit 
lm_cbcl_tot_1 <- with(imp.test_mids, lmer(cbcl_total_score ~ sumscore_inf_tot + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1|IDC))) #random intercept, to control for variation in IDC; so accounts for non-independence repeated measures 
summary(pool(lm_cbcl_tot_1))

lm_cbcl_tot_2 <- with(imp.test_mids, lmer(cbcl_total_score ~ sumscore_inf_tot + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1 + time||IDC))) #random slope for time which is uncorrelated to random intercept, to allow different slopes between groups (explanation: http://www.bristol.ac.uk/cmm/learning/videos/random-slopes.html)
summary(pool(lm_cbcl_tot_2))

anova(lm_cbcl_tot_1, lm_cbcl_tot_2) #likelihood ratio test to see if it improves model fit 

#Model building steps for fixed effects part
#'Non-linear effect of time?
lm_cbcl_tot_4 <- with(imp.test_mids, lmer(cbcl_total_score ~ sumscore_inf_tot + ns(time, df=3, B = c(0, 1)) + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1|IDC)), REML = F)
summary(pool(lm_cbcl_tot_4))

anova(lm_cbcl_tot_1, lm_cbcl_tot_4) #non-sign, so adding non-linear effect of time does not improve model fit so moving forward with linear model 

## after now obtaining final model fit we move forward with statistical analyses

#Mixed model loop for outcomes: CBCL tot, CBCL int, CBCL ext and exposures inf tot, tri1, tri2, tri3
outcome_vars <- c("cbcl_total_score", "cbcl_int_score", "cbcl_ext_score") #vector with outcomes

results_mixed_model <- data.frame() #create empty df 

#'Mixed model for total infection and every outcome 
for(x in outcome_vars) {
    f <- paste0(x, "~ sumscore_inf_tot + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1 + time||IDC)") #full model 
    g <- paste0(x, "~ sumscore_inf_tot*GENDER + time + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1 + time||IDC)") #sex interaction 
    h <- paste0(x, "~ sumscore_inf_tot*time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1 + time||IDC)") #time interaction
    #calculating beta, SE, confidence interval and pval for total model 
    bval1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,2]
    seval1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,3]
    lowerCI1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,7]
    upperCI1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,8]
    pval1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))), conf.int = T)[2,6]
    #calculating beta, SE, confidence interval and pval for sex interaction 
    bval2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[21,2]
    seval2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[21,3]
    lowerCI2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[21,7]
    upperCI2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[21,8]
    pval2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))), conf.int = T)[21,6]
    #calculating beta, SE, confidence interval and pval for time interaction 
    bval3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))),conf.int = T)[21,2]
    seval3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))),conf.int = T)[21,3]
    lowerCI3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))),conf.int = T)[21,7]
    upperCI3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))),conf.int = T)[21,8]
    pval3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))), conf.int = T)[21,6]
    #storing values in empty datafile  
    results_mixed_model[x,1] <- bval1
    results_mixed_model[x,2] <- seval1
    results_mixed_model[x,3] <- lowerCI1
    results_mixed_model[x,4] <- upperCI1
    results_mixed_model[x,5] <- pval1
    results_mixed_model[x,6] <- bval2
    results_mixed_model[x,7] <- seval2
    results_mixed_model[x,8] <- lowerCI2
    results_mixed_model[x,9] <- upperCI2
    results_mixed_model[x,10] <- pval2
    results_mixed_model[x,11] <- bval3
    results_mixed_model[x,12] <- seval3
    results_mixed_model[x,13] <- lowerCI3
    results_mixed_model[x,14] <- upperCI3
    results_mixed_model[x,15] <- pval3
    #assigning names to columns 
    colnames(results_mixed_model) <- c("bval1", "seval1", "lowerCI1", "upperCI1", "pval1", "bval2", "seval2", "lowerCI2", "upperCI2", "pval2", "bval3", "seval3", "lowerCI3", "upperCI3", "pval3")
    #saving results to excel 
    write.xlsx(results_mixed_model, "Results_total_infection.xlsx")
}

#'Mixed model for trimester specific infection score and every outcome 
results_mixed_model2 <- data.frame() #create empty df 

for(x in outcome_vars) {
  f <- paste0(x, "~ sumscore_inf_tri1 + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1 + time||IDC)") #full model 
  g <- paste0(x, "~ sumscore_inf_tri2 + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1 + time||IDC)") #sex interaction 
  h <- paste0(x, "~ sumscore_inf_tri3 + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1 + time||IDC)") #time interaction
  #calculating beta, SE, confidence interval and pval for tri1 
  bval1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,2]
  seval1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,3]
  lowerCI1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,7]
  upperCI1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,8]
  pval1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))), conf.int = T)[2,6]
  #calculating beta, SE, confidence interval and pval for tri2
  bval2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[2,2]
  seval2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[2,3]
  lowerCI2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[2,7]
  upperCI2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[2,8]
  pval2 <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))), conf.int = T)[2,6]
  #calculating beta, SE, confidence interval and pval for tri3 
  bval3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))),conf.int = T)[2,2]
  seval3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))),conf.int = T)[2,3]
  lowerCI3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))),conf.int = T)[2,7]
  upperCI3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))),conf.int = T)[2,8]
  pval3 <- summary(pool(with(imp.test_mids, lmer(as.formula(h)))), conf.int = T)[2,6]
  #storing values in empty datafile  
  results_mixed_model2[x,1] <- bval1
  results_mixed_model2[x,2] <- seval1
  results_mixed_model2[x,3] <- lowerCI1
  results_mixed_model2[x,4] <- upperCI1
  results_mixed_model2[x,5] <- pval1
  results_mixed_model2[x,6] <- bval2
  results_mixed_model2[x,7] <- seval2
  results_mixed_model2[x,8] <- lowerCI2
  results_mixed_model2[x,9] <- upperCI2
  results_mixed_model2[x,10] <- pval2
  results_mixed_model2[x,11] <- bval3
  results_mixed_model2[x,12] <- seval3
  results_mixed_model2[x,13] <- lowerCI3
  results_mixed_model2[x,14] <- upperCI3
  results_mixed_model2[x,15] <- pval3
  #assigning names to columns 
  colnames(results_mixed_model2) <- c("bval1", "seval1", "lowerCI1", "upperCI1", "pval1", "bval2", "seval2", "lowerCI2", "upperCI2", "pval2", "bval3", "seval3", "lowerCI3", "upperCI3", "pval3")
  #saving results to excel 
  write.xlsx(results_mixed_model2, "Results_trimester_infection.xlsx")
}

#'Mixed model for combined trimesters model and three way interaction and every outcome 
results_mixed_model3 <- data.frame() #create empty df

for(x in outcome_vars) {
  f <- paste0(x, "~ sumscore_inf_tri1 + sumscore_inf_tri2 + sumscore_inf_tri3 + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1 + time||IDC)") #full model 
  g <- paste0(x, "~ sumscore_inf_tot*time*GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + (1 + time||IDC)") #sex and time interaction 
  #calculating beta, SE, confidence interval and pval for combined trimester model
  bval_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,2]
  seval_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,3]
  lowerCI_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,7]
  upperCI_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,8]
  pval_tri1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))), conf.int = T)[2,6]
  bval_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[3,2]
  seval_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[3,3]
  lowerCI_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[3,7]
  upperCI_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[3,8]
  pval_tri2 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))), conf.int = T)[3,6]
  bval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[4,2]
  seval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[4,3]
  lowerCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[4,7]
  upperCI_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[4,8]
  pval_tri3 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))), conf.int = T)[4,6]
  #calculating beta, SE, confidence interval and pval for three way interaction model 
  bval_int <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[24,2]
  seval_int <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[24,3]
  lowerCI_int <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[24,7]
  upperCI_int <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))),conf.int = T)[24,8]
  pval_int <- summary(pool(with(imp.test_mids, lmer(as.formula(g)))), conf.int = T)[24,6]
  #storing values in empty datafile  
  results_mixed_model3[x,1] <- bval_tri1
  results_mixed_model3[x,2] <- seval_tri1
  results_mixed_model3[x,3] <- lowerCI_tri1
  results_mixed_model3[x,4] <- upperCI_tri1
  results_mixed_model3[x,5] <- pval_tri1
  results_mixed_model3[x,6] <- bval_tri2
  results_mixed_model3[x,7] <- seval_tri2
  results_mixed_model3[x,8] <- lowerCI_tri2
  results_mixed_model3[x,9] <- upperCI_tri2
  results_mixed_model3[x,10] <- pval_tri2
  results_mixed_model3[x,11] <- bval_tri3
  results_mixed_model3[x,12] <- seval_tri3
  results_mixed_model3[x,13] <- lowerCI_tri3
  results_mixed_model3[x,14] <- upperCI_tri3
  results_mixed_model3[x,15] <- pval_tri3
  results_mixed_model3[x,16] <- bval_int
  results_mixed_model3[x,17] <- seval_int
  results_mixed_model3[x,18] <- lowerCI_int
  results_mixed_model3[x,19] <- upperCI_int
  results_mixed_model3[x,20] <- pval_int
  #assigning names to columns 
  colnames(results_mixed_model3) <- c("bval_tri1", "seval_tri1", "lowerCI_tri1", "upperCI_tri1", "pval_tri1", "bval_tri2", "seval_tri2", "lowerCI_tri2", "upperCI_tri2", "pval_tri2","bval_tri3", "seval_tri3", "lowerCI_tri3", "upperCI_tri3", "pval_tri3", "bval_int", "seval_int", "lowerCI_int", "upperCI_int", "pval_int")
  #saving results to excel 
  write.xlsx(results_mixed_model3, "Results_sensitivity_analyses.xlsx")
}

##Multiple testing correction
pval <- c(4.32E-09, 8.43414E-09, 1.75688E-06, 0.885499, 0.854267, 0.945963, 0.552759, 0.817905, 0.325076, 0.000364046, 7.55723E-05, 0.012246, 3.25939E-07, 5.53041E-05, 3.70356E-06, 9.95989E-06, 2.64124E-06, 0.000328262, 0.071382025, 0.010933455, 0.359394407, 0.001259, 0.073033, 0.000903, 0.006162, 0.000869, 0.038513, 0.997832, 0.691802, 0.501697)
p.adjust(pval,method="BH")

###Figures----
setwd('set_path_to_figures_folder')

###FIGURE 3###
#for visualizaition purpose only we want to run regression per individual time point, but first we need to repeat the imputations but now with wide datasets, so with all individual time points
imp0 <- mice(df_final_long, maxit = 0, defaultMethod = c("norm", "logreg", "polyreg", "polr")) 

imp0$loggedEvents

#Imputation method matrix
meth <- imp0$method

#Variables that should not be imputed are set to ""
meth[c(1, 5:6)] <- "" #change column names 

#Predictor matrix
#' write a code to put outcome as predictor but not to be imputed (same was as in nihes project)
pred <- imp0$predictorMatrix

pred[, c(1, 5:6)] <- 0 #Variables that should not be used as predictors set to 0 

#Visit sequence
visSeq <- imp0$visitSequence

#Performing the imputation
#'30 iterations (generally recommended) with 30 datasets for reaching convergence
imp.test <- mice(df_final_long, method = meth, predictorMatrix = pred, visitSequence = visSeq, maxit = 30, m = 30, printFlag = TRUE, seed = 2023) #imputing long df

#repeat linear regressions with imp.test but now with all individual time points for figure purpose only 
outcomes <- c("cbcl_sum_sqrt_standardized.t1", "cbcl_sum_sqrt_standardized.t2", "cbcl_sum_sqrt_standardized.t3", "cbcl_sum_sqrt_standardized.t4", "cbcl_sum_sqrt_standardized.t5", "cbcl_int_sqrt_standardized.t1", "cbcl_int_sqrt_standardized.t2", "cbcl_int_sqrt_standardized.t3", "cbcl_int_sqrt_standardized.t4", "cbcl_int_sqrt_standardized.t5", "cbcl_ext_sqrt_standardized.t1", "cbcl_ext_sqrt_standardized.t2", "cbcl_ext_sqrt_standardized.t3", "cbcl_ext_sqrt_standardized.t4", "cbcl_ext_sqrt_standardized.t5")

results_crossect <- data.frame()

for(x in outcomes){
  f <- paste0(x, "~ sumscore_inf_tot + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101")
  bval <- summary(pool(with(imp.test, lm(as.formula(f)))),conf.int = T)[2,2]
  lowerCI <- summary(pool(with(imp.test, lm(as.formula(f)))),conf.int = T)[2,7]
  upperCI <- summary(pool(with(imp.test, lm(as.formula(f)))),conf.int = T)[2,8]
  pval <- summary(pool(with(imp.test, lm(as.formula(f)))),conf.int = T)[2,6]
  results_crossect [x,1] <- bval
  results_crossect [x,2] <- lowerCI
  results_crossect [x,3] <- upperCI
  results_crossect [x,4] <- pval
}

colnames(results_crossect) <- c("bval", "lowerCI", "upperCI", "pval")

write.xlsx(results_crossect, "Results_crosssect_imp_figure2.xlsx")

##Figure 3 
library(ggplot2)
library(dplyr)
library(cowplot)

#internalizing
dd1 <- data.frame(
  timepoint = c('Child mean age 1.5', 'Child mean age 3', 'Child mean age 6', 'Child mean age 10', 'Child mean age 14'),
  beta = c(0.0195, 0.0244, 0.0361, 0.0353, 0.0262),
  lower = c(0.0019, 0.0074, 0.0203, 0.0180, 0.0072),
  upper = c(0.0371, 0.0414, 0.0519, 0.0526, 0.0453)
)

dd1$behav <- 'Child internalizing symptoms'

#externalizing
dd2 <- data.frame(
  timepoint = c('Child mean age 1.5', 'Child mean age 3', 'Child mean age 6', 'Child mean age 10', 'Child mean age 14'),
  beta = c(0.017785778, 0.031665085, 0.024528723, 0.031101777,0.014740444),
  lower = c(-0.000243007,0.013576329, 0.008680487, 0.013739341, -0.003125753),
  upper = c(0.035814563, 0.049753841, 0.040376959, 0.048464212, 0.03260664)
)

dd2$behav <- 'Child externalizing symptoms'

#total
dd3 <- data.frame(
  timepoint = c('Child mean age 1.5', 'Child mean age 3', 'Child mean age 6', 'Child mean age 10', 'Child mean age 14'),
  beta = c(0.023141419, 0.033937386, 0.030806027, 0.037095699, 0.027003423),
  lower = c(0.005856, 0.016987413, 0.015156638, 0.020196641, 0.008978624),
  upper = c(0.040426815, 0.050887358, 0.046455417,0.053994756, 0.045028222)
)

dd3$behav <- 'Child total behavioral symptoms'

merge <- rbind(dd1, dd2, dd3)

#create panels 
intercepts <- data.frame(behav = unique(merge$behav), intercept = c(0, 0, 0))

merge$timepoint <- as.factor(merge$timepoint)

crossectional_fig <- ggplot(merge, 
       aes(x=factor(timepoint, level=c('Child mean age 1.5', 'Child mean age 3', 'Child mean age 6', 'Child mean age 10', 'Child mean age 14')), 
           y=beta, 
           ymin=lower, 
           ymax=upper)) + 
  #specify line position here
  geom_linerange(size=2,
                 position=position_dodge(width = 0.3), 
                 color = '#999999') +
  geom_hline(yintercept=1, lty=2) +
  #specify point position 
  geom_point(size=3, 
             shape=22, 
             colour="#003399", 
             fill = '#003399', 
             stroke = 0.5,
             position=position_dodge(width = 0.5)) + 
  #adjust theme and lab names 
  theme_bw() +
  scale_x_discrete(name="Time point") +
  scale_y_continuous(name="Beta coefficient", 
                     limits = c(-0.01, 0.1)) + 
  theme(legend.position = 'bottom') + 
  #add panels
  facet_wrap(~behav, ncol = 3) + 
  #change color of panels
  theme(strip.background = element_rect(fill="#99CCFF")) +
  theme(strip.text = element_text(colour = 'black', face = 'bold')) +
  #add mean regression line
  geom_hline(data = intercepts, 
             aes(yintercept = intercept), 
             linetype = "dashed") + 
  #change x labs to skew
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

### Supplementary Figure S3
## Supplementary Figure on trimester specific results = forest plot 
#df inf tri1 (individually adjusted)
df3 <- data.frame(
  outcomes = c('Total', 'Internalizing', 'Externalizing'),
  index = 1:3,
  beta = c(0.040, 0.042, 0.029),
  lower = c(0.018, 0.021, 0.006),
  upper = c(0.062, 0.063, 0.052)
)

df3$Exposure <- "Trimester 1"

#df inf tri2(individually adjusted)
df4 <- data.frame(
  outcomes = c('Total', 'Internalizing', 'Externalizing'),
  index = 1:3,
  beta = c(0.059, 0.045, 0.055),
  lower = c(0.037, 0.023, 0.032),
  upper = c(0.082, 0.066, 0.078)
)

df4$Exposure <- 'Trimester 2'

#df inf tri3(individually adjusted)
df5 <- data.frame(
  outcomes = c('Total', 'Internalizing', 'Externalizing'),
  index = 1:3,
  beta = c(0.051, 0.052, 0.041),
  lower = c(0.028, 0.030, 0.019),
  upper = c(0.074, 0.074, 0.064)
)

df5$Exposure <- 'Trimester 3'

#df inf tri1 (mutually adjusted)
df6 <- data.frame(
  outcomes = c('Total', 'Internalizing', 'Externalizing'),
  index = 1:3,
  beta = c(0.021, 0.028, 0.011),
  lower = c(-0.002, 0.006, -0.013),
  upper = c(0.044, 0.050, 0.035)
)

df6$Exposure <- 'Trimester 1'

#df inf tri2(mutually adjusted)
df7 <- data.frame(
  outcomes = c('Total', 'Internalizing', 'Externalizing'),
  index = 1:3,
  beta = c(0.041, 0.022, 0.043),
  lower = c(0.016, -0.002, 0.018),
  upper = c(0.066, 0.045, 0.068)
)

df7$Exposure <- 'Trimester 2'

#df inf tri3(mutually adjusted)
df8 <- data.frame(
  outcomes = c('Total', 'Internalizing', 'Externalizing'),
  index = 1:3,
  beta = c(0.033, 0.039, 0.025),
  lower = c(0.009, 0.016, 0.001),
  upper = c(0.057, 0.062,0.048)
)

df8$Exposure <- 'Trimester 3'

#combing df's
dd_ind<- rbind(df3, df4, df5) #individually adjusted
dd_mut <- rbind(df6, df7, df8) #mutually adjusted

dd_ind$outcomes <- as.factor(dd_ind$outcomes)
dd_mut$outcomes <- as.factor(dd_mut$outcomes)

dotCOLS = c("black", "black", "black", "black", "black", "black", "black", "black")
barCOLS = c("lightblue4", 'palegreen3', 'mistyrose2')

#creating  forest plot (one for total effect)
a <- ggplot(dd_ind, aes(x=outcomes, y=beta, ymin=lower, ymax=upper,col=Exposure,fill=Exposure)) + 
  #specify position here
  geom_linerange(size=4,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  #specify position here too
  geom_point(size=3, shape=22, colour="black", stroke = 0.5,position=position_dodge(width = 0.5)) + coord_flip() + theme_bw() +
  scale_x_discrete(name="Child psychiatric symptom outcome") +
  scale_y_continuous(name="Beta coefficient", limits = c(-0.020, 0.1)) + theme(legend.position = 'bottom') + theme(legend.text = element_text(size=11), axis.text = element_text(size=11), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_manual(values=dotCOLS)+
  scale_color_manual(values=barCOLS)

b <- ggplot(dd_mut, aes(x=outcomes, y=beta, ymin=lower, ymax=upper,col=Exposure,fill=Exposure)) + 
  #specify position here
  geom_linerange(size=4,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  #specify position here too
  geom_point(size=3, shape=22, colour="black", stroke = 0.5,position=position_dodge(width = 0.5)) + coord_flip() + theme_bw() +
  scale_x_discrete(name="Child psychiatric symptom outcome") +
  scale_y_continuous(name="Beta coefficient", limits = c(-0.020, 0.1)) + theme(legend.position = 'bottom') + theme(legend.text = element_text(size=11), axis.text = element_text(size=11), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_manual(values=dotCOLS)+
  scale_color_manual(values=barCOLS)

plot_grid(a, b, labels = c('A', 'B'))

###########Post-hoc anaysis###########
#Mixed model loop for outcomes: CBCL tot, CBCL int, CBCL ext and exposure inf tot 
setwd('set_path_to_results_folder') #set wd back to results 

outcome_vars <- c("cbcl_total_score", "cbcl_int_score", "cbcl_ext_score") #vector with outcomes

results_mixed_model_prs <- data.frame() #create empty df 

#'Mixed model for total infection and every outcome and now additionally adjusting for genetic confounding 
for(x in outcome_vars) {
  f <- paste0(x, "~ sumscore_inf_tot + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + scale(PRS_neurodev_resid)  + scale(PRS_moodpsych_resid) + scale(PRS_compuls_resid) + (1 + time||IDC)") #full model; main effect 
  #calculating beta, SE, confidence interval and pval for total model 
  bval1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,2]
  seval1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,3]
  lowerCI1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,7]
  upperCI1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))),conf.int = T)[2,8]
  pval1 <- summary(pool(with(imp.test_mids, lmer(as.formula(f)))), conf.int = T)[2,6]
  #storing values in empty datafile  
  results_mixed_model_prs[x,1] <- bval1
  results_mixed_model_prs[x,2] <- seval1
  results_mixed_model_prs[x,3] <- lowerCI1
  results_mixed_model_prs[x,4] <- upperCI1
  results_mixed_model_prs[x,5] <- pval1
  #assigning names to columns 
  colnames(results_mixed_model_prs) <- c("bval1", "seval1", "lowerCI1", "upperCI1", "pval1")
  #saving results to excel 
  write.xlsx(results_mixed_model_prs, "Results_total_infection_prs_correction.xlsx")
}

pval <- c(0.001142, 0.010336, 0.000844)
p.adjust(pval, method = 'fdr')

#\ 

#######################################################
#########EXTRA ANALYSIS FOR REVISIONS FOR JCPP#########
#######################################################
setwd('set_path_to_figures_folder')

###Prevalence plot infections---- 
library(RColorBrewer)

tri1_df <- data.frame(
  infection = as.factor(c("Upper respiratory infection", "Gastro-intestinal infection", "Flu", "Urinary tract infection", "Dermatitis", "Lower respiratory infection", "Herpes Zoster", "Eye infection", "STD")),
 prevalence = as.numeric(c(55.5, 17.4, 16.4, 2.6,2.6, 0.3, 0.1, 0.1, 0.1))
)
tri1_df$timing <- 'Trimester 1'

tri2_df <- data.frame(
  infection = as.factor(c("Upper respiratory infection", "Gastro-intestinal infection", "Flu", "Urinary tract infection", "Dermatitis", "Lower respiratory infection", "Herpes Zoster", "Eye infection", "STD")),
  prevalence = as.numeric(c(54.9, 15.3, 14.2, 2.3,2.7, 0.3, 0.2, 0.2, 0.2))
)
tri2_df$timing <- 'Trimester 2'

tri3_df <- data.frame(
  infection = as.factor(c("Upper respiratory infection", "Gastro-intestinal infection", "Flu", "Urinary tract infection", "Dermatitis", "Lower respiratory infection", "Herpes Zoster", "Eye infection", "STD")),
  prevalence = as.numeric(c(51.7, 15.9, 12.7, 3.5,1.7, 0.4, 0.2, 0.1, 0.1))
)
tri3_df$timing <- 'Trimester 3'

full_trimester_df <- rbind(tri1_df, tri2_df, tri3_df)
full_trimester_df$timing <- as.factor(full_trimester_df$timing)

figure <- ggplot(full_trimester_df, aes(infection, prevalence, fill = timing)) + geom_bar(stat = 'identity', position = 'dodge') + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.title = element_text(face = 'bold'), legend.position = c(0.009, 0.985), legend.justification = c(0,1), legend.box.background = element_rect(fill = 'darkgrey'),legend.box.margin=margin(1,1,1,1), legend.key = element_rect(fill = 'black'), legend.background = element_rect(fill = 'snow2', size = 0.5, linetype = 'solid')) + labs(x = "", y = 'Prevalence (%)') + scale_fill_brewer(palette = 'Pastel1') + guides(fill = guide_legend(title = 'Timing of infection')) + geom_text(aes(label = paste0(round(prevalence, 2), "%")), position = position_dodge(width = 0.9), vjust = -0.5, size = 2)
annotate_figure(figure, top = text_grob('Frequency plot: Prenatal infection types', color = 'black', face = 'bold', size =14))

###Additional sensitivity analyses---- 
setwd('set_path_to_data_folder')

imputed_df_revisions <- readRDS('imputed_df_after_revisions.rds')

imputed_df_revisions_long <- complete(imputed_df_revisions, include = T, action = "long")

setwd('set_path_to_results_folder')

## Sensitivity analysis 1: birth complications
## Sensitivity analysis 2: maternal illness 
## Sensitivity analysis 3: child illness 

outcome_vars <- c("cbcl_total_score", "cbcl_int_score", "cbcl_ext_score") #vector with outcomes

results_mixed_model_sens1 <- data.frame() #create empty df 
results_mixed_model_sens2 <- data.frame() #create empty df 
results_mixed_model_sens3 <- data.frame() #create empty df 

#'Mixed model for sens 1 and 3 
for(x in outcome_vars) {
  f <- paste0(x, "~ sumscore_inf_tot + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 +WEIGHT + GESTBIR + UMPI2 + UtRI2 + PLAWGHT + PLGF_g1+ APGAR5+ PH_NAV+sectio+promm+(1 + time||IDC)") #birth complications 
  g <- paste0(x, "~ sumscore_inf_tot + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + sumscore_childhood_infection + (1 + time||IDC)") #child chronic illness  
  #calculating beta, SE, confidence interval and pval for birth complications
  bval1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))),conf.int = T)[2,2]
  seval1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))),conf.int = T)[2,3]
  lowerCI1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))),conf.int = T)[2,7]
  upperCI1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))),conf.int = T)[2,8]
  pval1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))), conf.int = T)[2,6]
  #calculating beta, SE, confidence interval and pval for child illness
  bval2 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(g)))),conf.int = T)[2,2]
  seval2 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(g)))),conf.int = T)[2,3]
  lowerCI2 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(g)))),conf.int = T)[2,7]
  upperCI2 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(g)))),conf.int = T)[2,8]
  pval2 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(g)))), conf.int = T)[2,6]
  #storing values in empty datafile  
  results_mixed_model_sens1[x,1] <- bval1
  results_mixed_model_sens1[x,2] <- seval1
  results_mixed_model_sens1[x,3] <- lowerCI1
  results_mixed_model_sens1[x,4] <- upperCI1
  results_mixed_model_sens1[x,5] <- pval1
  results_mixed_model_sens3[x,1] <- bval2
  results_mixed_model_sens3[x,2] <- seval2
  results_mixed_model_sens3[x,3] <- lowerCI2
  results_mixed_model_sens3[x,4] <- upperCI2
  results_mixed_model_sens3[x,5] <- pval2
  #assigning names to columns 
  colnames(results_mixed_model_sens1) <- c("bval1", "seval1", "lowerCI1", "upperCI1", "pval1")
  colnames(results_mixed_model_sens3) <- c("bval2", "seval2", "lowerCI2", "upperCI2", "pval2")
  #saving results to excel 
  write.xlsx(results_mixed_model_sens1, "results_mixed_model_sens1.xlsx")
  write.xlsx(results_mixed_model_sens3, "results_mixed_model_sens3.xlsx")
}

#'Mixed model for sens 2, maternal chronic illness 
for(x in outcome_vars) { 
  #model formula 
  f <- paste0(x, "~ sumscore_inf_tot + time + GENDER + AGE_M_v2 + ETHNMv2 + EDUCM_3groups + EDUCP_3groups + INCOME + APM_IQ + WISC13_FSIQ + GSI + SMOKE_ALL + mdrink_updated + f1100101 + psymed_usage + inflam_med_usage +cortico_med_usage + preg_inflam_cond_score_bin + inflam_conditions_score_bin +(1 + time||IDC)") 
  #calculating beta, SE, confidence interval and pval for birth complications
  bval1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))),conf.int = T)[2,2]
  seval1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))),conf.int = T)[2,3]
  lowerCI1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))),conf.int = T)[2,7]
  upperCI1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))),conf.int = T)[2,8]
  pval1 <- summary(pool(with(imputed_df_revisions, lmer(as.formula(f)))), conf.int = T)[2,6]
  #storing values in empty datafile  
  results_mixed_model_sens2[x,1] <- bval1
  results_mixed_model_sens2[x,2] <- seval1
  results_mixed_model_sens2[x,3] <- lowerCI1
  results_mixed_model_sens2[x,4] <- upperCI1
  results_mixed_model_sens2[x,5] <- pval1
  #assigning names to columns 
  colnames(results_mixed_model_sens2) <- c("bval1", "seval1", "lowerCI1", "upperCI1", "pval1")
  #saving results to excel 
  write.xlsx(results_mixed_model_sens2, "results_mixed_model_sens2.xlsx")
}

#\ 

# apply fdr correction to the three new sensitivity analysis tables and save table

apply_fdr_adjustment_and_save <- function(df, column_index, file_name) {
  # Extract the values from the specified column
  val <- unlist(df[, column_index])
  
  # Perform FDR adjustment on the values
  df$fdr_pvalue <- p.adjust(val, method = 'fdr')
  
  # Save the modified dataframe to an Excel file
  file_name_with_extension <- paste0("fdr_", file_name, ".xlsx")
  write.xlsx(df, file_name_with_extension)
  
  return(df)
}

results_sens1_df <- apply_fdr_adjustment_and_save(results_mixed_model_sens1, 5, 'results_mixed_model_sens1')
results_sens2_df <- apply_fdr_adjustment_and_save(results_mixed_model_sens2, 5, 'results_mixed_model_sens2')
results_sens3_df <- apply_fdr_adjustment_and_save(results_mixed_model_sens3, 5, 'results_mixed_model_sens3')

# make a forest plot, figure 4, of main results + prs sensitivity analysis + 3 new sensitivity analyses after revisions 

# ' first create df with results 
df_results_main <- data.frame(
  outcome = c('Child total behavioral symptoms', 'Child internalizing symptoms', 'Child externalizing symptoms'),
  beta = c(0.032, 0.029, 0.026),
  lower = c(0.021, 0.019, 0.016),
  upper = c(0.042, 0.039, 0.037)
)
df_results_main$Model = 'Main model'

df_results_sens2 <- data.frame(
  outcome = c('Child total behavioral symptoms', 'Child internalizing symptoms', 'Child externalizing symptoms'),
  beta = c(0.02265388,0.02201210 , 0.01515914),
  lower = c(0.0081920461,0.0084594620,0.0005389758),
  upper = c(0.03711571,0.03556475,0.02977930)
)
df_results_sens2$Model = 'Additionally adjusted for chronic maternal illness'

df_results_sens1 <- data.frame(
  outcome = c('Child total behavioral symptoms', 'Child internalizing symptoms', 'Child externalizing symptoms'),
  beta = c(0.02604908 , 0.02458711 ,0.01878534),
  lower = c(0.011968480,0.011391474,0.004468637),
  upper = c(0.04012967,0.03778275, 0.03310204)
)
df_results_sens1$Model = 'Additionally adjusted for birth complications'

df_results_prs <- data.frame(
  outcome = c('Child total behavioral symptoms', 'Child internalizing symptoms', 'Child externalizing symptoms'),
  beta = c(0.031, 0.023, 0.033),
  lower = c(0.013, 0.005, 0.014),
  upper = c(0.051, 0.041, 0.053)
)
df_results_prs$Model = 'Additionally adjusted for child genetic confounding'

df_results_sens3 <- data.frame(
  outcome = c('Child total behavioral symptoms', 'Child internalizing symptoms', 'Child externalizing symptoms'),
  beta = c(0.02543672 ,0.02433388 ,0.01822335),
  lower = c(0.011428547,0.011204700,0.003963614),
  upper = c(0.03944489,0.03746307,0.03248309)
)
df_results_sens3$Model = 'Additionally adjusted for childhood infections'

df_forest_plot <- rbind(df_results_main, df_results_sens2,df_results_sens1, df_results_prs, df_results_sens3) #combine plot

#' choose colors of interest 
library(RColorBrewer)
dotCOLS = c("black", "black", "black", "black", "black", "black", "black", "black")
desired_color_palette <- brewer.pal(n = 8, name = "Pastel1") 

#' Define the order of levels for the Model column
model_order <- c('Main model', 
                 'Additionally adjusted for chronic maternal illness', 
                 'Additionally adjusted for birth complications', 
                 'Additionally adjusted for child genetic confounding', 
                 'Additionally adjusted for childhood infections')

#' Reorder the levels of the Model column
df_forest_plot$Model <- factor(df_forest_plot$Model, levels = model_order)

ggplot(df_forest_plot,
       aes(x=outcome, 
           y=beta, 
           ymin=lower, 
           ymax=upper,
           col=Model,
           fill=Model)) + 
  geom_linerange(size=4,
                 position=position_dodge(width = -0.5)) +  # Adjust the position dodge
  geom_point(size=3, 
             shape=22, 
             colour="black", 
             stroke = 0.5,
             position=position_dodge(width = -0.5)) +  # Adjust the position dodge
  coord_flip() +
  theme_bw() +
  scale_x_discrete(name="Child psychiatric symptom outcome") +
  scale_y_continuous(name="Beta coefficient", limits = c(-0.005, 0.06)) + 
  theme(panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.text = element_text(size=14), 
        legend.title = element_text(face = 'bold'), 
        legend.box.background = element_rect(fill = 'darkgrey'),
        legend.box.margin=margin(1,1,1,1), 
        legend.key = element_rect(fill = 'darkgrey'),
        legend.background = element_rect(fill = 'snow2', 
                                         size = 0.5, 
                                         linetype = 'solid'),
        axis.title = element_text(face = 'bold', size = 16),
        axis.text = element_text(size=14), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14)) + 
  scale_fill_manual(values=dotCOLS) +
  scale_color_manual(values=desired_color_palette)

#\ END OF SCRIPT 
