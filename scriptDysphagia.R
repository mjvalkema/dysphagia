# Dysphagia study
# M.J. Valkema April 2023

#############  Set up ############# 
rm(list=ls()) # clear global environment
library(dplyr) # function coalesce
library(lubridate) # function time differences
library(tableone)
library(ggplot2)
library(car)
library(tidyr) # function fill

###  plots
library("lattice")
library("nlme")
library("splines")

source('/Users/Maartje/Documents/Promotietraject/SANO/DataPrepSANOQoL.R') 
# The above script loads:
# data -> these are all available QoL questionnaires from SANO patients
# dataCastor -> baseline data from Castor using labels, to retrieve patients eligible for dysphagia study
# CastorIDs -> this is the link between Castor Participant.Ids and SANO patients
# SANOrecur -> for all Participant.Ids it lists the date of LRR or meta

# Export CREs within Repeating Data from Castor, use names (labels)
setwd(dir="/Users/Maartje/repos/dysphagia/data/Dysphagia_during_active_surveill_csv_export_20221221095425")
CREs <- read.csv(file="Dysphagia_during_active_surveill_CRE_export_20221221.csv", header=TRUE, sep = ";")

setwd(dir="/Users/Maartje/repos/dysphagia/")

############# Load baseline SANO data for this study ############# 
# Choose columns to be used from the SANO dataset in patients who underwent AS
names(SANOdata)
SANOvars <- c("SubjectKey", "SiteCode", "TreatmentChosen", "ARM", "PatientType", "age_diagnosis", "SEX", "histology", "grade",
              "tumlocation", "DEND", "cTstage", "cNstage","DENDOS_base")

BLdataSANO <- subset(SANOdata, c(TreatmentChosen == "AS" | TreatmentChosen == "SANO ALL" | SubjectKey == '5' | SubjectKey == '852'))
# include subjectnr 5 of surgery arm because patient underwent AS
# include subjectnr 852 of surgery arm because patient underwent AS at own request (had no-pass at CRE-2)

BLdataSANO <- BLdataSANO[, SANOvars]
BLdataSANO <- BLdataSANO[!duplicated(BLdataSANO), ]
BLdataSANO <- merge(BLdataSANO, SANOrecur, by = c("SubjectKey"), all.x = TRUE)

# Rename variables so it matches variable names in Castor database
names(BLdataSANO)[names(BLdataSANO) == "age_diagnosis"] <- "age"
names(BLdataSANO)[names(BLdataSANO) == "SEX"] <- "sex"
names(BLdataSANO)[names(BLdataSANO) == "grade"] <- "tumor_differentiation"
names(BLdataSANO)[names(BLdataSANO) == "tumlocation"] <- "tumor_location"
names(BLdataSANO)[names(BLdataSANO) == "DEND"] <- "end_ncrt"
names(BLdataSANO)[names(BLdataSANO) == "DENDOS_base"] <- "date_diagnosis"

# In Castor, only select SANO patients
dataCastor <- dataCastor[dataCastor$study == "SANO",]

# Dataprep Castor Baseline data
# Dates in correct format
vars <- names(dataCastor %>% dplyr::select(starts_with("D")))
drop <- c("diabetes") #specify all variables not dates
datevars <- c(vars[!(vars %in% drop)], "end_ncrt")
for (datevar in datevars) {
  dataCastor[[datevar]] <- as.Date(dataCastor[[datevar]], format = '%d-%m-%Y')
}

# Merge LRR, OK and meta dates to BLdataCastor
BLdataCastor <- merge(dataCastor, SANOrecur, by = c("SubjectKey", "Participant.Id"), all.x = TRUE) # only do once

# Enrich SANO with Castor columns
BLvars <- c("Participant.Id", "age", "sex", "histology", "tumor_differentiation", "tumor_location", "cTstage", "cNstage", "end_ncrt", "date_diagnosis", "LRR", "OK", "meta") # specify variables both in SANO and Castor
dataDys <- merge(BLdataSANO, BLdataCastor[BLdataCastor$study == "SANO",!(names(BLdataCastor) %in% BLvars)], by = 'SubjectKey', all.x = TRUE) 

dataDys <- subset(dataDys, c(AS_yn == "yes" & cCR != "no")) # Only use cCR patients who underwent active surveillance for the dysphagia study

dataDys$end_ncrt[dataDys$Participant.Id == "cCR_131"] <- "2018-12-11" # correct one nCRT date for this pt

dataDys$LRR_interval <- lubridate::interval(dataDys$end_ncrt, dataDys$LRR) %/% months(1)
dataDys$LRR16mo_yn <- ifelse(is.na(dataDys$LRR_interval), 0, ifelse(dataDys$LRR_interval < 17, 1, 0))
dataDys$LRRyn <- as.factor(ifelse(!is.na(dataDys$LRR), 1, 0))
dataDys$metayn <- as.factor(ifelse(!is.na(dataDys$meta), 1, 0))

# Exclude patients from Dysphagia study who have a LRR at CRE-3 (6 months after nCRT)
# Should be character to prevent NAs to be dropped if this stays numeric
dataDys$LRR_interval_character <- as.character(dataDys$LRR_interval)
dataDys$LRR_interval_character <- ifelse(is.na(dataDys$LRR_interval_character), "99", dataDys$LRR_interval_character) # assign 99 to patients LRR-free
unique(dataDys$LRR_interval_character)
dataDys <- dataDys[!(dataDys$LRR_interval_character %in% c("4", "5", "6")),]

dataDys$LRR_interval_character_months <- car::recode(as.factor(dataDys$LRR_interval_character),
                              "c('7', '8', '9', '10')='9';
                              c('11', '12', '13', '14')='12';
                              c('15', '16', '17', '18')='16';
                              c('19', '20', '21', '22')='20';
                              c('23', '24', '25', '26', '27')='24';
                              c('28', '29', '30', '31', '32', '33')='30';
                              c('34', '35', '36', '37', '38', '39')='36';
                              c('42', '45', '47', '48', '49')='48';
                              c('60')='60'")
dataDys$LRR_interval_character_months <- factor(dataDys$LRR_interval_character_months, ordered = TRUE, levels = c('0', '3', '6', '9', '12', '16', '20', '24', '30', '36', '48', '60'))

IDs_DysphagiaStudy <- unique(dataDys$Participant.Id)
length(IDs_DysphagiaStudy) # 131 patients included

# Rename values
dataDys$sex <- as.factor(car::recode(dataDys$sex, "c('Male')='male'; c('Female')='female'"))
dataDys$cTstage <- car::recode(dataDys$cTstage, "c('T2')='cT2'; c('T3')='cT3'; c('T4a')='cT4'; c('Tx')='cTx'")
dataDys$cNstage <- car::recode(dataDys$cNstage, "c('N0')='cN0'; c('N1')='cN1'; c('N2')='cN2'")
dataDys$histology <- car::recode(dataDys$histology, "c('Adenocarcinoma')='adenocarcinoma'; 
                   c('Squamous cell carcinoma')='squamous cell carcinoma';
                   c('adenosquamous', 'Other')='other'")
dataDys$tumor_differentiation <- car::recode(dataDys$tumor_differentiation, "c('(G1) Well differentiated', '(G2) Moderately differentiated', 'good', 'moderate')='good-moderate'; 
                   c('(G3) Poorly differentiated', '(G4) Undifferentiated')='poor'; 
                   c('Unknown', '(GX) Grade cannot be assessed')='unknown'")
dataDys$tumor_location <- car::recode(dataDys$tumor_location, "c('Middle esophagus', 'mid')='middle'; 
                   c('Distal', 'Distal esophagus', 'GEJ', 'Esophagogastric junction', 'distal')='distal-GEJ';
                   c('Proximal esophagus')='proximal'")
dataDys$circumference <- car::recode(dataDys$circumference, "c('##USER_MISSING_95##', '##USER_MISSING_99##', NA, '')='unknown';
                                     c('not done - exophytic')='N/A, exophytic lesion'")
dataDys$tumor_length_endoscopy <- replace(dataDys$tumor_length_endoscopy, which(dataDys$tumor_length_endoscopy < 0), NA) # convert negative values (missing) to NAs
dataDys$tumor_length_baselineCT <- replace(dataDys$tumor_length_baselineCT, which(dataDys$tumor_length_baselineCT < 0), NA) # convert negative values (missing) to NAs

write.csv(dataDys, "output/dataset_Dysphagiastudy.csv", row.names = FALSE, na = "")

#-----------
############# data descriptives baseline of dataDys #############
# This describes patients with cCR at CRE-2 and CRE-3, since cancer status has to be 0 at least one time-moment later after obtaining cCR

# Baseline table
baseVars <- c("age", "sex", "histology", "tumor_location", "tumor_differentiation", "cTstage", "cNstage", "circumference", "tumor_length_endoscopy", "tumor_length_baselineCT")
catVars <- c("sex", "histology", "tumor_location", "tumor_differentiation", "cTstage", "cNstage", "circumference")
tableBase <- CreateTableOne(vars = baseVars, factorVars = catVars, data = dataDys)
tableBase <- print(tableBase, nonnormal = c("age", "tumor_length_endoscopy", "tumor_length_baselineCT"), quote = FALSE, noSpaces = TRUE, digits=NULL)
write.csv(tableBase, "output/Table1.csv", row.names = TRUE, na = "")

############# 1. What are the dysphagia scores? #############
# Preparation of QoL dysphagia scores
# Use dataDys for patient set
dataDysQoL_IDs <- dataDys[dataDys$study == "SANO",]$SubjectKey 
length(dataDysQoL_IDs) # 131 SANO patients in the Dysphagia study who were asked to participate in QoL questionnaires

# Obtain the quality of life questionnaires for these patients by using "data" including all SANO QoL questionnaires
# Select patients who have participated in questionnaires until 16 months after nCRT
dataDysQoL <- subset(data, data$SubjectKey %in% dataDysQoL_IDs) # get QoL for these patients]
dataDysQoL <- subset(dataDysQoL, dataDysQoL$time_months <17) # only select patients who participated until 16 months
dataDysQoL <- subset(dataDysQoL, !is.na(dataDysQoL$dysphagiaScore)) # remove the questionnaire for patients who skipped a page filling in answers and have NA score for the dysphagia score
write.csv(dataDysQoL, "output/DataDysQoL.csv", row.names = FALSE, na = "")

# Describe number of patients in the set
length(unique(dataDysQoL[dataDysQoL$LRRstatus2 =="0" & dataDysQoL$timePoint_months != "during",]$SubjectKey)) # 114 number of patients after subselection, so who participated at least once between baseline and 16 months after nCRT
summary(dataDysQoL[dataDysQoL$time_months > 1 & dataDysQoL$LRRstatus2 =="0",]$dateQoL) # to get range of dates when questionnaires have been completed between 3 - 16 months after nCRT

# Exploration of the data
# plot dysphagia scores in patients in whom LRR occurs vs who remain cancer-free at the primary tumor site
xyplot(dysphagiaScore ~ time_months | LRR16mo_yn, group = SubjectKey, data = dataDysQoL[dataDysQoL$OKstatus == "0", ], 
       panel = function (x, y, ...) {
         #panel.xyplot(x, y, type = "l", col = 1, ...)
         panel.loess(x, y, col = 2, lwd = 2)
       }, xlab = "Time (months)", ylab = "Dysphagia score", ylim = c(0, 110), xlim = c(-2, 20))

# plot dysphagia scores in patients in patients who remain cancer-free at the primary tumor site
# only use scores until recurrence occurs at the next response evaluation
xyplot(dysphagiaScore ~ time_months, group = SubjectKey, data = dataDysQoL[dataDysQoL$LRRstatus2 == "0", ], 
       panel = function (x, y, ...) {
         #panel.xyplot(x, y, type = "l", col = 1, ...)
         panel.loess(x, y, col = 2, lwd = 2)
       }, xlab = "Time (months)", ylab = "Dysphagia score", ylim = c(0, 40), xlim = c(-2, 20))

# Analysis of dysphagia scores using questionnaires at which there is no residual tumor at that time and not within the next 4 months
dataDysQoL %>%
  filter(LRRstatus2 == 0) %>%
  group_by(timePoint_months, SubjectKey) %>%
  top_n(n=-1, wt = dateQoL) %>% #the first questionnaire was analyzed in case patients completed the questionnaire twice by mistake
  ungroup() %>%
  group_by(timePoint_months) %>%
  dplyr::summarize(n_pt_participating = n_distinct(SubjectKey),
            score_mean=mean(dysphagiaScore),
            score_sd=sd(dysphagiaScore),
            score_50=quantile(dysphagiaScore, probs = 0.5),
            score_25=quantile(dysphagiaScore, probs = 0.25),
            score_75=quantile(dysphagiaScore, probs = 0.75),
            score_5=quantile(dysphagiaScore, probs = 0.05),
            score_95=quantile(dysphagiaScore, probs = 0.95))

# Count number of eligible patients per time-point
# two patients have LRRstatus2 at 3 months after nCRT: this is because they developed recurrence at CRE-4 (time point 9 months) but their recurrence was at 7 months
dataDysQoL$LRRstatus2_num <- as.numeric(as.character(dataDysQoL$LRRstatus2))

dataDysQoL %>%
  group_by(SubjectKey) %>%
  arrange(-LRRstatus2_num, dateQoL, by_group = TRUE) %>%
  slice(1) %>% # take one row per patient with lowest date and highest LRRstatus
  filter(LRRstatus2_num == 1) %>% # filter out everything not LRRstatus 1
  group_by(timePoint_months) %>%
  dplyr::summarize(count=n(), count_check=n_distinct(SubjectKey))

# Describe patients with outlier dysphagia scores >35
dataDysQoL %>%
  filter(LRRstatus2 == 0) %>%
  filter(time_months > 1) %>%
  group_by(timePoint_months) %>%
  filter(dysphagiaScore >34) %>%
  dplyr::summarize(SubjectKey, dysphagiaScore, OG1, OG2, OG3, dateQoL)

# Check one patient
number <- "200"
dataDysQoL[dataDysQoL$SubjectKey == number,][names(dataDysQoL) %in% c("SubjectKey", "dysphagiaScore", "OG1", "OG2", "OG3", "timePoint_months", "time", "dateQoL", "LRRstatus1", "LRRstatus2", "OKstatus")][order(dataDysQoL[dataDysQoL$SubjectKey == number,][names(dataDysQoL) %in% c("SubjectKey", "dysphagiaScore", "OG1", "OG2", "OG3", "timePoint_months", "time", "dateQoL", "LRRstatus1", "LRRstatus2", "OKstatus")]$dateQoL),]
# check follow-up after 16 months after nCRT in the total QoL set (this dataframe is limited until 16 months after nCRT)

############# 2. Clinically relevant stenosis with or without dysphagia #############
CREdf <- CREs[CREs$Participant.Id %in% dataDys$Participant.Id,] # take CRE forms of patients within dataDys
CREdf <- CREdf[CREdf$CRE!="baseline",] # remove baseline reports

# Remove unused columns CREdf
names(CREdf)
drop_columns <- c("Participant.Status", "Repeating.Data.Creation.Date", "Repeating.data.Parent", "X" , "interval")
for (drop_column in drop_columns) {
  index <- grep(drop_column, colnames(CREdf))
  CREdf <- CREdf[,!(names(CREdf) %in% names(CREdf[index]))]
}

# Prepare dates in CREs
datevars <- c("date_CRE", "intervention_date")
for (datevar in datevars) {
  CREdf[[datevar]] <- as.Date(CREdf[[datevar]], format = '%d-%m-%Y')
}

# Add end of nCRT dates to CRE dataframe and calculate interval between CRE and interventions and end of nCRT
dfnCRT <- as.data.frame(dataDys %>% dplyr::select(Participant.Id, end_ncrt))
CREdf <- merge(CREdf, dfnCRT, by = "Participant.Id", all.x = T) # add end of nCRT to CRE data

CREdf$interval_CRE <- lubridate::interval(CREdf$end_ncrt, CREdf$date_CRE) %/% months(1)
CREdf$interval_intervention <- lubridate::interval(CREdf$end_ncrt, CREdf$intervention_date) %/% months(1)

# Add LRR date to CRE dataframe
CREdf <- merge(CREdf, SANOrecur, by = "Participant.Id", all.x = TRUE)

dfintervalLRR <- as.data.frame(dataDys %>% dplyr::select(Participant.Id, LRR_interval_character_months))
CREdf <- merge(CREdf, dfintervalLRR, by = "Participant.Id", all.x = TRUE)

# per CRE: compute whether there is LRR or not
CREdf$LRRstatus1 <- as.factor(ifelse(is.na(CREdf$LRR), 0, ifelse(CREdf$date_CRE >= CREdf$LRR, 1, 0))) # has LRR at that time moment
CREdf$LRRstatus2 <- as.factor(ifelse(is.na(CREdf$LRR), 0, ifelse(CREdf$date_CRE >= (CREdf$LRR %m-% period("4 month")), 1, 0))) # has LRR at the CRE or max 4 months later

# per CRE: compute whether there has been an intervention (at the particular or previous CRE) or not
CREdf$intervention_date_last <- CREdf$intervention_date # copy variable with intervention date since this variable will be mutated below
CREdf <- CREdf %>% group_by(Participant.Id) %>% arrange(date_CRE, .by_group = TRUE) %>% fill(intervention_date_last) %>% ungroup()

write.csv(CREdf, "output/CREs.csv", row.names = FALSE, na = "")

table(CREdf$intervention_dysphagia) # number of interventions

unique(CREdf[CREdf$intervention_dysphagia == "yes",]$Participant.Id) # unique patients who underwent intervention for dysphagia
# patient cCR_200 also patient with clinically relevant dysphagia, who did not undergo dilatation but received tube feeding for dysphagia


############# 3. What are the rates of stenosis after nCRT? #############
# In patients without locoregional recurrence at that particular CRE and the next 
# Excluding patients who (eventually) required an intervention for dysphagia 
CREdf_stenosis <- subset(CREdf, CREdf$Participant.Id %in% IDs_DysphagiaStudy) # get CREs for patients included in the study

# when did these patients have a LRR
# to get number of patients participating at each CRE being cancer-free
IDs_CREs <- unique(CREdf_stenosis$Participant.Id) 
intervals_CRE_LRR <- subset(dataDys, dataDys$Participant.Id %in% IDs_CREs)
table(intervals_CRE_LRR$LRR_interval_character_months) 

CREdf_stenosis <- CREdf_stenosis[!(CREdf_stenosis$interval_CRE > 18),] # excluding CREs after 18 months, corresponding with CRE at 16 months
CREdf_stenosis <- CREdf_stenosis[!(CREdf_stenosis$LRRstatus2 == "1"),] # excluding patients who had LRR at that CRE or the next

# Compare patients with stenosis to patients without stenosis (including patients with clinically relevant dyshpagia)
IDs_with_stenosis <- unique(CREdf_stenosis[CREdf_stenosis$stenosis == "yes",]$Participant.Id) # unique patients with stenosis
length(IDs_with_stenosis) # number of patients

dataDys$stenosis <- ifelse(dataDys$Participant.Id %in% IDs_with_stenosis, "yes", "no")

# Comparison baseline statistics in patients with and without ever a stenosis
dataDys$tumor_location_cat <- as.factor(car::recode(dataDys$tumor_location, "c('middle', 'proximal')='prox-mid'"))
dataDys$cT_cat <- as.factor(car::recode(dataDys$cTstage, "c('cT3', 'cT4')='cT3-4'"))
dataDys$circumference_cat <- as.factor(car::recode(dataDys$circumference, "c('>75%', '51-75%')='>50%'; c('0-25%', '26-50%')='0-50%'; c('N/A, exophytic lesion', 'unknown')='unknown'"))

# Table stenosis vs non-stenosis
baseVars <- c("sex", "histology", "tumor_location_cat", "tumor_differentiation", "cT_cat", "circumference_cat", "tumor_length_endoscopy")
catVars <- c("sex", "histology", "tumor_location_cat", "tumor_differentiation", "cT_cat", "circumference_cat")
tableStenosis <- CreateTableOne(vars = baseVars, factorVars = catVars, strata = "stenosis", data = dataDys)
tableStenosis <- print(tableStenosis, nonnormal = c("tumor_length_endoscopy"), exact = catVars, quote = FALSE, noSpaces = TRUE, digits=NULL)
write.csv(tableStenosis, "output/Stenosis.csv", row.names = TRUE, na = "")

# Check Fisher exact tests
fisher.test(rbind(c(80,20), c(24,7))) # sex
fisher.test(rbind(c(75,18), c(28, 8))) # histology
fisher.test(rbind(c(16,4), c(88, 23))) # location
fisher.test(rbind(c(73,14), c(29, 10))) # differentiation
fisher.test(rbind(c(21,6), c(76, 17))) # cT
fisher.test(rbind(c(41,14), c(41, 7))) # circumference

# Make table per time-moment: how many patients had clinically relevant dysphagia
# Manually add patients with clinically relevant dysphagia, censoring them after first intervention
CREdf_stenosis <- CREdf_stenosis[!(CREdf_stenosis$Participant.Id == "cCR_071"),] # excluding patients who eventually required an intervention for dysphagia

length(unique(CREdf_stenosis$Participant.Id)) # number of patients
table(CREdf_stenosis$stenosis) # number of stenosis
table(CREdf_stenosis$no_pass) # type of stenosis

# List number and type of stenosis per time-point
CREdf_stenosis$interval_CRE_months <- car::recode(as.factor(CREdf_stenosis$interval_CRE),
                              "c('1')='1';
                              c('2', '3', '4')='3';
                              c('5', '6', '7')='6';
                              c('8', '9', '10')='9';
                              c('11', '12', '13', '14')='12';
                              c('15', '16', '17', '18')='16';
                              c('19', '20', '21', '22')='20';
                              c('23', '24', '25', '26', '27')='24';
                              c('28', '29', '30', '31', '32', '33')='30';
                              c('34', '35', '36', '37', '38', '39')='36';
                              c('42', '45', '47', '48', '49')='48';
                              c('60')='60'")
CREdf_stenosis$interval_CRE_months <- factor(CREdf_stenosis$interval_CRE_months, ordered = TRUE, levels = c('0', '1', '3', '6', '9', '12', '16', '20', '24', '30', '36', '48', '60'))
write.csv(CREdf_stenosis, "output/CREdf_stenosis.csv", row.names = TRUE, na = "")

# Determine number of patients with a stenosis per time point
CREdf_stenosis %>%
  count(interval_CRE_months, stenosis)

# Number of patients per time point
CREdf_stenosis %>%
  group_by(interval_CRE_months) %>%
  dplyr::summarize(n_SubjectKeys = n_distinct(SubjectKey))

# Determine number of patients with stenosis per time point, with specification of stenosis type
CREdf_stenosis %>%
  count(interval_CRE_months, no_pass)

# Check if patients with stenosis are unique patients
CREdf_stenosis %>%
  filter(stenosis == "yes") %>%
  group_by(interval_CRE_months) %>%
  dplyr::summarize(n_SubjectKeysStenosis = n_distinct(SubjectKey))

# Calculate dysphagia scores for patients with a stenosis at a particular time point
# First define the patients (SubjectKeys) with a stenosis
uniquestenosis_3_easy <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "3" & CREdf_stenosis$no_pass == "no, relative stenosis, easy to pass",]$SubjectKey)
uniquestenosis_3_pressure <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "3" & CREdf_stenosis$no_pass == "no, relative stenosis, passable with pressure",]$SubjectKey)

uniquestenosis_6_easy <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "6" & CREdf_stenosis$no_pass == "no, relative stenosis, easy to pass",]$SubjectKey)
uniquestenosis_6_pressure <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "6" & CREdf_stenosis$no_pass == "no, relative stenosis, passable with pressure",]$SubjectKey)

uniquestenosis_9_easy <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "9" & CREdf_stenosis$no_pass == "no, relative stenosis, easy to pass",]$SubjectKey)
uniquestenosis_9_pressure <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "9" & CREdf_stenosis$no_pass == "no, relative stenosis, passable with pressure",]$SubjectKey)

uniquestenosis_12_easy <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "12" & CREdf_stenosis$no_pass == "no, relative stenosis, easy to pass",]$SubjectKey)
uniquestenosis_12_pressure <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "12" & CREdf_stenosis$no_pass == "no, relative stenosis, passable with pressure",]$SubjectKey)

uniquestenosis_16_easy <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "16" & CREdf_stenosis$no_pass == "no, relative stenosis, easy to pass",]$SubjectKey)
uniquestenosis_16_pressure <- unique(CREdf_stenosis[CREdf_stenosis$interval_CRE_months == "16" & CREdf_stenosis$no_pass == "no, relative stenosis, passable with pressure",]$SubjectKey)

# choose
months <- 16
ids <- uniquestenosis_16_easy
#ids <- uniquestenosis_3_pressure
dataDysQoL %>%
  filter(LRRstatus2 == 0) %>%
  filter(timePoint_months == months) %>%
  filter(SubjectKey %in% ids) %>%
  group_by(SubjectKey) %>% 
  top_n(n=-1, wt = dateQoL) %>% # the first questionnaire was analyzed in case patients completed the questionnaire twice by mistake
  ungroup() %>%
  dplyr::summarize(n_pt_participating = n_distinct(SubjectKey),
                   score_mean=mean(dysphagiaScore),
                   score_sd=sd(dysphagiaScore),
                   score_50=quantile(dysphagiaScore, probs = 0.5),
                   score_25=quantile(dysphagiaScore, probs = 0.25),
                   score_75=quantile(dysphagiaScore, probs = 0.75),
                   score_5=quantile(dysphagiaScore, probs = 0.05),
                   score_95=quantile(dysphagiaScore, probs = 0.95))
