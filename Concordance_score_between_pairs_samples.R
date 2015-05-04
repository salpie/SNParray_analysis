
##########################################################################################################################################
######################################################### SCRIPT TO CALCULATE SIMILARITY BETWEEN SAMPLES #################################
########################################################################################################################################## 
######### Salpie Nowinski ############
######### date: 25/03/2015 ###########
# calculates concordance between pairs of samples of tumour from same patient to identify which tumour subtypes are more closely related #
# combines both break points and absolute copy number before calculating the concordance between the samples #


#####pairs of samples - files_DCIS contains one, files_IDC contains the other########
files_DCIS <- (Sys.glob("/Users/salpienowinski/Documents/pleomorphic_oncoscanv3_100probes/new_pleo/segmentCN/*pLCIS*segmentCN*"))
files_IDC <- (Sys.glob("/Users/salpienowinski/Documents/pleomorphic_oncoscanv3_100probes/new_pleo/segmentCN/*pILC*segmentCN*"))


###################FILES FOR REFERENCE SIMILARITY SCORE##########################################################################

#files_DCIS <- (Sys.glob("/Users/salpienowinski/Documents/100_probes_for_similarity_score_reference/*segmentCN.txt"))
#files_IDC <- (Sys.glob("/Users/salpienowinski/Documents/vandna_cloe_100/mixed_samples_for_clonality_test/*DC*segmentCN.txt"))

##############################################     FOR SPECIFIC COPY NUMBER CHANGES     ##########################################
#loss <- list()
#loss_idc <- list()
#gain <- list()
#gain_idc <- list()
#cnloh <- list()
#cnloh_idc <- list()

########################################################################################################################################## 

####load data and plot data
#load("/Users/salpienowinski/Similarity_results_for_Cloe_Vandna.RData")
#barplot(thefinal_results*100, beside=TRUE, main="Similarity Score for DCIS, LCIS, and pLCIS with paired IDC", las=2, cex.names=0.38, ylab="percentage similarity (%)", col=c("gold", "darkgreen", "darkblue"))
#legend(x="topleft", legend=c("DCIS vs IDC", "LCIS vs IDC", "DCIS vs LCIS"), fill=c("gold", "darkgreen", "darkblue"))




missing <- list()
thefinal_sample_name <- matrix()
thefinal_score <- matrix()
idc_result <- as.numeric(list())

dcis_result <- as.numeric(list())

idc_result_chr <- as.numeric(list()))

dcis_result_chr <- as.numeric(list())
dcis_only <- as.numeric(list())


compareTo2 <- function(dcis, idc) {
  as.numeric(!((dcis > 2 & idc > 2) | (dcis < 2 & idc < 2) | (dcis == idc) )) 
}

compareToMinorCopy <- function(idc_cn, idc_mcn, dcis) {
  as.numeric(idc_cn > 1 & idc_mcn < 1 & dcis > 2)
}

remove_chry <- function(dataframe) {
  dataframe[!dataframe$Chromosome == "chrY", ]
}

bp_segment_leeway <- function(idc, End, Start) {
  (sum(as.numeric(idc$End -idc$Start)) / sum(idc$probes) ) * 10
}

ds <- data.frame(Chromosome=character(),
                  Start=integer(), 
                  End=integer(), 
                  Cn_idc=integer(),
                  mCn_idc=integer(),
                  Cn_dcis=integer(),
                  mCn_dcis=integer(),
                  sampleID=character(),
                  both=integer()
                  ) 
new <- data.frame(Chromosome=character(),
                  Start=integer(), 
                  End=integer(), 
                  Cn_idc=integer(),
                  mCn_idc=integer(),
                  Cn_dcis=integer(),
                  mCn_dcis=integer(),
                  sampleID=character(),
                  result=integer()
                  ) 

n <- nrow(unique(idc$sampleID))

clonality_score <- data.frame(sampleID=character(),
                  #Chromosome=numeric(),
                  #result=integer(),
                  final_score=numeric()
                  ) 

results_lcis <- data.frame(result_dcis=numeric())

library(data.table)



for (i in 1:length(files_IDC)) {

 
score_gain <- 0
score_loss <- 0
score_cnloh <- 0
score_amp <- 0
 
 #combine segment files from samples
 
 idc <- read.table(files_IDC[i], header=TRUE)
 dcis <- read.table(files_DCIS[i], header=TRUE)

 idc <- idc[!idc$Chromosome == "Y", ]
 dcis <- dcis[!dcis$Chromosome == "Y", ]

 idc$mCn[is.na(idc$mCn)] <- 1
 dcis$mCn[is.na(dcis$mCn)] <- 1

 probes_per_bp <- ((idc$End[1] - idc$Start[1]) / idc$probes[1]) * 1 ## can alter last digit to alter number of probes - to decrease or increase number distance between break points of each sample


dcis_We <- setDT(dcis)[, .ind:= cumsum(c(TRUE,Start[-1]!=End[-.N])),
        list(Cn, mCn)][,
       list(chr=Chromosome[1], start=Start[1], stop=End[.N]),
       list(Cn, mCn, .ind)][,.ind:=NULL][]


idc_We <- setDT(idc)[, .ind:= cumsum(c(TRUE,Start[-1]!=End[-.N])),
        list(Cn, mCn)][,
       list(chr=Chromosome[1], start=Start[1], stop=End[.N]),
       list(Cn, mCn, .ind)][,.ind:=NULL][]


###########All similarities#################
#idc$result <-  as.numeric(idc$Cn > 2 | idc$Cn < 2 | idc$Cn == 2 & idc$mCn == 0)
#dcis$result <- as.numeric(dcis$Cn > 2 | dcis$Cn < 2 | dcis$Cn == 2 & dcis$mCn == 0)


##############################################     FOR SPECIFIC COPY NUMBER CHANGES     ##########################################


#######GAIN

idc_We$result_gain <-  as.numeric(idc_We$Cn > 2 & idc_We$Cn < 5)

dcis_We$result_gain <-  as.numeric(dcis_We$Cn > 2)

#######LOSS

idc_We$result_loss <-  as.numeric(idc_We$Cn < 2)

dcis_We$result_loss <-  as.numeric(dcis_We$Cn < 2)

#######CNLOH

idc_We$result_cnloh <-  as.numeric(idc_We$Cn == 2 & idc_We$mCn == 0)

dcis_We$result_cnloh <-  as.numeric(dcis_We$Cn == 2 & dcis_We$mCn == 0)

######AMPLIFICATION

idc_We$result_amp <-  as.numeric(idc_We$Cn > 4)

dcis_We$result_amp <-  as.numeric(dcis_We$Cn > 4)

#print(dcis_We)
#print(idc_We)
##########################################################################################################################################

dcis_we <- setDT(dcis_We)[, .ind:= cumsum(c(TRUE,start[-1]!=stop[-.N])),
        list(result_gain, result_loss, result_cnloh)][,
       list(chr=chr[1], start=start[1], stop=stop[.N]),
       list(result_gain, result_loss, result_cnloh, result_amp, .ind)][,.ind:=NULL][]


idc_we <- setDT(idc_We)[, .ind:= cumsum(c(TRUE,start[-1]!=stop[-.N])),
        list(result_gain, result_loss, result_cnloh)][,
       list(chr=chr[1], start=start[1], stop=stop[.N]),
       list(result_gain, result_loss, result_cnloh, result_amp, .ind)][,.ind:=NULL][]

#print(dcis_we)
#print(idc_we)


idc_we$chr <- gsub("chr", "", idc_we$chr)
idc_we$chr <- gsub("X", "23", idc_we$chr)

dcis_we$chr <- gsub("chr", "", dcis_we$chr)
dcis_we$chr <- gsub("X", "23", dcis_we$chr)

idc_we <- idc_we[!idc_we$chr == "Y", ]
dcis_we <- dcis_we[!dcis_we$chr == "Y", ]

idc_all <- subset(idc_we, idc_we$result_gain > 0 | idc_we$result_loss > 0 | idc_we$result_cnloh > 0 | idc_we$result_amp > 0)
dcis_all <- subset(dcis_we, dcis_we$result_gain > 0 | dcis_we$result_loss > 0 | dcis_we$result_cnloh > 0 | dcis_we$result_amp > 0)

###############REMOVE CHR1 AND CHR16##############################

#dcis_no16 <- dcis_we[!dcis_we$chr == 16, ]
#dcis_no16_1 <- dcis_no16[!dcis_no16$chr == 1, ]

#idc_no16 <- idc_we[!idc_we$chr == 16, ]
#idc_no16_1 <- idc_no16[!idc_no16$chr == 1, ]


#dcis_we <- dcis_no16_1
#idc_we <- idc_no16_1

#dcis_we <- dcis_no16
#idc_we <- idc_no16

#############move chromosomes up##################################

#dcis_we$chr[dcis_we$chr == "17"] <- "16"
#dcis_we$chr[dcis_we$chr == "18"] <- "17"
#dcis_we$chr[dcis_we$chr == "19"] <- "18"
#dcis_we$chr[dcis_we$chr == "20"] <- "19"
#dcis_we$chr[dcis_we$chr == "21"] <- "20"
#dcis_we$chr[dcis_we$chr == "22"] <- "21"
#dcis_we$chr[dcis_we$chr == "23"] <- "22"

#idc_we$chr[idc_we$chr == "17"] <- "16"
#idc_we$chr[idc_we$chr == "18"] <- "17"
#idc_we$chr[idc_we$chr == "19"] <- "18"
#idc_we$chr[idc_we$chr == "20"] <- "19"
#idc_we$chr[idc_we$chr == "21"] <- "20"
#idc_we$chr[idc_we$chr == "22"] <- "21"
#idc_we$chr[idc_we$chr == "23"] <- "22"

##################################################################


for (y in 1:23) {
  idcc <- subset(idc_we, idc_we$chr == y)
  dciss <- subset(dcis_we, dcis_we$chr == y)       ####subset by chromosome
  gains_idcc <- subset(idcc, idcc$result_gain == 1)
  gains_dciss <- subset(dciss, dciss$result_gain == 1)
  gains_dciss$cal <- 0
  gains_idcc$cal <- 0
  print(y)
#####subset if gain
  if (nrow(gains_idcc) > 0 && nrow(gains_dciss) > 0) { 
  #  if (nrow(gains_idcc) == nrow(gains_dciss)) {     #####if both have gains on same chromosome
      for (f in 1:nrow(gains_idcc)) {
        for (y in 1:nrow(gains_dciss)) {  
          if (gains_dciss$start[y] <= (idcc$start[f] + probes_per_bp) | dciss$start[y] <= (idcc$start[f] - probes_per_bp) | dciss$start[y] == idcc$start[f] && dciss$stop[y] >= (idcc$stop[f] + probes_per_bp) | dciss$stop[y] >= (idcc$stop[f] - probes_per_bp) | dciss$stop[y] <= idcc$stop[f]) {
            if("cal" %in% colnames(gains_dciss)){
              if (gains_dciss$cal[y] == 1 | gains_idcc$cal[f] == 1){
                #print("already_matched")
              }
              else {
                  score_gain = score_gain + 1
                  gains_dciss$cal[y] <- 1
                  gains_idcc$cal[f] <- 1

              }
            } 
          }
        }
      }
    #}
  }
}
#print("FINAL ANSWER")
print(score_gain)

for (y in 1:23) {
  idcc <- subset(idc_we, idc_we$chr == y)
  dciss <- subset(dcis_we, dcis_we$chr == y)       ####subset by chromosome
  loss_idcc <- subset(idcc, idcc$result_loss == 1)
  loss_dciss <- subset(dciss, dciss$result_loss == 1)
  loss_dciss$cal <- 0
  loss_idcc$cal <- 0
#####subset if loss
  if (nrow(loss_idcc) > 0 && nrow(loss_dciss) > 0) { 
  #  if (nrow(loss_idcc) == nrow(loss_dciss)) {     #####if both have loss on same chromosome
      for (f in 1:nrow(loss_idcc)) {
        for (y in 1:nrow(loss_dciss)) {  
          if (loss_dciss$start[y] <= (idcc$start[f] + probes_per_bp) | dciss$start[y] <= (idcc$start[f] - probes_per_bp) | dciss$start[y] == idcc$start[f] && dciss$stop[y] >= (idcc$stop[f] + probes_per_bp) | dciss$stop[y] >= (idcc$stop[f] - probes_per_bp) | dciss$stop[y] <= idcc$stop[f]) {
            if("cal" %in% colnames(loss_dciss)){
              if (loss_dciss$cal[y] == 1 | loss_idcc$cal[f] == 1){
                #print("already_matched")
              }
              else {
                score_loss = score_loss + 1
                loss_dciss$cal[y] <- 1
                loss_idcc$cal[f] <- 1

            }
          } 
        }
      }
    }
  }
}


for (y in 1:23) {
  idcc <- subset(idc_we, idc_we$chr == y)
  dciss <- subset(dcis_we, dcis_we$chr == y)       ####subset by chromosome
  cnloh_idcc <- subset(idcc, idcc$result_cnloh == 1)
  cnloh_dciss <- subset(dciss, dciss$result_cnloh == 1)
  cnloh_dciss$cal <- 0
  cnloh_idcc$cal <- 0
#####subset if cnloh
  if (nrow(cnloh_idcc) > 0 && nrow(cnloh_dciss) > 0) { 
  #  if (nrow(cnloh_idcc) == nrow(cnloh_dciss)) {     #####if both have cnloh on same chromosome
      for (f in 1:nrow(cnloh_idcc)) {
        for (y in 1:nrow(cnloh_dciss)) {  
          if (cnloh_dciss$start[y] <= (idcc$start[f] + probes_per_bp) | dciss$start[y] <= (idcc$start[f] - probes_per_bp) | dciss$start[y] == idcc$start[f] && dciss$stop[y] >= (idcc$stop[f] + probes_per_bp) | dciss$stop[y] >= (idcc$stop[f] - probes_per_bp) | dciss$stop[y] <= idcc$stop[f]) {
            if("cal" %in% colnames(cnloh_dciss)){
              if (cnloh_dciss$cal[y] == 1 | cnloh_idcc$cal[f] == 1){
              }
              else {
                score_cnloh = score_cnloh + 1
                cnloh_dciss$cal[y] <- 1
                cnloh_idcc$cal[f] <- 1
            }
          } 
        }
      }
    }
  }
}



for (y in 1:23) {
  idcc <- subset(idc_we, idc_we$chr == y)
  dciss <- subset(dcis_we, dcis_we$chr == y)       ####subset by chromosome
  amp_idcc <- subset(idcc, idcc$result_amp == 1)
  amp_dciss <- subset(dciss, dciss$result_amp == 1)
  amp_dciss$cal <- 0
  amp_idcc$cal <- 0
#####subset if amp
  if (nrow(amp_idcc) > 0 && nrow(amp_dciss) > 0) { 
  #  if (nrow(amp_idcc) == nrow(amp_dciss)) {     #####if both have amp on same chromosome
      for (f in 1:nrow(amp_idcc)) {
        for (y in 1:nrow(amp_dciss)) {  
          if (amp_dciss$start[y] <= (idcc$start[f] + probes_per_bp) | dciss$start[y] <= (idcc$start[f] - probes_per_bp) | dciss$start[y] == idcc$start[f] && dciss$stop[y] >= (idcc$stop[f] + probes_per_bp) | dciss$stop[y] >= (idcc$stop[f] - probes_per_bp) | dciss$stop[y] <= idcc$stop[f]) {
            if("cal" %in% colnames(amp_dciss)){
              if (amp_dciss$cal[y] == 1 | amp_idcc$cal[f] == 1){
                #print("already_matched")
              }
              else {
                score_amp = score_amp + 1
                amp_dciss$cal[y] <- 1
                amp_idcc$cal[f] <- 1

            }
          } 
        }
      }
    }
  }
}




sim_dcis <- score_gain + score_loss + score_cnloh + score_amp
sim_idc <- score_gain + score_loss + score_cnloh + score_amp

#total_dcis_all <- subset(idc_we, idc_we$result_gain > 0 | idc_we$result_loss > 0 | idc_we$result_cnloh > 0 | idc_we$result_amp > 0)
#total_idc_all <- subset(dcis_we, dcis_we$result_gain > 0 | dcis_we$result_loss > 0 | dcis_we$result_cnloh > 0 | dcis_we$result_amp > 0)


total_idc <- sum(idc_we$result_gain) + sum(idc_we$result_loss) + sum(idc_we$result_cnloh) + sum(idc_we$result_amp)
total_dcis <- sum(dcis_we$result_gain) + sum(dcis_we$result_loss) + sum(dcis_we$result_cnloh) + sum(dcis_we$result_amp)


#total_dcis <- sum(total_dcis_all$score)
#total_idc <- sum(total_idc_all$score)

print(sim_dcis)
print(total_dcis)
print(total_idc)


####concordance#####  CALCULATE CONCORDANCE SCORE 
thefinal_score <- rbind(thefinal_score, (sim_dcis/(sim_idc + ((0.5)*((total_dcis-sim_dcis) + (total_idc - sim_idc))))))
log(thefinal_score)


}



#######################
###load results to plot
#######################


load("/Users/salpienowinski/Similarity_results_for_Cloe_Vandna.RData")
load("/Users/salpienowinski/Documents/100_probes_for_similarity_score_reference/reference_pairs_for_similarity_score.RData")

hist(all_reference, breaks=100, xlim=c(0, 1), col=rgb(0,0,1,1/4), xlab="Similarity Score", main="Similarity Score comparison of DCIS and IDC pairs and Reference unmatched pairs")
hist(results_$V1[c(1:3,5,7:14)], col=rgb(1,0,0,1/4), xlim=c(0,1), add=T, breaks=15)
legend(x="topright", legend=c("reference pairs", "matched DCIS-IDC pairs"), fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)))





hist(not_matched_most_samples, breaks=25, xlim=c(0,1), col=rgb(0,0,1,1/4), main="Clonality Score: Reference Pairs and Clonal Pairs")




barplot(final__, beside=TRUE, legend=c("DCIS", "LCIS"), main="Similarity Score for DCIS, LCIS, and pLCIS with paired IDC", las=2, cex.names=0.6, ylab="percentage similarity to idc (%)")





####################################
######load results for distance lm()
####################################

load("Similarity_score_and_distance.RData")

###dcis
plot(results_$V1, results_$dcis_distance, pch=18, main="distance between dcis and idc and similarity score", xlab="similarity score", ylab="distance (mm) to idc")  ####plot
abline(lm(results_$dcis_distance ~ results_$V1), col="red")  ###regression line

###lcis
plot(results_$V2, results_$lcis_distance, pch=18, main="distance between lcis and idc and similarity score", xlab="similarity score", ylab="distance (mm) to idc")  ####plot
abline(lm(results_$lcis_distance ~ results_$V2), col="red")



#### order your similarity score first for this and make sure you 
df <- as.data.frame(thefinal_score)
pvalue <- list()
similarity_score <- list()
t.test(all_reference, df)
for (i in 2:nrow(df)) {
  print(df[i])
  a <- t.test(all_reference, df[1:paste(i)])
  pvalue <- c(pvalue, a$p.value)
  similarity_score <- c(similarity_score, df[i])
}

FINAL_RESULT <- cbind(pvalue, similarity_score)
