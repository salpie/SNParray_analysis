#####################################################
############To produce Frequency Table in R##########
#################Salpie Nowinski#####################
###########Date: 25/03/15############################
# calculate similarities and differences of absolute copy number at break points between groups of samples ######################


library(dplyr)
library(data.table)
library(biomaRt)




####### Group ONE and TWO location
group_one <- (Sys.glob("*segmentCN.txt")) # path to *segmentCN.txt files produced from TAPS 2.0
group_two <- (Sys.glob("**segmentCN.txt")) # path to *segmentCN.txt files produced from TAPS 2.0


## functions for caluclating copy number changes

findChromothripsis <- function(Cn) {
  if (length(Cn[Cn >= 5]) > 2) { #if more than 2 segements in chromosome have copy number 5
    y1 <- as.numeric(Cn >= 7 & lead(Cn) %in% 1:2) #if there are at least 2 oscillation of from copy number 1/2 to 7 or more then chromothripsis then score 1
    as.numeric(Cn > 7 & lead(Cn) %in% 1:2) #then score 1
  }
  else {
    as.numeric(!(Cn))
  }
}

findAmplicon_single <- function(Cn) {
  if (length(Cn[Cn >= 5]) > 2) {
    y1 <- as.numeric(Cn > 7 & lead(Cn) %in% 1:2)
    as.numeric(Cn >= 7 & (lead(Cn) <= 4 | lag(Cn) <= 4) & !y1)     #if chromothripsis still look for amplicons
  }
  else {
    as.numeric(Cn >= 7 & (lead(Cn) <= 4 | lag(Cn) <= 4))
  }  
}

findAmplicon_stepwise <- function(Cn) {
    as.integer((Cn >= 5 & (lead(Cn) > 7 | lag(Cn) > 7) )| (Cn > 7 & lead(Cn) >= 5 & lag(Cn) >= 5)) #if amplification is 'stepwise' then score 1
}


findGain <- function(Cn) {
  as.numeric(Cn > 2)
}

findGain_not_amp <- function(Cn) {
  as.numeric(Cn > 2 & (!y1 | !y2 | !y3))
}

findLoss <- function(Cn) {  
  as.numeric(Cn < 2)
}

findCnLoh <- function(Cn, mCn) {
  as.numeric(Cn == 2 & mCn == 0)
}

findAmp <- function(Cn) { 
  as.numeric(Cn >= 5)
}

rm(y)
rm(x)
rm(d)
rm(r)
rm(i)


  chromothrip <- as.numeric()
  amplicon_single <- as.numeric()
  amplicon_step <- as.numeric()


group_one_first <- read.table(group_one[1], header=TRUE)
names <- names(group_one_first)

group_one_first$mCn[is.na(group_one_first$mCn)] <- 1

group_one_first$Chromosome <- gsub("chr", "", group_one_first$Chromosome)
group_one_first$Chromosome <- gsub("X", "23", group_one_first$Chromosome)
group_one_first <- group_one_first[!group_one_first$Chromosome == "Y", ]
group_one_first$Cn[group_one_first$Cn == 0 ] <- 1
ds_group_one <- group_one_first

for (i in 1:23) {
    ds_group_chr <- subset(ds_group_one, ds_group_one$Chromosome == i)
    chromothrip <- c(chromothrip, findChromothripsis(ds_group_chr$Cn))
    amplicon_single <- c(amplicon_single, findAmplicon_single(ds_group_chr$Cn))
    amplicon_step <- c(amplicon_step, findAmplicon_stepwise(ds_group_chr$Cn))

  }

  group_one_first$chromothrip <- chromothrip
  group_one_first$amplicon_single <- amplicon_single
  group_one_first$amplicon_step <- amplicon_step

group_one_first$mCn[is.na(group_one_first$mCn)] <- 1

group_one_first$loss <- findLoss(group_one_first$Cn)
group_one_first$gain <- findGain(group_one_first$Cn)
group_one_first$amp <- findAmp(group_one_first$Cn)
group_one_first$cnloh <- findCnLoh(group_one_first$Cn, group_one_first$mCn)

plcis <- group_one_first
ds_group_one <- group_one_first

ds_group_one$Start[ds_group_one$Chromosome == '22'][1] <- 16855618


for (d in 2:length(group_one)) {

  print(group_one[d])

  chromothrip <- as.numeric()
  amplicon_single <- as.numeric()
  amplicon_step <- as.numeric()

  group_one_file_ <- read.table(group_one[d], header=TRUE)
  colnames(group_one_file_) <- names
  group_one_file <- group_one_file_[,c(1,2,3,4,9,10)]

  group_one_file$Chromosome <- gsub("chr", "", group_one_file$Chromosome)
  group_one_file$Chromosome <- gsub("X", "23", group_one_file$Chromosome)
  group_one_file <- group_one_file[!group_one_file$Chromosome == "Y", ]

  group_one_file$Cn[group_one_file$Cn == 0 ] <- 1

  #Human cytosnp only 
  #group_one_file$Start[group_one_file$Chromosome == '22'][1] <- 16855618


###### for each chromosome find chromothripsis and amplicons
  for (r in 1:23) {
    group_one_chr <- subset(group_one_file, group_one_file$Chromosome == r)
    chromothrip <- c(chromothrip, findChromothripsis(group_one_chr$Cn))
    amplicon_single <- c(amplicon_single, findAmplicon_single(group_one_chr$Cn))
    amplicon_step <- c(amplicon_step, findAmplicon_stepwise(group_one_chr$Cn))

  }

  group_one_file$chromothrip <- chromothrip
  group_one_file$amplicon_single <- amplicon_single
  group_one_file$amplicon_step <- amplicon_step

  group_one_file$mCn[is.na(group_one_file$mCn)] <- 1

  group_one_file$loss <- findLoss(group_one_file$Cn)
  group_one_file$gain <- findGain(group_one_file$Cn)
  group_one_file$amp <- findAmp(group_one_file$Cn)
  group_one_file$cnloh <- findCnLoh(group_one_file$Cn, group_one_file$mCn)

  lst1 <- split(ds_group_one, ds_group_one$Chromosome)
  lst2 <- split(group_one_file, group_one_file$Chromosome)


  require(survival)


 # merge
 df <- do.call(rbind, mapply(FUN = function(x, y) {
 
   x$event <- y$event <- 0
   idc.spl <- survSplit(x, cut=y$End, start='Start', end='End', event='event')
   dcis.spl <- survSplit(y, cut=x$End, start='Start', end='End', event='event')
   mrg <- merge(idc.spl, dcis.spl, 
                by=c('Chromosome', 'Start', 'End'), 
                #all=TRUE, 
                suffixes = c("_first","_second"))
   mrg[c('Chromosome', 'Start', 'End', 'loss_first', 'loss_second', 'gain_first', 'gain_second', 'cnloh_first', 'cnloh_second', 'amp_first', 'amp_second', 'chromothrip_first', 'chromothrip_second', 'amplicon_single_first', 'amplicon_single_second', 'amplicon_step_first', 'amplicon_step_second')]
 },
 lst1, lst2, SIMPLIFY=FALSE))

df$loss <- df$loss_second + df$loss_first
 df$gain <- df$gain_second + df$gain_first
 df$cnloh <- df$cnloh_second + df$cnloh_first
 df$amp <- df$amp_second + df$amp_first
df$amplicon_step <- df$amplicon_step_first + df$amplicon_step_second
df$amplicon_single <- df$amplicon_single_first + df$amplicon_single_second
df$chromothrip <- df$chromothrip_first + df$chromothrip_second

ds_group_one <- df[,c(1:3,18:24)]

if (d == 27) {
  print("pLCIS_samples_done")
}

}
rm(y)
rm(x)
rm(d)
rm(r)
rm(i)


  chromothrip <- as.numeric()
  amplicon_single <- as.numeric()
  amplicon_step <- as.numeric()


group_two_first <- read.table(group_two[1], header=TRUE)
names <- names(group_two_first)

group_two_first$mCn[is.na(group_two_first$mCn)] <- 1

group_two_first$Chromosome <- gsub("chr", "", group_two_first$Chromosome)
group_two_first$Chromosome <- gsub("X", "23", group_two_first$Chromosome)
group_two_first <- group_two_first[!group_two_first$Chromosome == "Y", ]
group_two_first$Cn[ds_group_two$Cn == 0 ] <- 1

for (i in 1:23) {
    ds_group_chr <- subset(ds_group_two, ds_group_two$Chromosome == i)
    chromothrip <- c(chromothrip, findChromothripsis(ds_group_chr$Cn))
    amplicon_single <- c(amplicon_single, findAmplicon_single(ds_group_chr$Cn))
    amplicon_step <- c(amplicon_step, findAmplicon_stepwise(ds_group_chr$Cn))

  }

  group_two_first$chromothrip <- chromothrip
  group_two_first$amplicon_single <- amplicon_single
  group_two_first$amplicon_step <- amplicon_step

group_two_first$mCn[is.na(group_two_first$mCn)] <- 1

group_two_first$loss <- findLoss(group_two_first$Cn)
group_two_first$gain <- findGain(group_two_first$Cn)
group_two_first$amp <- findAmp(group_two_first$Cn)
group_two_first$cnloh <- findCnLoh(group_two_first$Cn, group_two_first$mCn)

plcis <- group_two_first
ds_group_two <- group_two_first
ds_group_two$Start[ds_group_two$Chromosome == '22'][1] <- 16855618


for (d in 2:length(group_two)) {

  print(group_two[d])

  chromothrip <- as.numeric()
  amplicon_single <- as.numeric()
  amplicon_step <- as.numeric()

  group_two_file_ <- read.table(group_two[d], header=TRUE)
  colnames(group_two_file_) <- names
  group_two_file <- group_two_file_[,c(1,2,3,4,9,10)]

  group_two_file$Chromosome <- gsub("chr", "", group_two_file$Chromosome)
  group_two_file$Chromosome <- gsub("X", "23", group_two_file$Chromosome)
  group_two_file <- group_two_file[!group_two_file$Chromosome == "Y", ]
  group_two_file$Cn[group_two_file$Cn == 0 ] <- 1

####################HUMAN CYTOSNP ONLY ############################
#group_two_file$Start[group_two_file$Chromosome == '22'][1] <- 16855618


###### for each chromosome find chromothripsis and amplicons
  for (r in 1:23) {
    group_two_chr <- subset(group_two_file, group_two_file$Chromosome == r)
    chromothrip <- c(chromothrip, findChromothripsis(group_two_chr$Cn))
    amplicon_single <- c(amplicon_single, findAmplicon_single(group_two_chr$Cn))
    amplicon_step <- c(amplicon_step, findAmplicon_stepwise(group_two_chr$Cn))

  }

  group_two_file$chromothrip <- chromothrip
  group_two_file$amplicon_single <- amplicon_single
  group_two_file$amplicon_step <- amplicon_step

  group_two_file$mCn[is.na(group_two_file$mCn)] <- 1

  group_two_file$loss <- findLoss(group_two_file$Cn)
  group_two_file$gain <- findGain(group_two_file$Cn)
  group_two_file$amp <- findAmp(group_two_file$Cn)
  group_two_file$cnloh <- findCnLoh(group_two_file$Cn, group_two_file$mCn)

  lst1 <- split(ds_group_two, ds_group_two$Chromosome)
  lst2 <- split(group_two_file, group_two_file$Chromosome)


  require(survival)


 # merge
 df <- do.call(rbind, mapply(FUN = function(x, y) {
 
   x$event <- y$event <- 0
   idc.spl <- survSplit(x, cut=y$End, start='Start', end='End', event='event')
   dcis.spl <- survSplit(y, cut=x$End, start='Start', end='End', event='event')
   mrg <- merge(idc.spl, dcis.spl, 
                by=c('Chromosome', 'Start', 'End'), 
                #all=TRUE, 
                suffixes = c("_first","_second"))
   mrg[c('Chromosome', 'Start', 'End', 'loss_first', 'loss_second', 'gain_first', 'gain_second', 'cnloh_first', 'cnloh_second', 'amp_first', 'amp_second', 'chromothrip_first', 'chromothrip_second', 'amplicon_single_first', 'amplicon_single_second', 'amplicon_step_first', 'amplicon_step_second')]
 },
 lst1, lst2, SIMPLIFY=FALSE))

df$loss <- df$loss_second + df$loss_first
 df$gain <- df$gain_second + df$gain_first
 df$cnloh <- df$cnloh_second + df$cnloh_first
 df$amp <- df$amp_second + df$amp_first
df$amplicon_step <- df$amplicon_step_first + df$amplicon_step_second
df$amplicon_single <- df$amplicon_single_first + df$amplicon_single_second
df$chromothrip <- df$chromothrip_first + df$chromothrip_second

ds_group_two <- df[,c(1:3,18:24)]

if (d == 28) {
  print("iLCIS_samples_done")
}

}


lst1 <- split(ds_group_one, ds_group_one$Chromosome)
lst2 <- split(ds_group_two, ds_group_two$Chromosome)


df <- do.call(rbind, mapply(FUN = function(x, y) {
 
   x$event <- y$event <- 0
   idc.spl <- survSplit(x, cut=y$End, start='Start', end='End', event='event')
   dcis.spl <- survSplit(y, cut=x$End, start='Start', end='End', event='event')
   mrg <- merge(idc.spl, dcis.spl, 
                by=c('Chromosome', 'Start', 'End'), 
                #all=TRUE, 
                suffixes = c("_group_one","_group_two"))
   mrg[c('Chromosome', 'Start', 'End', 'loss_group_one', 'loss_group_two', 'gain_group_one', 'gain_group_two', 'cnloh_group_one', 'cnloh_group_two', 'amp_group_one', 'amp_group_two', 'chromothrip_group_one', 'chromothrip_group_two', 'amplicon_single_group_one', 'amplicon_single_group_two', 'amplicon_step_group_one', 'amplicon_step_group_two')]
 },
 lst1, lst2, SIMPLIFY=FALSE))


df$cnloh_group_two[is.na(df$cnloh_group_two)] <- 0
df$cnloh_group_one[is.na(df$cnloh_group_one)] <- 0


df$group_one_loss <- length(group_one) - df$loss_group_one
df$group_two_loss <- length(group_two) - df$loss_group_two
df$group_one_gain <- length(group_one) - df$gain_group_one
df$group_two_gain <- length(group_two) - df$gain_group_two
df$group_one_cnloh <- length(group_one) - df$cnloh_group_one
df$group_two_cnloh <- length(group_two) - df$cnloh_group_two


df$pvalue_loss <- apply(as.matrix(df[,c(4,18,5,19)]), 1, function(x) 
         fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)


df$pvalue_gain <- apply(as.matrix(df[,c(6,20,7,21)]), 1, function(x) 
         fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)

df$pvalue_cnloh <- apply(as.matrix(df[,c(8,22,9,23)]), 1, function(x) 
         fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)


newdata <- df[order(df$Chromosome, df$Start),]


we <- setDT(newdata)[, .ind:= cumsum(c(TRUE,Start[-1]!=End[-.N])),
        list(loss_group_one, loss_group_two, gain_group_one, gain_group_two, cnloh_group_one, cnloh_group_two, amp_group_one, amp_group_two, chromothrip_group_one, chromothrip_group_two, amplicon_single_group_one, amplicon_single_group_two, amplicon_step_group_one, amplicon_step_group_two,pvalue_gain, pvalue_loss, pvalue_cnloh)][,
       list(chr=Chromosome[1], start=Start[1], stop=End[.N]),
       list(loss_group_one, loss_group_two, gain_group_one, gain_group_two, cnloh_group_one, cnloh_group_two, amp_group_one, amp_group_two, chromothrip_group_one, chromothrip_group_two, amplicon_single_group_one, amplicon_single_group_two, amplicon_step_group_one, amplicon_step_group_two,pvalue_gain, pvalue_loss, pvalue_cnloh, .ind)][,.ind:=NULL][]

write.table(we, file="NAME_OF_FILE.txt", sep="\t", quote=F)


for (i in 1:23) {
  grp1 <- subset(ds_group_one, ds_group_one$Chromosome == i)
  grp2 <- subset(ds_group_two, ds_group_two$Chromosome == i)

  print(min(grp1$Start))
  print(min(grp2$Start))
  print(max(grp1$Start))
  print(max(grp2$Start))

}

