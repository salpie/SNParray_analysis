

################################################################################################################################
########## similarities & differences between pairs ############################################################################
################################################################################################################################
################################################################################################################################
######################################  Salpie Nowinski     date:27/03/2015           ##########################################
################################################################################################################################
# calculate similarities and differences of absolute copy number at break points between pairs of samples ######################


library(plyr)
library(data.table)


###########  put in pairs into folders in same order
group_one <- (Sys.glob("*segmentCN.txt")) # add path to *segmentCN.txt files produced from TAPS 2.0
group_two <- (Sys.glob("*segmentCN.txt")) # add path to *segmentCN.txt files produced from TAPS 2.0


# functions to calculate copy number changes 

compareTo2Amp <- function(group_one, group_two) {
  as.numeric(group_one > 4 & group_two > 4)
}

compareTo2Gain <- function(group_one, group_two) {
  as.numeric(group_one > 2 & group_two > 2)
}

compareTo2Loss <- function(group_one, group_two) {
  as.numeric(group_one < 2 & group_two < 2)
}

compareTo2Gain_two <- function(group_one, group_two) {
  as.numeric(group_one < 3 & group_two > 2)
}

compareTo2Gain_one <- function(group_one, group_two) {
  as.numeric(group_one > 2 & group_two < 3)
}

compareTo2Loss_two <- function(group_one, group_two) {
  as.numeric(group_one > 1 & group_two < 2)
}

compareTo2Loss_one <- function(group_one, group_two) {
  as.numeric(group_one < 2 & group_two > 1)
}

compareTo2Gain_higher <- function(group_one, group_two) {
  as.numeric(group_one > 2 & group_two < 3)
}

compareToCNLOH_sim <- function(group_one, group_two, group_one_mcn, group_two_mcn){
	as.numeric(group_one == 2 & group_two == 3)
}

ishigher <- function(group_one, group_two) {
	as.numeric(group_one > 2 & group_two > group_one)
}


findAmp <- function(group_one, group_two) { 
  as.numeric(group_one >= 5 & group_two >=5)
}

findAmp_one <- function(group_one, group_two) {
	as.numeric(group_one >= 5 & group_two < 5)
}

findAmp_two <- function(group_one, group_two) {
	as.numeric(group_two >= 5 & group_one < 5)
}


group_one_file <- read.table(group_one[1], header=TRUE)
names <- colnames(group_one_file)

for (i in 1:length(group_one)) {
#for (i in 1:3) {
	print(group_one[i])

	group_one_file <- read.table(group_one[i], header=TRUE)
	colnames(group_one_file) <- names

	group_two_file <- read.table(group_two[i], header=TRUE)
	colnames(group_two_file) <- names
	lst1 <- split(group_one_file, group_one_file$Chromosome)
	lst2 <- split(group_two_file, group_two_file$Chromosome)

	require(survival)

	for (z in 1:length(lst1)) {
		if (nrow(lst1[[z]]) == 1 && nrow(lst2[[z]]) == 1) {
			if (lst1[[z]]$Start > lst2[[z]]$Start) {
				lst2[[z]]$Start <- lst1[[z]]$Start
			}
			else if (lst2[[z]]$Start > lst1[[z]]$Start){
				lst1[[z]]$Start <- lst2ww[[z]]$Start
			}
		}
	}


 # merge
 	df <- do.call(rbind, mapply(FUN = function(x, y) {
 
   		x$event <- y$event <- 0
   		idc.spl <- survSplit(x, cut=y$End, start='Start', end='End', event='event')
   		dcis.spl <- survSplit(y, cut=x$End, start='Start', end='End', event='event')
   		mrg <- merge(idc.spl, dcis.spl, 
                by=c('Chromosome', 'Start', 'End'), 
                #all=TRUE, 
                suffixes = c("_group_one","_group_two"))
   		mrg[c('Chromosome', 'Start', 'End', 'Cn_group_one', 'Cn_group_two', 'mCn_group_one', 'mCn_group_two')]
 		},
 	lst1, lst2, SIMPLIFY=FALSE))

 	df$loss_sim <- compareTo2Loss(df$Cn_group_one, df$Cn_group_two)
	df$gain_sim <- compareTo2Gain(df$Cn_group_one, df$Cn_group_two)
	df$x <- ifelse(df$gain_sim == 1,paste(group_one[i], group_two[i]), "")
	df$gain_samples_sim <- paste(df$gain_samples_sim,df$x)
  	df$x <- ifelse(df$loss_sim == 1,paste(group_one[i], group_two[i]), "")
	df$loss_samples_sim <- paste(df$loss_samples_sim,df$x)
  #df$x <- ifelse(df$cnloh == 1,paste(group_one[i]), "")
  #group_one_file$cnloh_samples <- paste(group_one_file$cnloh_samples,group_one_file$x)


	df$loss_grp1 <- compareTo2Loss_one(df$Cn_group_one, df$Cn_group_two)
	df$x <- ifelse(df$loss_grp1 == 1,paste(group_one[i]), "")
	df$loss_grp1_samples <- paste(df$loss_grp1_samples,df$x)


	df$gain_grp1 <- compareTo2Gain_one(df$Cn_group_one, df$Cn_group_two)
	df$x <- ifelse(df$gain_grp1 == 1,paste(group_one[i]), "")
	df$gain_grp1_samples <- paste(df$gain_grp1_samples,df$x)


	df$loss_grp2 <- compareTo2Loss_two(df$Cn_group_one, df$Cn_group_two)
	df$x <- ifelse(df$loss_grp2 == 1,paste(group_two[i]), "")
	df$loss_grp2_samples <- paste(df$loss_grp2_samples,df$x)


	df$gain_grp2 <- compareTo2Gain_two(df$Cn_group_one, df$Cn_group_two)
	df$x <- ifelse(df$gain_grp2 == 1,paste(group_two[i]), "")
	df$gain_grp2_samples <- paste(df$gain_grp2_samples,df$x)


	df$ishigher <- ishigher(df$Cn_group_one, df$Cn_group_two)
	
	df$cnloh_sim <- compareToCNLOH_sim(df$Cn_group_one, df$Cn_group_two, df$mCn_group_one, df$mCn_group_two)
	
	df$amp_sim <- findAmp(df$Cn_group_one, df$Cn_group_two)
	df$x <- ifelse(df$amp_sim == 1,paste(group_one[i], group_two[i]), "")
	df$amp_samples_sim <- paste(df$amp_samples_sim,df$x)
	
	df$amp_one <- findAmp_one(df$Cn_group_one, df$Cn_group_two)
	df$x <- ifelse(df$amp_one == 1,paste(group_one[i]), "")
	df$amp_samples_one <- paste(df$amp_samples_one,df$x)
	
	df$amp_two <- findAmp_two(df$Cn_group_one, df$Cn_group_two)
	df$x <- ifelse(df$amp_two == 1,paste(group_two[i]), "")
	df$amp_samples_two <- paste(df$amp_samples_two,df$x)

	#make dataframes of each pair

	if (i < 2){
		print("The first one is working!!!")
	}
	if (i == 1){
		please_work <- df
	}
	if (i > 1) {
	lst1 <- split(df, df$Chromosome)
	lst2 <- split(please_work, please_work$Chromosome)


for (w in 1:length(lst1)) {
		if (nrow(lst1[[w]]) == 1 && nrow(lst2[[w]]) == 1) {
			if (lst1[[w]]$Start > lst2[[w]]$Start) {
				lst1[[w]]$Start <- lst2[[w]]$Start
			}
			else if (lst2[[w]]$Start > lst1[[w]]$Start){
				lst2[[w]]$Start <- lst1[[w]]$Start
			}
		}
	}
	
print(nrow(please_work))

	please_work <- do.call(rbind, mapply(FUN = function(x, y) {
 
   x$event <- y$event <- 0
   idc.spl <- survSplit(x, cut=y$End, start='Start', end='End', event='event')
   dcis.spl <- survSplit(y, cut=x$End, start='Start', end='End', event='event')
   mrg <- merge(idc.spl, dcis.spl, 
                by=c('Chromosome', 'Start', 'End'), 
                #all=TRUE, 
                suffixes = c("_one","_two"))
   mrg[c('Chromosome', 'Start', 'End', 'loss_sim_one', 'loss_sim_two', 'gain_sim_one', 'gain_sim_two','loss_grp1_one', 'loss_grp2_one', 'gain_grp1_one', 'gain_grp2_one', 'loss_grp1_two', 'loss_grp2_two', 'gain_grp1_two', 'gain_grp2_two', 'gain_grp1_samples_one', 'gain_grp1_samples_two', 'gain_grp2_samples_one', 'gain_grp2_samples_two',  'loss_grp1_samples_one', 'loss_grp1_samples_two', 'loss_grp2_samples_one', 'loss_grp2_samples_two', 'loss_samples_sim_one', 'gain_samples_sim_one', 'loss_samples_sim_two', 'gain_samples_sim_two', 'ishigher_one', 'ishigher_two', 'cnloh_sim_one', 'cnloh_sim_two', 'amp_sim_one', 'amp_sim_two', 'amp_samples_sim_one', 'amp_samples_sim_two', 'amp_one_one', 'amp_one_two', 'amp_two_one', 'amp_two_two', 'amp_samples_one_one', 'amp_samples_one_two', 'amp_samples_two_one', 'amp_samples_two_two')]
 },
 lst1, lst2, SIMPLIFY=FALSE))

	please_work$loss_sim <- please_work$loss_sim_one + please_work$loss_sim_two
 	please_work$gain_sim <- please_work$gain_sim_one + please_work$gain_sim_two
 	please_work$loss_grp1 <- please_work$loss_grp1_one + please_work$loss_grp1_two
 	please_work$loss_grp1_samples <- paste(please_work$loss_grp1_samples, please_work$loss_grp1_samples_one, please_work$loss_grp1_samples_two)
 	please_work$gain_grp1 <- please_work$gain_grp1_one + please_work$gain_grp1_two
 	please_work$gain_grp1_samples <- paste(please_work$gain_grp1_samples, please_work$gain_grp1_samples_one, please_work$gain_grp1_samples_two)
 	please_work$loss_grp2 <- please_work$loss_grp2_one + please_work$loss_grp2_two
 	please_work$loss_grp2_samples <- paste(please_work$loss_grp2_samples, please_work$loss_grp2_samples_one, please_work$loss_grp2_samples_two)
 	please_work$gain_grp2 <- please_work$gain_grp2_one + please_work$gain_grp2_two
 	please_work$gain_grp2_samples <- paste(please_work$gain_grp2_samples, please_work$gain_grp2_samples_one, please_work$gain_grp2_samples_two)

 	please_work$loss_samples_sim <- paste(please_work$loss_samples_sim, please_work$loss_samples_sim_one, please_work$loss_samples_sim_two)
 	please_work$gain_samples_sim <- paste(please_work$gain_samples_sim, please_work$gain_samples_sim_one, please_work$gain_samples_sim_two)
 	please_work$ishigher <- please_work$ishigher_one + please_work$ishigher_two
 	please_work$cnloh_sim <- please_work$cnloh_sim_one + please_work$cnloh_sim_two
 	please_work$amp_sim <- please_work$amp_sim_one + please_work$amp_sim_two
 	please_work$amp_samples_sim <- paste(please_work$amp_samples_sim, please_work$amp_samples_sim_one, please_work$amp_samples_sim_two)
 	please_work$amp_one <- please_work$amp_one_one + please_work$amp_one_two
 	please_work$amp_samples_one <- paste(please_work$amp_samples_one, please_work$amp_samples_one_one, please_work$amp_samples_one_two)
 	please_work$amp_two <- please_work$amp_two_one + please_work$amp_two_two
 	please_work$amp_samples_two <- paste(please_work$amp_samples_two, please_work$amp_samples_two_one, please_work$amp_samples_two_two)

 	please_work <- please_work[,c(1,2,3,44:63)]

	}
}


newdata <- please_work[order(please_work$Chromosome, please_work$Start),]


newdata$group_one_loss <- newdata$loss_sim + newdata$loss_grp1
newdata$group_one_loss_final <- length(group_one) - newdata$group_one_loss
newdata$group_two_loss <- newdata$loss_sim + newdata$loss_grp2
newdata$group_two_loss_final <- length(group_one) - newdata$group_two_loss
newdata$group_one_gain <- newdata$gain_sim + newdata$gain_grp1
newdata$group_one_gain_final <- length(group_one) - newdata$group_one_gain
newdata$group_two_gain <- newdata$gain_sim + newdata$gain_grp2
newdata$group_two_gain_final <- length(group_one) - newdata$group_two_gain

#calculate fisher's exact p-value

newdata$pvalue_loss <- apply(as.matrix(newdata[,c(24,25,26,27)]), 1, function(x) 
         fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)


newdata$pvalue_gain <- apply(as.matrix(newdata[,c(28,29,30,31)]), 1, function(x) 
         fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)



############### if you want a raw txt file use the following change "file" location to where you want it saved and name it whatever you want

write.table(newdata, file="/Users/salpienowinski/Documents/Tables_for_Paper/RAW TABLES/Paired_cILC_invcLCIS_raw.txt", sep="\t", quote=F)


############## if you want a compact txt files use the following change "file" location to where you want it saved and name it whatever you want

we_raw <- setDT(newdata)[, .ind:= cumsum(c(TRUE,Start[-1]!=End[-.N])),
        list(loss_sim, gain_sim, loss_grp1, gain_grp1, loss_grp2, gain_grp2, loss_samples_sim, gain_samples_sim, pvalue_loss, pvalue_gain, ishigher)][,
       list(chr=Chromosome[1], start=Start[1], stop=End[.N]),
       list(loss_sim, gain_sim, loss_grp1, gain_grp1, loss_grp2, gain_grp2, loss_samples_sim, gain_samples_sim, pvalue_loss, pvalue_gain, ishigher, .ind)][,.ind:=NULL][]

write.table(we_raw, file="/Users/salpienowinski/Documents/Tables_for_Paper/RAW TABLES/Paired_cILC_invcLCIS_raw.txt", sep="\t", quote=F)
