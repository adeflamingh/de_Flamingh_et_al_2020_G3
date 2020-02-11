# Rx_identifier for ancient elephants
# This script was modified from Mittnik et al 2016 (original author: Chuan-Chao Wang)
# Based on the ratio of X chromosome-derived shotgun sequencing data to the autosomal coverage to establish the probability of an XX or XY karyotype for ancient samples.
# example:samtools view -q 30 -b sampleID.fa.merged.bam > SampleID_q30.bam
        # samtools index SampleID_q30.bam
        # samtools idxstats SampleID_q30.bam > SampleID.idxstats
        # Rscript Rx_identifier.r SampleID > SampleID.Rx
		

args=(commandArgs(TRUE))
PREFIX=as.character(args[1])

idxstats<-read.table(paste(PREFIX,'.idxstats',sep=''),header=F,nrows=28,row.names=1) #lox.afr4 has 27 autosomal chr and chr X, therefore nrows = 27 + 1 = 28 - modify this value according to the chr number of the taxon of interest
c1 <- c(as.numeric(idxstats[,1])) #this column in the idxstats file indicates the number of reference reads
c2 <- c(as.numeric(idxstats[,2])) #this column in the idxstats file indicates the number of mapped reads
total_ref <- sum(c1)
total_map <- sum(c2)
  
LM <- lm(c1~c2)
summary(LM)  

#calculate the normalized ratio of each chr:
Rt1 <- (idxstats[1,2]/total_map)/(idxstats[1,1]/total_ref)
Rt2 <- (idxstats[2,2]/total_map)/(idxstats[2,1]/total_ref)
Rt3 <- (idxstats[3,2]/total_map)/(idxstats[3,1]/total_ref)
Rt4 <- (idxstats[4,2]/total_map)/(idxstats[4,1]/total_ref)
Rt5 <- (idxstats[5,2]/total_map)/(idxstats[5,1]/total_ref)
Rt6 <- (idxstats[6,2]/total_map)/(idxstats[6,1]/total_ref)
Rt7 <- (idxstats[7,2]/total_map)/(idxstats[7,1]/total_ref)
Rt8 <- (idxstats[8,2]/total_map)/(idxstats[8,1]/total_ref)
Rt9 <- (idxstats[9,2]/total_map)/(idxstats[9,1]/total_ref)
Rt10 <- (idxstats[10,2]/total_map)/(idxstats[10,1]/total_ref)
Rt11 <- (idxstats[11,2]/total_map)/(idxstats[11,1]/total_ref)
Rt12 <- (idxstats[12,2]/total_map)/(idxstats[12,1]/total_ref)
Rt13 <- (idxstats[13,2]/total_map)/(idxstats[13,1]/total_ref)
Rt14 <- (idxstats[14,2]/total_map)/(idxstats[14,1]/total_ref)
Rt15 <- (idxstats[15,2]/total_map)/(idxstats[15,1]/total_ref)
Rt16 <- (idxstats[16,2]/total_map)/(idxstats[16,1]/total_ref)
Rt17 <- (idxstats[17,2]/total_map)/(idxstats[17,1]/total_ref)
Rt18 <- (idxstats[18,2]/total_map)/(idxstats[18,1]/total_ref)
Rt19 <- (idxstats[19,2]/total_map)/(idxstats[19,1]/total_ref)
Rt20 <- (idxstats[20,2]/total_map)/(idxstats[20,1]/total_ref)
Rt21 <- (idxstats[21,2]/total_map)/(idxstats[21,1]/total_ref)
Rt22 <- (idxstats[22,2]/total_map)/(idxstats[22,1]/total_ref)
Rt23 <- (idxstats[23,2]/total_map)/(idxstats[23,1]/total_ref)
Rt24 <- (idxstats[24,2]/total_map)/(idxstats[24,1]/total_ref)
Rt25 <- (idxstats[25,2]/total_map)/(idxstats[25,1]/total_ref)
Rt26 <- (idxstats[26,2]/total_map)/(idxstats[26,1]/total_ref)
Rt27 <- (idxstats[27,2]/total_map)/(idxstats[27,1]/total_ref)
Rt28 <- (idxstats[28,2]/total_map)/(idxstats[28,1]/total_ref)

#calculate averaged normalized ratio of the X chr
tot <- c(Rt28/Rt1,Rt28/Rt2,Rt28/Rt3,Rt28/Rt4,Rt28/Rt5,Rt28/Rt6,Rt28/Rt7,Rt28/Rt8,Rt28/Rt9,Rt28/Rt10,Rt28/Rt11,Rt28/Rt12,Rt28/Rt13,Rt28/Rt14,Rt28/Rt15,Rt28/Rt16,Rt28/Rt17,Rt28/Rt18,Rt28/Rt19,Rt28/Rt20,Rt28/Rt21,Rt28/Rt22,Rt28/Rt23,Rt28/Rt24,Rt28/Rt25,Rt28/Rt26,Rt28/Rt27)
Rx <- mean(tot)
cat("Rx :",Rx,"\n")
confinterval <- 1.96*(sd(tot)/sqrt(27))
CI1 <- Rx-confinterval
CI2 <- Rx+confinterval
cat("95% CI :",CI1, CI2,"\n")

if (CI1 > 0.8) {print ("Sex assignment:The sample should be assigned as Female")
} else if (CI2 < 0.6) {print ("Sex assignment:The sample should be assigned as Male")
} else if (CI1 > 0.6 & CI2 > 0.8) {print ("Sex assignment:The sample is consistent with XX but not XY")
} else if (CI1 < 0.6 & CI2 < 0.8) {print ("Sex assignment:The sample is consistent with XY but not XX")
} else print ("Sex assignment:The sample could not be assigned")

print ("***It is important to realize that the assignment is invalid, if there is no correlation between the number of reference reads and that of the mapped reads***")


