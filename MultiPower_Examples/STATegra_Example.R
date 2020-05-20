############################################################################
########             Power study for STATegra data types            ########
########             Nature Communications                          ########
############################################################################

## By Sonia Tarazona
## 05-Oct-2017
## Last modified: 28-nov-2018


options(stringsAsFactors = FALSE)

# Loading MultiPower functions
source("~/Dropbox/STATegra_STC/FiguresOfMerit/MultiPower/R/MultiOmicsPower10.R")


# Set up
myfigdir = "~/Dropbox/STATegra_STC/FiguresOfMerit/powerResults/"

miscolores = c("red2", "darkslategray3", "azure4", "sienna2",
               "dodgerblue4", "darkolivegreen4")
omicas = c("RNA-seq", "miRNA-seq", "ChIP-seq", "DNase-seq",  "Metabolomics", "Proteomics")
names(miscolores) = omicas




# Required data -----------------------------------------------------------
statdata = statdesign = vector("list")

## RNA-seq (LOG2-transformed)
statdata$RNAseq = read.delim("~/Dropbox/STATegra_STC/STATegraData/RNAseq/STATegra.RNAseq.CQN.Combat.Annotated.positive_2014_09.csv",
                             row.names = 1, as.is = TRUE, sep = ",")[,-1]
statdata$RNAseq = statdata$RNAseq[,c(grep("_Ik_0H", colnames(statdata$RNAseq)),
                                     grep("_Ik_24H", colnames(statdata$RNAseq)))]
statdesign$RNAseq = rep(c("0h", "24h"), each = 3)
min(statdata$RNAseq) # 1.14
statdata$RNAseq = statdata$RNAseq - min(statdata$RNAseq) # min = 0


## miRNA-seq (LOG2-transformed)
statdata$miRNAseq = read.delim("~/Dropbox/STATegra_STC/STATegraData/miRNA/miRNAseq_FinalData_NOtecrep.txt", 
                               row.names = 1, as.is = TRUE)
statdata$miRNAseq = statdata$miRNAseq[,c(grep("I.0.", colnames(statdata$miRNAseq)),
                                         grep("I.24.", colnames(statdata$miRNAseq)))]
statdesign$miRNAseq = rep(c("0h", "24h"), each = 3)
min(statdata$miRNAseq) # 0.08775962
statdata$miRNAseq = statdata$miRNAseq - min(statdata$miRNAseq)  # min = 0



## ChIP-seq 
statdata$ChIPseq = read.delim("~/Dropbox/STATegra_STC/STATegraData/ChIPseq/STAT_CPM_ChIPseq_FILT.txt",
                              row.names = 1)
statdata$ChIPseq = log2(statdata$ChIPseq + 1)
statdesign$ChIPseq = rep(c("0h", "24h"), each = 2)
min(statdata$ChIPseq) # 0


## DNase-seq (LOG2-transformed) 
statdata$DNaseSeq = read.delim("~/Dropbox/STATegra_STC/STATegraData/DNaseSeq/STAT_DNaseSeq_homer_RPKM_TMM_ARSyN.txt", 
                               row.names = 1, as.is = TRUE)
statdata$DNaseSeq = statdata$DNaseSeq[,c(grep("Ik0h_", colnames(statdata$DNaseSeq)),
                                         grep("Ik24h_", colnames(statdata$DNaseSeq)))]
statdesign$DNaseSeq = statdesign$RNAseq
min(statdata$DNaseSeq) # 1
statdata$DNaseSeq = statdata$DNaseSeq - min(statdata$DNaseSeq)  # min = 0 


## Metabolomics 
statdata$metabolomics = read.delim("~/Dropbox/STATegra_STC/STATegraData/metabolomics/oldData_original.txt", 
                                   as.is = TRUE, header = TRUE, row.names = 1)
statdata$metabolomics = statdata$metabolomics[,c(grep("I.0.", colnames(statdata$metabolomics)),
                                                 grep("I.24.", colnames(statdata$metabolomics)))]
statdata$metabolomics = log2(statdata$metabolomics)
statdesign$metabolomics = statdesign$RNAseq
min(statdata$metabolomics)  # 8.82e-05


## Proteomics (LOG2-transformed) 
statdata$proteomics = read.delim("~/Dropbox/STATegra_STC/STATegraData/proteomics/proteomicsImputedTMM.txt", 
                                 as.is = TRUE, header = TRUE, row.names = 1)
statdata$proteomics = statdata$proteomics[,c(grep("IKA_0h", colnames(statdata$proteomics)),
                                             grep("IKA_24h", colnames(statdata$proteomics)))]
statdesign$proteomics = statdesign$RNAseq
min(statdata$proteomics) # 1.97
statdata$proteomics = statdata$proteomics - min(statdata$proteomics)


sapply(statdata, dim)
#      RNAseq miRNAseq  ChIPseq DNaseSeq metabolomics proteomics
# [1,]  12762      469   23875    52788           60       1077
# [2,]      6        6       4        6            6          6

names(statdata) = names(statdesign) = omicas



# Global MultiPower results -----------------------------------------------

par(mfrow = c(1,2))
type1 = rep(2,6)
p1omics = c(0.4, 0.2, 0.2, 0.2, 0.6, 0.2)


statResultsEQ = MultiPower(data = statdata, groups = statdesign, type = type1, d0 = 0.8, p1 = p1omics,
                         omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
                         max.size = 200, omicCol = miscolores, dispPerc = 75)

statResultsNE = MultiPower(data = statdata, groups = statdesign, type = type1, d0 = 0.8, p1 = p1omics,
                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = FALSE,
                           max.size = 200, omicCol = miscolores, dispPerc = 75)

mycost = c(1, 1.3, 1.5, 1.6, 1, 1)
statResultsNEcost = MultiPower(data = statdata, groups = statdesign, type = type1, d0 = 0.8, p1 = p1omics,
                               omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = mycost, equalSize = FALSE,
                               max.size = 200, omicCol = miscolores, dispPerc = 75)


save(statResultsEQ, statResultsNE, statResultsNEcost,
     file = "STATegraResultsMultiPower6.RData")




# Post-analysis -----------------------------------------------------------

STATpostEQ = postMultiPower(data = statdata, groups = statdesign, optResults = statResultsEQ, 
                            max.size = 4, omicCol = miscolores) 

STATpostNE = postMultiPower(data = statdata, groups = statdesign, optResults = statResultsNE, 
                            max.size = 5, omicCol = miscolores) 

STATpostNEcost = postMultiPower(data = statdata, groups = statdesign, optResults = statResultsNEcost,
                                max.size = 3, omicCol = miscolores) 


save(statResultsEQ, statResultsNE, statResultsNEcost,
     STATpostEQ, STATpostNE, 
     file = "STATegraResultsMultiPower6.RData")

load("STATegraResultsMultiPower6.RData", verbose = TRUE)

par(mfrow = c(1,2))
postPowerPlot(postResults = STATpostEQ, equalSize = TRUE, omicCol = miscolores, max.size = 3)
postPowerPlot(postResults = STATpostNE, equalSize = FALSE, omicCol = miscolores, max.size = 5)
postPowerPlot(postResults = STATpostNEcost, equalSize = FALSE, omicCol = miscolores, max.size = 5)


powerPlot(statResultsNEcost$parameters, statResultsNEcost$optimalSampleSize, omicCol = miscolores)




# Generating figures for the paper ----------------------------------------

setwd("~/Dropbox/STATegra_STC/FiguresOfMerit/powerResults")
load("STATegraResultsMultiPower6.RData")

pdf(file = "plots/stategraEQ.pdf", width = 3.5*2, height = 3.5*2)
par(mfrow = c(2,2))
powerPlot(parameters = statResultsEQ$parameters, optimalSampleSize = statResultsEQ$optimalSampleSize, 
          omicCol = miscolores)
postPowerPlot(postResults = STATpostEQ, equalSize = TRUE, omicCol = miscolores, max.size = 4)
dev.off()


pdf(file = "plots/stategraNE.pdf", width = 4.5*2, height = 4.5)
par(mfrow = c(1,2))
powerPlot(parameters = statResultsNE$parameters, optimalSampleSize = statResultsNE$optimalSampleSize, omicCol = miscolores)
dev.off()

pdf(file = "plots/stategraNEcost.pdf", width = 4.5*2, height = 4.5)
par(mfrow = c(1,2))
powerPlot(parameters = statResultsNEcost$parameters, optimalSampleSize = statResultsNEcost$optimalSampleSize, omicCol = miscolores)
dev.off()





# Figures for estimated parameters per omic -------------------------------------

## Number of features, %DE features

par(mar = c(1,7,2,2), mfrow = c(1,2))

bb = barplot(sapply(statResultsEQ$parameters, function (x) x$m)/1000, las = 1,
             col = miscolores, horiz = TRUE, beside = TRUE, log = "x", 
             border = miscolores, main = "# features", axes = FALSE)
text(0.5, bb, sapply(statResultsEQ$parameters, function (x) x$m), adj = 1, 
     col = 1, las = 2)

par(mar = c(1,2,2,2))
bb = barplot(100*sapply(statResultsEQ$parameters, function (x) x$p1), 
             las = 1, col = miscolores, horiz = TRUE, names.arg = "", 
             border = miscolores, main = "% DE features", axes = FALSE)
text(10, bb, paste0(round(100*sapply(statResultsEQ$parameters, function (x) x$p1),0), "%"), 
     adj = 1, col = "white", las = 2)



## Dispersion

par(mar = c(7,4,2,1))#, mfrow = c(1,2))
boxplot(lapply(statResultsEQ$parameters, function (x) x$allDispersions), col = miscolores,
        las = 2, main = "Pooled standard deviation", log = "y", ylab = "log-scale (dispersion)")




