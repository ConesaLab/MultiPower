############################################################################
########             Power study for TCGA-GBM                       ########
########             Nature Communications                          ########
############################################################################

## By Sonia Tarazona
## 14-Dic-2017
## Last modified: 28-nov-2018


options(stringsAsFactors = FALSE)

# Loading MultiPower functions
source("~/Dropbox/STATegra_STC/FiguresOfMerit/MultiPower/R/MultiOmicsPower10.R")



myfigdir = "~/Dropbox/STATegra_STC/FiguresOfMerit/powerResults/"

miscolores = c("red2", "orchid4", "dodgerblue3", "darkolivegreen4")
omicas = c("expression", "methylation", "miRNAs", "proteomics")
names(miscolores) = omicas



# Required data -----------------------------------------------------------

setwd("~/Dropbox/STATegra_STC/FiguresOfMerit/powerResults/")
load("dataTCGA.RData", verbose = TRUE)
# tcgadata
# tcgadesign
# DEresults


# Figures for estimated parameters per omic -------------------------------------

## Number of features, %DE features

par(mar = c(1,7,2,2), mfrow = c(1,2))
bb = barplot(sapply(tcgadata, nrow)/1000, las = 1,
             col = miscolores, horiz = TRUE, beside = TRUE, log = "x", 
             border = miscolores, main = "# features", axes = FALSE)
text(200, bb, sapply(tcgadata, nrow), adj = 1, 
     col = c(1, "white", rep(1,3)), las = 2)

par(mar = c(1,2,2,2))

bb = barplot(round(100*sapply(DEresults, function (x) sum(x[,"adj.P.Val"] < 0.05))/sapply(tcgadata, nrow),2), 
             las = 1, col = miscolores, horiz = TRUE, names.arg = "", 
             border = miscolores, main = "% DE features", axes = FALSE)
text(10, bb, paste0(round(100*sapply(DEresults, function (x) sum(x[,"adj.P.Val"] < 0.05))/sapply(tcgadata, nrow),0), "%"), 
     adj = 1, col = rep("white",5), las = 2)





# Global MultiPower results -----------------------------------------------

mytype = rep(2,4)
p1omics = c(0.64, 0.14, 0.56, 0.44)

tcgadesign = lapply(tcgadesign, function (x) x$type)

## d0 = 0.8
par(mfrow = c(1,2))
tcgaResultsEQ = MultiPower(data = tcgadata, groups = tcgadesign, type = mytype, d0 = 0.8, p1 = p1omics,
                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
                           max.size = 300, omicCol = miscolores)
tcgaResultsEQ$summary[,c(2:4,6:9,14:16)]


save(tcgaResultsEQ, tcgaResultsNE, file = "TCGAresultsMultiPower5.RData")





# Post-analysis -----------------------------------------------------------

TCGApostEQ = postMultiPower(data = tcgadata, groups = tcgadesign, optResults = tcgaResultsNE, 
                            max.size = 20, omicCol = miscolores) 



# Generating figures for the paper ----------------------------------------

setwd("~/Dropbox/STATegra_STC/FiguresOfMerit/powerResults")
load("TCGAresultsMultiPower5.RData", verbose = TRUE)


## Supplementary Figure
par(mfrow = c(1,2))
powerPlot(parameters = tcgaResultsEQ$parameters, optimalSampleSize = tcgaResultsEQ$optimalSampleSize, 
          omicCol = miscolores)





# Filtering low variability features in methylation -----------------------

table(tcgadesign$methylation)
dim(tcgadata$methylation)

head(DEresults$methylation)
plot(density(abs(DEresults$methylation$logFC)))

min(abs(DEresults$methylation$logFC))
min(abs(DEresults$methylation$logFC[DEresults$methylation$adj.P.Val < 0.05]))

sum(abs(DEresults$methylation$logFC) < 0.05) # 77020  (20%)
dim(DEresults$methylation)
# 384349
sum(DEresults$methylation$adj.P.Val < 0.05)
# 53064

## New proportion of DE methylation
53064 / (384349 - 77020) # = 17%

tcgadata2 = tcgadata
tcgadata2$methylation = tcgadata2$methylation[abs(DEresults$methylation$logFC) >= 0.05,]

sapply(tcgadata2, dim)




# Global MultiPower results AFTER filtering -----------------------------------------------

p1omics = c(0.64, 0.17, 0.56, 0.44)

## d0 = 0.8
par(mfrow = c(1,2))
tcgaResultsEQ2 = MultiPower(data = tcgadata2, groups = tcgadesign, type = mytype, d0 = 0.8, p1 = p1omics,
                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
                           max.size = 300, omicCol = miscolores)
tcgaResultsEQ2$summary[,c(2:4,6:9,14:16)]


save(tcgaResultsEQ, tcgaResultsNE, tcgaResultsEQ2, file = "TCGAresultsMultiPower5.RData")





# Post-analysis -----------------------------------------------------------
load("dataTCGA.RData", verbose = TRUE)
load("TCGAresultsMultiPower5.RData", verbose = TRUE)

TCGApostEQ = postMultiPower(data = tcgadata, groups = tcgadesign, optResults = tcgaResultsEQ, 
                            max.size = 30, omicCol = miscolores) 







# Answers to reviewers ---------------------------

p1omics
# 0.64 0.17 0.56 0.44


## Selecting n samples and computing DE --> do this R=5 times

myN = c(22, 20, 15, 10, 5) 
resum = NULL

for (r in 1:5) {
  
  dades = tcgadata
  disseny = tcgadesign
  
  print(paste0("r = ", r))
  
  for (n in myN) {
    print(paste0("n = ", n))
    
    for (i in 1:length(tcgadata)) {
      pro = which(disseny[[i]] == "Proneural")
      mes = which(disseny[[i]] == "Mesenchymal")
      selec = c(sample(pro, n), sample(mes, n))
      dades[[i]] = dades[[i]][,selec]
      disseny[[i]] = disseny[[i]][selec]
    }
    print(sapply(dades, ncol))
    
    DEres = lapply(1:4, function (i) apply(dades[[i]], 1, function (x) t.test(x ~ disseny[[i]])$p.value))
    DEresFDR = lapply(DEres, p.adjust, method = "fdr")
    tmp = sapply(DEresFDR, function (x) round(sum(x < 0.05)/length(x),2))
    resum = rbind(resum, c(r, n, tmp))
  }
  
}

colnames(resum) = c("r", "n", names(tcgadata))


head(resum)
resum[resum[,"n"] == 22,]
resum[resum[,"n"] == 20,]
resum[resum[,"n"] == 15,]
resum[resum[,"n"] == 10,]
resum[resum[,"n"] == 5,]

mediana = aggregate(resum[,3:6], by = list("n" = resum[,"n"]), median)
mediana

mitjana = aggregate(resum[,3:6], by = list("n" = resum[,"n"]), mean)

minim = aggregate(resum[,3:6], by = list("n" = resum[,"n"]), min)
minim

maxim = aggregate(resum[,3:6], by = list("n" = resum[,"n"]), max)
maxim

png(filename = "validationTCGA.png", pointsize = 22)
matplot(mitjana[,"n"], mitjana[,-1]*100, type = "o", col = miscolores, pch = 16, lwd = 4, ylim = c(0,30),
        xlab = "Sample size", ylab = "Mean % DE features", main = "TCGA glioblastoma data")
legend("topleft", omicas, col = miscolores, pch = 16, bty = "n")
dev.off()

