######################################################################################
###### MULTIPOWER
###### Optimization model to maximize power of multi-omics integration models
######################################################################################


## By Sonia Tarazona and David Gómez-Cabrero
## 05-Oct-2017
## Last modified: -Nov-2018
# New alpha computations (see getPower function)


require(RnaSeqSampleSize)

#### PACKAGES READ
# install.packages("slam")
# install.packages("lpmodeler")
require(lpmodeler)
require(Rsymphony)
# 1) Install SYMPHONY:  (now 5.6.16, first MultiPower version was with 5.6.10)
# cd
# svn checkout https://projects.coin-or.org/svn/SYMPHONY/releases/5.6.16 SYMPHONY-5.6.16
# cd SYMPHONY-5.6.16
# ./configure
# make
# make install
# 2) Install other components:
# sudo apt-get install coinor-libcgl-dev coinor-libclp-dev coinor-libcoinutils-dev coinor-libosi-dev
# sudo apt-get install coinor-libsymphony-dev
# sudo apt-get install autotools-dev
# 3) Intall package in R (previously downloaded from CRAN):
# install.packages("Rsymphony_0.1-26.tar.gz", repos = NULL)

# https://cran.r-project.org/src/contrib/Rsymphony_0.1-26.tar.gz
# https://projects.coin-or.org/SYMPHONY
#Download: wget http://www.coin-or.org/download/source/SYMPHONY/SYMPHONY-5.6.6.tgz


# require(boot)

geomean = function (x) exp(mean(log(x)))

cohen.h = function (p) abs(2*asin(sqrt(p[1])) - 2*asin(sqrt(p[2]))) # p is a vector with two components


# Estimating parameters needed for power calculation ----------------------
## Two-group comparison
paramEst = function (data, groups, type = 1, d0 = 0.8, p1 = 0.2, dispPerc = 75) {
  
  # type = 1 (counts), 2 (gaussian), 3 (binary variables: 0/1 or FALSE/TRUE)
  
  dispPerc = dispPerc/100
  
  # Number of features and DE features
  M = nrow(data)
  m1 = round(M*p1, 0)
  
  # Sample size per group
  nGroup = table(groups)
  
  # Sequencing depth correction for count data
  if (type == 1) {
    seqdepth = colSums(data)
    med = median(seqdepth)
    normD = seqdepth/med
    normDgroup = tapply(normD, INDEX = groups, geomean)
    myW = normDgroup[2]/normDgroup[1]
    data  = apply(data, 2, function (x) med*x/sum(x))
  } 
  
  
  # Mean per group (or p1 and p2 for binary data)
  meanPerGroup = t(apply(data, 1, tapply, INDEX = groups, mean, na.rm = TRUE))
  
  if (type != 3) {
    
    # Standard deviation per group
    sdPerGroup = t(apply(data, 1, tapply, INDEX = groups, sd, na.rm = TRUE))
    sdPerGroup = sdPerGroup[,names(nGroup)]
    
    # Pooled Standard Deviation
    SDpooled = sqrt((nGroup[1]*sdPerGroup[,1]^2 + nGroup[2]*sdPerGroup[,2]^2)/(sum(nGroup)-2))
    
    # Cohen's d per feature
    deltaPerFeat = abs(meanPerGroup[,1] - meanPerGroup[,2])
    d = deltaPerFeat/SDpooled
  } else {
    d = apply(meanPerGroup, 1, cohen.h)
  }
  
  
  # DEfeat
  DEfeat = which(d >= d0)
  if (length(DEfeat) < 10) stop("Less than 10 DE features for d=d0. Please, decrease d0 value. \n")
  if (length(DEfeat)/M > 0.9) stop("More than 90% of the features are DE for d=d0. Please, increase d0 value. \n")
  
  meanPerGroup = meanPerGroup[DEfeat,]
  if (type != 3) {
    sdPerGroup = sdPerGroup[DEfeat,]
    SDpooled = SDpooled[DEfeat]
    deltaPerFeat = deltaPerFeat[DEfeat]
    
    # Intervals P70-P80
    qqq = quantile(SDpooled, probs = c(max(0, dispPerc-0.05), min(1, dispPerc+0.05)), na.rm = TRUE)
    seleSDpooled = intersect(which(SDpooled <= qqq[2]), which(SDpooled >= qqq[1]))
  }
  
  
  if (type == 1) {   ## COUNT DATA
    cat("Parameters are to be estimated for count data \n")
    
    if(min(data, na.rm = TRUE) < 0) stop("Negative values were found. Are you sure these are count data?\n")
    
    # Fold-change estimation
    allFC = apply(meanPerGroup, 1, function (x) max(0.0000001, x[2]) / max(x[1], 0.0000001))
    abslogFC =  abs(log2(allFC))
    # minFC = quantile(allFC[seleSDpooled], probs = 1-dispPerc, na.rm = TRUE)
    minFC = quantile(abslogFC, probs = 1-dispPerc, na.rm = TRUE)
    minFC = which.min(abs(abslogFC - minFC))
    minFC = allFC[minFC]
    
    # Average counts
    # meanCounts = quantile(as.numeric(meanPerGroup[seleSDpooled,]), probs = 1-dispPerc, na.rm = TRUE)
    meanCounts = quantile(meanPerGroup[,1], probs = 1-dispPerc, na.rm = TRUE)*normDgroup[1]
    
    # Dispersion estimation
    # allPhis = (as.numeric(sdPerGroup)^2 - as.numeric(meanPerGroup)) / as.numeric(meanPerGroup)^2
    allPhis = ((sdPerGroup)^2 - meanPerGroup) / meanPerGroup^2
    allPhis = apply(allPhis, 1, max)
    maxPhi = quantile(allPhis, probs = dispPerc, na.rm = TRUE)
    # var75 = quantile(as.numeric(sdPerGroup), probs = dispPerc, na.rm = TRUE)^2
    # maxPhi = max(0.000001, (var75 - meanCounts) / meanCounts^2) 

    # Estimated parameter for count data
    myparameters = list("type" = type, "allDispersions" = allPhis, "dispersion" = maxPhi,
                        "p1" = p1, "d" = d0, "delta" = NA, "minFC" = minFC, "meanCounts" = meanCounts,
                        "m" = M, "m1" = m1, "alld" = d, "dispPerc" = dispPerc, "w" = myW)
    
  } 
  
  if (type == 2) {   ## NORMAL DATA
    cat("Parameters are to be estimated for normally distributed data \n")
    
    # Dispersion
    sdValue = quantile(SDpooled, probs = dispPerc, na.rm = TRUE)
    
    # Delta estimation
    myDelta = quantile(deltaPerFeat[seleSDpooled], probs = 1-dispPerc, na.rm = TRUE)
    
    # Estimated parameters for normal data
    myparameters = list("type" = type, "allDispersions" = SDpooled, "dispersion" = sdValue, 
                        "p1" = p1, "d" = d0, "delta" = round(myDelta,4), "minFC" = NA, "meanCounts" = NA,
                        "m" = M, "m1" = m1, "alld" = d, "dispPerc" = dispPerc, "w" = NA)
  }
  
  if (type == 3) {  ## BINARY DATA
    cat("Parameters are to be estimated for binary data \n")
    
    # Cohen's h for DE feature
    h = d[DEfeat]
    qqq = quantile(h, probs = 1-dispPerc, na.rm = TRUE)
    qqq = which.min(abs(h - qqq))
    prop = meanPerGroup[qqq,]
    
    # Estimated parameters for normal data
    myparameters = list("type" = type, "allDispersions" = NA, "dispersion" = NA, 
                        "p1" = p1, "d" = d0, "delta" = prop, "minFC" = NA, "meanCounts" = NA,
                        "m" = M, "m1" = m1, "alld" = d, "dispPerc" = dispPerc, "w" = NA)
  }
  
  return(myparameters)
  
}



# Computing power or sample size given the rest of parameters -------------

getPower = function (parameters, power = NULL, n = NULL, fdr = 0.05, alpha = 0.05) {
  
  if (is.null(fdr)) fdr = alpha*(parameters$m - parameters$m1)/(parameters$m1 + alpha*(parameters$m- parameters$m1))
  
  if (is.null(power)) { # Compute power for given n
    if (is.null(n)) stop("Please, indicate a value for either power or n arguments. \n")
    
    if (parameters$type == 1) { # COUNT DATA
      potencia = est_power(n = n, w = parameters$w, rho = parameters$minFC, lambda0 = parameters$meanCounts, 
                           phi0 = parameters$dispersion, f = fdr, m = parameters$m, m1 = parameters$m1)
      
    } 
    
    if (parameters$type == 2) {  # NORMAL DATA
      r1 = parameters$m1
      myAlpha = min(r1 * fdr / ((parameters$m - parameters$m1) * (1 - fdr)), 0.05)  ## NEW
      potencia = power.t.test(n = n, delta = parameters$delta, sd = parameters$dispersion, 
                              sig.level = myAlpha, type = "two.sample", alternative = "two.sided")$power
      # print(potencia); print(parameters$dispersion); print(n); print(parameters$delta); print(myAlpha)
    }
    
    
    if (parameters$type == 3) {  # BINARY DATA
      r1 = parameters$m1
      myAlpha = min(r1 * fdr / ((parameters$m - parameters$m1) * (1 - fdr)), 0.05)  ## NEW
      potencia = power.prop.test(n = n, p1 = parameters$delta[1], p2 = parameters$delta[2], 
                                 sig.level = myAlpha, alternative = "two.sided")$power
    }
    
    
    return(potencia)
    
  } else {  # Compute n for given power
    
    if (parameters$type == 1) { # COUNT DATA
      tamany = sample_size(power = power, m = parameters$m, m1 = parameters$m1, f = fdr,
                           k = 1, w = parameters$w, rho = parameters$minFC, lambda0 = parameters$meanCounts, 
                           phi0 = parameters$dispersion)
      
    } 
    
    if (parameters$type == 2) {  # NORMAL DATA
      # Alpha estimation
      # r1 = power * parameters$m1
      r1 = parameters$m1
      myAlpha = min(r1 * fdr / ((parameters$m - parameters$m1) * (1 - fdr)), 0.05)  ## NEW
      
      tamany = try(power.t.test(power = power, delta = parameters$delta, sd = parameters$dispersion, 
                            sig.level = myAlpha, type = "two.sample", alternative = "two.sided")$n, silent = TRUE)
      if (inherits(tamany, "try-error")) tamany = 2
    }
    
    
    if (parameters$type == 3) {  # BINARY DATA
      
      r1 = parameters$m1
      myAlpha = min(r1 * fdr / ((parameters$m - parameters$m1) * (1 - fdr)), 0.05)  ## NEW
      
      tamany = try(power.prop.test(power = power, p1 = parameters$delta[1], p2 = parameters$delta[2],
                                   sig.level = myAlpha, alternative = "two.sided")$n, silent = TRUE)
      if (inherits(tamany, "try-error")) tamany = 2
    }
    
    return(max(ceiling(tamany), 2))
    
  }
  
}





# Optimal sample size -----------------------------------------------------

optimalRep = function (parameters, omicPower = 0.6, averagePower = 0.85, fdr = 0.05, alpha = 0.05, cost = 1, 
                       equalSize = TRUE, max.size = 200) {
  
  omics = names(parameters)
  
  if (length(omicPower) == 1) omicPower = rep(omicPower, length(omics))
  names(omicPower) = omics
  
  if (equalSize) {  ## Same sample size for all omics
    
    # Compute n for each omic
    n1 = sapply(omics, function (oo) getPower(parameters[[oo]], power = omicPower[oo], n = NULL, fdr = fdr, alpha = alpha))
    names(n1) = omics
    
    n1max = max(n1, 2, na.rm = TRUE)
    
    allPowers = sapply(omics, function (oo) getPower(parameters[[oo]], power = NULL, n = n1max, fdr = fdr, alpha = alpha))
    n2 = n1max
    if (n2 > max.size) stop("Maximum size allowed has been exceed. 
                            Please, increase max.size parameter to get the optimal sample size. \n")
    
    # Compute n to satisfy global power
    while(sum(allPowers)/length(omics) < averagePower) {
      n2 = n2 + 1
      if (n2 > max.size) stop("Maximum size allowed has been exceed. 
                            Please, increase max.size parameter to get the optimal sample size. \n")
      allPowers = sapply(omics, function (oo) getPower(parameters[[oo]], power = NULL, n = n2, fdr = fdr, alpha = alpha))
    }
    
    return(list("n0" = n1, "n" = n2, "finalPower" = allPowers, "fdr" = fdr, "alpha" = alpha,
                "omicPower" = omicPower, "averagePower" = averagePower, "cost" = cost))
  
  } else {   ## Different sample size for each omic

    sss = optiSSnotEqual(parameters, fdr, alpha, cost, max.size, omicPower, averagePower)
    n2 = as.numeric(sss$summary[,"SampleSize"])
    allPowers = as.numeric(sss$summary[,"Power"])
    names(allPowers) = names(n2) = sss$summary[,"Omic"]
    
    return(list("n0" = NA, "n" = n2, "finalPower" = allPowers, "fdr" = fdr, "alpha" = alpha,
                "omicPower" = omicPower, "averagePower" = averagePower, "cost" = cost))
    
  }
  
}







# Summary of results ------------------------------------------------------

powerSummary = function(parameters, optimalSampleSize) {
  tabla = data.frame("omic" = names(parameters), "type" = sapply(parameters, function (x) x$type),
                     "numFeat" = sapply(parameters, function (x) x$m),
                     "DEperc" = sapply(parameters, function (x) x$p1),
                     "CohenD" = sapply(parameters, function (x) x$d),
                     "delta" = sapply(parameters, function (x) paste(x$delta, collapse = "-")),
                     "minFC" = sapply(parameters, function (x) x$minFC),
                     "meanCounts" = sapply(parameters, function (x) x$meanCounts),
                     "dispersion" = sapply(parameters, function (x) round(x$dispersion,4)),
                     # "FDR" = optimalSampleSize$fdr,
                     "minPower" = optimalSampleSize$omicPower,
                     "averPower" = optimalSampleSize$averagePower,
                     "cost" = optimalSampleSize$cost,
                     "minSampleSize" = optimalSampleSize$n0,
                     "optSampleSize" = optimalSampleSize$n,
                     "power" = round(optimalSampleSize$finalPower,4))
  print(tabla)
  return(tabla)
}






# Plots for power study ---------------------------------------------------

powerPlot = function(parameters, optimalSampleSize, omicCol = NULL) {
  
  if (is.null(omicCol)) {
    if (length(parameters) > 12) {
      stop("Too many omics to be plotted. Please, select a lower number of omics to plot. \n")
    }
    omicCol = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)]
    omicCol = omicCol[1:length(parameters)]
  } 
  
  omicShape = 1:length(parameters)
  names(omicCol) = names(omicShape) = names(parameters)
  
  
  ## Power versus Sample Size
  
  # Sample Sizes
  nmax = max(optimalSampleSize$n)
  xMin = 2
  xMax = round(max(nmax+20, (3*nmax - xMin)/2),0)
  xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)), optimalSampleSize$n)))
  
  # Powers
  yValues = matrix(NA, ncol = length(parameters), nrow = length(xValues))
  rownames(yValues) = xValues
  colnames(yValues) = names(parameters)
  
  for (i in 1:nrow(yValues)) {
    for (j in 1:ncol(yValues)) {
      yValues[i,j] = getPower(parameters[[j]], power = NULL, n = xValues[i], fdr = optimalSampleSize$fdr, alpha = optimalSampleSize$alpha)
    }
  }
  
  matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", ylab = "Statistical power",
          main = "Power vs Sample Size", col = omicCol, lty = omicShape)
  optiSS = optimalSampleSize$n
  if (length(optiSS) == 1) optiSS = rep(optiSS, length(parameters))
  points(optiSS, diag(yValues[as.character(optiSS),]), pch = 15, col = omicCol, cex = 1.2)
  legend("bottomright", names(parameters), lwd = 2, col = omicCol, lty = omicShape, bty = "n")
  
  
  ## Power vs Dispersion
  
  # Quantiles of dispersion
  xValues = seq(0,1,0.05)
  
  # Powers
  yValues2 = matrix(NA, ncol = length(parameters), nrow = length(xValues))
  rownames(yValues2) = xValues
  colnames(yValues2) = names(parameters)
  
  parameters2 = parameters
  
  optiSS = optimalSampleSize$n
  if (length(optiSS) == 1) optiSS = rep(optiSS, length(parameters))
  
  for (i in 1:nrow(yValues2)) {
    for (j in 1:ncol(yValues2)) {
      if (parameters2[[j]]$type != 3) {  ## no dispersion plot for binary variables
        parameters2[[j]]$dispersion = quantile(parameters[[j]]$allDispersions, probs = xValues[i], na.rm = TRUE)
        if (parameters2[[j]]$dispersion < 0) parameters2[[j]]$dispersion = NA
        if (!is.na(parameters2[[j]]$dispersion)) {
          # print(paste(i,j))
          yValues2[i,j] = getPower(parameters2[[j]], power = NULL, n = optiSS[j], fdr = optimalSampleSize$fdr, alpha = optimalSampleSize$alpha)
        } else {
          yValues2[i,j] = NA
        }
      }
    }
  }
  
  if (!all(is.na(yValues2))) {
    matplot(xValues*100, yValues2, type = "l", lwd = 2, xlab = "Dispersion percentiles", ylab = "Statistical power",
            main = "Power vs Dispersion", col = omicCol, lty = omicShape)
    points(rep(parameters2[[1]]$dispPerc*100, length(parameters)), 
           as.numeric(yValues2[as.character(c(parameters2[[1]]$dispPerc)),]), 
           pch = 15, col = omicCol, cex = 1.2)
    legend("bottomleft", names(parameters), lwd = 2, col = omicCol, lty = omicShape, bty = "n")
  }
  
  ## Data to plot
  return(list("PowerVsSsampleSize" = yValues,
              "PowerVsDispersion" = yValues2))
  
}




PowerDispersionPlot = function(n = 5, parameters, fdr = 0.05, alpha = 0.05, omicCol = NULL) {
  
  if (length(n) < length(parameters)) n = rep(n[1], length(parameters))
  
  if (is.null(omicCol)) {
    if (length(parameters) > 12) {
      stop("Too many omics to be plotted. Please, select a lower number of omics to plot. \n")
    }
    omicCol = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)]
    omicCol = omicCol[1:length(parameters)]
  } 
  
  omicShape = 1:length(parameters)
  names(omicCol) = names(omicShape) = names(parameters)
  
  ## Power vs Dispersion
  
  # Quantiles of dispersion
  xValues = seq(0,1,0.05)
  
  # Powers
  yValues2 = matrix(NA, ncol = length(parameters), nrow = length(xValues))
  rownames(yValues2) = xValues
  colnames(yValues2) = names(parameters)
  
  parameters2 = parameters
  
  optiSS = n

  for (i in 1:nrow(yValues2)) {
    for (j in 1:ncol(yValues2)) {
      if (parameters2[[j]]$type != 3) {
        parameters2[[j]]$dispersion = quantile(parameters[[j]]$allDispersions, probs = xValues[i], na.rm = TRUE)
        if (!is.na(parameters2[[j]]$dispersion)) {
          yValues2[i,j] = getPower(parameters2[[j]], power = NULL, n = optiSS[j], fdr = fdr, alpha = alpha)
        } else {
          yValues2[i,j] = NA
        }
      }
    }
  }
  
  if (!all(is.na(yValues2))) {
    matplot(xValues*100, yValues2, type = "l", lwd = 2, xlab = "Dispersion percentiles", ylab = "Statistical power",
            main = "Power vs Dispersion", col = omicCol, lty = omicShape)
    legend("bottomleft", names(parameters), lwd = 2, col = omicCol, lty = omicShape, bty = "n")
  }
  
  ## Data to plot
  return(yValues2)
  
}





# Computing optimal sample size when it is not equal for all omics ----------------------------------------------------------

optiSSnotEqual = function (parameters, fdr = 0.05, alpha = 0.05, cost = 1, max.size = 100, 
                           omicPower = 0.6, averagePower = 0.8) {
  
  ##### GENERATION OF MATRICES FOR THE PROBLEM
  
  K = length(parameters)  # number of omics
  
  # if (length(delta.vector) < K) delta.vector = rep(delta.vector[1], K)
  # if (length(sd.vector) < K) sd.vector = rep(sd.vector[1], K)
  if (length(cost) < K) cost = rep(cost[1], K)
  if (length(max.size) < K) max.size = rep(max.size[1], K)
  if (length(omicPower) < K) omicPower = rep(omicPower[1], K)
  # if (length(omic.vars) < K) omic.vars = rep(omic.vars[1], K)  ## number of features per omic -> 1 if no multiple testing correction is to be considered
  
  num.var = sum(max.size) - length(max.size)  ## sample size = 1 is not considered
  # coeffs of objective function; coeffs power; constraint sum(vars) = 1, to have only 1 sample size, coeffs average power 
  my.a = my.A2 = my.A3 = my.A4 = NULL    
  
  for (k in 1:K) {
    coef.power = coef1 = rep(0, num.var)
    
    for (i in 2:max.size[k]) {
      my.a = c(my.a, cost[k]*i*2)  # coefficients of objective function
      
      # power of each (omic, sample size)
      my.power = getPower(parameters[[k]], power = NULL, n = i, fdr = fdr, alpha = alpha)
      
      # coeff average power
      my.A4 = c(my.A4, my.power)
      
      if (k == 1) {
        coef.power[i-1] = my.power
        coef1[i-1] = 1
      }  
      
      if (k > 1) {
        coef.power[i-k+sum(max.size[1:(k-1)])] = my.power      
        coef1[i-k+sum(max.size[1:(k-1)])] = 1
      } 
            
    }
    my.A2 = rbind(my.A2, coef.power)
    my.A3 = rbind(my.A3, coef1)
  }  
  

  ##### DEFINITION OF THE PROBLEM
  
  prob<-newProblem(max=F)
  
  #### Add variables
  #my.a
  for(i in 1:length(my.a))
  {
    prob <- addVariable(prob, "B", my.a[i])
  }
  
  #### Add inequalities
  #A2 = my.A2
  #b2 = min.power
  
  for(i in 1:nrow(my.A2))
  {
    prob <- addConstraint(prob, sense=">=", rhs = omicPower[i],coefs=my.A2[i,])
  }
  
  # constraint average power
  prob <- addConstraint(prob, sense=">=", rhs = averagePower*K, coefs = my.A4)
  
  
  #### Add equalities
  #A3 = my.A3
  b3 = rep(1,nrow(my.A3))
  
  for(i in 1:nrow(my.A3))
  {
    prob <- addConstraint(prob, sense="==", rhs=b3[i],coefs=my.A3[i,])
  }
  
  #### Check consistency
  # print(prob)
  # print(checkDims(prob))
  
  
  #### SOLUTION
  mysol = mipSolve(prob)
  
  myvars = unlist(sapply(1:K, function (k) paste(names(parameters)[k], 2:max.size[k], sep = "=")))
  mysolution = myvars[which(mysol$solution == 1)]
  mysolution = as.data.frame(do.call("rbind", strsplit(mysolution, "=")), stringsAsFactors = FALSE)
  colnames(mysolution) = c("Omic", "SampleSize")
  mysolution = data.frame(mysolution, "OmicCost" = cost*as.numeric(mysolution[,2])*2,
                          "Power" = diag(my.A2[,mysol$solution == 1]))
  
  mysol$summary = mysolution

  return(mysol)
  
}






# Wrapper function: MULTIPOWER --------------------------------------------

MultiPower = function(data, groups, type, d0 = 0.8, p1 = 0.2, omicPower = 0.6, averagePower = 0.85, dispPerc = 75,
                      fdr = 0.05, alpha = 0.05,cost = 1, equalSize = TRUE, max.size = 200, omicCol = NULL, powerPlots = TRUE) {
  
  if (length(p1) == 1) p1 = rep(p1, length(data))
  if (length(d0) == 1) d0 = rep(d0, length(data))
  
  parameters = lapply(1:length(data), function (i) {
    cat(paste0("Estimating parameters for omic: ", names(data)[i], " \n")) 
    paramEst(data[[i]], groups[[i]], type[i], d0[i], p1[i], dispPerc = dispPerc)})
  names(parameters) = names(data)
  
  cat("Computing optimal sample size... \n")
  optimalSampleSize = optimalRep(parameters, omicPower, averagePower, fdr, alpha, cost, equalSize, max.size)
  
  resum = powerSummary(parameters, optimalSampleSize)
  
  if (powerPlots) {
    cat("Generating power plots... \n")
    data2plot = powerPlot(parameters, optimalSampleSize, omicCol)
  } else { data2plot = NULL }
  
  return(list("parameters" = parameters,
              "optimalSampleSize" = optimalSampleSize,
              "summary" = resum,
              "data2plot" = data2plot))
  
}




# postMultiPower ----------------------------------------------------------

postMultiPower = function(data, groups, optResults, max.size = 5, omicCol = NULL, dispPerc = 75) {
  
  p1 = optResults$summary$DEperc
  type = optResults$summary$type
  d0 = unique(optResults$summary$CohenD)
  if (length(d0) > 1) stop("Effect size (d0) must be the same for all the omics in order to use this function.")
  dmax = min(sapply(optResults$parameters, function (x) {quantile(na.omit(x$alld), probs = 0.9)}), na.rm = TRUE)
  
  if (d0 >= dmax) stop("It is not possible to increase Cohen's d to reduce sample size.")
  
  equalSS = TRUE
  if (sum(is.na(optResults$summary$minSampleSize))  == nrow(optResults$summary)) equalSS = FALSE 
  
  lasD = unique(round(seq(d0, dmax, (dmax-d0)/5),2))
  
  mySize = myPower =  matrix(NA, nrow = length(data), ncol = length(lasD))
  colnames(mySize) = colnames(myPower) = paste0("d=", lasD)
  rownames(mySize) = rownames(myPower) = names(data)
  
  mySize[,1] = optResults$summary$optSampleSize
  myPower[,1] = optResults$summary$power
  
  for (j in 2:length(lasD)) {
    
    print(paste("Cohen's d =", lasD[j]))
    parameters = lapply(1:length(data), function (i) paramEst(data[[i]], groups[[i]], type[i], d0 = lasD[j], p1[i], dispPerc = dispPerc))
    names(parameters) = names(data)
    
    tmp = optimalRep(parameters, omicPower = optResults$summary$minPower,
                     averagePower = optResults$summary$averPower[1], 
                     fdr = optResults$summary$FDR[1], cost = optResults$summary$cost, 
                     equalSize = equalSS, max.size = max(optResults$summary$optSampleSize, na.rm = TRUE))
    
    mySize[,j] = tmp$n
    myPower[,j] = tmp$finalPower
  }
  
  # Post-results
  myresult = list("SampleSize" = mySize, "Power" = myPower, "d" = lasD)
  
  # Plot
  postPowerPlot(postResults = myresult, equalSize = equalSS, omicCol = omicCol, max.size = max.size)
  
  return(myresult)
              
}




# postPowerPlot -----------------------------------------------------------

postPowerPlot = function(postResults, equalSize, omicCol = NULL, max.size = 10) {
  
  if (is.null(omicCol)) {
    omicCol = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)]
    omicCol = omicCol[1:nrow(postResults$SampleSize)]
  } 
  names(omicCol) = rownames(postResults$SampleSize)
  
  if (min(postResults$SampleSize) > max.size) {
    cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
    max.size = min(postResults$SampleSize)
    cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
  }
  
  if (equalSize) {
    
    fff = approxfun(postResults$SampleSize[1,], postResults$d)
    myD = round(fff(max.size),2)
    
    plot(postResults$d, postResults$SampleSize[1,], type = "l", lwd = 2, xlab = "Cohen's d", 
         ylab = "Number of replicates", main = "Sample size vs Cohen's d", 
         ylim = c(2, max(postResults$SampleSize)))
    arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size, lty = 2, col = 2)
    arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2, col = 2)
    text(min(postResults$d)+0.5, max(postResults$SampleSize, na.rm = TRUE), paste0("Cohen's d = ", myD), col = 2)
    
  } else {
    
    fff = sapply(1:nrow(postResults$SampleSize), function (i) { fff = approxfun(postResults$SampleSize[i,], postResults$d)
                                                                return(fff(max.size)) })
    myD = round(max(fff, na.rm = TRUE),2)
    
    matplot(postResults$d, t(postResults$SampleSize), type = "l", lwd = 2, col = omicCol, lty = 1, 
            xlab = "Cohen's d", ylab = "Number of replicates", main = "Sample size vs Cohen's d",
            ylim = c(2, max(postResults$SampleSize)))
    legend("topright", rownames(postResults$SampleSize), lwd = 2, col = omicCol, bty = "n")
    arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size, lty = 2, col = 1)
    arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2, col = 1)
    text(min(postResults$d)+0.5, max(postResults$SampleSize, na.rm = TRUE), paste0("Cohen's d = ", myD), col = 1)
  }
  
  mypar = par()
  suppressWarnings(par(xpd = TRUE, mar = c(6.2,4,3,0.8)))
  barplot(t(postResults$Power[,c(1,min(which(postResults$d >= myD)))]), col = rep(omicCol, each = 2), 
          beside = TRUE, las = 2, ylab = "Statistical power", ylim = c(0,1),
          border = rep(omicCol, each = 2), density = rep(c(30,100), nrow(postResults$Power)))
  legend(x = 0.5, y = 1.2, c("Optimal SS", "User's SS"), col = 1, density = c(30,100), ncol = 2, bty = "n")
  suppressWarnings(par(mypar))
  
}




