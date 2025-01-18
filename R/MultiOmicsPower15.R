######################################################################################
###### MULTIPOWER
###### Optimization model to maximize power of multi-omics integration models
######################################################################################


## By Sonia Tarazona and David Gomez-Cabrero
## 05-Oct-2017
## Last modified: March-2023
# 




#### PACKAGES READ
# install.packages("FDRsampsize")
require(FDRsampsize)
require(lpSolve)


# install.packages("slam")
# install.packages("lpmodeler")
# require(lpmodeler)
# require(Rsymphony)
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



# Auxiliary functions -----------------------------------------------------

geomean = function (x) exp(mean(log(x)))

cohen.h = function (p) abs(2*asin(sqrt(p[1])) - 2*asin(sqrt(p[2]))) # p is a vector with two components

power.binary = function(n, sig.level=0.05, p1_p2, p1) {
  if (length(p1_p2) != length(p1)) stop("p1_p2 and p1 must be of the same length")
  p2 = p1 - p1_p2
  mypower = sapply(1:length(p1), function (i) power.prop.test(n=n, p1 = p1[i], p2 = p2[i], sig.level = sig.level)$power)
  return(mypower)
}




# Estimating parameters needed for power calculation ----------------------

## Two-group comparison
paramEst = function (data, groups, type = 1) {
  
  # type = 1 (counts), 2 (gaussian), 3 (binary variables: 0/1 or FALSE/TRUE)
  
  # Sample size per group
  nGroup = table(groups)
  
  sd0 = apply(data, 1, sd)
  sd0 = which(sd0 == 0)
  if (length(sd0) > 0) {
    print(paste0(length(sd0), " constant features are to be removed from the analysis."))
    data = data[-sd0,]
  }

  # Number of features
  M = nrow(data)
  
  # Sequencing depth correction for count data
  if (type == 1) {
    seqdepth = colSums(data)
    med = median(seqdepth)
    data  = apply(data, 2, function (x) med*x/sum(x))
  }


  # Mean counts per group 
  meanPerGroup = t(apply(data, 1, tapply, INDEX = groups, mean, na.rm = TRUE))   

  if (type != 3) {

    # Standard deviation per group
    sdPerGroup = t(apply(data, 1, tapply, INDEX = groups, sd, na.rm = TRUE))
    sdPerGroup = sdPerGroup[,names(nGroup)]

    # Pooled Standard Deviation
    SDpooled = sqrt((nGroup[1]*sdPerGroup[,1]^2 + nGroup[2]*sdPerGroup[,2]^2)/(sum(as.numeric(nGroup))-2))

    # Cohen's d per feature
    deltaPerFeat = abs(meanPerGroup[,1] - meanPerGroup[,2])
    d = deltaPerFeat/SDpooled
  } else {
    d = apply(meanPerGroup, 1, cohen.h)
  }
  

  if (type == 1) {   ## COUNT DATA
    cat("Parameters are to be estimated for count data \n")

    if(min(data, na.rm = TRUE) < 0) stop("Negative values were found. Are you sure these are count data?\n")

    # Fold-change estimation
    allFC = log2(apply(meanPerGroup, 1, function (x) max(0.0000001, x[2]) / max(x[1], 0.0000001)))
    
    # Average counts
    mu = rowMeans(data, na.rm = TRUE)

    # CV 
    myCV = SDpooled/mu
    
    # Estimated parameters for count data
    myparameters = list("type" = type, "logFC" = allFC, "pooledSD" = SDpooled, "CV" = myCV,
                        "delta" = deltaPerFeat, "mu" = mu,"m" = M, "d" = d, "nGroup" = nGroup)

  }

  if (type == 2) {   ## NORMAL DATA
    cat("Parameters are to be estimated for normally distributed data \n")

    # Estimated parameters for normal data
    myparameters = list("type" = type, "delta" = deltaPerFeat, "pooledSD" = SDpooled, 
                        "m" = M, "d" = d, "nGroup" = nGroup)
  }

  if (type == 3) {  ## BINARY DATA
    cat("Parameters are to be estimated for binary data \n")

    # Estimated parameters for normal data
    myparameters = list("type" = type, "p1_p2" = meanPerGroup[,1]-meanPerGroup[,2], "p1" = meanPerGroup[,1],
                        "m" = M, "d" = d, "nGroup" = nGroup)
  }

  return(myparameters)

}



# Computing power or sample size given the rest of parameters -------------

getPower = function (parameters, power = NULL, n = NULL, fdr = 0.05, 
                     null.effect = 0, max.n = 500) {

  # Compute power for given n
  
  if (is.null(power)) { 
    if (is.null(n)) stop("Please, indicate a value for either power or n arguments. \n")

    if (parameters$type == 1) { # COUNT DATA (Negative Binomial)
      potencia = fdr.power(fdr = fdr, n = n, pow.func = power.hart, eff.size = parameters$logFC, null.effect = null.effect, 
                           mu = parameters$mu, sig = parameters$CV) 

    }

    if (parameters$type == 2) {  # NORMAL DATA
      potencia = fdr.power(fdr = fdr, n = n, pow.func = power.twosampt, eff.size = parameters$delta, null.effect = null.effect,
                           sigma = parameters$pooledSD)
    }


    if (parameters$type == 3) {  # BINARY DATA
      potencia = fdr.power(fdr = fdr, n = n, pow.func = power.binary, eff.size = parameters$p1_p2, null.effect = null.effect,
                           p1 = parameters$p1)
    }


    return(potencia)   } else {  # Compute n for given power
      
      if (parameters$type == 1) { # COUNT DATA
        
        tamany = fdr.sampsize(fdr = fdr, ave.pow = power, pow.func = power.hart, eff.size = parameters$logFC, 
                              null.effect = null.effect, mu = parameters$mu, sig = parameters$CV,
                              max.n = max.n, min.n = 2)$n

    }

    if (parameters$type == 2) {  # NORMAL DATA
      
      tamany = fdr.sampsize(fdr = fdr, ave.pow = power, pow.func = power.twosampt, eff.size = parameters$delta, 
                            null.effect = null.effect, sigma = parameters$pooledSD,
                            max.n = max.n, min.n = 2)$n
    }


    if (parameters$type == 3) {  # BINARY DATA

      tamany = fdr.sampsize(fdr = fdr, ave.pow = power, pow.func = power.binary, eff.size = parameters$p1_p2, 
                            null.effect = null.effect, p1 = parameters$p1, 
                            max.n = max.n, min.n = 2)$n

    }
      
      return(max(ceiling(tamany), 2))
      
    }
}



# Optimal sample size -----------------------------------------------------

optimalRep = function (parameters, omicPower = 0.6, averagePower = 0.85, fdr = 0.05, cost = 1,
                       equalSize = TRUE, max.size = 200, null.effect = 0) {

  omics = names(parameters)

  if (length(omicPower) == 1) omicPower = rep(omicPower, length(omics))
  names(omicPower) = omics

  if (equalSize) {  ## Same sample size for all omics

    # Compute n for each omic
    n1 = sapply(omics, function (oo) getPower(parameters[[oo]], 
                                              power = omicPower[oo], n = NULL, 
                                              fdr = fdr, max.n = max.size,
                                              null.effect = null.effect))
    names(n1) = omics

    n1max = max(n1, 2, na.rm = TRUE)

    allPowers = sapply(omics, function (oo) getPower(parameters[[oo]], 
                                                     power = NULL, n = n1max, 
                                                     fdr = fdr))
    n2 = n1max
    if (n2 > max.size) stop("Maximum size allowed has been exceed.
                            Please, increase max.size parameter to get the optimal sample size. \n")

    # Compute n to satisfy global power
    while(sum(allPowers)/length(omics) < averagePower) {     
      n2 = n2 + 1
      if (n2 > max.size) stop("Maximum size allowed has been exceed.
                            Please, increase max.size parameter to get the optimal sample size. \n")
      allPowers = sapply(omics, function (oo) getPower(parameters[[oo]], power = NULL, 
                                                       n = n2, fdr = fdr, 
                                                       null.effect = null.effect))
    }

    return(list("n0" = n1, "n" = n2, "finalPower" = allPowers, "fdr" = fdr, 
                "omicPower" = omicPower, "averagePower" = averagePower, "cost" = cost))

  } else {   ## Different sample size for each omic

    sss = optiSSnotEqual(parameters, fdr, cost, max.size, omicPower, averagePower, 
                         null.effect)
    n2 = as.numeric(sss[,"SampleSize"])
    allPowers = as.numeric(sss[,"Power"])
    names(allPowers) = names(n2) = sss[,"Omic"]

    return(list("n0" = NA, "n" = n2, "finalPower" = allPowers, "fdr" = fdr, 
                "omicPower" = omicPower, "averagePower" = averagePower, "cost" = cost))

  }

}







# Summary of results ------------------------------------------------------

powerSummary = function(parameters, optimalSampleSize) {
  
  tabla = data.frame("omic" = names(parameters), "type" = sapply(parameters, function (x) x$type),
                     "numFeat" = sapply(parameters, function (x) x$m),
                     "minCohenD" = sapply(parameters, function (x) round(min(x$d, na.rm = TRUE),2)),
                     "maxCohenD" = sapply(parameters, function (x) round(max(x$d, na.rm = TRUE),2)),
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
  ngroup = unique(as.numeric(sapply(parameters, function (x) x$nGroup)))
  xMin = 2
  xMax = round(max(nmax+20, (3*nmax - xMin)/2),0)
  xValues = c(round(seq(xMin, xMax, (xMax - xMin)/10)), optimalSampleSize$n)
  xValues = sort(unique(c(xValues, ngroup)))

  # Powers
  yValues = matrix(NA, ncol = length(parameters), nrow = length(xValues))
  rownames(yValues) = xValues
  colnames(yValues) = names(parameters)

  for (i in 1:nrow(yValues)) {
    for (j in 1:ncol(yValues)) {
      yValues[i,j] = getPower(parameters[[j]], power = NULL, n = xValues[i], fdr = optimalSampleSize$fdr)   ### null.effect
    }
  }

  matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", ylab = "Statistical power",
          main = "Power vs Sample Size", col = omicCol, lty = omicShape)
  optiSS = optimalSampleSize$n
  if (length(optiSS) == 1) optiSS = rep(optiSS, length(parameters))
  points(optiSS, diag(yValues[as.character(optiSS),]), pch = 15, col = omicCol, cex = 1.2)
  legend("bottomright", names(parameters), lwd = 2, col = omicCol, lty = omicShape, bty = "n")


  ## Power vs Effect Size

  # Quantiles of effect size
  xValues = seq(0,0.75,0.05)  # max P75 

  # Powers
  yValues2 = matrix(NA, ncol = length(parameters), nrow = length(xValues))
  rownames(yValues2) = xValues
  colnames(yValues2) = names(parameters)

  parameters2 = parameters

  optiSS = optimalSampleSize$n
  if (length(optiSS) == 1) optiSS = rep(optiSS, length(parameters))
  
  percentiles = lapply(parameters, function (x) {quantile(x$d, probs = xValues, na.rm = TRUE)})

  for (i in 1:nrow(yValues2)) {
    for (j in 1:ncol(yValues2)) {
      selefeat = which(parameters2[[j]]$d >= percentiles[[j]][i])
      parameters2[[j]]$d = parameters2[[j]]$d[selefeat]
      
      if (parameters2[[j]]$type == 1) {
        parameters2[[j]]$logFC = parameters2[[j]]$logFC[selefeat]
        parameters2[[j]]$pooledSD = parameters2[[j]]$pooledSD[selefeat]
        parameters2[[j]]$CV = parameters2[[j]]$CV[selefeat]
        parameters2[[j]]$delta = parameters2[[j]]$delta[selefeat]
        parameters2[[j]]$mu = parameters2[[j]]$mu[selefeat]
      }
      
      if (parameters2[[j]]$type == 2) {
        parameters2[[j]]$pooledSD = parameters2[[j]]$pooledSD[selefeat]
        parameters2[[j]]$delta = parameters2[[j]]$delta[selefeat]
      }
      
      if (parameters2[[j]]$type == 3) {
        parameters2[[j]]$p1_p2 = parameters2[[j]]$p1_p2[selefeat]
        parameters2[[j]]$p1 = parameters2[[j]]$p1[selefeat]
      }
      
      yValues2[i,j] = getPower(parameters2[[j]], power = NULL, n = optiSS[j], fdr = optimalSampleSize$fdr)
        
    }
  }

  if (!all(is.na(yValues2))) {
    matplot(xValues*100, yValues2, type = "l", lwd = 2, xlab = "Percentiles for effect size cutoff", ylab = "Statistical power",
            main = "Power vs Effect size", col = omicCol, lty = omicShape)
    points(rep(0, length(parameters)), as.numeric(yValues2[1,]),
           pch = 15, col = omicCol, cex = 1.2)
    legend("bottomright", names(parameters), lwd = 2, col = omicCol, lty = omicShape, bty = "n")
  }

  ## Data to plot
  return(list("PowerVsSsampleSize" = yValues,
              "PowerVsEffectSize" = yValues2))

}






# Computing optimal sample size when it is not equal for all omics ----------------------------------------------------------

optiSSnotEqual = function (parameters, fdr = 0.05, cost = 1, max.size = 100,
                           omicPower = 0.6, averagePower = 0.8, null.effect = 0) {

  ##### GENERATION OF MATRICES FOR THE PROBLEM

  K = length(parameters)  # number of omics

  if (length(cost) < K) cost = rep(cost[1], K)
  if (length(max.size) < K) max.size = rep(max.size[1], K)
  if (length(omicPower) < K) omicPower = rep(omicPower[1], K)

  num.var = sum(max.size) - length(max.size)  ## sample size = 1 is not considered
  # coeffs power;  
  # constraint sum(vars) = 1 (to have only 1 sample size); coeffs average power
  myC = NULL  # coeffs of objective function
  # A1: coeffs of power per omic
  # A3: sum(Zij) = 1 for each omic i
  A1 = A3 = matrix(0, nrow = K, ncol = num.var) 
  # A2: coeffs for average power for all omics
  A2 = NULL
  
  for (k in 1:K) {
    # coef.power --> A1, A2
    # coef1 --> A3
    coef.power = coef1 = rep(0, num.var)

    for (i in 2:max.size[k]) {
      myC = c(myC, cost[k]*i*2)  # coefficients of objective function

      # power of each (omic, sample size)
      my.power = getPower(parameters[[k]], power = NULL, n = i, fdr = fdr,
                          null.effect = null.effect, max.n = max.size)

      # coeff average power
      A2 = c(A2, my.power)

      if (k == 1) {   
        coef.power[i-1] = my.power
        coef1[i-1] = 1
      }

      if (k > 1) {
        coef.power[i-k+sum(max.size[1:(k-1)])] = my.power
        coef1[i-k+sum(max.size[1:(k-1)])] = 1
      }

    }
    A1[k,] = coef.power
    A3[k,] = coef1
  }
  
  A2 = matrix(A2, nrow = 1, byrow = TRUE)
  
  myA = rbind(A1, A2, A3)
  
  mydir = c(rep(">=", K+1), rep("=", K))
  
  myRHS = c(omicPower, K*averagePower, rep(1,K))


  #### SOLUTION OF THE PROBLEM
  mysol = lp(direction = "min", objective.in = myC, const.mat = myA,
             const.dir = mydir, const.rhs = myRHS, all.int = TRUE, all.bin = TRUE)
  
  if (mysol$status == 2) {
    
    stop("No feasible solution was found for these requirements.")
    
  } else {
    
    myvars = unlist(sapply(1:K, function (k) paste(names(parameters)[k], 
                                                   2:max.size[k], sep = "=")))
    mysolution = myvars[which(mysol$solution == 1)]
    mysolution = as.data.frame(do.call("rbind", strsplit(mysolution, "=")), 
                               stringsAsFactors = FALSE)
    colnames(mysolution) = c("Omic", "SampleSize")
    mysolution = data.frame(mysolution, "OmicCost" = cost*as.numeric(mysolution[,2])*2,
                            "Power" = diag(A1[,mysol$solution == 1]))

    return(mysolution)
    
  }

}






# Wrapper function: MULTIPOWER --------------------------------------------

MultiPower = function(data, groups, type, omicPower = 0.6, averagePower = 0.85, 
                      fdr = 0.05, cost = 1, equalSize = TRUE, max.size = 200, omicCol = NULL,
                      powerPlots = TRUE, null.effect = 0) {

  parameters = lapply(1:length(data), function (i) {
    cat(paste0("Estimating parameters for omic: ", names(data)[i], " \n"))
    paramEst(data[[i]], groups[[i]], type[i])})
  names(parameters) = names(data)

  cat("Computing optimal sample size... \n")
  optimalSampleSize = optimalRep(parameters, omicPower, averagePower, fdr, cost, 
                                 equalSize, max.size = max.size, null.effect)

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

postMultiPower = function(optResults, max.size = 5, omicCol = NULL) {

  omics = names(optResults$parameters)
  
  equalSS = TRUE
  if (sum(is.na(optResults$summary$minSampleSize))  == nrow(optResults$summary)) equalSS = FALSE
  
  maxD = round(min(optResults$summary$maxCohenD),1)
  maxD = min(c(maxD, 4))
  lasD = seq(0,maxD-0.1,0.1) 
  
  mySize = myPower = myM =  matrix(NA, ncol = length(omics), nrow = length(lasD))
  colnames(mySize) = colnames(myPower) = colnames(myM) = omics
  rownames(mySize) = rownames(myPower) = rownames(myM) = paste0("d=", lasD)

  mySize[1,] = optResults$summary$optSampleSize
  myPower[1,] = optResults$summary$power
  myM[1,] = sapply(optResults$parameters, function (x) x$m)
  
  parameters2 = optResults$parameters
  
  for (i in 2:nrow(mySize)) {
    for (j in 1:ncol(mySize)) {
      selefeat = which(parameters2[[j]]$d >= lasD[i])
      parameters2[[j]]$d = parameters2[[j]]$d[selefeat]
      
      if (parameters2[[j]]$type == 1) {
        parameters2[[j]]$logFC = parameters2[[j]]$logFC[selefeat]
        parameters2[[j]]$pooledSD = parameters2[[j]]$pooledSD[selefeat]
        parameters2[[j]]$CV = parameters2[[j]]$CV[selefeat]
        parameters2[[j]]$delta = parameters2[[j]]$delta[selefeat]
        parameters2[[j]]$mu = parameters2[[j]]$mu[selefeat]
      }
      
      if (parameters2[[j]]$type == 2) {
        parameters2[[j]]$pooledSD = parameters2[[j]]$pooledSD[selefeat]
        parameters2[[j]]$delta = parameters2[[j]]$delta[selefeat]
      }
      
      if (parameters2[[j]]$type == 3) {
        parameters2[[j]]$p1_p2 = parameters2[[j]]$p1_p2[selefeat]
        parameters2[[j]]$p1 = parameters2[[j]]$p1[selefeat]
      }
      
    }
    
    tmp = optimalRep(parameters2, omicPower = optResults$summary$minPower,
                     averagePower = optResults$summary$averPower[1],
                     fdr = optResults$optimalSampleSize$fdr, cost = optResults$summary$cost,
                     equalSize = equalSS, max.size = max(optResults$summary$optSampleSize, na.rm = TRUE))
    
    mySize[i,] = tmp$n
    myPower[i,] = tmp$finalPower
    myM[i,] = sapply(parameters2, function (x) length(x$d))
  }

  # Post-results
  myresult = list("SampleSize" = mySize, "Power" = myPower, "NumFeat" = myM, "d" = lasD)

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
  names(omicCol) = colnames(postResults$SampleSize)

  if (min(postResults$SampleSize) > max.size) {
    cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
    max.size = min(postResults$SampleSize)
    cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
  }

  if (equalSize) {

    fff = approxfun(postResults$SampleSize[,1], postResults$d)
    myD = round(fff(max.size),1)
    
    cat(paste0("For having a sample size of ", max.size, " and maintain the desired power, you need to remove features with Cohen's d below ", myD, ". \n"))
    cat("The number of remaining features in each omic is: \n")
    print(postResults$NumFeat[paste0("d=",myD),])

    plot(postResults$d, postResults$SampleSize[,1], type = "l", lwd = 2, xlab = "Cohen's d cutoff",
         ylab = "Number of replicates", main = "Sample size vs Cohen's d",
         ylim = c(2, max(postResults$SampleSize)))
    arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size, lty = 2, col = 2)
    arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2, col = 2)
    text(min(postResults$d)+1, max(postResults$SampleSize, na.rm = TRUE)-2, paste0("Cohen's d = ", myD), col = 2)

  } else {

    fff = sapply(1:ncol(postResults$SampleSize), function (i) { fff = approxfun(postResults$SampleSize[,i], postResults$d)
                                                                return(fff(max.size)) })
    myD = round(max(fff, na.rm = TRUE),2)

    matplot(postResults$d, postResults$SampleSize, type = "l", lwd = 2, col = omicCol, lty = 1,
            xlab = "Cohen's d", ylab = "Number of replicates", main = "Sample size vs Cohen's d",
            ylim = c(2, max(postResults$SampleSize)))
    legend("topright", colnames(postResults$SampleSize), lwd = 2, col = omicCol, bty = "n")
    arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size, lty = 2, col = 1)
    arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2, col = 1)
    text(mean(postResults$d), max(postResults$SampleSize, na.rm = TRUE), 
         paste0("Cohen's d = ", myD), col = 1, adj = 0.5)
  }

  mypar = par()
  suppressWarnings(par(xpd = TRUE, mar = c(6.2,4,3,0.8)))
  barplot(postResults$Power[c(1,min(which(postResults$d >= myD))),], col = rep(omicCol, each = 2),
          beside = TRUE, las = 2, ylab = "Statistical power", ylim = c(0,1),
          border = rep(omicCol, each = 2), density = rep(c(30,100), ncol(postResults$Power)))
  legend(x = 0.5, y = 1.2, c("Optimal SS", "User's SS"), col = 1, density = c(30,100), ncol = 2, bty = "n")
  suppressWarnings(par(mypar))

}





# Wrapper function MultiGroupPower ---------------------------------------------------------

MultiGroupPower = function(data, groups, type, comparisons = NULL,
                           omicPower = 0.6, averagePower = 0.85, 
                           fdr = 0.05, cost = 1, equalSize = TRUE, max.size = 200, omicCol = NULL,
                           powerPlots = FALSE, summaryPlot = TRUE) {
  
  grupsComuns = Reduce(intersect, groups)
  
  if (is.null(comparisons)) {  # Generating all possible comparisons
    
    comparisons = combn(grupsComuns, m = 2)
    
  } else {  # Checking that required comparisons are possible for all omics
    
    grupsComparats = unique(as.vector(comparisons))
    
    if (!setequal(grupsComparats, grupsComuns)) stop("Groups to be compared are not available for all omics.")
  }
  
  nomsCompa = apply(comparisons, 2, paste, collapse = "_")
  
  output = vector("list", length = length(nomsCompa)); names(output) = nomsCompa
  
  for (i in 1:length(output)) {
    
    cat(nomsCompa[i], sep = "\n")
    
    quines = lapply(1:length(data), function (j) which(is.element(groups[[j]], comparisons[,i])))
    data2 = lapply(1:length(data), function (j) data[[j]][,quines[[j]]]); names(data2) = names(data)
    groups2 = lapply(1:length(data), function (j) groups[[j]][quines[[j]]]); names(groups2) = names(groups)
    
    output[[i]] = MultiPower(data = data2, groups = groups2, type, omicPower, averagePower,
                             fdr, cost, equalSize, max.size, omicCol, powerPlots)
  }
  
  output$GlobalSummary = output[[1]]$summary
  
  output$GlobalSummary[,"minSampleSize"] = apply(sapply(output[-length(output)], function(x) x$summary[,"minSampleSize"]), 1,
                                                 function (y) {
                                                   if (length(unique(y)) == 1) return(unique(y))
                                                   if (length(unique(y)) != 1) return(paste(range(y, na.rm = T), collapse = "-"))})
  
  output$GlobalSummary[,"optSampleSize"] = apply(sapply(output[-length(output)], function(x) x$summary[,"optSampleSize"]),
                                                 1, max,  na.rm = TRUE)
  
  tmpSS = sapply(output[-length(output)], function(x) x$summary[,"optSampleSize"])
  tmpPower = sapply(output[-length(output)], function(x) x$summary[,"power"])
  output$GlobalSummary[,"power"] = sapply(1:nrow(tmpSS), function (i) tmpPower[i,which.max(tmpSS[i,])])
  
  cat("=============================== \n")
  cat("Global summary \n")
  cat("=============================== \n")
  print(output$GlobalSummary)
  
  
  if (summaryPlot) MultiCompaPlot(multiOutput = output, omicCol = omicCol, equalSize = equalSize, 
                                  legendLoc = "bottom")
  
  return(output)
  
}



# Plot for multiple comparisons -------------------------------------------

MultiCompaPlot = function(multiOutput, omicCol = NULL, equalSize, legendLoc = "bottomright") {

  if (is.null(omicCol)) {
    if (length(multiOutput[[1]]$parameters) > 12) {
      stop("Too many omics to be plotted. \n")
    }
    omicCol = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)]
    omicCol = omicCol[1:length(multiOutput[[1]]$parameters)]
  }

  omicShape = 1:length(multiOutput[[1]]$parameters)
  names(omicCol) = names(omicShape) = names(multiOutput[[1]]$parameters)


  ## Comparisons versus Sample Size
  if (equalSize) {
    apintar = sapply(multiOutput[-length(multiOutput)], function (x) x$summary[,"minSampleSize"])
    rownames(apintar) = names(multiOutput[[1]]$parameters)
    laopt = multiOutput$GlobalSummary[1,"optSampleSize"]
    optComp = sapply(multiOutput[-length(multiOutput)], function (x) x$summary[1,"optSampleSize"])

    matplot(1:ncol(apintar), t(apintar), type = "l", lwd = 2, xlab = "Comparisons", ylab = "Sample Size",
            main = "Sample Size per Comparison", col = omicCol, lty = omicShape, xaxt='n')
    axis(side = 1, at=1:ncol(apintar), labels=colnames(apintar))
    abline(h = laopt, lty = 1, lwd = 4)
    points(1:ncol(apintar), optComp, pch = 15, col = 1, cex = 1.3)
    legend(legendLoc, c(rownames(apintar), "  ", "Comparison Optimal SS", "Global Optimal SS"), lwd = 2, col = c(omicCol,"white", 1,1),
           lty = c(omicShape,1,NA,1), bty = "n", pch = c(rep(NA, nrow(apintar)+1), 15, NA), title = "Min SS per comparison", cex = 0.8)

  } else {
    apintar = sapply(multiOutput[-length(multiOutput)], function (x) x$summary[,"optSampleSize"])
    rownames(apintar) = names(multiOutput[[1]]$parameters)
    laOpt = multiOutput$GlobalSummary[,"optSampleSize"]
    apintar = data.frame(apintar, "Optimal" = laOpt, check.names = FALSE)

    matplot(1:ncol(apintar), t(apintar), type = "l", lwd = 2, xlab = "Comparisons", ylab = "Sample Size",
            main = "Sample Size per Comparison", col = omicCol, lty = omicShape, xaxt='n')
    axis(side = 1, at=1:ncol(apintar), labels=colnames(apintar))
    points(rep(ncol(apintar), nrow(apintar)), laOpt, pch = 15, col = omicCol, cex = 1.3)
    legend(legendLoc, rownames(apintar), lwd = 2, col = omicCol, cex = 0.8,
           lty = omicShape, bty = "n")

  }


  ## Comparisons vs Power
  apintar = sapply(multiOutput[-length(multiOutput)], function (x) x$summary[,"power"])
  rownames(apintar) = names(multiOutput[[1]]$parameters)
  optPow = multiOutput$GlobalSummary[,"power"]
  apintar = data.frame(apintar, "Optimal" = optPow, check.names = FALSE)

  matplot(1:ncol(apintar), t(apintar), type = "l", lwd = 2, xlab = "Comparisons", ylab = "Statistical Power",
          main = "Power per Comparison at Optimal SS", col = omicCol, lty = omicShape, xaxt='n')
  axis(side = 1, at=1:ncol(apintar), labels=colnames(apintar))
  points(rep(ncol(apintar), nrow(apintar)), optPow, pch = 15, col = omicCol, cex = 1.3)
  legend(legendLoc, rownames(apintar), lwd = 2, col = omicCol,
         lty = omicShape, bty = "n", cex = 0.8)

}





# Filtering data with Cohen's d cutoff ------------------------------------

CohenFilter = function (data, d, parameters) {
  
  if (length(d) == 1) d = rep(d, length(data))
  if (length(d) != length(data)) { 
    stop("Please, provide a single value for d or as many values a omic data types in data") 
  }
  
  dataF = data
  
  for (i in 1:length(data)) {
    quitar = which(parameters[[i]]$d < d[i])
    dataF[[i]] = dataF[[i]][-quitar,]
  }
  
  return(dataF)
  
}
