###########################################################
############     FoMPredictiveFunctions.R     #############
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)
# Last update: October/09/2019

# This script is to gather together the functions used for the MultiML predictive analysis 

###########################################################
# Functions:

permutations <- function(n){
  # Internal checkpoint function to create permutations to evaluate 
  # if all samples included are the same, and if they are in the 
  # same order so they can be compared across all Omics datasets
  
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

###########################################################

g_legend <- function(a.gplot){
  # Function to break ggplots legends and reconstruct them somewhere else
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
###########################################################

CombinationsFunction <- function (OmicsData, AnalysisType="complete", levels=NULL) {
  # Function to create the combinations of datasets to be studied
  # Input :
  # OmicsData: The data as list of predictor datasets including also the response (Y) matrix as the last element of the list
  # AnalysisType: The type of analysis, it can be 
  #               "exploratory", where all predictors vs Y along with each dataset of predictors vs Y are evaluated
  #               "complete", where all combinations are created or
  #               "detailed", where the user can choose by levels indicating if two datasets vs Y or three vs Y is required
  # levels: A vector indicating the levels to be analyzed if Analysis Type is detailed, if 2, then two predictors datasets will be contrasted against Y and so on 
  
  # Output:
  # A list of two lists:
  # Combinations: A list of the combinations required by the user to be created
  # Data: All required lists of datasets to be used as input for the error rate analysis
  
  vectito<-names(OmicsData)[1:length(OmicsData)-1]
  Combinatory<-list()
  combita<-combitachica<-NULL
  
  for (factores in 1:length (vectito)) {
    namecito<-paste(factores,"Omics", sep = " ")
    #print(namecito)
    Secuencia<-seq (1:factores)
    combita<-combn(x=vectito, m=factores)
    
    for (om in 1:ncol(combita)){
      namecitoChico<-paste(namecito,om, sep="")
      combitachica<-combita[,om]
      Combinatory[namecitoChico]<-list(combitachica)
      
    }
  }
  Combinations<-list()
  if (AnalysisType=="exploratory") {
    Combinations<-Combinatory[c(1:length(OmicsData)-1,(length(Combinatory)))]
    print(paste("Performing the analysis of",length(Combinations),"out of",length(Combinatory) ,"possible combinations", sep=" " ) )
  }
  if (AnalysisType=="complete") {
    Combinations<-Combinatory
    print(paste("Performing the analysis with all",length(Combinatory) ,"possible combinations", sep=" " ) )
  }
  if (AnalysisType=="detailed") {
    if (is.null(levels) ) {
      stop("If AnalysisType = detailed, you must indicate the levels")
    }
    SelList<-list()
    for (lev in 1:length(levels)) {
      levs=levels[lev]
      for (All in 1:length(Combinatory)) {
        if(length(Combinatory[[All]])==levs) {
          neiminho<-names(Combinatory[All])
          SelList<-Combinatory[[All]]
          Combinations[neiminho]<-list(SelList)
        }
      }
    }
    print(paste("Performing the analysis of",length(Combinations),"out of",length(Combinatory) ,"possible combinations", sep=" " ) )
  }
  ListofOmics<-OmicList<-list()
  for (branch in 1:length(Combinations)) {
    piece<-Combinations[[branch]]
    NameOmicsData<-paste0(piece, collapse = "_" )
    OmicsData2<-OmicsData[c(piece,names(OmicsData[length(OmicsData)]) )]
    OmicList<-list(OmicsData2)
    ListofOmics[NameOmicsData]<-OmicList
  }
  res<-list(Combinations=Combinations, Data=ListofOmics)
  return(res)
} 

###########################################################
SlurmFunction <- function (X) {
  # Function that will be included in the cluster for the required analysis
  # Input: 
  # X: A list of datasets created in "CombinationsFunction" 
  # comps: Number of componets to be calculated after each iteration through "ClassificationErrorRate" function
  # CrosVal: Type of cross validation to be applied in "ClassificationErrorRate" function, Leave-One-Out (LOOCV) or ten fold change (TenF) 
  # ticks: Number of segments (groups) of samples to evaluate.
  # iterations: Number of iterations in which error rate will be calculated
  
  ER_Calculator(Predictors=X, Response=length(X),Previous_CER=NULL,Ticks=NULL,WhichTicks=NULL,Function=Random.Forest.MP,Comps=10,crosval = "LOOCV",ERIterations=2,LassoIterations=2,TheoreticER=0.02)
}

###########################################################
Slurm_Creator<-function (Function, Parameters, jobname = NA, nodes = 2, cpus_per_node = 2, time=30,
                         pkgs = rev(.packages()) ) 
{
  # Function to create the files to be uploaded into a cluster to calculate the results of all particular combinations of predictor datasets.
  # Input: 
  # Function: An object class SlurmFunction that will calculate the error rates
  # Parameters: All required lists of datasets to be used as input for the error rate calculation
  # jobname: The desired name of the Slurm job. If NA, a random name of the form "slr####" will be assigned
  # nodes: The maximum number of nodes to be used in the cluster.
  # cpus_per_node: The number of CPUs per node in the cluster. It determines how many processes are run in parallel per node
  # pkgs: A character vector containing the names of packages that must be loaded on each cluster node. By default, it includes all packages loaded by the user when slurm_apply is called.
  # time: Time in days to run the data in the cluster
  
  # Output:
  # Four elements to be uploaded into cluster:
  # Function.RDS: A .RDS file with the Slurm Function that will be calculated
  # Parameters.RDS: A .RDS file with all lists of datasets that will be used for the error rate calculation
  # SlurmCreator_Run.R: A .R file with the parameters to run the job
  # SlurmCreator_submit_sh: A .sh file that indicates all parameters that the cluster requires to perform the job

  if (!is.function(Function)) {
    stop("first argument to slurm_apply should be a function")
  }
  
  if (!is.numeric(nodes) || length(nodes) != 1) {
    stop("nodes should be a single number")
  }
  if (!is.numeric(cpus_per_node) || length(cpus_per_node) != 
      1) {
    stop("cpus_per_node should be a single number")
  }
  tmpdir <- paste0("rslurm_", jobname)
  dir.create(tmpdir, showWarnings = FALSE)
  saveRDS(Parameters, file = file.path(tmpdir, "Parameters.RDS"))
  saveRDS(Function, file = file.path(tmpdir, "Function.RDS"))
  if (length(Parameters) < cpus_per_node * nodes) {
    nchunk <- cpus_per_node
  }
  else {
    nchunk <- ceiling(length(Parameters)/nodes)
  }
  nodes <- ceiling(length(Parameters)/nchunk)
  template_r <- readLines("templates/Slurm_run_R.txt")
  script_r <- whisker::whisker.render(template_r, list(pkgs = pkgs, 
                                                       nchunk = nchunk, cpus_per_node = cpus_per_node))
  writeLines(script_r, file.path(tmpdir, "SlurmCreator_Run.R"))
  template_sh <- readLines("templates/submit_sh.txt")
  
  rscript_path <- file.path(R.home("bin"), "Rscript")
  script_sh <- whisker::whisker.render(template_sh, list(max_node = nodes - 1, 
                                                         jobname = jobname,
                                                         time=time,
                                                         rscript = rscript_path))
  writeLines(script_sh, file.path(tmpdir, "SlurmCreator_submit_sh"))
  slurm_job(jobname, nodes)
}

###########################################################
ClassificationErrorRate<- function (Predictors, Response=length(Predictors),Function=Random.Forest.MP,Comps=10,crosval = "LOOCV",Ticks=10,WhichTicks=NULL,ERIterations=15,LassoIterations=15,ErrorRateType="ER",...) {
  # Function to evaluate through Error rate the predictive capability of the Multipower package
  # Input: 
  # Predictors: A list of different Omics Datasets and the Response matrix
  # Response: A number indicating the response matrix included in Predictors
  # Comps: Number of componets to be calculated after each iteration
  # crosval: Type of cross validation to be applied, Leave-One-Out (LOOCV) or ten fold change (TenF)
  # Ticks: Number of segments (groups) of samples to evaluate.
  # ERIterations: Number of iterations in which ER will be calculated
  # LassoIterations: Number of iterations of the Lasso selection per each Error rate analysis
  
  # Output:
  # Omics: A vector of the evaluated Omics
  # A list of two lists:
  # Minimums: A list of the minimum value of error rate, balanced (BER) or not (ER) obtained per each ten
  #           component analysis. This is measured through three distance metrics to evaluate the 
  #           classification performance of the model. Maximal distance (max.dist), distance to 
  #           centroids (centroids) or Mahalanobis distance (Mahalanobis)
  #           Thus, each table contains the results per each iteration at different subsets of samples
  
  # CompWinner: A list of the number of components in which the minimum value of error rate, 
  #             balanced (BER) or not (ER) was obtainde per each iteration. This is measured through 
  #             the three mentioned distance metrics to evaluate the classification performance 
  #             of the model. Thus, each table contains the components per each iteration at different 
  #             subsets of samples
  
  ##########
  
  components=Comps
  
  if (! is.list(Predictors) ) {
    stop("\\nOmics dataset must be a list with at least two elements")
  }
  ###########
  # Y
  Y<-t(as.matrix(Predictors[[Response]], drop=FALSE))
  if (is.character(Y[,1])) {
    Ycita<-transform(Y, Type = as.numeric(Type))
    Y<-Ycita
    rm(Ycita)
  }
  if (ncol(Y) != 1) {
    stop("\\nResponse must be a single variable")
  }
  if (any(is.na(Y))) {
    stop("\\nResponse must not contain missing values")
  }
  if (is.null(colnames(Y))) {
    colnames(Y) = "Y"
  }
  if (is.null(rownames(Y))) {
    rownames(Y) = 1:n
  }
  
  # Step1: Match the sample size
  LosIndivs<- rownames(Y)
  for (i in 1:length(Predictors)){
    LosIndivs = intersect(LosIndivs, colnames(Predictors[[i]]))
  }
  #print(paste("This analysis will be performed with",length(LosIndivs),"samples, since those are the ones repeated in all layers"))
  
  NewList<-Predictors
  
  for (i in 1:length(NewList)) {
    NewList[[i]]<-NewList[[i]][,colnames(NewList[[i]]) %in% LosIndivs]
  }
  
  ###########
  # Step2: Match the order
  LosColnames<-colnames(NewList[[1]])
  
  for (i in 1:length(NewList)) {
    NewList[[i]]<-NewList[[i]][,sort(LosColnames) ]
  }
  
  TestdeMatchTable<-unique(matrix(permutations(length(NewList)),ncol=2))
  
  for (i in 1:nrow(TestdeMatchTable)) {
    a=TestdeMatchTable[i,1]
    b=TestdeMatchTable[i,2]
    
    if (all(colnames(NewList[[a]])==colnames(NewList[[b]]))){
      #print(paste("Columns of lists",a,"and",b,"are equally sorted"))
    } else{
      #print("The colnames are not equally sorted")
    }
  }
  if (crosval=="LOOCV"){
    valid="loo"
  }
  if (crosval=="TenF"){
    valid="Mfold"
  }
  
  # Y
  Y<-t(as.matrix(NewList[[Response]], drop=FALSE))
  if (is.character(Y[,1])) {
    Ycita<-transform(Y, Type = as.numeric(Type))
    Y<-Ycita
    rm(Ycita)
  }
  
  #########
  Omics<-data.frame(Omics=names(NewList[-Response]))
  
  MinN<-round(length(table(as.factor(Y[,1]))) * 2, digits = 0)
  MaxN<-nrow(Y)
  Ngroups<-length(table(as.factor(Y[,1])))
  if(!is.null(Ticks)) {
    vectTicks<-round(seq(from=MinN,to = MaxN,length.out = Ticks),digits = 0)
  } else {
    if (is.null(WhichTicks)) {
      print("Please insert which ticks you want to calculate")
    } else {
    vectTicks<-WhichTicks
    }
  }
  muestra<-seq(1:length(vectTicks))
  
  muestrita<-muestritaLASSO<-NULL
  Muestra<-muestrasas<-list()
  MuestraLASSO<-muestrasasLASSO<-list()
  Minimums<-minimossas<-list()
  winnerSas<-CompWinner<-list()
  
  for (i in 1:length(vectTicks)) {
    
    #print(paste("Tick",i))
    print(paste("Processing your data... Please wait"))
    
    Mins<-ComponentWinner<-preMins<-NULL
    Minsbigrfc<-NULL
    TableMins<-NULL
    
    for (iter in 1:ERIterations){
      iter=iter
      #print(paste("performing iteration ",iter," of tick ",i,sep=""))
      #print(paste("Tick ",i,": ", vectTicks[i], " samples", sep = ""))
      namecito<-paste(vectTicks[i], " samples", sep = "")
      resample <- TRUE
      index <- rownames(Y)
      fun <- function(x) sample(x, round((vectTicks[i])/Ngroups,digits=0), replace = resample)
      a <- aggregate(index, by = list(group = Y[,1]), FUN = fun )
      a<-a[,-1]
      aLASSO<-aggregate(index, by = list(group = Y[,1]), FUN = fun )
      aLASSO<-aLASSO[,-1]
      Premuestrita<-as.vector(unlist(a))
      PremuestritaLASSO<-as.vector(unlist(aLASSO))
      if(vectTicks[i]>=length(Premuestrita)){
        muestrita<-c(Premuestrita,sample(rownames(Y), size=vectTicks[i]-length(Premuestrita), replace = TRUE))
        muestritaLASSO<-c(PremuestritaLASSO,sample(rownames(Y), size=vectTicks[i]-length(PremuestritaLASSO), replace = TRUE))
      } else{
        muestrita<-Premuestrita[1:(length(Premuestrita)-(length(Premuestrita)-vectTicks[i]))]
        muestritaLASSO<-PremuestritaLASSO[1:(length(PremuestritaLASSO)-(length(PremuestritaLASSO)-vectTicks[i]))]
      }
      muestrasas<-list(muestrita)
      muestrasasLASSO<-list(muestritaLASSO)
      Muestra[namecito]<-muestrasas
      MuestraLASSO[namecito]<-muestrasasLASSO
      #####
      Xchica2<-NewList
      Xchica2LASSO<-NewList
      
      # Matching the order again
      LasMuestras<-Muestra[[i]]
      LasMuestrasLASSO<-MuestraLASSO[[i]]
      # Sorting the tables
      for (long in 1:length(Xchica2)) {
        Xchica2[[long]]<-Xchica2[[long]][,sort(LasMuestras) ]
      }
      for (long2 in 1:length(Xchica2LASSO)) {
        Xchica2LASSO[[long2]]<-Xchica2LASSO[[long2]][,sort(LasMuestrasLASSO) ]
      }
      
      TestdeMatchTable<-unique(matrix(permutations(length(Xchica2)),ncol=2))
      for (lulu in 1:nrow(TestdeMatchTable)) {
        cola=TestdeMatchTable[lulu,1]
        colb=TestdeMatchTable[lulu,2]
        
        if (all(colnames(Xchica2[[cola]])==colnames(Xchica2[[colb]]))){
          #print(paste("Columns of lists",cola,"and",colb,"are equally sorted"))
        } else{
          print("Columns of lists",cola,"and",colb,"are NOT equally sorted")
        }
      }
      
      TestdeMatchTableLASSO<-unique(matrix(permutations(length(Xchica2LASSO)),ncol=2))
      for (lulu in 1:nrow(TestdeMatchTableLASSO)) {
        cola=TestdeMatchTableLASSO[lulu,1]
        colb=TestdeMatchTableLASSO[lulu,2]
        
        if (all(colnames(Xchica2LASSO[[cola]])==colnames(Xchica2LASSO[[colb]]))){
          #print(paste("Columns of lists",cola,"and",colb,"are equally sorted"))
        } else{
          print("Columns of lists",cola,"and",colb,"are NOT equally sorted")
        }
      }
      #
      
      for (ii in 1:length(Xchica2)) {
        Xchica2[[ii]]<-Xchica2[[ii]][,Muestra[[i]] %in% colnames(Xchica2[[ii]])]
      }
      for (ii in 1:length(Xchica2LASSO)) {
        Xchica2LASSO[[ii]]<-Xchica2LASSO[[ii]][,MuestraLASSO[[i]] %in% colnames(Xchica2LASSO[[ii]])]
      }
      
      ####
      Xchica<-Xchica2[-Response]
      Ychica<-as.matrix(Xchica2[[Response]])
      rm(Xchica2)
      if (is.character(Ychica[1,])) {
        Ycita<-as.numeric(as.factor(c(Ychica)))
        Ychica<-rbind(Ychica,as.numeric(Ycita))
      }
      for (lili in 1:length(Xchica)) {
        Xchica[[lili]]<-t(Xchica[[lili]])
      }
      
      XchicaLASSO<-Xchica2LASSO[-Response]
      YchicaLASSO<-as.matrix(Xchica2LASSO[[Response]])
      rm(Xchica2LASSO)
      if (is.character(YchicaLASSO[1,])) {
        Ycita2<-as.numeric(as.factor(c(YchicaLASSO)))
        YchicaLASSO<-rbind(YchicaLASSO,as.numeric(Ycita2))
      }
      for (lili in 1:length(XchicaLASSO)) {
        XchicaLASSO[[lili]]<-t(XchicaLASSO[[lili]])
      }
      
      PreRes<-operator (X=Xchica,Y=Ychica,XLasso=XchicaLASSO, YLasso=YchicaLASSO, LassoIterations=LassoIterations, FUN=Function,Comps = components)
      
      if (length(PreRes)>1){
        Mins<-as.matrix(cbind(Mins,PreRes$preMins))
        ComponentWinner<-cbind(ComponentWinner,PreRes$ComponentWinner2)
        class(ComponentWinner)<-"numeric"  
      } else {
        Mins<-cbind(Mins,PreRes)  
      }
    } # Closing Iteration
    nameSample<-paste(vectTicks[i], " samples", sep = "")
    if (dim(Mins)[1]==1){
      colnames(Mins)<-paste0("iter",seq(1:ERIterations));rownames(Mins)<-"OOB_ER"
      
    }
    minimossas<-list(Mins)
    Minimums[nameSample]<-minimossas
    winnerSas<-list(ComponentWinner)
    CompWinner[nameSample]<-winnerSas
  }
  ##### END of for (i in 1:length(vectTicks))
  ListofMins<-Minimums
  SDsas<-StdDev<-list()
  for(lista in 1:length(ListofMins)){
    Table<-ListofMins[[lista]]
    SDtita<-SDTota<-NULL
    for (row in 1:nrow(Table)){
      SDtita<-sd(Table[row,]) /sqrt(ncol(Table)) 
      SDTota<-rbind(SDTota,SDtita)
    }
    rownames(SDTota)<-rownames(ListofMins[[1]])
    nameList<-names(ListofMins[lista])
    SDsas<-list(SDTota)
    StdDev[nameList]<-SDsas
    
  }
  ###
  PreSDstable<-data.frame(matrix(data=NA,nrow=nrow(StdDev[[1]]),ncol = length(StdDev)))
  rownames(PreSDstable)<-rownames(StdDev[[1]])
  colnames(PreSDstable)<-t(data.frame(strsplit(names(StdDev),split=" "))[1,])[,1]
  
  for (iiiii in 1: length(StdDev)) {
    PreSDstable[,iiiii]<-rowMeans(StdDev[[iiiii]])
  }
  
  ###
  Premeanstable<-data.frame(matrix(data=NA,nrow=nrow(ListofMins[[1]]),ncol = length(ListofMins)))
  rownames(Premeanstable)<-rownames(ListofMins[[1]])
  colnames(Premeanstable)<-t(data.frame(strsplit(names(ListofMins),split=" "))[1,])[,1]
  
  for (i in 1: length(ListofMins)) {
    Premeanstable[,i]<-rowMeans(ListofMins[[i]])
  }
  
  ListofComps<-CompWinner
  if (is.null(ListofComps[[1]]) ){
    PreCompstable<-NULL
  } else {
    ToNumListofComps2<-PreCompstabletita<-PreCompstable<-NULL
    for (i in 1: length(ListofComps)) {
      ToNumListofComps2<-ListofComps[[i]]
      class(ToNumListofComps2) <- "numeric"
      for (row in 1:nrow(ToNumListofComps2)){
        PreCompstabletita<-rowMedians(ToNumListofComps2)
      }
      PreCompstable<-cbind(PreCompstable,PreCompstabletita)
    }
    rownames(PreCompstable)<-rownames(ListofComps[[1]])
    colnames(PreCompstable)<-t(data.frame(strsplit(names(ListofComps),split=" "))[1,])[,1]
    PreCompstable<-as.data.frame(PreCompstable)
  }
  
  
  MeansTable2<-Premeanstable
  SDstable2<-PreSDstable
  Compstable2<-PreCompstable
  
  MeansTable2$Category <- rownames(MeansTable2)
  SDstable2$Category <- rownames(SDstable2)
  Compstable2$Category <- rownames(Compstable2)
  meltedMeans<-melt(MeansTable2, id.vars="Category")
  meltedSDs<-melt(SDstable2, id.vars="Category")
  colnames(meltedSDs)[3]<-"StdDev"
  if (is.null(Compstable2) ){
    meltedComps<-NULL
    allmelted<-cbind(meltedMeans,StdDev=meltedSDs$StdDev)
    allmelted$ErrorRate<-t(data.frame(strsplit(allmelted$Category, split="_"))[2,])[,1]
    allmelted$Metric<-t(data.frame(strsplit(allmelted$Category, split="_"))[1,])[,1]  
  } else {
    meltedComps<-melt(Compstable2, id.vars="Category")
    colnames(meltedComps)[3]<-"Comp"
    allmelted<-cbind(meltedMeans,StdDev=meltedSDs$StdDev,Comp=meltedComps$Comp)
    allmelted$ErrorRate<-t(data.frame(strsplit(allmelted$Category, split="_"))[2,])[,1]
    allmelted$Metric<-t(data.frame(strsplit(allmelted$Category, split="_"))[1,])[,1]
  }
  
  if (!is.null(Ticks)) { 
    seqdeTicks<-as.numeric(as.character(unique(allmelted$variable)))
    MinN<-min(as.numeric(as.character(unique(allmelted$variable))))
    MaxN<-max(as.numeric(as.character(unique(allmelted$variable))))
    vectdeTicks<-round(seq(from=MinN,to = MaxN,length.out = Ticks),digits = 0)
    ElAusente<-vectdeTicks[!(vectdeTicks %in% seqdeTicks)]
    
    Elreemplazito<-Elreemplazo<-NULL
    for (ab in 1:length(ElAusente)){
      Elreemplazito<-seqdeTicks[which.min(abs(seqdeTicks - ElAusente[ab])) ]
      Elreemplazo<-c(Elreemplazo,Elreemplazito)
    }
    while (length(Elreemplazo[duplicated(Elreemplazo)])>0) { 
      Elnuevoreemplazo<-seqdeTicks[seqdeTicks>Elreemplazo[duplicated(Elreemplazo)] & seqdeTicks< Elreemplazo[which(duplicated(Elreemplazo))+1] ]
      Elreemplazo<-sort(c(unique(Elreemplazo),Elnuevoreemplazo ))  
    }
    
    vectdeTicks2<-vectdeTicks[(vectdeTicks %in% seqdeTicks)]
    NvectdeTicks<-sort(c(vectdeTicks2,Elreemplazo))
    allmelted<-allmelted[allmelted$variable %in% NvectdeTicks,]
  }
  
  if (!is.null(WhichTicks)) {
    NvectdeTicks<-sort(c(WhichTicks))
  }
    
  
  allmelted2<-allmelted
  ##################################################
  Metrica<-unique(allmelted2$Metric)
  errate<-unique(allmelted2$ErrorRate)
  allmelted3<-allmeltedsas<-list()
  CItita<-CItota<-CI<-NULL
  CL<-c(90,95,99)
  Z<-c(1.64,1.96, 2.58)
  CItita<-CItota<-NULL
  
  if (length(Metrica)==1){
    allmelted3<-allmelted2[allmelted2$Metric==Metrica,]
    for (zeta in 1:length(Z)) {
      namecito=paste("CI_",CL[zeta], sep="")
      for (al in 1:dim(allmelted3)[1]) {
        M=allmelted3$value[al]
        SM=sqrt( (allmelted3$StdDev[al])^2/as.numeric(as.character(allmelted3$variable[al]) ) )
        lower<-round(M-(Z[zeta]* SM),digits=3)
        upper<-round(M+(Z[zeta]* SM),digits=3)
        CItita<-cbind(namecito,paste("[",lower," , ",upper,"]", sep="") )
        CItota<-rbind(CItota,CItita)  
      }
      allmelted3[namecito]<-CItota[,2]
      CItita<-CItota<-NULL
    }
    
  } else {
    for(me in 1:length(Metrica) ) {
      tab<-allmelted2[allmelted2$Metric==Metrica[me],]
      
      for(eler in 1:length(errate)) {
        namecito<-paste(Metrica[me],errate[eler],sep="_")
        tab2<-tab[tab$ErrorRate==errate[eler],]
        allmeltedsas<-list(tab2)
        allmelted3[namecito]<-allmeltedsas
      }
    }
    for (uu in 1:length(allmelted3)){
      for (zeta in 1:length(Z)) {
        namecito=paste("CI_",CL[zeta], sep="")
        for (al in 1:dim(allmelted3[[uu]])[1]) {
          M=allmelted3[[uu]]$value[al]
          SM=sqrt( (allmelted3[[uu]]$StdDev[al])^2/as.numeric(as.character(allmelted3[[uu]]$variable[al]) ) )
          lower<-round(M-(Z[zeta]* SM),digits=3)
          upper<-round(M+(Z[zeta]* SM),digits=3)
          CItita<-cbind(namecito,paste("[",lower," , ",upper,"]", sep="") )
          CItota<-rbind(CItota,CItita)  
        }
        allmelted3[[uu]][namecito]<-CItota[,2]
        CItita<-CItota<-NULL
      }
    }
  }
  
  allmelted3
  #######
  
  if (is.null(winnerSas[[1]]) ) {
    result<-list(TestedTicks=NvectdeTicks,Omics=Omics,Minimums=Minimums,
                 TablebyTick=allmelted3)
  } else {
  result<-list(TestedTicks=NvectdeTicks,Omics=Omics,Minimums=Minimums,CompWinner=CompWinner ,
               TablebyTick=allmelted3)
  }
  return(result)
}


###########################################################
operator<-function (X,Y, XLasso,YLasso,FUN,LassoIterations,...) {
  args<-list(...)
  Comps<-args$Comps
  iter<-args$iter
  FUN (X,Y,XLasso,YLasso,LassoIterations,...)
}

###########################################################
PLSDA.MP<-function (X,Y,XLasso, YLasso, LassoIterations,...) {
  args<-list(...)
  Comps<-args$Comps
  iter<-args$iter
  
    if (length(X) >1){
    summary(X)
    ### LassoSelection ###
    SelVarstita<-SelVars<-NULL
    for (lst in 1:length(X)){
      X[[lst]]<-scale(X[[lst]],center=TRUE, scale = FALSE)
      SelVars<-LassoSelection (X=XLasso[[lst]],Y=as.factor(YLasso[2,]), IterationsforVarSelections=LassoIterations)
      X[[lst]]<-X[[lst]][,colnames(X[[lst]]) %in% SelVars$SelectedVariables]
    }
    ######################
    res2<- block.plsda(X,as.vector(Y[2,]), ncomp=Comps)  
    perf.plsda <- try (perf(res2, validation = valid, folds = 5, 
                            progressBar = FALSE, auc = TRUE, nrepeat = 1) ,silent=TRUE) 
    ##
    #TableMeans<-table<-NULL
    if (class(perf.plsda) == "try-error") {
      preMins<-matrix(1,nrow=6,ncol=1)
      colnames(preMins)<-iter;rownames(preMins)<-c("max.dist_ER","max.dist_BER","centroids.dist_ER","centroids.dist_BER","mahalanobis.dist_ER","mahalanobis.dist_BER")
      ComponentWinner2<-matrix(1,nrow=6,ncol=1)
      colnames(ComponentWinner2)<-iter;rownames(ComponentWinner2)<-c("max.dist_ER","max.dist_BER","centroids.dist_ER","centroids.dist_BER","mahalanobis.dist_ER","mahalanobis.dist_BER")
    } else{
      preMeans<-perf.plsda$WeightedVote.error.rate
      
      TableMeans<-table<-NULL
      for (p in 1: length (preMeans)){
        namecito2<-names(preMeans[p])
        table<-preMeans[[p]][grep(pattern = "Overall",x = rownames(preMeans[[p]])),]
        tipodeER<-gsub("Overall.", "", rownames(table))
        rownames(table)<-paste(namecito2,tipodeER, sep="_")
        TableMeans<-rbind(TableMeans,table)
      }
      
      preMins<-as.matrix(apply(TableMeans, 1, FUN=min))
      preComponentWinner<-data.frame(ComponentWinner=colnames(TableMeans)[apply(TableMeans,1,which.min)])
      ComponentWinner2<-t(data.frame(strsplit(as.vector(preComponentWinner[,1]),split=" "))[2,])
      rownames(ComponentWinner2)<-rownames(preMins);colnames(ComponentWinner2)<-iter
    }
    Mins<-as.matrix(cbind(Mins,preMins))
    ComponentWinner<-cbind(ComponentWinner,ComponentWinner2)
    ComponentWinner2<-ComponentWinner
    class(ComponentWinner2)<-"numeric" 
    
    ################################################
  } else {
    ### LassoSelection ###
    SelVars<-LassoSelection (X=XLasso[[1]],Y=as.factor(YLasso[2,]), IterationsforVarSelections=LassoIterations)
    X[[1]]<-X[[1]][,colnames(X[[1]]) %in% SelVars$SelectedVariables]
    ######################
    CompsNEW<-min( (nrow(X[[1]])-1), Comps )
    res2<- plsda(X[[1]],as.factor(Y[2,]), ncomp=CompsNEW)  
    perf.plsda <- try (perf(res2, validation = valid, folds = 5, 
                            progressBar = FALSE, auc = TRUE, nrepeat = 1) ,silent=TRUE) 
    
    if (class(perf.plsda)[[1]] == "try-error") {
      preMins<-matrix(1,nrow=6,ncol=1)
      colnames(preMins)<-iter;rownames(preMins)<-c("max.dist_ER","max.dist_BER","centroids.dist_ER","centroids.dist_BER","mahalanobis.dist_ER","mahalanobis.dist_BER")
      ComponentWinner2<-matrix(1,nrow=6,ncol=1)
      colnames(ComponentWinner2)<-iter;rownames(ComponentWinner2)<-c("max.dist_ER","max.dist_BER","centroids.dist_ER","centroids.dist_BER","mahalanobis.dist_ER","mahalanobis.dist_BER")
    } else {
      preMeans<-perf.plsda$error.rate
      
      TableMeans<-table<-NULL
      for (p in 1: length (preMeans)){
        namecito2<-names(preMeans[p])
        if (namecito2=="overall"){
          namecito2<-"ER"
        }
        table<-t(preMeans[[p]])
        rownames(table)<-paste(rownames(table),namecito2, sep="_")
        TableMeans<-rbind(TableMeans,table)
      }
      preMins<-as.matrix(apply(TableMeans, 1, FUN=min))
      preComponentWinner<-data.frame(ComponentWinner=colnames(TableMeans)[apply(TableMeans,1,which.min)])
      ComponentWinner2<-t(data.frame(strsplit(as.vector(preComponentWinner[,1]),split=" "))[2,])
      rownames(ComponentWinner2)<-rownames(preMins);colnames(ComponentWinner2)<-iter
    }
    
    
  } # Closing separation between block.plsda and plsda
  
  res=list(preMins=preMins,ComponentWinner2=ComponentWinner2 )
  return (res)
} 

###########################################################
Random.Forest.MP<-function (X,Y,XLasso,YLasso, LassoIterations,...) {
  args<-list(...)
  Comps<-args$Comps
  iter<-args$iter
  Xchicaflat2Lasso<-do.call(cbind.data.frame, XLasso)
  Xchicaflat2Lasso<-scale(Xchicaflat2Lasso,center=TRUE, scale = FALSE)
  names(Xchicaflat2Lasso) <- make.names(names(Xchicaflat2Lasso))
  YYLasso<-as.data.frame(YLasso[1,]);colnames(YYLasso)<-"Y"
  XchicaflatLasso<-merge(YYLasso,Xchicaflat2Lasso,by=0 );rownames(XchicaflatLasso)<-XchicaflatLasso[,1];XchicaflatLasso<-XchicaflatLasso[,-1]
  names(XchicaflatLasso) <- make.names(names(XchicaflatLasso))
  rm(Xchicaflat2Lasso,YYLasso)
  
  Xchicaflat2<-do.call(cbind.data.frame, X); dim(Xchicaflat2)
  Xchicaflat2<-scale(Xchicaflat2,center=TRUE, scale = FALSE)
  names(Xchicaflat2) <- make.names(names(Xchicaflat2))
  YY<-as.data.frame(Y[1,]);colnames(YY)<-"Y"
  Xchicaflat<-merge(YY,Xchicaflat2,by=0 );rownames(Xchicaflat)<-Xchicaflat[,1];Xchicaflat<-Xchicaflat[,-1]
  names(Xchicaflat) <- make.names(names(Xchicaflat))
  rm(Xchicaflat2,YY)
  
  ### LassoSelection ###
  SelVars<-LassoSelection (X=XchicaflatLasso[,-1],Y=XchicaflatLasso[,1], IterationsforVarSelections=LassoIterations)
  SelVars$SelectedVariables<-c("Y",SelVars$SelectedVariables)
  Xchicaflat<-Xchicaflat[,colnames(Xchicaflat) %in% SelVars$SelectedVariables]#; dim(Xchicaflat)
  set.seed(1)
  model <- randomForest(Y ~., data=Xchicaflat)
  
  preMins<-model$err.rate[500,1]
  rm(model)
  return(preMins)
}  
###########################################################

LassoSelection<- function (X,Y, IterationsforVarSelections=15, ...){
  # Input:
  # X: This is the table of predictor variables where n is individuals and p is variables
  # Y: This is the vector of response variables
  # IterationsforVarSelections: Number of iterations of the variable selection step
  ##################
  require ("glmnet")
  
  X<-as.matrix(X)
  alpha_val=1
  SelVars<- SelectedVariables<-NULL
  SelecVariablesalliterations<-NULL 
  CoefficientsVarstita<-CoefficientsVars<-NULL
  for (i in 1:IterationsforVarSelections) {
    #print (paste ("LASSO Variable Selection iteration ",i, sep=" "))
    
    repeat {
      cfit<-try(cv.glmnet(as.matrix(X),as.vector(Y), 
                          standardize=TRUE, family="multinomial", 
                          alpha = alpha_val, grouped = FALSE,
                          type.measure = "mae"),silent=TRUE ) 
      
      if (class(cfit) == "try-error") {
        X=rbind(X,X)
        Y=c(Y,Y)
      } else { 
        
        for (ls in 1:length(coef(cfit, s = "lambda.min"))) {
          CoefficientsVarstita<-as.matrix(coef(cfit, s = "lambda.min")[[ls]])
          CoefficientsVarstita<-CoefficientsVarstita[-1,,drop=FALSE]
          CoefficientsVars<-rbind(CoefficientsVars,CoefficientsVarstita)
        }
        SelVarsPositions<-which(CoefficientsVars != 0)
        SelVarssmall<-CoefficientsVars[which(CoefficientsVars != 0),,drop=FALSE]
        SelVarssmalllista<- unique(rownames(SelVarssmall))
        SelVars<-c(SelVars,SelVarssmalllista)
        SelVars<-unique(SelVars)
        break 
      }
    }
    SelecVariablesalliterations<-c(SelecVariablesalliterations,SelVars)  
    SelectedVariables<-unique(c(SelecVariablesalliterations))
  }
  res<-list(SelectedVariables=SelectedVariables)
  return(res)
}
###########################################################
ErrorRateplot<-function (x, ErrorRateType="BER",MetricsType="max.dist",DoF=NULL,Projection=FALSE,Spline=TRUE, TheoreticER=NULL,ConfInt=0.95) {
  # Input:
  # x: list of tables with the error rates results and the number of components per tick 
  # ErrorRateType: Character to indicate which error rate to plot. It can be error rate "ER", 
  #                balanced error rate "BER" or "Both"
  # MetricsType: Character to indicate the metrics to be plotted. It can be Maximal 
  #              distance "max.dist", distance to centroids "centroids.dist", Mahalanobis 
  #              distance "mahalanobis.dist" or "All"
  # DoF: It can be either NULL to indicate that the degrees of freedom of the Spline model are ticks-1, 
  #             or a value so the model is created with a different degree of freedom. The chosen DoF affects the samples required so the user must be careful.
  # Projection: Character to indicate if the user needs for the Spline model to be projected until 
  #             the minimum error rate Min_ER (TRUE) or not (FALSE).
  # Spline: Character to indicate if the user needs to plot the Spline model (TRUE) or not (FALSE).
  # TheoreticER: Character to indicate the minimum value of ER in order to calculate the adequate number of samples.
  # ConfInt: Character to indicate the confidence interval for the calculation of margin of error and plot of the confidence area of the Spline model. The user can use 0.90, 0.95 or 0.99.
  
  # Output:
  # A Linechart with standard deviation. If the plot represent just one error type and one metrics
  # The number of components will appear indicating that was the best number of component.
  suppressWarnings({ 
  Layers<-x$Omics
  rownames(Layers)<-NULL
  ListofMins<-x$Minimums
  ticks<-DoF
  
  if(is.null(ticks)){
    DegreesOfFreedom<-length(ListofMins)-1
  } else {
    DegreesOfFreedom<-ticks-1
  }
  if (MetricsType=="All" || ErrorRateType=="Both"){
    allmelted3<-do.call(rbind.data.frame, x$TablebyTick)
    if (MetricsType!="All" ){
      allmelted3<-allmelted3[allmelted3$Metric==MetricsType,]
    }
    if (ErrorRateType!="Both" ){
      allmelted3<-allmelted3[allmelted3$ErrorRate==ErrorRateType,]
    }
  } else{
    selector<-paste(MetricsType,ErrorRateType,sep="_")  
  }
  
  if (dim(ListofMins[[1]])[1]==1){
    allmelted3<-x$TablebyTick
  
  } else{
    allmelted3<-x$TablebyTick[[selector]]
  }
  
  if (ConfInt==0.90){
    df<-allmelted3[,c(7:9)]
    df<-lapply (df, function (x) gsub("\\[|\\]","",x))
    df2<-as.data.frame(strsplit(df$CI_90, split=","))
    lowers<-as.numeric(as.character(unname(unlist(df2[1,]))))
    uppers<-as.numeric(as.character(unname(unlist(df2[2,]))))
  } 
  if (ConfInt==0.95) {
    df<-allmelted3[,c(7:9)]
    df<-lapply (df, function (x) gsub("\\[|\\]","",x))
    df2<-as.data.frame(strsplit(df$CI_95, split=","))
    lowers<-as.numeric(as.character(unname(unlist(df2[1,]))))
    uppers<-as.numeric(as.character(unname(unlist(df2[2,]))))
  }
  if(ConfInt==0.99){
    df<-allmelted3[,c(7:9)]
    df<-lapply (df, function (x) gsub("\\[|\\]","",x))
    df2<-as.data.frame(strsplit(df$CI_99, split=","))
    lowers<-as.numeric(as.character(unname(unlist(df2[1,]))))
    uppers<-as.numeric(as.character(unname(unlist(df2[2,]))))
  }
  
  allmelted3$lower<-lowers
  allmelted3$upper<-uppers
  
  allmeltedmodel<-allmelted3
  colnames(allmeltedmodel)[colnames(allmeltedmodel)=="SEM"] <- "StdDev"
  colnames(allmeltedmodel)[colnames(allmeltedmodel)=="Samples"] <- "variable"
  colnames(allmeltedmodel)[colnames(allmeltedmodel)=="ER_value"] <- "value"
  colnames(allmeltedmodel)[colnames(allmeltedmodel)=="ER_type"] <- "ErrorRate"
  colnames(allmelted3)[colnames(allmelted3)=="SEM"] <- "StdDev"
  colnames(allmelted3)[colnames(allmelted3)=="Samples"] <- "variable"
  colnames(allmelted3)[colnames(allmelted3)=="ER_value"] <- "value"
  colnames(allmelted3)[colnames(allmelted3)=="ER_type"] <- "ErrorRate"
  
  ro=dim(allmeltedmodel)[1]
  
  if (allmeltedmodel$value[ro-1]<allmeltedmodel$value[ro]) {
        allmeltedmodel$value[ro]<-allmeltedmodel$value[ro-1]-abs(allmeltedmodel$StdDev[ro-1])
  } else{}
  
  
  #
  if (is.null(TheoreticER) ) {
    MinimumError=min(allmelted3$value)
    SamplesEvaluated=SamplesRequired=as.numeric(as.character(allmelted3$variable[allmelted3$value==min(allmelted3$value)]))
  } else{
    MinimumError=TheoreticER 
    
  }
  
  if (length(unique(allmelted3$Category))>1) {
    p<-ggplot(allmelted3, aes(x=variable, y= value, group=Category)) +
          geom_line(aes(linetype=ErrorRate, color=Metric), size=1) + 
          ylim(NA,1.1) +
          ylab("Classification Error Rate") +
          xlab("Number of Samples") +
          theme(axis.title=element_text(face="bold",size="14"),
                axis.text.x = element_text(face="bold", size=18),
                axis.text.y = element_text(face="bold", size=18),
                plot.title = element_text(size = "16", face = "bold")
          ) +
          ggtitle(paste(length(unique(allmelted2$variable))-1, "segments", sep= " ")) +
          scale_linetype_manual(values=c("solid", "twodash"))+
          geom_errorbar(aes(ymin = value - StdDev,
                            ymax = value + StdDev, color=Metric))
        legend <- g_legend(p)
        grid.newpage()
        vp1 <- viewport(width = 0.75, height = 1, x = 0.375, y = .5)
        vpleg <- viewport(width = 0.25, height = 0.5, x = 0.85, y = 0.75)
        subvp <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.25)
        print(p + theme(legend.position = "none"), vp = vp1)
        upViewport(0)
        pushViewport(vpleg)
        grid.draw(legend)
        upViewport(0)
        pushViewport(subvp)
        my_table <- tableGrob(Layers) 
        grid.draw(my_table)
  } else {
    if (Projection==TRUE & Spline==FALSE){
          Spline=TRUE
          print("Since Projection==TRUE, the Spline will be plotted")
        }
    if (Projection==FALSE & Spline==FALSE){
      p<-ggplot() + 
            geom_line(aes(x=as.numeric(as.character(variable)), y=value, group=Category,
                          linetype=ErrorRate, color=Metric), allmelted3, size=1) +
            geom_errorbar(aes(x=as.numeric(as.character(variable)), group=Category,
                              ymin = value - StdDev, ymax = value + StdDev, color=Metric), allmelted3) +
            ylim(NA,1.1) +
            ylab("Classification Error Rate") +
            xlab("Number of Samples") +
            theme(axis.title=element_text(face="bold",size="14"),
                  axis.text.x = element_text(face="bold", size=18),
                  axis.text.y = element_text(face="bold", size=18),
                  plot.title = element_text(size = "16", face = "bold")
            ) +
            ggtitle(paste(length(unique(allmelted3$variable))-1, "segments", sep= " ")) +
            scale_linetype_manual(values=c("solid", "twodash"))
          
          legend <- g_legend(p)
          grid.newpage()
          vp1 <- viewport(width = 0.75, height = 1, x = 0.375, y = .5)
          vpleg <- viewport(width = 0.25, height = 0.5, x = 0.85, y = 0.75)
          subvp <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.25)
          print(p + theme(legend.position = "none"), vp = vp1)
          upViewport(0)
          pushViewport(vpleg)
          grid.draw(legend)
          upViewport(0)
          pushViewport(subvp)
          my_table <- tableGrob(Layers) 
          grid.draw(my_table)
        }
    if (Projection==FALSE & Spline==TRUE){
      MinimumError=min(allmelted3$value)
      fm1 <- lm(value ~ bs(as.numeric(as.character(variable)), degree = 1,df = DegreesOfFreedom), data = allmeltedmodel)
      DESVESTprediction<-lm(StdDev ~ bs(as.numeric(as.character(variable)), degree = 1,df = DegreesOfFreedom), data = allmeltedmodel)
      
      ht01 <- seq(min(as.numeric(as.character(allmeltedmodel$variable))),2000, length.out = 2000)
      prediction<-as.data.frame(cbind(predict(fm1, data.frame(variable = ht01)),predict(DESVESTprediction, data.frame(variable = ht01))   )  );colnames(prediction)<-c("Prediction","STDpred")
      ht01Table<-cbind(ht01,as.data.frame(prediction));colnames(ht01Table)<-c("ht","ErrorPred", "STDpred")
      
      ht <- seq(min(as.numeric(as.character(allmeltedmodel$variable))),max(as.numeric(as.character(allmeltedmodel$variable))), length.out = 200)
      htTable<-cbind(ht,as.data.frame(predict(fm1, data.frame(variable = ht))));colnames(htTable)<-c("ht","ErrorPred")
      
      NSampleMaxTable<-as.data.frame(ht01Table[ht01Table[,2]>=MinimumError,])
      NSampleMax<-1+round(NSampleMaxTable$ht[dim(NSampleMaxTable)[1]])
      NSampleMaxPosition<-as.numeric(rownames(NSampleMaxTable)[dim(NSampleMaxTable)[1]])
      
      NSampleMaxTable$EplusStd<-NSampleMaxTable$ErrorPred - NSampleMaxTable$STDpred
      Theplusminustable<-as.data.frame(NSampleMaxTable[NSampleMaxTable[,4]>=MinimumError,])
      Theplusminus<-round(Theplusminustable$ht[dim(Theplusminustable)[1]])
      TheplusminusPosition<-as.numeric(rownames(Theplusminustable)[dim(Theplusminustable)[1]])
      
      SamplesEvaluated<-as.numeric(as.character(allmelted3$variable[dim(allmelted3)[1]]))
      SamplesRequired<-round(NSampleMaxTable$ht[NSampleMaxPosition])
      MOE<-paste("±",round(NSampleMax-Theplusminus))
      
      p<-ggplot() + 
        geom_line(aes(x=as.numeric(as.character(variable)), y=value, group=Category,
                      linetype=ErrorRate, color=Metric), allmelted3, size=1) +
        geom_errorbar(aes(x=as.numeric(as.character(variable)), group=Category,
                          ymin = lower, ymax = upper, color=Metric), allmelted3) +
        geom_smooth(method="loess", span=0.2, aes(x=ht, y=ErrorPred,color="SPline"), data=NSampleMaxTable, size=1.3 ) +
        ylim(NA,1.1) +
        ylab("Classification Error Rate") +
        xlab("Number of Samples") +
        theme(axis.title=element_text(face="bold",size="14"),
              axis.text.x = element_text(face="bold", size=18),
              axis.text.y = element_text(face="bold", size=18),
              plot.title = element_text(size = "16", face = "bold")
        ) +
        ggtitle(paste(length(unique(allmelted3$variable))-1, "segments", sep= " ")) +
        scale_linetype_manual(values=c("solid", "twodash"))
      
      legend <- g_legend(p)
      grid.newpage()
      vp1 <- viewport(width = 0.75, height = 1, x = 0.375, y = .5)
      vpleg <- viewport(width = 0.25, height = 0.5, x = 0.85, y = 0.75)
      subvp <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.5)
      subvp2 <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.25)
      print(p + theme(legend.position = "none"), vp = vp1)
      upViewport(0)
      pushViewport(vpleg)
      grid.draw(legend)
      upViewport(0)
      pushViewport(subvp)
      rownames(Layers)<-NULL
      my_table <- tableGrob(Layers) 
      grid.draw(my_table)
      
      
      TERtable<-rbind(ERtarget=round(MinimumError,digits = 3),
                      PSS=SamplesRequired,
                      MOE= MOE) 
      TERtableGrob <- tableGrob(TERtable) 
      upViewport(0)
      pushViewport(subvp2)
      grid.draw(TERtableGrob)
    }
    
    if (Projection==TRUE & Spline==TRUE){
      fm1 <- lm(value ~ bs(as.numeric(as.character(variable)), degree = 1,df = DegreesOfFreedom), data = allmeltedmodel)
      DESVESTprediction<-lm(StdDev ~ bs(as.numeric(as.character(variable)), degree = 1,df = DegreesOfFreedom), data = allmeltedmodel)
      
      ht01 <- seq(min(as.numeric(as.character(allmeltedmodel$variable))),2000, length.out = 2000)
      prediction<-as.data.frame(cbind(predict(fm1, data.frame(variable = ht01)),predict(DESVESTprediction, data.frame(variable = ht01))   )  );colnames(prediction)<-c("Prediction","STDpred")
      ht01Table<-cbind(ht01,as.data.frame(prediction));colnames(ht01Table)<-c("ht","ErrorPred", "STDpred")
      
      ht <- seq(min(as.numeric(as.character(allmeltedmodel$variable))),max(as.numeric(as.character(allmeltedmodel$variable))), length.out = 200)
      htTable<-cbind(ht,as.data.frame(predict(fm1, data.frame(variable = ht))));colnames(htTable)<-c("ht","ErrorPred")
      
      NSampleMaxTable<-as.data.frame(ht01Table[ht01Table[,2]>=MinimumError,])
      rownames(NSampleMaxTable)<-seq(1:dim(NSampleMaxTable)[1])
      NSampleMax<-1+round(NSampleMaxTable$ht[dim(NSampleMaxTable)[1]])
      NSampleMaxPosition<-as.numeric(rownames(NSampleMaxTable)[dim(NSampleMaxTable)[1]])

      NSampleMaxTable$EplusStd<-NSampleMaxTable$ErrorPred - NSampleMaxTable$STDpred
      Theplusminustable<-as.data.frame(NSampleMaxTable[NSampleMaxTable[,4]>=MinimumError,])
      Theplusminus<-round(Theplusminustable$ht[dim(Theplusminustable)[1]])
      TheplusminusPosition<-as.numeric(rownames(Theplusminustable)[dim(Theplusminustable)[1]])

      SamplesEvaluated<-as.numeric(as.character(allmelted3$variable[dim(allmelted3)[1]]))
      SamplesRequired<-round(NSampleMaxTable$ht[NSampleMaxPosition])
      MOE<-paste("±",round(NSampleMax-Theplusminus))
  
        p<-ggplot() + 
          geom_line(aes(x=as.numeric(as.character(variable)), y=value, group=Category,
                        linetype=ErrorRate, color=Metric), allmelted3, size=1) +
          geom_errorbar(aes(x=as.numeric(as.character(variable)), group=Category,
                            ymin = lower, ymax = upper, color=Metric), allmelted3) +
          
          
          geom_smooth(method="loess", span=0.2, aes(x=ht, y=ErrorPred,color="SPline"), data=NSampleMaxTable, size=1.3 ) +
          ylim(NA,1.1) +
          ylab("Classification Error Rate") +
          xlab("Number of Samples") +
          theme(axis.title=element_text(face="bold",size="14"),
                axis.text.x = element_text(face="bold", size=18),
                axis.text.y = element_text(face="bold", size=18),
                plot.title = element_text(size = "16", face = "bold")
          ) +
          ggtitle(paste(length(unique(allmelted3$variable))-1, "segments", sep= " ")) +
          scale_linetype_manual(values=c("solid", "twodash"))
        
        legend <- g_legend(p)
        grid.newpage()
        vp1 <- viewport(width = 0.75, height = 1, x = 0.375, y = .5)
        vpleg <- viewport(width = 0.25, height = 0.5, x = 0.85, y = 0.75)
        subvp <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.5)
        subvp2 <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.25)
        print(p + theme(legend.position = "none"), vp = vp1)
        upViewport(0)
        pushViewport(vpleg)
        grid.draw(legend)
        upViewport(0)
        pushViewport(subvp)
        rownames(Layers)<-NULL
        my_table <- tableGrob(Layers) 
        grid.draw(my_table)
        
        TERtable<-rbind(ERtarget=round(MinimumError,digits = 3),
                        PSS=SamplesRequired,
                        MOE= MOE) 
        TERtableGrob <- tableGrob(TERtable) 
        upViewport(0)
        pushViewport(subvp2)
        grid.draw(TERtableGrob)
        
    }
    res<-TERtable
    return(res)
  }
  })
}

###########################################################
Comparative_ERPlot<-function (L, ErrorRateType = "ER", MetricsType ="mahalanobis.dist"){
  # Input:
  # L: list of tables obtained from "UntilAllmelted" analysis 
  
  # Output:
  # A plot of number of samples vs Classification Error Rate in order to compare and evaluate different analyzes
  
  LosRowNames<-Ltransp1<-Ltransp2<-Ltransp3<-Ltransp4<-NULL
  Ltransp5<-Ltransp<-list()
  
  if (length(unique(L[[1]]$TablebyTick$Category))==1) {
    LosRowNames<-Ltransp1<-Ltransp2<-Ltransp3<-Ltransp4<-NULL
    Ltransp5<-Ltransp<-list()
    for (ls in 1:length(L)){
      Namecito<-names(L[ls])
      Ltransp1<-data.frame(L[[ls]]$TablebyTick)
      Ltransp1$Omics<-rep(Namecito, dim(Ltransp1)[1])
      Ltransp[[Namecito]]<-Ltransp1
    }
  
  TableAll <- do.call(rbind,Ltransp)
  rownames(TableAll)<-seq(1:dim(TableAll)[1])
  
    
  Metrica<-t(data.frame(cbind(unique(TableAll$Metric),unique(TableAll$ErrorRate))))
  rownames(Metrica)<-NULL
  p<-ggplot() + 
      geom_line(aes(x=as.numeric(as.character(variable)), y=value, group=Omics,colour = Omics), data = TableAll, size=1) +
      geom_errorbar(aes(x=as.numeric(as.character(variable)),group=Omics ,
                        ymin = value - StdDev, ymax = value + StdDev, color=Omics), TableAll) +
      ylab("Classification Error Rate") +
      xlab("Number of Samples") +
      theme(axis.title=element_text(face="bold",size="14"),
            axis.text.x = element_text(face="bold", size=18),
            axis.text.y = element_text(face="bold", size=18),
            plot.title = element_text(size = "16", face = "bold")
      ) +
      ggtitle(paste("Comparative Classification Plot "))
    legend <- g_legend(p)
    grid.newpage()
    vp1 <- viewport(width = 0.75, height = 1, x = 0.375, y = .5)
    vpleg <- viewport(width = 0.25, height = 0.5, x = 0.85, y = 0.75)
    subvp <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.25)
    print(p + theme(legend.position = "none"), vp = vp1)
    upViewport(0)
    pushViewport(vpleg)
    grid.draw(legend)
    upViewport(0)
    pushViewport(subvp)
    my_table <- tableGrob(Metrica) 
    grid.draw(my_table)  
  } else {
    LosRowNames<-Ltransp1<-Ltransp2<-Ltransp3<-Ltransp4<-NULL
    Ltransp5<-Ltransp<-list()
    selector<-paste(MetricsType,ErrorRateType, sep="_")
    loscolnames<-colnames(L[[1]]$TablebyTick[[1]])
    for (ls in 1:length(L)){
      Namecito<-names(L[ls])
      Ltransp1<-data.frame(L[[ls]]$TablebyTick [ names(L[[ls]]$TablebyTick)==selector ])
      colnames(Ltransp1)<-loscolnames
      Ltransp1$Omics<-rep(Namecito, dim(Ltransp1)[1])
      Ltransp[[Namecito]]<-Ltransp1
    }
    
    TableAll <- do.call(rbind,Ltransp)
    rownames(TableAll)<-seq(1:dim(TableAll)[1])
    
    Metrica<-t(data.frame(cbind(unique(TableAll$Metric),unique(TableAll$ErrorRate))))
    rownames(Metrica)<-NULL
    p<-ggplot() + 
      geom_line(aes(x=as.numeric(as.character(variable)), y=value, group=Omics,colour = Omics), data = TableAll, size=1) +
      geom_errorbar(aes(x=as.numeric(as.character(variable)),group=Omics ,
                        ymin = value - StdDev, ymax = value + StdDev, color=Omics), TableAll) +
      ylab("Classification Error Rate") +
      xlab("Number of Samples") +
      theme(axis.title=element_text(face="bold",size="14"),
            axis.text.x = element_text(face="bold", size=18),
            axis.text.y = element_text(face="bold", size=18),
            plot.title = element_text(size = "16", face = "bold")
      ) +
      ggtitle(paste("Comparative Classification Plot "))
    
    legend <- g_legend(p)
    grid.newpage()
    vp1 <- viewport(width = 0.75, height = 1, x = 0.375, y = .5)
    vpleg <- viewport(width = 0.25, height = 0.5, x = 0.85, y = 0.75)
    subvp <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.25)
    print(p + theme(legend.position = "none"), vp = vp1)
    upViewport(0)
    pushViewport(vpleg)
    grid.draw(legend)
    upViewport(0)
    pushViewport(subvp)
    my_table <- tableGrob(Metrica) 
    grid.draw(my_table)  
  }
}

###########################################################
RequiredtimeTest<-function (Predictors, Response, Comps = 10, Function,crosval = "LOOCV",Ticks = 20,ERIterations = 20,LassoIterations=50,cpus_per_node=1, ...) {
  # Function to calculate the time required for a particular MultiML analysis. 
  # Input: 
  # Predictors: A list of the different Omics Datasets and the Response matrix.
  # Response: A number indicating the position of the response matrix included in "Predictors" object.
  # Comps: Number of componets to be calculated after each iteration. Just applicable to PLSDA approach.
  # Function: Modular function used to calculate the error rate. It can be PLSDA.MP or Random.Forest.MP to indicate the approach to be used.
  # crosval: Type of cross validation to be applied, Leave-One-Out (LOOCV) or ten fold change (TenF).
  # Ticks: Number of segments (groups) of samples to evaluate.
  # ERIterations: Number of iterations in which the error rate (ER) will be calculated.
  # LassoIterations: Number of iterations of the Lasso selection per each error rate analysis.
  # cpus_per_node: Number of CPUs that will be used in the analysis.
  # ...: Arguments to be passed to methods.
  
  # Output:
  # Summary: A table indicating the conditions of the analysis, established by the user
  # EstimatedTime: A table of the estimated time that the process will last showed in different metrics (seconds, minutes, hours, or days).

  components=Comps
  print("Calculating the estimated time required for the classification analysis")
  Omics<-data.frame(Omics=names(Predictors[-Response]))
  Omics2<-paste(Omics$Omics,collapse="_")
  TestSumm<-data.frame(Omics=Omics2,CrossValidation=crosval,Ticks = Ticks,ERIterations=ERIterations,LassoIterations=LassoIterations,CPUs=cpus_per_node)
  rep<-length(Predictors)-1
  NewList<-Predictors
  NewList<-NewList[-Response]
  if (length(NewList)>1){
    tabe<-data.frame(sapply(NewList,dim))
    wti<-as.numeric(tabe[2,][apply(tabe[2,],1,which.min)])
    
  } else {
    tabe<-data.frame(sapply(NewList,dim))
    wti<-as.numeric(tabe[2,])
    }
  TPLSDA<-t(data.frame(summary(system.time(capture.output(ProtsPLSDA<-ClassificationErrorRate(Predictors=Predictors,Response=Response,Comps = components,
                                                                                              Function=Function,
                                                                                              crosval = crosval,Ticks = NULL,WhichTicks = wti,ERIterations = 1,LassoIterations = 1)  )))))
  TEval<-(2 + TPLSDA[,3]*ERIterations*LassoIterations*Ticks*rep)/cpus_per_node
  
  template<-c("seconds","minutes","hours","days")
  Result<-c(round(TEval, digits=3),round(TEval/60, digits=3),round(TEval/3600, digits=3),
                round(TEval/86400, digits=3))
  Result2<-as.data.frame(cbind(Result,template))
      
      res<-list(Summary=TestSumm,EstimatedTime=Result2)
      
      return(res)
      
}

###########################################################
ER_Calculator<-function(Predictors, Response,Previous_CER=NULL,Ticks=5,WhichTicks=NULL,Function=Random.Forest.MP,Comps=10,crosval = "LOOCV",ERIterations=2,LassoIterations=2,TheoreticER=NULL,ErrorRateType="ER") {
  
  # Function to evaluate through Error rate the predictive capability of the Multipower package. 
  # It encompasses "ClassificationErrorRate" and it is the fusion between the mentioned and "ER_Adder"
  # Input: 
  # Predictors: A list of different Omics Datasets and the Response matrix
  # Response: A number indicating the response matrix included in "Predictors" object
  # Previous_CER: A previous result of class "ClassificationErrorRate" (CER) to be fusioned to the one that is going to be calculated
  # Ticks: Number of segments (groups) of samples to evaluate. If NULL, the calculation is made on its own considering the TheoreticER
  # WhichTicks: Vector of numbers of ticks the user wants to analyze. If NULL, a random selection between the minimum and maximum number of samples is calculated. It results useful in posterior rounds of analysis.
  # Function: Modular function used to calculate the the error rate. It can be PLSDA.MP or Random.Forest.MP to indicate the approach to be used.
  # Comps: Number of componets to be calculated after each iteration. Just applicable to PLSDA approach.
  # crosval: Type of cross validation to be applied, Leave-One-Out (LOOCV) or ten fold change (TenF)
  # ERIterations: Number of iterations in which the error rate (ER) will be calculated.
  # LassoIterations: Number of iterations of the Lasso selection per each error rate analysis.
  
  # Output:
  # TestedTicks: A vector of the number of evaluated number of samples
  # Omics: A vector of the evaluated Omics
  # A list of two lists:
  # Minimums: A list of the minimum value of error rate, balanced (BER) or not (ER) obtained per each ten
  #           component analysis. This is measured through three distance metrics to evaluate the 
  #           classification performance of the model. Maximal distance (max.dist), distance to 
  #           centroids (centroids) or Mahalanobis distance (Mahalanobis)
  #           Thus, each table contains the results per each iteration at different subsets of samples
  
  # CompWinner: A list of the number of components in which the minimum value of error rate, 
  #             balanced (BER) or not (ER) was obtainde per each iteration. This is measured through 
  #             the three mentioned distance metrics to evaluate the classification performance 
  #             of the model. Thus, each table contains the components per each iteration at different 
  #             subsets of samples
  suppressWarnings({ 
  
  if (! is.list(Predictors) ) {
    stop("\\nOmics dataset must be a list with at least two elements")
  }
  ###########
  # Y
  Y<-t(as.matrix(Predictors[[Response]], drop=FALSE))
  if (is.character(Y[,1])) {
    Ycita<-transform(Y, Type = as.numeric(Type))
    Y<-Ycita
    rm(Ycita)
  }
  if (ncol(Y) != 1) {
    stop("\\nResponse must be a single variable")
  }
  if (any(is.na(Y))) {
    stop("\\nResponse must not contain missing values")
  }
  if (is.null(colnames(Y))) {
    colnames(Y) = "Y"
  }
  if (is.null(rownames(Y))) {
    rownames(Y) = 1:n
  }
  # Step1: Match the sample size
  LosIndivs<- rownames(Y)
  for (i in 1:length(Predictors)){
    LosIndivs = intersect(LosIndivs, colnames(Predictors[[i]]))
  }
  
  print(paste("This analysis will be performed with",length(LosIndivs),"samples, since those are the ones repeated in all layers"))
  
  NewList<-Predictors
  
  for (i in 1:length(NewList)) {
    NewList[[i]]<-NewList[[i]][,colnames(NewList[[i]]) %in% LosIndivs]
  }
  
  ###########
  # Step2: Match the order
  LosColnames<-colnames(NewList[[1]])
  
  for (i in 1:length(NewList)) {
    NewList[[i]]<-NewList[[i]][,sort(LosColnames) ]
  }
  
  TestdeMatchTable<-unique(matrix(permutations(length(NewList)),ncol=2))
  
  for (i in 1:nrow(TestdeMatchTable)) {
    a=TestdeMatchTable[i,1]
    b=TestdeMatchTable[i,2]
    
    if (all(colnames(NewList[[a]])==colnames(NewList[[b]]))){
      #print(paste("Columns of lists",a,"and",b,"are equally sorted"))
    } else{
      #print("The colnames are not equally sorted")
    }
  }
  if (crosval=="LOOCV"){
    valid="loo"
  }
  if (crosval=="TenF"){
    valid="Mfold"
  }
  
  # Y
  Y<-t(as.matrix(NewList[[Response]], drop=FALSE))
  if (is.character(Y[,1])) {
    Ycita<-transform(Y, Type = as.numeric(Type))
    Y<-Ycita
    rm(Ycita)
  }
  
  #########
  Omics<-data.frame(Omics=names(NewList[-Response]))
  
  MinN<-round(length(table(as.factor(Y[,1]))) * 2, digits = 0)
  MaxN<-nrow(Y)
  Ngroups<-length(table(as.factor(Y[,1])))
  
  if(!is.null(Ticks)) {
    vectTicks<-round(seq(from=MinN,to = MaxN,length.out = Ticks),digits = 0)
    setdeTicks = 1
    DegreesOfFreedom<-Ticks-1
  } else {
    if (is.null(WhichTicks)) {
      setdeTicks<-seq(from=5,to = MaxN,by=2)
    } else {
      vectTicks<-WhichTicks
      setdeTicks<-1
    }
  }
  if (length(setdeTicks)==1) {
    tab<-  ClassificationErrorRate (Predictors=Predictors, Response=length(Predictors),
                                    Comps=Comps,crosval = crosval,ERIterations=ERIterations,
                                    LassoIterations=LassoIterations,Ticks=NULL, 
                                    WhichTicks=vectTicks,Function = Function)

    if (is.null(Previous_CER)){
      
      Previous_CER<-tab
      allmeltedmodel<-Previous_CER$TablebyTick
    
    } else{
        ttdticks<-sort(c(tab$TestedTicks,Previous_CER$TestedTicks))
      
      if (length(tab$TablebyTick)==6) {
        tbyticks2<-NULL
        tbyticks<-list()
        for (tod in 1:length(tab$TablebyTick)) {
          namecito<-names(tab$TablebyTick[tod])
          tbyticks2<-rbind(tab$TablebyTick[[tod]],Previous_CER$TablebyTick[[tod]])
          tbyticks2<-tbyticks2[order(as.numeric(as.character(tbyticks2$variable))),]
          rownames(tbyticks2)<-seq(1:dim(tbyticks2)[1])
          tbyticks[namecito]<-list(tbyticks2)
        }
        orden<-as.numeric(as.character(tbyticks[[1]]$variable))
        orden2<-paste(orden,"samples", sep=" ")
        orden2
      } else {
        tbyticks<-rbind(tab$TablebyTick,Previous_CER$TablebyTick)
        tbyticks<-tbyticks[order(as.numeric(as.character(tbyticks$variable))),]
        rownames(tbyticks)<-seq(1:dim(tbyticks)[1])
        
        orden<-as.numeric(as.character(tbyticks$variable))
        orden2<-paste(orden,"samples", sep=" ")
        
      } 
      ######
      om<-as.data.frame(tab$Omics, drop=F)
      Omics<-om
      ######
      lstdeMins<-append(tab$Minimums,Previous_CER$Minimums)
      lstdeMins<-lstdeMins[match(orden2,names(lstdeMins))]
      
      
      if (!is.null(tab$CompWinner)) {
        lstdewinners<-append(tab$Minimums,Previous_CER$Minimums)
        lstdewinners<-lstdewinners[match(orden2,names(lstdewinners))]
      } else {}
      
      if (is.null(tab$CompWinner)) {
        Previous_CER<-list(TestedTicks=ttdticks,Omics=Omics,Minimums=lstdeMins,TablebyTick=tbyticks)
      } else {
        Previous_CER<-list(TestedTicks=ttdticks,Omics=Omics,Minimums=lstdeMins,
                           CompWinner=lstdewinners,TablebyTick=tbyticks)
      }
      
      allmeltedmodel<-Previous_CER$TablebyTick
    
    }
    
    ######
  } else {
    Rankk<-ranksas<-list()
    for (sdt in 1:length(setdeTicks)){
      numTicks<-setdeTicks[sdt]
      vectTicks<-round(seq(from=MinN,to = MaxN,length.out = numTicks),digits = 0)
      howmany<-length(vectTicks)
      
      if (is.null(Previous_CER) ){
        already=0
      } else {
        already<-length(Previous_CER$TestedTicks)
        howmany<-already+2
        vectTicks<-round(seq(from=MinN,to = MaxN,length.out = howmany),digits = 0)
        vectTicks<-c(setdiff(as.vector(vectTicks),as.vector(Previous_CER$TestedTicks) ))
        vectTicks<-sort(sample(vectTicks,size = (howmany-already), replace = FALSE ))
      }
      ######
      tab<-  ClassificationErrorRate (Predictors=Predictors, Response=length(Predictors),
                                      Comps=Comps,crosval = crosval,ERIterations=ERIterations,
                                      LassoIterations=LassoIterations,Ticks=NULL, 
                                      WhichTicks=vectTicks,Function = Function)
      
      if (is.data.frame(tab$TablebyTick)){
      colnames(tab$TablebyTick)[colnames(tab$TablebyTick)=="StdDev"] <- "SEM"
      colnames(tab$TablebyTick)[colnames(tab$TablebyTick)=="variable"] <- "Samples"
      colnames(tab$TablebyTick)[colnames(tab$TablebyTick)=="value"] <- "ER_value"
      colnames(tab$TablebyTick)[colnames(tab$TablebyTick)=="ErrorRate"] <- "ER_type"
      } else {}
      
      if (is.null(Previous_CER)){
        
        Previous_CER<-tab
        allmeltedmodel<-Previous_CER$TablebyTick
      
      } else{
        ttdticks<-sort(c(tab$TestedTicks,Previous_CER$TestedTicks))
      
        if (!is.data.frame(tab$TablebyTick)) {
          tbyticks2<-NULL
          tbyticks<-list()
          colnames <- c("Category","Samples","ER_value","SEM","Comp","ER_type","Metric","CI_90","CI_95","CI_99") 
            for (i in seq_along(tab$TablebyTick)){
                colnames(tab$TablebyTick[[i]]) <- colnames
            }
          for (tod in 1:length(tab$TablebyTick)) {
            namecito<-names(tab$TablebyTick[tod])
            tbyticks2<-rbind(tab$TablebyTick[[tod]],Previous_CER$TablebyTick[[tod]])
            tbyticks2<-tbyticks2[order(as.numeric(as.character(tbyticks2$Samples))),]
            rownames(tbyticks2)<-seq(1:dim(tbyticks2)[1])
            tbyticks[namecito]<-list(tbyticks2)
          }
          orden<-as.numeric(as.character(tbyticks[[1]]$Samples))
          orden2<-paste(orden,"samples", sep=" ")
          orden2
        } else {
          tbyticks<-rbind(tab$TablebyTick,Previous_CER$TablebyTick)
          tbyticks<-tbyticks[order(as.numeric(as.character(tbyticks$Samples))),]
          rownames(tbyticks)<-seq(1:dim(tbyticks)[1])
          
          orden<-as.numeric(as.character(tbyticks$Samples))
          orden2<-paste(orden,"samples", sep=" ")
        } 
        ######
        om<-as.data.frame(tab$Omics, drop=F)
        Omics<-om
        ######
        lstdeMins<-append(tab$Minimums,Previous_CER$Minimums)
        lstdeMins<-lstdeMins[match(orden2,names(lstdeMins))]
      
        if (!is.null(tab$CompWinner)) {
          lstdewinners<-append(tab$Minimums,Previous_CER$Minimums)
          lstdewinners<-lstdewinners[match(orden2,names(lstdewinners))]
        } else {}
        
       #####
        if (is.null(tab$CompWinner)) {
          Previous_CER<-list(TestedTicks=ttdticks,Omics=Omics,Minimums=lstdeMins,TablebyTick=tbyticks)
        } else {
          Previous_CER<-list(TestedTicks=ttdticks,Omics=Omics,Minimums=lstdeMins,
                             CompWinner=lstdewinners,TablebyTick=tbyticks)
        }
        
        allmeltedmodel<-Previous_CER$TablebyTick
        
        
      } # Close of Previous_CER
      
      if (is.null(TheoreticER) ) {
        MinimumError=min(tbyticks$value)
        SamplesEvaluated=SamplesRequired=as.numeric(as.character(
                                        Previous_CER$TablebyTick$mahalanobis.dist_ER$variable[Previous_CER$TablebyTick$mahalanobis.dist_ER$value==min(Previous_CER$TablebyTick$mahalanobis.dist_ER$value)]  ))
      } else{
        MinimumError=TheoreticER 
      }
      

      if (!is.data.frame(allmeltedmodel)){
        allmeltedmodel<-allmeltedmodel$mahalanobis.dist_ER
        ro=dim(allmeltedmodel)[1]
        colnames(allmeltedmodel)[colnames(allmeltedmodel)=="StdDev"] <- "SEM"
        colnames(allmeltedmodel)[colnames(allmeltedmodel)=="variable"] <- "Samples"
        colnames(allmeltedmodel)[colnames(allmeltedmodel)=="value"] <- "ER_value"
        colnames(allmeltedmodel)[colnames(allmeltedmodel)=="ErrorRate"] <- "ER_type"
      } else {
        ro=dim(allmeltedmodel)[1]
        } 
      
      #####
      if (allmeltedmodel$ER_value[ro-1]<allmeltedmodel$ER_value[ro]) {
        allmeltedmodel$ER_value[ro]<-allmeltedmodel$ER_value[ro-1]-abs(allmeltedmodel$SEM[ro-1])
      } else{}
      
      #####
      
      namecito<-paste(dim(allmeltedmodel)[1],"ticks", sep=" ")
      DegreesOfFreedom<-length(Previous_CER$Minimums)-1
      
      fm1 <- lm(ER_value ~ bs(as.numeric(as.character(Samples)), degree = 1,df = DegreesOfFreedom), data = allmeltedmodel)
      DESVESTprediction<-lm(SEM ~ bs(as.numeric(as.character(Samples)), degree = 1,df = DegreesOfFreedom), data = allmeltedmodel)
      ht01 <- seq(min(as.numeric(as.character(allmeltedmodel$Samples))),2000, length.out = 2000)
      prediction<-as.data.frame(cbind(predict(fm1, data.frame(Samples = ht01)),predict(DESVESTprediction, data.frame(Samples = ht01))   )  );colnames(prediction)<-c("Prediction","STDpred")
      ht01Table<-cbind(ht01,as.data.frame(prediction));colnames(ht01Table)<-c("ht","ErrorPred", "STDpred")
      
      NSampleMaxTable<-as.data.frame(ht01Table[ht01Table[,2]>=MinimumError,])
      NSampleMax<-1+round(NSampleMaxTable$ht[dim(NSampleMaxTable)[1]])
      NSampleMaxPosition<-as.numeric(rownames(NSampleMaxTable)[dim(NSampleMaxTable)[1]])
      
      NSampleMaxTable$EplusStd<-NSampleMaxTable$ErrorPred - NSampleMaxTable$STDpred
      Theplusminustable<-as.data.frame(NSampleMaxTable[NSampleMaxTable[,4]>=MinimumError,])
      Theplusminus<-round(Theplusminustable$ht[dim(Theplusminustable)[1]])
      TheplusminusPosition<-as.numeric(rownames(Theplusminustable)[dim(Theplusminustable)[1]])
      
      SamplesEvaluated<-as.numeric(as.character(allmeltedmodel$variable[dim(allmeltedmodel)[1]]))
      SamplesRequired<-round(NSampleMaxTable$ht[NSampleMaxPosition])
      MOE<-paste("±",round(NSampleMax-Theplusminus))
      
      TERtable<-rbind(ERtarget=round(MinimumError,digits = 3),
                      PSS=SamplesRequired,
                      MOE= MOE) 
      ####################
      rango<-seq(SamplesRequired-round(NSampleMax-Theplusminus),SamplesRequired+round(NSampleMax-Theplusminus), by=1)
      ranksas<-list(rango)
      Rankk[namecito]<-ranksas
      
      # STOP protocol
      if (length(Rankk)>2) {
        testedfrags<-testedfragsSas<-list()
        for (lngt in 1:length(Rankk)) {
          namecito<-names(Rankk[lngt])
          if (length(Rankk[[lngt]])<20) {
            testedfragsSas<- list(Rankk[[lngt]])
            testedfrags[namecito]<-testedfragsSas
            
          } else {}
        }
        if (length(testedfrags)>2){
          fin<-length(testedfrags)
          intrsct<-intersect(intersect(testedfrags[[fin]],testedfrags[[fin-1]]),testedfrags[[fin-2]])
          if(length(intrsct)>=3){
            if (is.data.frame(tbyticks)){
            colnames(tbyticks)[colnames(tbyticks)=="StdDev"] <- "SEM"
            colnames(tbyticks)[colnames(tbyticks)=="variable"] <- "Samples"
            colnames(tbyticks)[colnames(tbyticks)=="value"] <- "ER_value"
            colnames(tbyticks)[colnames(tbyticks)=="ErrorRate"] <- "ER_type"
            } else {
            colnames <- c("Category","Samples","ER_value","SEM","Comp","ER_type","Metric","CI_90","CI_95","CI_99") 
              for (i in seq_along(tbyticks)){
                  colnames(tbyticks[[i]]) <- colnames
              }
            }
            
            if (is.null(tab$CompWinner)) {
              res<-list(TestedTicks=ttdticks,Omics=Omics,Minimums=lstdeMins,Prediction_table=TERtable,TablebyTick=tbyticks)
            } else {
              res<-list(TestedTicks=ttdticks,Omics=Omics,Minimums=lstdeMins,Prediction_table=TERtable,
                        CompWinner=lstdewinners,TablebyTick=tbyticks)
            }
            return(res)
            break
          }else{}
        } else {}
      }
    } 
  }
  })
}

# END 

