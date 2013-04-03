#####################################################################
#
# \file script.R
#
# Copyright (c) 2009-2010, Danny Arends, with modificationd from: Ronny Joosen
# last modified Apr, 2013
# first written Nov, 2009
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: Seed Quality QTL mapping pipeline with MQM in R
#
######################################################################
set.seed(1336)

# Convert Traits to Z-scores
# Z = (x-mean(x))/sd(x)
# http://nl.wikipedia.org/wiki/Z-score
traitstoZscores <- function(traits){ invisible(apply(traits, 2, function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})) }

# Get third column from a matrix
getThird <- function(x){ x[,3] }

# Get the first column 'chr' from a matrix
getChr <- function(x){ x[,1] }

# Removes outliers with a Z-scores above threshold
# Returns: List of [[1]] trait matrix with outliers set to NA and [[2]] Number of values removed
removeoutliers <- function(traits,cutoff=3,logf,verbose=FALSE){
  zscores <- traitstoZscores(traits)
  cnt <- 0
  removed <- NULL
  for(x in 1:ncol(zscores)){
    n <- which(abs(zscores[,x])>3)
    cnt <- cnt + length(n)
    traits[n,x] = NA
    if(verbose){
      cat("Removed",length(n),"out of",nrow(traits),"as outliers\n")
    }else{
      cat("Removed",length(n),"out of",nrow(traits),"as outliers\n",file=logf,append=T)
    }
    removed <- c(removed,length(n))
  }
  cat("Summary: Removed",cnt,"out of",ncol(traits)*nrow(traits),"as outliers\n",file=logf,append=T)
  list(traits,removed)
}

# Removeoutliers function on a cross object
# Applies the: removeoutliers funtion to the cross object
removeoutliersfromcross <- function(cross,cutoff=3,logf){
  traits <- pull.pheno(cross)
  newtraits <- removeoutliers(traits,cutoff,logf=logf)
  cross$pheno <- newtraits[[1]]
  list(cross,newtraits[[2]])
}

# MQMmodel to SIF file
# Create a MQM model from the SIF data
mqmmodelsasnetwork <- function(cross,result){
  if(is.null(cross)) stop("No cross object. Please supply a valid cross object.")
  if(!any(class(result)=="mqmmulti")) stop("No mqmmulti object. Please supply a valid mqmmulti object.")

  models <- lapply(FUN=mqmgetmodel,result)
  namez <- colnames(pull.pheno(cross))
  cat(file="QTLnetwork.sif","",append=F)
  cat(file="QTLnodes.sif","",append=F)
  for(x in 1:length(models)){
    cat(file="QTLnodes.sif",namez[x],"Trait\t0\n",sep="\t",append=TRUE)
    for(y in 1:length(models[[x]][[2]])){
      cat(file="QTLnodes.sif",models[[x]][[2]][y],"Marker",models[[x]][[4]][y],"\n",sep="\t",append=TRUE)
      cat(file="QTLnetwork.sif",namez[x],"QTLeffect",models[[x]][[2]][y],result[[x]][models[[x]][[2]][y],3],"\n",sep="\t",append=TRUE)
    }
  }
}


# Get multi chromosome QTL support intervals
# Adaptation of the original work by K. Broman supportint
mqmsupportint <- function(result, marker){
  num <- which(rownames(result)==marker)
  result <- result[result[num,1]==result[,1],]
  num <- which(rownames(result)==marker)
  if(result[num,3] > 1.5){
    drop <- 1.5
  }else{
    drop <- 0.75 * result[num,3]
  }
  temp <- which(result[,3] < result[num,3]-drop)
  min <- temp[which(temp < num)[length(which(temp < num))]]
  max <- temp[which(temp > num)[1]]
  if(is.na(min&&1)){ min <- 1 }
  if(is.na(max&&1)){ max <- nrow(result) }
  list(result[min,],result[max,])
}

# Plot all circleplots after a multi trait MQM mapping
doplotcircles <- function(cross,result,save=FALSE,spacing=100,plotFUN,extension,w,h){
  n <- 1
  traitnames <- colnames(pull.pheno(cross))
  while(n <= length(result)){
    if(save) plotFUN(file=paste(traitnames[n],"/circleplot.",extension,sep=""),w=w,h=h)
    mqmplot.circle(cross,result,spacing=spacing,highlight=n)
    if(save) dev.off()
    n <- n+1
  }
}

# Calculate the reciproce
reciproce <- function(x){ return(1/x) }

# Does nothing but return the original values
donothing <- function(x){ return(x) }

# Do normalization on a cross object (return the same crossobject with traits normalized)
# Normalization methods are stored in: cross$transformations
# Methods tried: donothing, probit, log, sqrt, reciproce
DoNormalization <- function(cross,logf,threshold=0.05,transformations=NULL,verbose=TRUE){
  originalfile <- cross
  newcross <- cross
  methods <- c(donothing,probit,log,sqrt,reciproce)
  chosen <- NULL
  for(x in 1:nphe(cross)){
    cnt <- 1
    suc <- FALSE
    if(!is.null(transformations) && length(transformations)==ncol(pull.pheno(cross)) && transformations[x]!=0){
      cross$pheno[[x]] <- methods[transformations[x]](originalfile$pheno[[x]])
      chosen <- c(chosen,transformations[x])
    }else{
      for(m in methods){
        if(!any(originalfile$pheno[[x]]<0,na.rm=T)){
          cross$pheno[[x]] <- m(originalfile$pheno[[x]])
        }else{
          if(verbose) cat("[script] Method ",cnt," skipped for phenotype ",x,"\n",file=logf,append=T)
        }
        if(sum(is.na(cross$pheno[[x]])) > sum(is.na(originalfile$pheno[[x]])) ){
          if(verbose) cat("[script] NA produced skipping ",cnt," for phenotype ",x,"\n",file=logf,append=T)
        }else if(any(is.infinite(cross$pheno[[x]]))){
          if(verbose) cat("[script] Infinity produced skipping ",cnt," for phenotype ",x,"\n",file=logf,append=T)
        }else{
          if(var(cross$pheno[[x]],na.rm=T)==0){
            if(verbose) cat("[script] No variation, skipping ",cnt," for phenotype ",x,"\n",file=logf,append=T)
          }else{
            if(!mqmtestnormal(cross, pheno.col = x,significance=threshold)){
              if(verbose) cat("[script] Phenotype ",x," NOT normal after ",cnt,"\n",file=logf,append=T)
            }else{
              if(verbose) cat("[script] Phenotype ",x," normal after ",cnt," (1=No transformation, 2=probit, 3=log, 4=sqrt, 5=reciprocal)\n",file=logf,append=T)
              chosen <- c(chosen,cnt)
              newcross$pheno[[x]] <- cross$pheno[[x]]
              suc <- TRUE
              break;
            }
          }
        }
        cnt <- cnt+1
      }
      if(!suc){
        chosen <- c(chosen,0)
        newcross$pheno[[x]] <- originalfile$pheno[[x]]
      }
    }
  }
  newcross$transformations <- chosen
  newcross
}

#Main routine:
#1.0 Read the configFILE
#1.1 Check if the required parameters are supplied 
#2.0 Remove Outliers
#2.1 Creating a normalized cross object
#3.0 Main analysis loop
#3.1 PLot of basic statistics, Genetic Map, Distributions, Etc
#4.0 QTL by MQM
#4.1 QTL by scanone and scan.two
#4.2 QTL plots
#5.0 Based on MQM results -> Circleplots, Heatmaps, Effect plots, Grouping
#5.1 PCA analysis on the first 3 PCs
#5.2 Output of SIF network images
DoAnalysisOn <- function(filename,directory){
  if(is.null(filename)){stop("Supply a configfile")}
  if(is.null(directory)){stop("Specify the full path to the directory with configfile")}
  if(!file.exists(directory)) stop(paste("No such directory:",directory))
  setwd(directory)
  if(!file.exists(filename)) stop(paste("No such file:",filename))
  log <- file("log.txt", "w")
  configfile <- scan(file=filename,what="raw",sep="\n",quiet = TRUE)
  stamp <- file.info(filename)$mtime
  required <- c("memorylimitN", "maindir", "verbose", "datafileraw", "crosstype", "resultfile", "tablesep", "imageoutput", "imagewidthI",
                "imageheightI", "stepsizeN", "windowsizeN", "cofactorfile","setalfaN", "npermutationsI", "ncoresI", "interactionstrenghtN", "plothistogram",
                "ploteffects", "plotinteractions", "plotcircles", "circlespacingI", "plotheatmap", "heatmapcolors", "heatmapbreaks","plotmap",
                "plotclusteredheatmap", "plotclustergroups", "clustergroupcutoffI", "normalitythresholdN", "normalize")
  cat("[script] log started at ",Sys.time(),"\n",file=log)
  for(line in configfile){
    line <- gsub(" ","",line)
    ta = strsplit(line,"=",fixed=TRUE)
    if(strsplit(line,"")[[1]][1] != "#" && !is.na(ta[[1]][2])) {
      temp = ta[[1]][2]
      if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") { temp = as.integer(temp) }
      if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") { temp = as.numeric(temp) }
      assign(ta[[1]][1],temp)
      cat("> assigned to ",ta[[1]][1]," value: ",temp,"\n",file=log,append=T)  
    }
  }
  cat("[script] Configuration parsed\n",file=log,append=T)
  #Check if the user supplied all
  Error <- 0
  for(x in 1:length(required)){
    if(exists(required[x])){
    if(is.null(get(required[x]))){
      Error <- 1
      cat("[script] Parameter: ",required[x]," missing in configfile\n",file=log,append=T)
      stop("Parameter: ",required[x]," missing in configfile")
    }
    }else{
      Error <- 1
      cat("[script] Parameter: ",required[x]," missing in configfile\n",file=log,append=T)
      stop("Parameter: ",required[x]," missing in configfile")  
    }
  }
  if(!Error){
    mycolorrange <- strsplit(heatmapcolors,",")[[1]]
    mycolorbreaks <- as.numeric(strsplit(heatmapbreaks,",")[[1]])
    if(length(mycolorbreaks)-1!=length(mycolorrange)){
      cat("There should be 1 more break then colors defined in the configfile\nUsing default colors and range")
      mycolorrange <- c("darkblue","blue","lightblue","yellow","orange","red")
      mycolorbreaks <- c(-100,-10,-3,0,3,10,100)
    }else{
      cat("[script] Colors seem owkay\n",file=log,append=T)
    }
    cat("[script] No errors\n",file=log,append=T)
    cat(" -- Remember: Cancelation is always possible [Esc], your result will be saved.\n")
    cat("Phase 1/3 - Data loading and pre-processing.\n")
    memory.limit(memorylimitN)          # Set maximum memory for R
    step.size <- stepsizeN
    window.size <- windowsizeN
    nbootstrap <- npermutationsI
    compare <- mapmethodI               # Set to compare mqm with scanone = 1, no comparison = 0
    setalfa <- setalfaN
    cat("[script] Changing to working area\n",file=log,append=T)
    setwd(maindir)
    require(qtl,warn.conflicts = FALSE,quietly = TRUE)           # R/QTL library for QTL mapping
    require(snow,warn.conflicts = FALSE,quietly = TRUE)          # SNOW for multicore analysis
    require(VGAM,warn.conflicts = FALSE,quietly = TRUE)          # VGAM for probit and other transformations
    qtlversion()
    cat("[script] Trying to read in the raw data\n",file=log,append=T)
    mycrossraw <- read.cross(format="csv",file=datafileraw,genotypes=c("AA","BB"))
    cat("[script] Gonna remove outliers using a basic Z-score transformation on the raw data\n",file=log,append=T)
    output <- removeoutliersfromcross(mycrossraw,3,logf=log)
    mycrossraw <- output[[1]]
    cat("[script] Gonna start normalization\n",file=log,append=T)
    if(!exists("datafilenormalized")){
      if(normalize=="on"){
        cat("[script] No normalized file, doing normalization\n",file=log,append=T)
        mycross <- DoNormalization(mycrossraw,log,threshold=normalitythresholdN)
        normalizations <- mycross$transformations
        cat("[script] Going to write out the normalized data file\n",file=log,append=T)
        write.cross.csv(mycross,file="normalized_cross")
      }else{
        cat("[script] No normalized file, NOT doing normalization\n",file=log,append=T)
        mycross <- mycrossraw
      }
    }else{
      cat("[script] Normalized file, using that\n",file=log,append=T)
      mycross <- read.cross(format="csv",file=datafilenormalized,genotypes=c("AA","BB"))
      normalizations <- NULL
    }

    class(mycross)[1] <- crosstype
    class(mycrossraw)[1] <- crosstype

    # Print a summary of the loaded data
    summary(mycross)
    cat("[script] Filling genotype data\n",file=log,append=T)
    # Fill holes in genotype data
    mycross <- fill.geno(mycross)
    mycrossraw <- fill.geno(mycrossraw)
    
    #Remove any phenotypes with no varaition
    toremove <- NULL
    toremove <- (which(apply(mycrossraw$pheno,2,var,na.rm=T)==0))
    if(!is.na(toremove&&1)){
      cat("[script] Removing phenotype(s) without variation\n",file=log,append=T)
      cat("[script] Removing phenotype(s) without variation\n")
      mycrossraw$pheno <- mycrossraw$pheno[,-toremove]
      mycross$pheno <- mycross$pheno[,-toremove]
    }
    
    cat("[script] Sim genotype data\n",file=log,append=T)
    # Simulate genotypes given observed marker data
    mycrossraw <- sim.geno(mycrossraw,n.draws=32,step=stepsizeN)
    mycross <- calc.genoprob(mycross,step=stepsizeN)
    
    cat("[script] Starting cofactors\n",file=log,append=T)
    # Try to find if we use cofactors in the analysis
    cof <- read.table(file=cofactorfile,header=F,sep="\n")
    if(is.numeric(cof[1,])){
      cat("[script] Numeric cofactor file\n",file=log,append=T)
      cof <- as.numeric(unlist(cof))
      if(sum(nmar(mycross))!= length(cof)){
        cat("[error] Cofactor file doesn't match in length\n",file=log,append=T)
        stop("Cofactor file doesn't match in length found:",length(cof)," needed:",sum(nmar(mycross)),"\n")
      }
    }else{
      cat("[script] Text cofactor file\n",file=log,append=T)
      cofactors <- as.character(unlist(cof))
      toset <- NULL
      for(x in cofactors){
        toset <- c(toset,find.markerindex(mycross,x))
      }
      cat("WARNING removing",length(which(is.na(toset))),"Cofactors, they couldn't be matched to the cross object")
      toset <- toset[-which(is.na(toset))]
      cof <- mqmsetcofactors(mycross,cofactors = toset)
    }
    cat("[script] Cofactors:",sum(cof),"\n",file=log,append=T)
    cat("[script] Done with the cofactor file\n",file=log,append=T)
    
    # Initialize the plot functions
    plotFUN <- jpeg
    verbose <- TRUE
    if(imageoutput=="jpg"){ plotFUN <- jpeg; extension <- ".jpg" }
    if(imageoutput=="bmp"){ plotFUN <- bmp;  extension <- ".bmp" }
    if(imageoutput=="png"){ plotFUN <- png;  extension <- ".png" }
    if(imageoutput=="pdf"){ plotFUN <- pdf;  extension <- ".pdf" }
    if(imageoutput=="eps"){ plotFUN <- postscript; extension <- ".eps" }    
    if(verbose=="false") verbose <- FALSE

    cat("[script] Starting analysis\n",file=log,append=T)    
    mqmallres <- generateall(cross=mycross,crossraw=mycrossraw,cof=cof,
                             workdir=maindir,nbootstrap=nbootstrap,compare=compare,
                             setalfa=setalfa,verbose=verbose,imagewidthI=imagewidthI,
                             imageheightI=imageheightI,ncoresI=ncoresI,
                             interactionstrenghtN=interactionstrenghtN,stepsizeN=stepsizeN,
                             windowsizeN=windowsizeN,ploteffects=ploteffects,plothistogram=plothistogram,
                             plotinteractions=plotinteractions,plotcircles=plotcircles,plotheatmap=plotheatmap,
                             plotclusteredheatmap=plotclusteredheatmap,plotclustergroups=plotclustergroups,
                             plotmap=plotmap,plotcorrelations=plotcorrelations, circlespacingI=circlespacingI, 
                             clustergroupcutoffI=clustergroupcutoffI, normalizations=normalizations, outliers=output[[2]],
                             plotFUN=plotFUN,extension=extension,log=log,stamp=stamp,tablesep=tablesep,mycolorrange=mycolorrange,mycolorbreaks=mycolorbreaks)
  }
  invisible(mqmallres)
}

########################
# generateall
########################
# Loop function to analyse all traits in mycross / mycrossraw
########################
generateall <- function(cross,crossraw,cof,workdir,nbootstrap = 0,compare = 0, setalfa = 0.05,
                        verbose=FALSE,imagewidthI=800,imageheightI=600,ncoresI=2,
                        interactionstrenghtN=2,stepsizeN=5,windowsizeN=20,ploteffects="on",
                        plothistogram="on",plotinteractions="on",plotcircles="on",plotheatmap="on",
                        plotclusteredheatmap="on",plotclustergroups="on", plotmap="on",plotcorrelations="on",
                        circlespacingI=100,clustergroupcutoffI=10,normalizations,outliers,plotFUN=jpeg,extension=".jpg",log=NULL,stamp="",tablesep=";",mycolorrange,mycolorbreaks){
	mqmmulti <- vector(mode = "list", length = nphe(cross))
	if(plotmap=="on"){
    cat("[script] Going to do recombinationmap plots\n",file=log,append=T)
    filename = paste("recombinationmap",extension,sep="")
    plotFUN(file = filename, width = imagewidthI, height = imageheightI)
      crossraw <- est.rf(crossraw)
      plot.rf(crossraw)
    dev.off()
  }
  if(plotcorrelations=="on"){
      cat("[script] Going to do trait correlation plot\n",file=log,append=T)
      cat(" -- Trait correlation plot\n")
      plotFUN(file = paste("TraitCorrelations",extension,sep=""), width = imagewidthI, height = imageheightI)
      image(x=1:nphe(cross),y=1:nphe(cross),z=cor(pull.pheno(cross),use="pairwise.complete.obs"),xlab="Traits",ylab="Traits",main="Trait correlation structure")
    dev.off()
  }
  #Loop through all the phenotypes scanning them
  cat("Phase 2/3 - Creating QTL statistics.\n")
  traitnames <- colnames(pull.pheno(cross))
  cat("",file="interactions.sif")
  cat("... Cleared interactions\n")
  for(x in 1:nphe(cross)){
    waitForEscape()
		cat("|",traitnames[x]," [",x,"/",nphe(cross),"]")
    cat("[script] Phenotype",traitnames[x]," [",x,"/",nphe(cross),"]\n",file=log,append=T)
		setwd(workdir)
    dir.create(paste(traitnames[x]),showWarnings=FALSE)
    setwd(paste(workdir,"/",traitnames[x],sep=""))
    reloaded <- FALSE
    if(file.exists("s.log")){
      if(readLines("s.log",1)==stamp){
        cat("... Match timestamp 'configuration'\n")
        if(file.exists(paste("mqmdata.Rdata",sep=""))){
          load(paste("mqmdata.Rdata",sep="")) #LOAD mqmres that we saved
          mqmmulti[[x]] <- mqmres
          cat("... Skipping\n")
          reloaded <- TRUE
        }
      }else{
        cat("... Mismatch timestamp 'configuration'\n")
      }
    }else{
       cat("... Missing timestamp file\n")
    }
    if(!reloaded){
      tryCatch(
      mqmmulti[[x]] <- plotalleffects(cross=cross,crossraw,cof,
        pheno.col=x,workdir=paste(workdir,"/",traitnames[x],sep=""), nbootstrap=nbootstrap, 
        compare=compare,setalfa=setalfa,verbose=verbose,imagewidthI=imagewidthI,
        imageheightI=imageheightI,ncoresI=ncoresI,stepsizeN=stepsizeN,windowsizeN=windowsizeN,ploteffects=ploteffects,plothistogram=plothistogram,
        plotinteractions=plothistogram,normalizations=normalizations,outliers=outliers,plotFUN=plotFUN,extension=extension,log=log,tablesep=tablesep,maindir=workdir)
      , error = function(e){cat("MAYOR error with trait:",x,": ",e[[1]],"\nTrying to continue\n")})
      cat(as.character(stamp),"\n",sep="",file="s.log")
    }
	}
	class(mqmmulti) <- c(class(mqmmulti),"mqmmulti")
	setwd(workdir)
  cat("Phase 3/3 - Creating summary statistics.\n")
  allqtlmatrix <- matrix(unlist(lapply(mqmmulti,getThird)),nrow(mqmmulti[[1]]),nphe(cross))
  if(plotcorrelations=="on"){
    cat(" -- Correlation plots, please wait...\n")
    cat("[script] Going to do QTL correlation plot\n",file=log,append=T)
    plotFUN(file = paste("QTLCorrelations",extension,sep=""), width = imagewidthI, height = imageheightI)
      image(1:nphe(cross),1:nphe(cross),cor(allqtlmatrix,use="pairwise.complete.obs"),xlab="Traits",ylab="Traits",main="QTL correlation structure")
    dev.off()
  }
  if(plotcircles=="on"){
    cat(" -- Circle plots, please wait...\n")
    cat("[script] Going to do overview circle plot\n",file=log,append=T)
    filename = paste("CircleAllTraits",extension,sep="")
    plotFUN(file = filename, width = imagewidthI, height = imageheightI)
      mqmplot.circle(crossraw,mqmmulti,spac=circlespacingI,interactstrength=interactionstrenghtN)
    dev.off()
    cat("[script] Per trait circle plots\n",file=log,append=T)    
    doplotcircles(crossraw,mqmmulti,save=T,spac=circlespacingI, plotFUN, extension, imagewidthI, imageheightI)
    cat("[script] Per trait circleplots done\n",file=log,append=T)    
  }
  #Do a heatmap
  if(plotheatmap=="on"){
    cat(" -- Heat maps, please wait...\n")
    cat("[script] Going to do Heatmaps\n",file=log,append=T)
    filename = paste("Heatmap",extension,sep="")
    plotFUN(file = filename, width = imagewidthI, height = imageheightI)
      LODprofiles <- mqmplot.heatmap(crossraw,mqmmulti,col=mycolorrange,breaks=mycolorbreaks)
    dev.off()
    filename=paste("QTL_profiles",".csv",sep="")
    write.table(LODprofiles,file=filename, sep="\t")    
  }

  #Do a clustering on that heatmap
  if(plotclusteredheatmap=="on"){
    cat("[script] Going to do Clustering on the heatmap\n",file=log,append=T)
    filename = paste("HeatmapCluster",extension,sep="")
    plotFUN(file = filename, width = imagewidthI, height = imageheightI)
      cresults <- mqmplot.clusteredheatmap(crossraw,mqmmulti,col=mycolorrange,breaks=mycolorbreaks)
    dev.off()
    cat("[script] Going to do Clustering on the heatmap\n",file=log,append=T)
    filename = paste("HeatmapCluster_tree",extension,sep="")
    plotFUN(file = filename, width = imagewidthI, height = imageheightI)
       plot(cresults$Rowv)
    dev.off()
  }

  #Extract the groups (25 is a cutoff, use :  plot(cresults$Rowv) to get an idea what a good cutoff would be for the tree
  if(plotclustergroups=="on"){
    cat(" -- Cluster analysis.\n")
    dir.create("GroupProfiles",showWarnings=FALSE)
    cat("[script] Groups from data based on QTLs\n",file=log,append=T)
    groups <- groupclusteredheatmap(cross,cresults,clustergroupcutoffI)
    cnt <- 1
    temp <- lapply(mqmmulti, getThird)
    qtlmatrix <- do.call("rbind", temp)

    for(group in groups){
      cat("[script] Mean QTL profile calculation for group:",group,"\n",file=log,append=T)
      filename = paste("GroupProfiles/QTLsCluster",cnt,extension,sep="")
      cat(paste("Cluster",cnt),"\t",names(pull.pheno(cross))[group],"\n",file="QTLclusterDescr.txt",append=TRUE)
      plotFUN(file = filename, width = imagewidthI, height = imageheightI)
        mqmplot.multitrait(cross,type="lines",result = mqmmulti,group=group,meanpro="mean")
      dev.off()
      cnt <- cnt + 1
      #todo: pca en dan mqm/scan.one
      if(length(group) > nind(cross)){
        cat("[script] Principale component profile calculation for group: ",group,"size:",length(group),"\n",file=log,append=T)
        tryCatch(
          doPrincipalComponents(cross,group, cnt,cof,setalfa,stepsizeN,windowsizeN, plotFUN, imagewidthI, imageheightI,extension, verbose,log)
        , error = function(e){cat("PCA error: ",e[[1]],"\nMoving on\n")})
      }else{
        cat("[script] Principale component profile calculation for group: ",group,"size:",length(group)," skipped\n",file=log,append=T)
      }
    }
  }
  #if(plotnetwork=="on"){
    cat("[script] Creating a cytoscape network\n",file=log,append=T)
    mqmmodelsasnetwork(cross,mqmmulti)
  #}
  cat("[script] Done\n",file=log,append=T)
  cat("Analysis complete\n")
	mqmmulti
}

########################
# waitForEscape
########################
waitForEscape <- function(){
  if(runif(1) < 0.3){
    if(runif(1) < 0.05){
      cat(" -- Remember: Cancelation is always possible [Esc].\n")
    }else{
      cat(" -- Remember: Results are saved when stopping the analysis [Esc].\n")
    }
  }else{
    #Small easter egg for people running the script a lot
    if(runif(1) < 0.05){
      cat(" -- Remember: http://en.wikipedia.org/wiki/Repetitive_strain_injury\n")
    }
  }
}

########################
# doPrincipalComponents
########################
doPrincipalComponents <- function(cross, group, cnt,cof,setalfa,stepsizeN,windowsizeN, plotFUN, imagewidthI, imageheightI, extension, verbose,log){
  cat("[script] PCA\n",file=log,append=T)
  first3comps <- prcomp(na.exclude(t(cross$pheno[,group])),scale=T)[[2]][,1:3]
  tempcross <- cross
  tempcross$pheno <- first3comps
  cat("[script] PCA done",dim(first3comps)[1],dim(first3comps)[2],"=",dim(cross$pheno)[1],dim(cross$pheno)[2],"\n",file=log,append=T)
  plotFUN(file = paste("PcompsCluster",cnt,extension,sep=""), width = imagewidthI, height = imageheightI)
    plot(scanone(fill.geno(tempcross),pheno.col=1, model="normal", method="hk"),
         scanone(fill.geno(tempcross),pheno.col=2, model="normal", method="hk"),
         scanone(fill.geno(tempcross),pheno.col=3, model="normal", method="hk"),col=c("red","green","blue"),lwd=2)
    legend("topright",c("pc1","pc2","pc3"),col=c("red","green","blue"))
  dev.off()
  cat("[script] PCA MQM scanning done\n",file=log,append=T)
  first3comps
}

########################
# plotalleffects
########################
plotalleffects <- function(cross,crossraw,cof,pheno.col=1,workdir,nbootstrap = 0,compare = 0, setalfa = 0.05, 
                           verbose=FALSE, imagewidthI=800,imageheightI=600,ncoresI=2,stepsizeN=5,windowsizeN=20,
                           ploteffects="on",plothistogram="on",plotinteractions="on",normalizations,outliers,plotFUN=jpeg,extension=".jpg",log=NULL,tablesep=";",maindir){
  setwd(workdir)
  traitname <- colnames(pull.pheno(cross))[pheno.col]
	if(compare > 0){
		#scanone method for interval mapping without cofactors
    cat("[script] Scanone because we compare back\n",file=log,append=T)
    cat(" -- Scanone\n")
		capture.output(scanone.hk <- scanone(cross=cross, pheno.col=pheno.col, model="normal", method="hk"),file=log,append=T)
		scanonesum <- summary(scanone.hk)
    cat("[script] Scanone Done\n",file=log,append=T)
	}
  
  #perform effectscan on non transformed data, plot and save as .csv
  if(ploteffects=="on"){
    cat("[script] Effects on RAW datafile\n",file=log,append=T)
    filename=paste("QTL_effectscan",pheno.col,extension,sep="")
    plotFUN(file = filename, width = imagewidthI, height = imageheightI)
    eff <- effectscan(crossraw,pheno.col=pheno.col,chr=1:nchr(crossraw),get.se=T,add.legend=F,ylab="Trait value",main=colnames(crossraw$pheno)[pheno.col])
    dev.off()
    #eff <- eff[which(rownames(eff)==rownames(mqmres)),]
    #
    filename=paste("QTL_eff",pheno.col,".csv",sep="")
    write.table(eff,file=filename,sep="\t")
    cat("[script] Effects on RAW datafile\n",file=log,append=T)
  }
  
	#mqm method with automatic cofactor selection via backward selection and ML
  cat("[script] MQM starting\n",file=log,append=T)
  cat(" -- Multiple QTL Mapping\n")
	capture.output(mqmres <- mqmscan(fill.geno(cross),cofactors=cof,pheno.col=pheno.col,cofactor.significance=setalfa, step.size=stepsizeN, window.size=windowsizeN,verbose=verbose),file=log,append=T)
  filename = paste("QTL_trait",pheno.col,extension,sep="")
  
  plotFUN(file = filename, width = imagewidthI, height = imageheightI)
    mqmplot.singletrait(mqmres)
  dev.off()
  
  cat("[script] MQM Done\n",file=log,append=T)
	mqmsum<-summary(mqmres)
  #Scale and save data for MQM before we do anything else
  if(ploteffects=="on"){
    mqmres <- mqmres[which(rownames(mqmres) %in% rownames(eff)),]
  }
	write.table(mqmres,file=paste("MQM_output",pheno.col,".csv",sep=""), sep="\t")
  cat("[script] Results saved\n",file=log,append=T)
	if(compare > 0){
    cat("[script] Scantwo with cofactors\n",file=log,append=T)
    cat(" -- Scantwo plot\n")
		#scantwo method with selected cofactors to identify putative interactions
		if(!is.null(attr(mqmres,"mqmmodel"))){
			capture.output(scantwo.hk <- scantwo(cross=cross, pheno.col=pheno.col, method="hk", addcovar = attr(mqmres,"mqmmodel")[1][[1]][,,1]),file=log,append=T)
		}else{
			capture.output(scantwo.hk <- scantwo(cross=cross, pheno.col=pheno.col, method="hk"),file=log,append=T)
		}
		sumscantwo <- summary(scantwo.hk)
		filename=paste("QTL_scantwo",pheno.col,extension,sep="")    
		plotFUN(file = filename, width = imagewidthI, height = imageheightI)
      plot(scantwo.hk)
		dev.off()	 
    cat("[script] Scantwo done\n",file=log,append=T)    
	}
  #save QTL data for scanone n scantwo methods to .csv file 
	if(compare > 0){
		write.table(scanone.hk,file=paste("SCANONE_output",pheno.col,".csv",sep=""),sep="\t")
		write.table(sumscantwo,file=paste("SCANTWO_output",pheno.col,".csv",sep=""), sep="\t")
    cat("[script] Comparison results saved\n",file=log,append=T)

    filename=paste("QTL_compare",pheno.col,extension,sep="")
    plotFUN(file = filename, width = imagewidthI, height = imageheightI)
		plot(mqmres,scanone.hk,incl.markers=TRUE, alternate.chrid=TRUE,lwd=c(2,1),lty=c(1,2),col=c("black","green"))
    legend("topright",lwd=c(2,1),lty=c(1,2),col=c("black","green"),c("mqmres","scanone"))
    dev.off()
    cat("[script] Comparison plot saved\n",file=log,append=T)
	}else{
		if(verbose) cat("No comparison\n")
  }
	
  #plot histogram of phenotype distribution from transformed data and non transformed data
  if(plothistogram=="on"){
    filename=paste("PhenotypeHistogram",pheno.col,extension,sep="")
    plotFUN(file = filename, width = imagewidthI, height = imageheightI)
    hist(cross$pheno[[pheno.col]],freq=F,breaks=25)
    dev.off()	
    
    filename=paste("PhenotypeHistogramraw",pheno.col,extension,sep="")
    plotFUN(file = filename, width = imagewidthI, height = imageheightI)
    hist(crossraw$pheno[[pheno.col]],freq=F,breaks=25)
    dev.off()
    cat("[script] Histograms done\n",file=log,append=T)
  }
  #bootstrap analysis per trait
	if(nbootstrap > 0){
    cat(" --Trait permutation, please wait...\n")
    cat("[script] Going to permutate trait\n",file=log,append=T)
		filename=paste("QTL_bootstrap",pheno.col,extension,sep="")
		plotFUN(file = filename, width = imagewidthI, height = imageheightI)
		bresults <- mqmpermutation(fill.geno(cross),mqmscan,pheno.col=pheno.col, cofactors=cof, plot=T, n.perm=nbootstrap,batchsize=25,multicore=TRUE,n.clusters=ncoresI,verbose=T)
		dev.off()
	
		permObject <- mqmprocesspermutation(bresults)
		filename=paste("QTL_perm",pheno.col,".csv",sep="")
		write.table(summary(permObject),file=filename,sep="\t")  
    cat("[script] Bootstrap done\n",file=log,append=T)
    cat(" --Trait permutation, DONE...\n")
	}
  #Create Model output 
  if(pheno.col==1){
    cat("ID", tablesep, "Trait", tablesep, "QTL", tablesep, "NAvalues", tablesep, "Outliers", tablesep, "Normalization", tablesep,file="../QTLsummary.csv",append=T,sep="")   
    cat("TopMarker",tablesep,"Chromosome",tablesep,"LOD",tablesep,"Direction",tablesep,"Interval_Start",tablesep,"Peak",tablesep,"Interval_End\n",file="../QTLsummary.csv",append=T,sep="")
  }
  if(!is.null(attr(mqmres,"mqmmodel"))){
    model <- summary(mqmgetmodel(mqmres))
    for(x in 1:nrow(model)){
      cat(paste(pheno.col,".",x,sep=""),tablesep,traitname,tablesep,file="../QTLsummary.csv",append=T,sep="")
      cat(paste("QTL",x,sep=""),tablesep,sum(is.na(cross$pheno[,pheno.col])),tablesep,outliers[pheno.col],file="../QTLsummary.csv",append=T,sep="")
      if(!is.null(normalizations)){
        cat(tablesep,normalizations[pheno.col],file="../QTLsummary.csv",append=T,sep="")
      }else{
        cat(tablesep,"NA",file="../QTLsummary.csv",append=T,sep="")
      }      
      effRatio <- eff[as.character(model[x,1]),3]/abs(eff[as.character(model[x,1]),3])
      cat(tablesep,as.character(model[x,1]),tablesep,as.character(model[x,2]),tablesep,mqmres[as.character(model[x,1]),3],tablesep, effRatio,file="../QTLsummary.csv",append=T,sep="")
      cat(tablesep,as.character(mqmsupportint(mqmres,model[x,1])[[1]][1,2]),tablesep,as.character(model[x,3]),tablesep,as.character(mqmsupportint(mqmres,model[x,1])[[2]][1,2]),"\n",file="../QTLsummary.csv",append=T,sep="")
    }
  }else{
    cat(paste(pheno.col,".",0,sep=""),tablesep,traitname,tablesep,sum(is.na(cross$pheno[,pheno.col])),tablesep,outliers[pheno.col],file="../QTLsummary.csv",append=T)
    if(!is.null(normalizations)){
      cat(tablesep,normalizations[pheno.col],file="../QTLsummary.csv",append=T)
    }else{
      cat(tablesep,"NA",file="../QTLsummary.csv",append=T)
    }
    cat(rep(paste(tablesep,"NA"),7),"\n",file="../QTLsummary.csv",append=T)
  }
 
  #Create effect and interaction plots for all significant QTL
  if(plotinteractions=="on"){
    cat(" -- Interaction plots\n")
    if(!is.null(attr(mqmres,"mqmmodel"))){
      cat("[script]Interaction plots between significant cofactors\n",file=log,append=T)
      markers <- NULL
      for(x in 1:length(attr(mqmres,"mqmmodel")$chr)){
        sumres <- summary(mqmres)[x,c(1,2,3)]
        chr <- as.numeric(attr(mqmres,"mqmmodel")$chr[x])
        pos <- as.numeric(attr(mqmres,"mqmmodel")$pos[x])
        markers <- c(markers,find.marker(cross,chr,pos))
      }
      for(x in 1:length(markers)){
        for(y in x:length(markers)){
          name1 <- markers[x]
          name2 <- markers[y]
          if(name1!=name2){
            filename=paste("Interaction_",name1,"_",name2,extension,sep="")
            plotFUN(file = filename, width = imagewidthI, height = imageheightI)
            eff <- effectplot(crossraw,pheno.col=pheno.col,mname1=name1,mname2=name2)
            dev.off()

            setwd(maindir)
            changeA <- (eff$Means[1,2]-eff$Means[1,1])
            changeB <- (eff$Means[2,2]-eff$Means[2,1])
            if(!is.na(changeA) && !is.na(changeB) && !any(is.na(eff$SEs))){
              if(abs(abs(changeA)-abs(changeB)) > 2*mean(eff$SEs)){
                cat(name1,traitname,name2,"\n",file="interactions.sif",append=T)
              }
            }
            setwd(workdir)
          }else{
            filename=paste("Effect",name1,extension,sep="")
            plotFUN(file = filename, width = imagewidthI, height = imageheightI)
            plot.pxg(crossraw,pheno.col=pheno.col,marker=name1)
            dev.off()
          }
        }
      }
    }else{ cat("[script] No significant cofactors\n",file=log,append=T) }
    cat("[script] Interaction plots finished\n",file=log,append=T)
  }
  cat("[script] Returning result\n",file=log,append=T)
  save(mqmres,file="mqmdata.Rdata") #Save the file as last, so if we redo we know that if its created the plots are also there
	mqmres
}

#res <- DoAnalysisOn("config.txt","D:/GBIC/Joosen/Publication")

