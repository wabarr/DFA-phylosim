if(!require(phytools)) install.packages("phytools")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(reshape2)) install.packages("reshape2")
if(!require(MASS)) install.packages("MASS")
if(!require(caper)) install.packages("caper")
if(!require(plyr)) install.packages("plyr")
if(!require(stringr)) install.packages("stringr")

###########define functions to be used########
###########define functions to be used########
###########define functions to be used########

geomean<-function(df,varz){
  #takes a dataframe of measurements as argument, returns vector of geomeans for each row
  #for the selected measurements only
  df<-df[,varz]
  apply(df,MARGIN=1,FUN=function(x) prod(x)^(1/length(x)))
}

getSimulationParametersFromFileName<-function(filename){
  parameters=list()
  parameters[["filename"]]<-filename
  #identify the source dataset
  ifelse(str_detect(filename,"lowR"),parameters[["datafile"]]<-"~/Dropbox/WAB Dissertation/Chapter 2 - Methods/dataForSims_fixedS_lowR.txt",
         ifelse(str_detect(filename,"0R"),parameters[["datafile"]]<-"~/Dropbox/WAB Dissertation/Chapter 2 - Methods/dataForSims_fixedS_0R.txt",
                ifelse(str_detect(filename,"highR"),parameters[["datafile"]]<-"~/Dropbox/WAB Dissertation/Chapter 2 - Methods/dataForSims_fixedS_highR.txt",
                       ifelse(str_detect(filename,"NOPHYLOSIGNAL"),parameters[["datafile"]]<-"~/Dropbox/WAB Dissertation/Chapter 2 - Methods/dataForSims_fixedS_NOPHYLOSIGNAL.txt","NA"))))
  
  #identify whether or not crossvalidation is needed
  parameters[["doCrossValidataion"]]<-ifelse(str_detect(filename,"CROSSVALIDATED"),TRUE,FALSE)
  
  
  #identify weather or not to control for body size
  parameters[["doCorrectForBodySize"]]<-ifelse(str_detect(filename,"DontCorrect"),FALSE,
                                               ifelse(str_detect(filename,"Correct"),TRUE,NA))
  
  parameters[["randomizeHabs"]]<-ifelse(str_detect(filename,"RANDOMIZEDHABITATS"),TRUE,FALSE)
  
  return(parameters)
}

simulateData<-function(r,s) {
  #r is correlation between body mass and x var
  #s is scaling coef between 0 and 1 indicating relative strength of function versus phylo s=1 = all phylo
  habitatMultipliers<-seq(from=-var(log(habs$bmKilos))/5,to=var(log(habs$bmKilos))/5,length.out=4)[habs$Hab]
  as.numeric(r*log(habs$bmKilos)+s*sqrt(1-r^2)*fastBM(nsim=1,myTree,sig2=mean(pic(log(habs$bmKilos),multi2di(myTree))^2))+(1-s)*habitatMultipliers)
}

doDFA<-function(dataframe,correctForBodySize=TRUE,crossValidate=TRUE){
  #function takes a dataframe, assumes any vars with X in the name are the important variables
  varNames<-colnames(dataframe)[grep("X",colnames(dataframe))]
  if(correctForBodySize == TRUE){#correct by geomean if flagged to do so
    dataframe$geomean<-geomean(dataframe,varNames)
    dataframe[,varNames]<-dataframe[,varNames]/dataframe$geomean
  }
  
  if(crossValidate==TRUE){model<-lda(formula(paste("Hab ~ ", paste(varNames,collapse=" + "))),data= dataframe,CV=TRUE)
                          preds<-model$class
                          regPerCorr<-sum(preds==habs$Hab)/nrow(habs)#percent correct
  }
  if(crossValidate==FALSE){model<-lda(formula(paste("Hab ~ ", paste(varNames,collapse=" + "))),data= dataframe)
                           preds<-predict(model, dataframe[, varNames])
                           regPerCorr<-sum(preds$class==habs$Hab)/nrow(habs)#percent correct
                           
  }
  
  WilksPval<-summary(manova(as.matrix(dataframe[,varNames])~dataframe$Hab),test="Wilks")$stats[1,6]
  PillaiPval<-summary(manova(as.matrix(dataframe[,varNames])~dataframe$Hab),test="Pillai")$stats[1,6]
  HotellingPval<-summary(manova(as.matrix(dataframe[,varNames])~dataframe$Hab),test="Hotelling")$stats[1,6]
  RoyPval<-summary(manova(as.matrix(dataframe[,varNames])~dataframe$Hab),test="Roy")$stats[1,6]
  return(data.frame(regPerCorr=regPerCorr,WilksPval=WilksPval,PillaiPval=PillaiPval,HotellingPval=HotellingPval,RoyPval=RoyPval))
}

doPGLS<-function(dataframe,correctForBodySize=TRUE){
  varNames<-colnames(dataframe)[grep("X",colnames(dataframe))]
  geostring<-""
  if(correctForBodySize == TRUE){#correct by geomean if flagged to do so
    dataframe$geomean<-geomean(dataframe,varNames)
    dataframe[,varNames]<-dataframe[,varNames]/dataframe$geomean
    #geostring<-"+ log(geomean)"
  }
  comp<-comparative.data(myTree,dataframe,names.col='taxon')
  lapply(varNames,FUN=function(x){pgls(formula(paste(x, " ~ Hab ",geostring,sep="")),data=comp,lambda="ML")})
}


###########end function defs##################
###########end function defs##################
###########end function defs##################

#do some prepatory work like reading database for habitats etc
#that will get reused for each analysis category

read.tree("~/Dropbox/Public/ruminants.phy")->myTree
#get habitat data from my database

SQLstring = "SELECT DISTINCT PantheriaBodyMass.bodyMassGrams / 1000.0 as bmKilos, 
  Habitat.habitat as Hab, 
  taxonomy.Fernandez_Vrba_2005_Name, 
  taxonomy.BinindaEmonds_2008_Name 
  
  FROM
  ((PantheriaBodyMass LEFT OUTER JOIN taxonomy ON PantheriaBodyMass.taxon_id = taxonomy.id)
  LEFT OUTER JOIN Habitat ON Habitat.taxon_id = taxonomy.id);"

###get hab results from the database
habs<-read.table(pipe(paste("sqlite3 -header ~/Dropbox/bovidecomorph/bovidecomorph.db", "'",SQLstring,"'")),header=T,sep="|")
rm("SQLstring")
myTree<-drop.tip(myTree,myTree$tip.labe[!myTree$tip.labe %in% habs$Fernandez_Vrba_2005_Name])


############PLOTS################
############PLOTS################
############PLOTS################
#   #plot the phylogeny used
#   row.names(habs)<-habs$Fernandez_Vrba_2005_Name
#   habsForPlot<-habs[myTree$tip.labe,]
#   
#   colz<-c("grey20","grey45","grey70","grey90")
#   lineColz<-rep(NA,length(myTree$edge))
#   lineColz[myTree$edge[,2]<=nrow(habsForPlot)]<-colz[habsForPlot$Hab] #update lineColz vector only for tip edges (which by definition connect to a node number <= nrow(habsForPlot))
#   lineColz[is.na(lineColz)]<-1
#   
#   ltyz=c(1,3,3,1)
#   lineTypes<-rep(NA,length(myTree$edge))
#   lineTypes[myTree$edge[,2]<=nrow(habsForPlot)]<-ltyz[habsForPlot$Hab] #update ltyz vector only for tip edges (which by definition connect to a node number <= nrow(habsForPlot))
#   lineTypes[is.na(lineTypes)]<-1

#uncomment to save as tiff
#tiff("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/phylogenyUSED_4habitats.tiff",width=10,height=10,units="in",res=600)
#par(cex=1)
#plot(myTree,edge.width=5,edge.col=lineColz,cex=1.1,show.tip.label=T,edge.lty=lineTypes,label.offset=.15)
#par(xpd=NA)#setting xpd=NA shuts off truncation in the margins
#legend(3,-2,horiz=T,legend=levels(habs$Hab),lwd=5,lty=ltyz,cex=1,seg.len=2.5,col=colz,title.col=1,title="Habitat")
#dev.off()
#end phylogeny plot
############END PLOTS################
############END PLOTS################
############END PLOTS################



############ANALYSIS################
############ANALYSIS################
############ANALYSIS################


#load the list of all the desired filenames
filesToDo<-as.character(read.table("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/simFilesToWrite.txt",header=F)[,1])
#filesToDo<-filesToDo[c(1:2,18,20:21)]

fileCounter<-0
# loop through each file and do the analysis
lapply(filesToDo,FUN=function(theFile)
  {
  fileCounter<-fileCounter + 1
  parameters<-getSimulationParametersFromFileName(theFile)
  
  myData<-source(parameters[["datafile"]])[[1]]
  
  # #code to resimulate chars if need be
  # myData<-lapply(1:10000,FUN=function(x){
  #   r<-0
  #   s<-1
  #   measurements<-simulateData(r,s)
  #   list(r=r,s=s,measurements=measurements)
  # })
  # dump("myData","/Users/andrewbarr/Dropbox/WAB Dissertation/Chapter 2 - Methods/dataForSims_fixedS_0R.txt")
  
  
  nSims<-2000
  
  results<-lapply(1:nSims,FUN=function(DFAcounter){
    
    vars<-sample(x=1:length(myData),size=sample(7:12,1),replace=FALSE)
    
    #make a dataframe from these indices
    myDataFrame<-data.frame(lapply(vars,FUN=function(x){
      myData[[x]]$measurements
    }))
    names(myDataFrame)<-paste("X",1:length(vars),sep="")
    
    #if there are negative numbers, add in an arithmetic scaling factor
    #because negative numbers will break the geomean calculation
    if(sum(as.vector(myDataFrame)<0)>0){myDataFrame<-myDataFrame + 100}
    
    #add in Hab and taxon columns and randomize habitats if called for
    randomizeHabs<-parameters$randomizeHabs
    ifelse(randomizeHabs,
           yes=myDataFrame<-data.frame(myDataFrame,taxon=habs$Fernandez_Vrba_2005_Name,Hab=sample(levels(habs$Hab),nrow(myDataFrame),replace=TRUE)),
           no=myDataFrame<-data.frame(myDataFrame,taxon=habs$Fernandez_Vrba_2005_Name,Hab=habs$Hab))
    
    
    
    #extract the r and s value for each variable
    longResults<-rbind.fill(lapply(vars,FUN=function(x){
      data.frame(r=myData[[x]]$r,s=myData[[x]]$s,measurementID=x)
    }))
    doCrossValidataion<-parameters[["doCrossValidataion"]]
    doCorrectForBodySize<-parameters[["doCorrectForBodySize"]]
    
    #do the DFA
    df<-tryCatch(doDFA(myDataFrame,correctForBodySize=doCorrectForBodySize,crossValidate=doCrossValidataion),error=function(e) return(NA))
    
    if(doCrossValidataion==FALSE){#only do this stuff if we are not crossvalidating.  There is no need to re-run the PGLS stuff if we just need the cross-val DFA stuff 
      if(sum(is.na(df))==0){
        
        pglss<-tryCatch(doPGLS(myDataFrame,correctForBodySize=doCorrectForBodySize),error=function(e) return(rep(NA,length(vars))))
        
        if(sum(is.na(pglss))==0){
          monotonic<-unlist(lapply(pglss,FUN=function(y){paste(order(summary(y)$coef[2:4,1]),collapse="") %in% c("321","123")}))
          longResults$pvalues<-unlist(lapply(pglss,FUN=function(y){summary(y)$coef[4,4]}))
          longResults$overallPvalues<-unlist(lapply(pglss,FUN=function(z){pf(summary(z)$fstatistic[1], summary(z)$fstatistic[2], summary(z)$fstatistic[3], lower.tail = FALSE)}))
          longResults$monotonic<-monotonic
          
        }
      }
    }
    
    if(doCrossValidataion==TRUE){#put in placeholders for PGLS stuff if we are skipping PGLS
      longResults$pvalues<-rep(-999,nrow(longResults))
      longResults$overallPvalues<-rep(-999,nrow(longResults))
      longResults$monotonic<-rep(FALSE,nrow(longResults))
    }
    
    longResults$dfaSuccessRate<-rep(df$regPerCorr,nrow(longResults))
    longResults$wilkesLambda<-rep(df$Wilks,nrow(longResults))        
    longResults$PillaiPval<-rep(df$Pillai,nrow(longResults))
    longResults$HotellingPval<-rep(df$Hotell,nrow(longResults))
    longResults$RoyPval<-rep(df$Roy,nrow(longResults))
    longResults$dfaID<-DFAcounter
    longResults$nvars<-rep(nrow(longResults),nrow(longResults))
    #print(paste(DFAcounter," successful iterations in file # ",fileCounter,"/",length(filesToDo)))
    return(longResults)
  })
  results<-rbind.fill(results)
  write.table(results,parameters[["filename"]],sep="\t",row.names=FALSE)
  #print(paste0("Finished file # ",fileCounter,"/",length(filesToDo)))
})
############END ANALYSIS################
############END ANALYSIS################
############END ANALYSIS################