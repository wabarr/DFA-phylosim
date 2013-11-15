require(reshape2)
require(ggplot2)
require(plyr)
require(stringr)
require(scales)
require(xtable)
require(gridExtra)
require(parallel)

#environmental variable for number of cores to use when doing parallel
Sys.setenv(MC_CORES=4)

nice<-function(x,places=2){round(x,places)}

#theme_set(theme_bw(20) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))

setwd("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/")
inputFiles<-list.files(path=".",pattern="BrownianMotionSimResults")



#whole analysis in a function
doTheWholeShebang<-function(theFile){
  allRESULTS<-list()
  rez<-lapply(inputFiles[[theFile]],FUN=function(x){read.table(x,header=T,sep="\t")})
  #print(inputFiles[[theFile]])
  
  #parse the file names to add a column called simGroup
  lapply(1:length(rez),FUN=function(x){
    rez[[x]]$simGroup<<-rep(gsub(pattern=".txt",replacement="",paste0(str_split(string=inputFiles[[theFile]],pattern="_")[[1]][3:4],collapse="-")),nrow(rez[[x]]))
    rez[[x]]$isCrossvalidated<<-length(grep("CROSSVALIDATED",paste0(str_split(string=inputFiles[[theFile]],pattern="_")[[1]])))>0
    rez[[x]]$isRandomizedHabs<<-length(grep("RANDOMIZEDHABITATS",paste0(str_split(string=inputFiles[[theFile]],pattern="_")[[1]])))>0
  })
  
  rez<-do.call(rbind,rez)
  #get rid of NAs
  rez<-subset(rez,!is.na(dfaSuccessRate))
  
  #create a new binning variable for nvars with only three levels, for creating cleaner plots
  #rez$binnedNvars<-factor(cut(rez$nvars,3))->recoded
  #levels(rez$binnedNvars)<-c("7-8 variables","8-9 variables","10-11 variables")
  
  #calculate the average s value for each simulation, add to the dataframe then delete the temp variable
  #TEMPaverageS<-data.frame(averageS=tapply(rez$s,rez$dfaID,FUN=mean),dfaID=as.numeric(names(tapply(rez$s,rez$dfaID,FUN=mean))))
  #rez<-merge(TEMPaverageS,rez)
  #rm(TEMPaverageS)
  
  
  #calculate the average r value for each simulation, add to the dataframe then delete the temp variable
  #TEMPaverageR<-data.frame(averageR=tapply(rez$r,rez$dfaID,FUN=mean),dfaID=as.numeric(names(tapply(rez$r,rez$dfaID,FUN=mean))))
  #rez<-merge(TEMPaverageR,rez)
  #rm(TEMPaverageR)
  
  #calculate the average p value for each simulation, add to the dataframe then delete the temp variable
  #TEMPaverageP<-data.frame(averageP=tapply(rez$pvalues,rez$dfaID,FUN=mean),dfaID=as.numeric(names(tapply(rez$pvalues,rez$dfaID,FUN=mean))))
  #rez<-merge(TEMPaverageP,rez)
  #rm(TEMPaverageP)
  
  #rez$sCategory<-cut(rez$s,4)
  #rez$rCategory<-cut(rez$r,4)
  
  #create column that uniquely identifies each DFA~simGroup combination
  rez$dfaID_simGroup<-paste(rez$dfaID,rez$simGroup)
  
  rez$bonnferroniP<-unlist(tapply(rez$overall,INDEX=rez$dfaID_simGroup, FUN=function(x) p.adjust(p=x,method="bonferroni")))
  rez$HolmBonnferroniP<-unlist(tapply(rez$overall,INDEX=rez$dfaID_simGroup, FUN=function(x) p.adjust(p=x,method="holm")))
  rez$fdrP<-unlist(tapply(rez$overall,INDEX=rez$dfaID_simGroup, FUN=function(x) p.adjust(p=x,method="BH")))
  #logical vector of PGLSs that are significant overall and at Bonneferronin corrected post-hoc level
  rez$sigs<-(rez$overall<.05)# + (rez$p<.0083333333333)==2
  rez$sigsBonferroni<-(rez$bonnferroniP<.05)
  rez$sigsHolmes<-(rez$HolmBonnferroni<.05)
  rez$sigFDR<-(rez$fdrP<.05)
  
  #create a shortForm dataframe that includes only columns with apply to a single dfa subset
  shortFormVarnames<-c("dfaID_simGroup","dfaSuccessRate","nvars","simGroup","wilkesLambda","PillaiPval","HotellingPval","RoyPval","isCrossvalidated","isRandomizedHabs")
  shortForm<-ddply(rez[,shortFormVarnames], "dfaID_simGroup", .fun=function(x) x[1,])
  
  splitSimGroup<-str_split(shortForm$simGroup[1:2],"-")
  shortForm$BM_correlation<-factor(unlist(lapply(str_split(shortForm$simGroup,"-"),FUN=function(x) x[1])),levels=c("NOPHYLOSIGNAL","0R","lowR","highR"),ordered=TRUE)
  levels(shortForm$BM_correlation)<-c("r = 0 no phylo","r = 0 phylo","low r","high r")
  
  
  shortForm$BMCorrectionOrNot<-factor(unlist(lapply(str_split(shortForm$simGroup,"-"),FUN=function(x) x[2])))
  levels(shortForm$BMCorrectionOrNot)<-c("Size Corrected","Raw Data")
  #shortForm$sCategory<-cut(shortForm$averageS,3)
  #levels(shortForm$sCategory)<-c("low S","medium S","high S")
  #shortForm$wilkes<-tapply(rez$wilkes,rez$dfaID_simGroup,FUN=unique)
  #shortForm$pillai<-tapply(rez$Pillai,rez$dfaID_simGroup,FUN=unique)
  #shortForm$hotelling<-tapply(rez$Hotelling,rez$dfaID_simGroup,FUN=unique)
  #shortForm$roy<-tapply(rez$Roy,rez$dfaID_simGroup,FUN=unique)
  #annotateDF<-ddply(rez,"binnedNvars",.fun=function(x){paste("Mean = ",round(mean(x$dfaSuccessRate),2)*100,"%",sep="")})
  #annotateDF2<-ddply(rez,"binnedNvars",.fun=function(x){paste("95th percentile = ",round(quantile(x$dfaSuccessRate,.95),2)*100,"%",sep="")})
  
  
  simGroupSummaries_all<-ddply(shortForm,"simGroup",.fun=function(x){
    data.frame(
      #BMCorrectionOrNot=unique(x$BMCorrectionOrNot),
      isCrossvalidated=unique(x$isCrossvalidated),
      isRandomizedHabs=unique(x$isRandomizedHabs),
      BM_correlation=unique(x$BM_correlation),
      proportionOfDFAsSignificant=sum(x$wilkesLambda<.05)/nrow(x),
      meanDFASuccessRate=mean(x$dfaSuccessRate*100,na.rm=TRUE),
      sdDFASuccessRate=sd(x$dfaSuccessRate*100,na.rm=TRUE)
    )
  })
  
  
  #   simGroupSummaries<-ddply(subset(shortForm,wilkesLambda<.05),"simGroup",.fun=function(x){
  #     data.frame(
  #       #BMCorrectionOrNot=unique(x$BMCorrectionOrNot),
  #       BM_correlation=unique(x$BM_correlation),
  #       isCrossvalidated=unique(x$isCrossvalidated),
  #       isRandomizedHabs=unique(x$isRandomizedHabs),
  #       proportionOfDFAsSignificant=sum(x$wilkesLambda<.05)/nrow(x),
  #       meanDFASuccessRate=mean(x$dfaSuccessRate*100,na.rm=TRUE),
  #       sdDFASuccessRate=sd(x$dfaSuccessRate*100,na.rm=TRUE)
  #     )
  #   })
  #   
  countSigPGLS<-as.data.frame(tapply(rez$sigs,INDEX=rez$simGroup,FUN=function(x) sum(x,na.rm=T) / (length(x) - sum(is.na(x)))))
  countSigPGLSBonferroni<-as.data.frame(tapply(rez$sigsBonferroni,INDEX=rez$simGroup,FUN=function(x) sum(x,na.rm=T) / (length(x) - sum(is.na(x)))))
  countSigPGLSHolm<-as.data.frame(tapply(rez$sigsHolmes,INDEX=rez$simGroup,FUN=function(x) sum(x,na.rm=T) / (length(x) - sum(is.na(x)))))
  countSigPGLSfdr<-as.data.frame(tapply(rez$sigFDR,INDEX=rez$simGroup,FUN=function(x) sum(x,na.rm=T) / (length(x) - sum(is.na(x)))))
  
  names(countSigPGLSBonferroni)<-"countSigPGLSBonferroni"
  names(countSigPGLS)<-"countSigPGLS"
  names(countSigPGLSHolm)<-"countSigPGLSHolm"
  names(countSigPGLSfdr)<-"countSigPGLSfdr"
  
  simGroupSummaries_all$countSigPGLS<-countSigPGLS
  simGroupSummaries_all$countSigPGLSBonferroni<-countSigPGLSBonferroni
  simGroupSummaries_all$countSigPGLSHolm<-countSigPGLSHolm
  simGroupSummaries_all$countSigPGLSfdr<-countSigPGLSfdr
  #print(xtable(simGroupSummaries_all),type="html")
  allRESULTS[[theFile]]<-simGroupSummaries_all
  # #count of how many significant chars per subset
  # countTable<-table(tapply(rez$sigs,rez$dfaID,FUN=sum))
  # countSigPerSample<-qplot(y=as.numeric(countTable),x=1:length(countTable)-1,size=I(5),xlab="Number of Significant Characters",ylab="Number of Samples") + theme_bw(20)
  # #ggsave("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/countSigPerSample.pdf",countSigPerSample)
  
  
  #splitSimGroup<-str_split(rez$simGroup[1:2],"-")
  #rez$BM_correlation<-factor(unlist(mclapply(str_split(rez$simGroup,"-"),FUN=function(x) x[1])),levels=c("0R","lowR","highR"),ordered=TRUE)
  #levels(rez$BM_correlation)<-c("r = 0","low r","high r")
  #distributionR<- qplot(x=r,data=rez,fill=I('grey'),color=I("black")) + facet_grid(facets=BM_correlation~.) + ylab("proportion") + scale_y_continuous(breaks=c(0,37000/2,37000),labels=c("0","0.5","1"))
  # ggsave("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/distributionRAcrossCharacterSets.pdf",distributionR)
  
  
  #linear model investigating which factors influence the success rate the most
  #mod1<-lm(dfaSuccessRate ~ nvars , data=shortForm)
  #mod2<-lm(dfaSuccessRate ~ nvars + BM_correlation , data=shortForm)
  #mod3<-lm(dfaSuccessRate ~ nvars + BM_correlation + BMCorrectionOrNot, data=shortForm)
  #mod4<-lm(dfaSuccessRate ~ nvars * BM_correlation * BMCorrectionOrNot, data=shortForm)
  
  #anova(mod3,mod4)
  
  #anova(lm(dfaSuccessRate ~ nvars + BM_correlation + BMCorrectionOrNot + nvars:BM_correlation + BM_correlation:BMCorrectionOrNot , data=shortForm))
  
  #write.table(x=shortForm,file="~/Dropbox/WAB Dissertation/Chapter 2 - Methods/shortForm.txt",sep="\t",row.names=FALSE)
  
  #change these values from dataframes to numerics
  allRESULTS<-lapply(allRESULTS,FUN=function(df){
    df$countSigPGLS<-as.numeric(df$countSigPGLS)
    df$countSigPGLSBonferroni<-as.numeric(df$countSigPGLSBonferroni)
    df$countSigPGLSHolm<-as.numeric(df$countSigPGLSHolm)
    df$countSigPGLSfdr<-as.numeric(df$countSigPGLSfdr)
    return(df)
  })
  
  #bind it up
  allRESULTS<-do.call(rbind,allRESULTS)
  allRESULTS$isCrossvalidated<-as.character(allRESULTS$isCrossvalidated)
  #toWrite<-print(xtable(allRESULTS,digits=4,),type="html")
  #writeLines(toWrite,con=paste0("~/Desktop/",str_split(inputFiles[theFile],".txt")[[1]][1],".html"))
  write.table(x=allRESULTS,file=paste0("~/Desktop/",str_split(inputFiles[theFile],".txt")[[1]][1],".txt"),row.names=FALSE)
  
}

#debug(doTheWholeShebang)
#doTheWholeShebang(2)
mclapply(1:length(inputFiles),mc.cores=4,FUN=doTheWholeShebang)

summarizeResults<-TRUE
if(summarizeResults) {
  temporaryFILES<-list.files("~/Desktop/","BrownianMotionSim",full.names=TRUE)
  toWrite<-mclapply(temporaryFILES,FUN=function(x) read.table(x,header=TRUE))
  toWrite<-do.call(rbind,toWrite)
  write.table(toWrite,file="SimulationResults.txt",row.names=FALSE,sep="\t")
  #deletes all the temporary files
  lapply(temporaryFILES,unlink)
}

rez<-read.table("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/SimulationResults.txt",header=TRUE,sep="\t")
####TABLE 3
ddply(rez,
      .variables=.(simGroup,isRandomizedHabs),
      .fun=summarize,
          ##additional args to summarize 
          Percent_DFA_significant=round(mean(proportionOfDFAsSignificant)*100,2),
          IsCrossValidated = paste(isCrossvalidated,collapse=","),
          mean_percent_correct_classifications=paste(round(meanDFASuccessRate,2),collapse=","),
          PGLS_Type_I_error_rate=round(countSigPGLS[2]*100,2),
          PGLS_Type_I_error_rate_FDR=round(countSigPGLSfdr[2]*100,2)
      )


####################################PLOTS################################
####################################PLOTS################################
####################################PLOTS################################

# theFiles<-list.files("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/","BrownianMotionSim",full.names=TRUE)
# theData<-mclapply(theFiles,function(x) {
#                                         TEMP<-read.table(x,header=TRUE,sep="\t")
#                                         TEMP$simGroup<-rep(gsub(pattern=".txt",replacement="",paste0(str_split(string=x,pattern="_")[[1]][3:4],collapse="-")),nrow(TEMP))
#                                         TEMP$isCrossvalidated<-length(grep("CROSSVALIDATED",paste0(str_split(string=x,pattern="_")[[1]])))>0
#                                         TEMP$isRandomizedHabs<-length(grep("RANDOMIZEDHABITATS",paste0(str_split(string=x,pattern="_")[[1]])))>0
#                                         TEMP$BMCorrectionOrNot<-str_split(TEMP$simGroup[1],pattern="-")[[1]][2]
#                                         TEMP$BM_correlation<-str_split(TEMP$simGroup[1],pattern="-")[[1]][1]
#                                         
#                                         if(str_split(TEMP$simGroup[1],pattern="-")[[1]][2]=="CorrectBodySize" && TEMP$isRandomizedHabs==FALSE) {return(TEMP)}
#                                         })
# 
# theData<-do.call(rbind,theData)
# theData$isCrossvalidated<-factor(theData$isCrossvalidated)
# levels(theData$isCrossvalidated)<-c("Resubstitution","Crossvalidation")
# 
# theData$BM_correlation<-factor(theData$BM_correlation)
# levels(theData$BM_correlation)<-c("r = 0 (phylo)  ","high r   ","low r   ","r = 0 (no phylo)")
# 
# 
# boxplot_SuccessByR_all<-qplot(nvars,dfaSuccessRate*100,data=subset(theData,BMCorrectionOrNot="CorrectBodySize"),geom="boxplot",fill=isCrossvalidated,group=interaction(nvars,isCrossvalidated),xlab="# of characters",ylab="Classification Success Rate (%)",facets=~BM_correlation,outlier.size=0) + 
#   scale_fill_grey(start = 0.9, end = 0.5) + 
#   labs(fill="Validation",title="DFA success by number of characters\nwith size-corrected data, actual habitats") + 
#   guides(fill = guide_legend(keywidth = 1, keyheight = 3)) + 
#   theme_bw(20) + 
#   theme(legend.position="bottom",panel.grid.minor=element_blank(),panel.grid.major=element_line(size=.7)) +
#   geom_hline(y=25,lty=2)
# 
# setEPS()
# postscript("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/boxplotsSuccessByR.eps",width=14,height=9.15)
# print(boxplot_SuccessByR_all)
# dev.off()
#   boxplot_SuccessByR<-qplot(nvars,dfaSuccessRate,data=subset(shortForm, wilkes<.05),geom="boxplot",group=interaction(nvars,BMCorrectionOrNot),xlab="# of characters",ylab="Classification Success Rate",fill=BMCorrectionOrNot,facets=~BM_correlation,outlier.size=0) + 
#     scale_y_continuous(labels=percent) + 
#     scale_fill_grey(start = 0.9, end = 0.5) + 
#     labs(fill="Size Correction") + 
#     guides(fill = guide_legend(keywidth = 1, keyheight = 3)) + 
#     theme_bw(20) + 
#     theme(legend.position="bottom",panel.grid.minor=element_blank(),panel.grid.major=element_line(size=.7))
#   #ggsave("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/boxplotsSuccessByR.pdf",boxplot_SuccessByR,height=9.15,width=14,units="in")
#   


# densityDFA_all<-qplot(x=dfaSuccessRate*100,data=shortForm,lty=BM_correlation,geom="density",alpha=I(.35),fill=BM_correlation, facets=BMCorrectionOrNot~.,main="Fig 2A - All DFAs",xlab="Success Rate") + 
#   scale_linetype_manual(values=c("dotted", "solid","longdash","dotdash")) +
#   scale_x_continuous(limits=c(30,90)) + 
#   theme(legend.title = element_blank())
# ggsave(densityDFA_all,"~/Dropbox/WAB Dissertation/Chapter 2 - Methods/densityDFA.pdf",width=14,height=11.15)
#   densityDFA<-qplot(x=dfaSuccessRate*100,data=subset(shortForm,wilkes<.05),lty=BM_correlation,geom="density",alpha=I(.35),fill=BM_correlation, facets=BMCorrectionOrNot~.,main="Fig 2B - Only Significant DFAs",xlab="Success Rate") + 
#     scale_linetype_manual(values=c("dotted", "solid","longdash","dotdash")) +
#     scale_x_continuous(limits=c(30,90)) + 
#     theme(legend.title = element_blank())
#   pdf(file="~/Dropbox/WAB Dissertation/Chapter 2 - Methods/densityDFA.pdf",width=14,height=11.15)
#   grid.arrange(densityDFA_all,densityDFA,ncol=1)
#   dev.off()
#   

rez$Dataset<-sapply(str_split(string=rez$simGroup,pattern="-"), FUN=function(x) x[[1]])
rez$Dataset<-factor(rez$Dataset)
levels(rez$Dataset)<-c("r = 0 (phylo)  ","high r   ","low r   ","r = 0 (no phylo)")

rez$isCrossvalidated<-factor(rez$isCrossvalidated)
levels(rez$isCrossvalidated)<-c("Resubstitution","Crossvalidation")

rez$isRandomizedHabs<-factor(rez$isRandomizedHabs)
levels(rez$isRandomizedHabs)<-c("Actual Habitats","Randomized Habitats")

DFAresults<-ggplot(data=rez,aes(x=proportionOfDFAsSignificant*100,y=meanDFASuccessRate,color=Dataset,shape=Dataset)) + 
    geom_point(size=7) + 
    facet_grid(facets=isRandomizedHabs~isCrossvalidated) + 
    theme_bw(20) + 
    scale_color_grey(start=.8,end=.2) + 
    scale_shape_manual(values=16:19) +
    labs(x="% of DFAs significant",y="Mean Success Rate (%)",title="Summary of DFA Results") + 
    theme(legend.position="bottom")
    ggsave("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/DFA_results_points.eps",DFAresults,width=7,height=7,units="in")
    
####################################ENDPLOTS################################
####################################ENDPLOTS################################
####################################ENDPLOTS################################

##Example of how to test for PhyloSignal in categorical variable
require(geiger)

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
habitatDF <- read.table(pipe(paste("sqlite3 -header ~/Dropbox/bovidecomorph/bovidecomorph.db", "'",SQLstring,"'")),header=T,sep="|")
habs <- habitatDF$Hab
names(habs) <- habitatDF$Fernandez_Vrba_2005_Name

myTree<-drop.tip(myTree,myTree$tip.labe[!myTree$tip.labe %in% habitatDF$Fernandez_Vrba_2005_Name])

zeroTree <- transform(myTree,"lambda",0)


####with Actual hab - single parameter model
lambdaActualData<-fitDiscrete(multi2di(myTree),habs, model="ER")
lambdaNoPhyloSignal<-fitDiscrete(multi2di(zeroTree),habs, model="ER")
pchisq(2 * (lambdaActualData$opt$lnL - lambdaNoPhyloSignal$opt$lnL),1,lower.tail=FALSE)

#####with random habs
randomHabs <- sample(unique(habs),length(habs),replace=TRUE)
names(randomHabs) <- names(habs)
lambdaActualData<-fitDiscrete(multi2di(myTree),randomHabs, model="ER")
lambdaNoPhyloSignal<-fitDiscrete(multi2di(zeroTree),randomHabs, model="ER")
pchisq(2 * (lambdaActualData$opt$lnL - lambdaNoPhyloSignal$opt$lnL),1,lower.tail=FALSE)


####with Actual hab - symetric model
ALTlambdaActualData<-fitDiscrete(multi2di(myTree),habs, model="SYM")
ALTlambdaNoPhyloSignal<-fitDiscrete(multi2di(zeroTree),habs,model="SYM")
pchisq(2 * (ALTlambdaActualData$opt$lnL - ALTlambdaNoPhyloSignal$opt$lnL),1,lower.tail=FALSE)

####with Actual hab - meristic model
ALT2lambdaActualData<-fitDiscrete(multi2di(myTree),habs,model="meristic")
ALT2lambdaNoPhyloSignal<-fitDiscrete(multi2di(zeroTree),habs,model="meristic")
pchisq(2 * (ALT2lambdaActualData$opt$lnL - ALT2lambdaNoPhyloSignal$opt$lnL),1,lower.tail=FALSE)