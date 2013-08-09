require(reshape2)
require(ggplot2)
require(plyr)
require(stringr)
require(scales)
require(xtable)
require(gridExtra)
require(parallel)

nice<-function(x,places=2){round(x,places)}

theme_set(theme_bw(20) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))

setwd("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/")
inputFiles<-list.files(path=".",pattern="RANDOMIZEDHABITATS")


allRESULTS<-list()

lapply(1:length(inputFiles),FUN=function(theFile){
  rez<-lapply(inputFiles[theFile],FUN=function(x){read.table(x,header=T,sep="\t")})
  
  #parse the file names to add a column called simGroup
  lapply(1:length(rez),FUN=function(x){
    rez[[x]]$simGroup<<-rep(gsub(pattern=".txt",replacement="",paste0(str_split(string=inputFiles[theFile],pattern="_")[[1]][3:4],collapse="-")),nrow(rez[[x]]))
    rez[[x]]$isCrossvalidated<<-length(grep("CROSSVALIDATED",paste0(str_split(string=inputFiles[theFile],pattern="_")[[1]])))>0
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
  shortFormVarnames<-c("dfaID_simGroup","dfaSuccessRate","nvars","simGroup","wilkesLambda","PillaiPval","HotellingPval","RoyPval","isCrossvalidated")
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
  
  ####################################PLOTS################################
  ####################################PLOTS################################
  ####################################PLOTS################################
  # boxplot_SuccessByR_all<-qplot(nvars,dfaSuccessRate,data=shortForm,geom="boxplot",group=interaction(nvars,BMCorrectionOrNot),xlab="# of characters",ylab="Classification Success Rate",fill=BMCorrectionOrNot,facets=~BM_correlation,outlier.size=0) + 
  #   scale_y_continuous(labels=percent) + 
  #   scale_fill_grey(start = 0.9, end = 0.5) + 
  #   labs(fill="Size Correction") + 
  #   guides(fill = guide_legend(keywidth = 1, keyheight = 3)) + 
  #   theme_bw(20) + 
  #   theme(legend.position="bottom",panel.grid.minor=element_blank(),panel.grid.major=element_line(size=.7))
  # 
  # 
  # boxplot_SuccessByR<-qplot(nvars,dfaSuccessRate,data=subset(shortForm, wilkes>.05),geom="boxplot",group=interaction(nvars,BMCorrectionOrNot),xlab="# of characters",ylab="Classification Success Rate",fill=BMCorrectionOrNot,facets=~BM_correlation,outlier.size=0) + 
  #   scale_y_continuous(labels=percent) + 
  #   scale_fill_grey(start = 0.9, end = 0.5) + 
  #   labs(fill="Size Correction") + 
  #   guides(fill = guide_legend(keywidth = 1, keyheight = 3)) + 
  #   theme_bw(20) + 
  #   theme(legend.position="bottom",panel.grid.minor=element_blank(),panel.grid.major=element_line(size=.7))
  # #ggsave("~/Dropbox/WAB Dissertation/Chapter 2 - Methods/boxplotsSuccessByR.pdf",boxplot_SuccessByR,height=9.15,width=14,units="in")
  # 
  # 
  # 
  # densityDFA_all<-qplot(x=dfaSuccessRate*100,data=shortForm,lty=BM_correlation,geom="density",alpha=I(.35),fill=BM_correlation, facets=BMCorrectionOrNot~.,main="Fig 2A - All DFAs",xlab="Success Rate") + 
  #   scale_linetype_manual(values=c("dotted", "solid","longdash","dotdash")) +
  #   scale_x_continuous(limits=c(30,90)) + 
  #   theme(legend.title = element_blank())
  # 
  # densityDFA<-qplot(x=dfaSuccessRate*100,data=subset(shortForm,wilkes<.05),lty=BM_correlation,geom="density",alpha=I(.35),fill=BM_correlation, facets=BMCorrectionOrNot~.,main="Fig 2B - Only Significant DFAs",xlab="Success Rate") + 
  #   scale_linetype_manual(values=c("dotted", "solid","longdash","dotdash")) +
  #   scale_x_continuous(limits=c(30,90)) + 
  #   theme(legend.title = element_blank())
  # pdf(file="~/Dropbox/WAB Dissertation/Chapter 2 - Methods/densityDFA.pdf",width=14,height=11.15)
  # grid.arrange(densityDFA_all,densityDFA,ncol=1)
  # dev.off()
  
  ####################################ENDPLOTS################################
  ####################################ENDPLOTS################################
  ####################################ENDPLOTS################################
  
  simGroupSummaries_all<-ddply(shortForm,"simGroup",.fun=function(x){
    data.frame(
      #BMCorrectionOrNot=unique(x$BMCorrectionOrNot),
      isCrossvalidated=unique(x$isCrossvalidated),
      BM_correlation=unique(x$BM_correlation),
      proportionOfDFAsSignificant=sum(x$wilkesLambda<.05)/nrow(x),
      meanDFASuccessRate=mean(x$dfaSuccessRate*100,na.rm=TRUE),
      sdDFASuccessRate=sd(x$dfaSuccessRate*100,na.rm=TRUE)
    )
  })
  
  
  simGroupSummaries<-ddply(subset(shortForm,wilkesLambda<.05),"simGroup",.fun=function(x){
    data.frame(
      #BMCorrectionOrNot=unique(x$BMCorrectionOrNot),
      BM_correlation=unique(x$BM_correlation),
      isCrossvalidated=unique(x$isCrossvalidated),
      proportionOfDFAsSignificant=sum(x$wilkesLambda<.05)/nrow(x),
      meanDFASuccessRate=mean(x$dfaSuccessRate*100,na.rm=TRUE),
      sdDFASuccessRate=sd(x$dfaSuccessRate*100,na.rm=TRUE)
    )
  })
  
  countSigPGLS<-as.data.frame(tapply(rez$sigs,INDEX=rez$simGroup,FUN=function(x) sum(x)/length(x)))
  countSigPGLSBonferroni<-as.data.frame(tapply(rez$sigsBonferroni,INDEX=rez$simGroup,FUN=function(x) sum(x)/length(x)))
  countSigPGLSHolm<-as.data.frame(tapply(rez$sigsHolmes,INDEX=rez$simGroup,FUN=function(x) sum(x)/length(x)))
  countSigPGLSfdr<-as.data.frame(tapply(rez$sigFDR,INDEX=rez$simGroup,FUN=function(x) sum(x)/length(x)))
  
  names(countSigPGLSBonferroni)<-"countSigPGLSBonferroni"
  names(countSigPGLS)<-"countSigPGLS"
  names(countSigPGLSHolm)<-"countSigPGLSHolm"
  names(countSigPGLSfdr)<-"countSigPGLSfdr"
  
  simGroupSummaries_all$countSigPGLS<-countSigPGLS
  simGroupSummaries_all$countSigPGLSBonferroni<-countSigPGLSBonferroni
  simGroupSummaries_all$countSigPGLSHolm<-countSigPGLSHolm
  simGroupSummaries_all$countSigPGLSfdr<-countSigPGLSfdr
  #print(xtable(simGroupSummaries_all),type="html")
  allRESULTS[[theFile]]<<-simGroupSummaries_all
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
  
})

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
toWrite<-print(xtable(allRESULTS,digits=4,),type="html")
writeLines(toWrite,con="~/Desktop/randomizedResults.html")
