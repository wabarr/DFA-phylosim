
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