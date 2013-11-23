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


####with Actual hab - Equal Rates Model (ER)
ERActualTree<-fitDiscrete(multi2di(myTree),habs, model="ER")
ERStarPhylogeny<-fitDiscrete(multi2di(zeroTree),habs, model="ER")
pchisq(2 * (ERActualTree$opt$lnL - ERStarPhylogeny$opt$lnL),1,lower.tail=FALSE)

#####with random habs  - Equal Rates Model (ER)
randomHabs <- sample(unique(habs),length(habs),replace=TRUE)
names(randomHabs) <- names(habs)
ERActualTree<-fitDiscrete(multi2di(myTree),randomHabs, model="ER")
ERStarPhylogeny<-fitDiscrete(multi2di(zeroTree),randomHabs, model="ER")
pchisq(2 * (ERActualTree$opt$lnL - ERStarPhylogeny$opt$lnL),1,lower.tail=FALSE)


####with Actual hab  - symetric rates model (SYM)
SYMActualTree<-fitDiscrete(multi2di(myTree),habs, model="SYM")
SYMStarPhylogeny<-fitDiscrete(multi2di(zeroTree),habs,model="SYM")
pchisq(2 * (SYMActualTree$opt$lnL - SYMStarPhylogeny$opt$lnL),1,lower.tail=FALSE)

#####with random habs  - symetric rates model (SYM)
randomHabs <- sample(unique(habs),length(habs),replace=TRUE)
names(randomHabs) <- names(habs)
SYMActualTree<-fitDiscrete(multi2di(myTree),randomHabs, model="SYM")
SYMStarPhylogeny<-fitDiscrete(multi2di(zeroTree),randomHabs,model="SYM")
pchisq(2 * (SYMActualTree$opt$lnL - SYMStarPhylogeny$opt$lnL),1,lower.tail=FALSE)

####with Actual hab - all rates different model (ARD)
ARDActualTree<-fitDiscrete(multi2di(myTree),habs,model="meristic")
ARDStarPhylogeny<-fitDiscrete(multi2di(zeroTree),habs,model="meristic")
pchisq(2 * (ARDActualTree$opt$lnL - ARDStarPhylogeny$opt$lnL),1,lower.tail=FALSE)

####with random habs - all rates different model (ARD)
randomHabs <- sample(unique(habs),length(habs),replace=TRUE)
names(randomHabs) <- names(habs)
ARDActualTree<-fitDiscrete(multi2di(myTree),randomHabs,model="meristic")
ARDStarPhylogeny<-fitDiscrete(multi2di(zeroTree),randomHabs,model="meristic")
pchisq(2 * (ARDActualTree$opt$lnL - ARDStarPhylogeny$opt$lnL),1,lower.tail=FALSE)

####with Actual hab - meristic model (ordered states)
MeristicActualTree<-fitDiscrete(multi2di(myTree),habs,model="meristic")
MeristicStarPhylogeny<-fitDiscrete(multi2di(zeroTree),habs,model="meristic")
pchisq(2 * (MeristicActualTree$opt$lnL - MeristicStarPhylogeny$opt$lnL),1,lower.tail=FALSE)

#####with random habs  - meristic model (ordered states)
randomHabs <- sample(unique(habs),length(habs),replace=TRUE)
names(randomHabs) <- names(habs)
MeristicActualTree<-fitDiscrete(multi2di(myTree),randomHabs,model="meristic")
MeristicStarPhylogeny<-fitDiscrete(multi2di(zeroTree),randomHabs,model="meristic")
pchisq(2 * (MeristicActualTree$opt$lnL - MeristicStarPhylogeny$opt$lnL),1,lower.tail=FALSE)