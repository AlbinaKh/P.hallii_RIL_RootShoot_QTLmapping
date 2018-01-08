
#This is an example of the workflow for QTL analysis in R/qtl for shoot and root traits measured in the P.hallii RIL mapping population



#loading rQTL library
library(qtl)


#read in data 
cross<-read.cross("csv", dir="",file="Table_S2.csv", genotypes=c("AA","BB"))
#Seperating markers that are in the same location
cross<- jittermap(cross)
#Convert a cross to type "riself" (RIL by selfing)
class(cross)[1]<-"riself"
#Imputation of genotype data
cross<-calc.genoprob(cross, step=2,
                     error.prob=0.0001,
                     map.function="kosambi")


# select the trait (Example on Shoot biomass (SQRT_SHMASS))

phe<-"SQRT_SHMASS"
plot(cross, pheno.col=phe)



#Scantwo analysis to calaculate penalties for stepwise QTL

SHMASS1000perms<-scantwo(cross, pheno.col=phe, method="hk", model="normal", n.perm=1000)

save(SHMASS1000perms, file="SHMASS1000perms.RData")


print(pens<-calc.penalties(SHMASS1000perms, alpha=0.05))
#Stepwise QTL using penalties from scantwo permutations

SHMASSstepout<-stepwiseqtl(cross, pheno.col=phe, 
                            method="hk", model="normal", penalties=pens, max.qtl=6,
                            refine.locations=TRUE, keeptrace=TRUE, keeplodprofile=TRUE)

plotLodProfile(SHMASSstepout)

title (main="Shoot Biomass, Scantwo, stepwiseqtl")
plot(SHMASSstepout)
plotModel(SHMASSstepout)

summary(fitqtl(cross, qtl=SHMASSstepout, formula=formula(SHMASSstepout), pheno.col=phe,
               method="hk", model="normal", get.ests=T, dropone=T))


#Visualizing the effects for epistatic interactions

effectplot(cross, pheno.col=phe, mname1="05@58.6",mname2= "05@136.0" )


#Calculate confidence interval for QTLs


stats<-qtlStats(cross, pheno.col = phe, mod = SHMASSstepout)

print(stats)
  
# Visualizing QTL on the MAP

cis2<-stats[c("phenotype","LOD","chr","pos","lowposition","highposition")]

cis<-cis2[complete.cases(cis2),c("phenotype","LOD","chr",
                                 "pos","lowposition","highposition")]

with(cis, segmentsOnMap(cross, phe = phenotype, chr = chr, l = lowposition, h = highposition,
                        peaklod = LOD, peakcM = pos,  showPeaks = TRUE,
                        chrBuffer = c(0.1,0.2) , tick.width=0.05,lwd=3,
                        leg.inset=.55, legendCex=0.6, legendPosition="topleft"))





#Loop for Visualizing QTL on the MAP (include all traits)

#loading qtlTools library
library(qtlTools)

# Selecting the list of traits (here named as X1, X2, X3, X4 ...)
phes <- c("X1",
          "X2",
          "X3",
          "X4")

#loading previously saved scantwo permutation files for each trait

perm.files = c("X1perms.RData", 
               "X2perms.RData",
               "X3perms.RData",
               "X4perms.RData")

for(i in perm.files) load(i)

# add all perm objects in same order as "phes"

perms = list(  X1perms,
               X2perms,
               X3perms,
               X4perms)

              
#Loop for Stepwise QTL using penalties from scantwo permutations (include all previously selected traits)

stats.list = lapply(1:length(phes), function(x){
  print(phes[x])
  print(pens<-calc.penalties(perms[[x]], alpha=0.05))
  mod<-stepwiseqtl(cross, pheno.col=phes[x], 
                   method="hk", model="normal", penalties=pens, max.qtl=6,
                   additive.only=FALSE, refine.locations=TRUE,
                   keeptrace=TRUE, keeplodprofile=TRUE)
  stats<-qtlStats(cross, pheno.col = phes[x], mod = mod)
  return(stats)
})

stats.df<-do.call(rbind, stats.list)

print(stats.df)

cis2<-stats.df[c("phenotype","LOD","chr","pos","lowposition","highposition")]

cis<-cis2[complete.cases(cis2),c("phenotype","LOD","chr",
                                 "pos","lowposition","highposition")]

print(cis)


#Visualize the positions of QTL on a genetic map.
with(cis, segmentsOnMap(cross, phe = phenotype, chr = chr, l = lowposition, h = highposition,
                        peaklod = LOD, peakcM = pos,  showPeaks = TRUE,
                        chrBuffer = c(0.1,0.2) , tick.width=0.05,lwd=3,
                        leg.inset=.55, legendCex=0.6, legendPosition="topleft"))




