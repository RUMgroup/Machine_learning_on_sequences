#Table S3 from this paper lists 108 organisms, 60 that signal to a greater or lesser extent, 48 that don't.

#The modelling challenge is to work out which genes might be responsible for signalling, i.e. which genes are shared by the individuals that express the signal compared to strains taht do not express the signal.

#This contains as many of the organisms as I could find with complete genomes:

setwd("~/Desktop/")
load(file="geneContentJarozPlus.RData")

#not so impressive
table(genomes$name,genomes$GAR_inducer)
genomes[,c('name', 'GAR_inducer')]

#create a table of gene family occurrences for these organisms using seqinr
library(seqinr)
# You can run the following loop for all 28 strains. Here, for tme reasons, we only do it for one:
for(i in 1){
  #open and close each time, possibly non-ideal, give it a few tries each time
  cat(i, "\n")
  for(j in 1:10){
    #accessing Hogenome
    try(a <- choosebank("hogenom", timeout=10))
    if(class(a)!="try-error") break
  }
  if(class(a) =="try-error") closebank(); next
  #retrieve what is actually in the database belonging to the strain
  nm <- genomes$name[i]
  hg <- try(seqinr::query(listname="hgList", query=paste("sp", nm, sep="=")))
  if(class(hg)=="try-error") closebank(); next
  #now, let's try to download the real data; potentially slow so try just a few
#   kw <- try(getKeyword(hg$req))
  kw <- try(getKeyword(hg$req[1:10]))
  if(class(kw)=="try-error") closebank(); next
  closebank()
  #looking for hogenom families in keywords
  mch <- sapply(kw, pmatch, x = "HOG")
  #there is one element, that has not a complete proteome:
  #mch
  fams <- mapply(FUN = function(x,y) if (is.na(x)) "NoFamily" else y[x], y=kw, x=mch)
  tf <- table(fams)
  existing <- match(names(tf),names(genomes))
# modify existing data
#   genomes[i,existing[!is.na(existing)]] <- tf[!is.na(existing)]
  #where you don't get a match for the existing families of previous strains, add a new column for the additional gene family
  extras <- data.frame(matrix(0, nrow = nrow(genomes),ncol = sum(is.na(existing))))
  names(extras)  <- names(tf)[is.na(existing)]
  extras[i,]  <- tf[is.na(existing)]
  genomes <- cbind(genomes, extras)
}

#end up with huge table:
dim(genomes)
#there are >20k gene families but only 28 observations of the signal
#challenging problem!

n <- 17
#look at distribution of occurrence
nAny <- apply(genomes[,n:ncol(genomes)], 2, function(x) sum(x>0))
table(nAny)
#many only appear in one, only 82 are universal

#remove all where 1 level (only 44 of those)
nUniq <- apply(genomes, 2, function(x) length(unique(x)))
genomes2 <- genomes
genomes2 <- genomes2[,nUniq>1]

#also remove where have two levels, just separating one strain, more than 11k of these, i.e. down to ~11k
nKeep <- apply(genomes2, 2, function(x) {tx <- table(x); if(length(tx)==2 & min(tx)==1) FALSE else TRUE})
genomes2 <- genomes2[,nKeep]

#However, this still contains redundancy (different gene families with the same patterns of occurrence), so get rid of that:
genomes4 <- cbind(genomes2[,1:(n-1)], unique(as.matrix(genomes2[,n:ncol(genomes2)]), MARGIN=2))

#becoming a slightly more realistic problem (??)

#and still contains some co-linearity (but not a lot), which could be removed
cr4 <- cor(genomes4[,n:ncol(genomes4)])
pc4 <- which(cr4==1, arr.ind=TRUE)
pc4 <- pc4[pc4[,1]<pc4[,2],]

#Decide to send the lone 'weak' inducer as positive or negative (can later test if answers are robust to this choice)
genomes4$gia <- genomes4$GAR_inducer
levels(genomes4$gia) <- c("Negative", "Strong", "Negative")

#see if this can be worked with:

library(Boruta)
rf1 <- randomForest( genomes4[,n:(ncol(genomes4)-1)], genomes4$gia)

#you can also save more information on the explanatory variables ...
rf1 <- randomForest( genomes4[,n:(ncol(genomes4)-1)], genomes4$gia, importance=TRUE)
#... and plot it
varImpPlot(rf1)
#could be worse

#try random variable selection (testing relative to 'shadow' variables)- actually run with a lot higher maxRuns (I used 10k)

set.seed(311)
#which are significantly important
boA <- Boruta(genomes4[,n:(ncol(genomes4)-4)], genomes4$gia, maxRuns=11)
boA
#reasonable start - comes out with between 10 and 20 gene families, but really doesn't seem stable with different seeds.

#Can we use some evolutionary information to focus down more clearly?
#plot out the correlations among organisms
crts <- cor(t(genomes4[,n:(ncol(genomes4)-4)]), method="spearman")
rownames(crts) <- paste(genomes4$name, genomes4$Strain,sep="_") -> colnames(crts)
library(ggplot2)
library(reshape2)
crtsm <- melt(crts)
ggplot(crtsm, aes(Var1, Var2, fill = value)) + geom_tile()+scale_fill_gradient(low = "blue",  high = "red")

genomes$name
dif <- which(genomes4[4,] !=genomes4[5,])
genomes5 <- genomes4[,c(1:(n-1),dif[dif>=n]) ]
#comes down to 62

rf2 <- randomForest( genomes5[,n:(ncol(genomes5)-1)], genomes5$gia, importance = TRUE)
varImpPlot(rf2)
set.seed(7784)
boB <- Boruta(genomes5[,n:(ncol(genomes5)-4)], genomes4$gia, maxRuns=10000)
plot(boB)

# really want to
impE <- getSelectedAttributes(boB)

rfE <- randomForest(formula(paste("gia", paste(impE, collapse="+"), sep="~")),data=genomes4, importance=TRUE)

library(GGally)
ggpairs(genomes5[,c("GAR_inducer", impE)])

for(i in 1:1){
  df <- data.frame(genomes5$name, genomes5$Strain,genomes5$GAR_inducer, predict(rfE),genomes5[,impE[i]])
  names(df) <- c("name", "Strain", "GAR_inducer", "predicted", impE[i])
  print(df)
}

#so what ARE these things?? (try 4 first then 2)
library(seqinr)
choosebank("hogenom")
for(i in 2){
  cat(impE[i],"\n")
  cat("E. coli MG1655, Signaller", "\n")
  qa <- paste("K=", impE[i], " AND SP=", "ESCHERICHIA COLI STR. K-12 SUBSTR. MG1655", sep="")
  hga <- query(listname="hgList", query=qa)
  print(getKeyword(hga))
  cat("E. coli W3110, non-Signaller", "\n")
  qp <- paste("K=", impE[i], " AND SP=", "ESCHERICHIA COLI STR. K-12 SUBSTR. W3110", sep="")
  hgp <- query(listname="hgList", query=qp)
  print(getKeyword(hgp))
}

#pull out sequences (could look at at http://doua.prabi.fr/databases/hogenom/home.php)
seqSig <- getSequence(hga)
seqNon <- getSequence(hgp)
closebank()

seqq <- c(seqSig,seqNon)
nm <- c(sapply(hga$req,paste),sapply(hgp$req,paste))
write.fasta(seqq,names=nm,file="sequences.fasta")

#could attempt to align within R using tools from ape and/or phyloch, however seems a little flakey
#leave R and look for outlier online http://www.ebi.ac.uk/goldman-srv/webprank/
#check out GatR, e.g. at http://ecoliwiki.net or http://www.xbase.ac.uk/colibase/
