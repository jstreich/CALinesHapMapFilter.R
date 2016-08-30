#change working directorysetwd("../pel2b")
setwd("/home/bvzlab/")

# read in numberic raw file
geno <- read.table("P9_10.txt",header=T,na.strings=".")

str(geno)

#strip marker info
g <- geno[,-(1:11)] 

# KDM HACKING
g2 = g
g2 = sapply(g, function (x) as.numeric(as.character(x)))
image(g2)
g = g2

#format names for data sheet matching
names.mat <- matrix(unlist(strsplit(names(g),split="_")),nrow=5)
names(g) <- names.mat[1,]
#image raw data
jpeg(file="bd6.jpg",width = 1024, height = 768)
image(1:nrow(g),1:ncol(g),as.matrix(g),xlab = "SNPs", ylab = "samples", 
       main = paste(c("SNPs","samples"),dim(g)) )
dev.off()
## polyploids look obvious in last plate and end of plate 3 ~110 samples

# identify the number of genotype calls by marker
samplecalls <- apply(g,1, function(x) sum(!is.na(x)));
#look at distribution and try to determine a threshold for markers
hist(samplecalls,breaks=ncol(g))
samp.thres <- 75
abline(v=samp.thres,lty=2)

#test a couple thresholds
table(samplecalls > samp.thres)
#FALSE  TRUE 
#33376 34970 
g.snp <- g[samplecalls > samp.thres,]

# look at counts across samples
snpsPerSample<-apply(g.snp,2,function(x) sum(!is.na(x)));
pdf("snpsSampleCALines.pdf")
hist(snpsPerSample,breaks = 40)
#test some sample count thresholds
snpCallThresh <- 1600
abline(v = snpCallThresh,lty= 2)
dev.off()
table(snpsPerSample > snpCallThresh)
names(g.snp)[snpsPerSample < snpCallThresh]

names(g.snp)[snpsPerSample < snpCallThresh]
 #[1] "X10bx2s"   "X2sx5m"    "ALM.1"     "BDTR.13o"  "BRO1.1"    "BUR.7"     "BdTR.8c"   #"CFW.3"     "CFWA.1"    "GUL1.3"    "Kah.3x10b" "Koz.4x1f"  "LAL.4"     "LAL5.2out" #"MEN.2029" 
#[16] "MRY.8"     "OSB.2"    


#FALSE  TRUE 
#    5    91 

names(g.snp)[snpsPerSample < snpCallThresh]
[1] "X2sx5m"   "BDTR.13o" "GUL1.3"   "MRY.8"    "OSB.2"   

g.minor <- g.snp[,snpsPerSample > snpCallThresh]

names.mat <- names.mat[,snpsPerSample > snpCallThresh]

# rare varients
# 1 het
rare.var <- rowSums(g.minor,na.rm=T) == 1
table(rare.var)
g.minor <- g.minor[!rare.var,]
# 1 hom
rare.var <- (rowSums(g.minor==2,na.rm=T) == 1) & (rowSums(g.minor,na.rm=T) == 2)
table(rare.var)
g.minor <- g.minor[!rare.var,]

# paralogs?
image(1:nrow(g.minor),1:ncol(g.minor),g.minor==1)
# paralogs/duplicate markers 
het.count <- rowSums(g.minor==1,na.rm=T)
hist(het.count,breaks=ncol(g.minor))
paraThresh <- 40
abline(v = paraThresh, lty = 2)
table(het.count < paraThresh)
#FALSE  TRUE 
#18644 16134 
g.minor <- g.minor[het.count < paraThresh,]
# how many total calls?
table(is.na(g.minor))
image(as.matrix(g.minor))

# minor allele frequency distribution
hist(rowSums(g.minor,na.rm=T)/2)

test <- t(as.matrix(g))
test[test==0] <- "c"
test[test==1] <- NA
test[test==2] <- "t"

require.package(ape)
library(ape)
hap.genotypes.bin <- as.DNAbin(test)
## "raw" is the proportion of sites that differ between each pair of sequences.
hap.gene.dist <- dist.dna(hap.genotypes.bin, model = "raw" , pairwise.deletion=T)

hist(hap.gene.dist,breaks=100)
image(1:nrow(test),1:nrow(test),as.matrix(hap.gene.dist))
write.csv(as.matrix(hap.gene.dist),file = "bd6geneticDist.csv")

pdf(file = "CaliforniaBrachyKinshipMatrix.pdf", paper = "a4", width = 11, height = 15)
image(1:nrow(test),1:nrow(test),as.matrix(hap.gene.dist))
dev.off()


pdf(file = "dendrogramcalibrachy.pdf",paper = "a4",width = 11, height = 15)
plot(nj(hap.gene.dist),cex=0.3)
dev.off()

pdf(file = "caliBrachy_no_hang.pdf",paper = "a4",width = 11, height = 15)
plot(hclust(hap.gene.dist),cex=0.3)
dev.off()

swil.1 <- cmdscale(hap.gene.dist, k = 8)

swil.1_df <- data.frame(swil.1)
tot.var <- sum(sapply(swil.1_df, var))
per.var <- round((sapply(swil.1_df, var)/tot.var)*100,2)

pdf(file = "caliPCoA.pdf",paper = "a4r",width = 11, height = 15)
plot(swil.1[,2]~swil.1[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Brachy PCoA" )
text(I(swil.1[,2] + .0006)~I(swil.1[,1]), 
                            labels =rownames(swil.1) ,
                            cex=0.4,col = names.mat[3,] )
#legend(500,400,paste("plate",c(2,3,4,6)),col = c(2,3,4,6),lty=1)
dev.off()

pdf(file = "caliBrachyPCoAsm.pdf",paper = "a4r",width = 11, height = 15)
plot(swil.1[,2]~swil.1[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Brachy PCoA",type="n" )
text(I(swil.1[,2] + .0006)~I(swil.1[,1]), 
                            labels =rownames(swil.1) ,
                            cex=0.6,col = names.mat[3,] )
legend(500,400,paste("plate",c(2,3,4,6)),col = c(2,3,4,6),lty=1)
dev.off()


### The other way to calculate distance
snp.cor.dist <- as.dist(1-cor(g.minor,use = "pairwise.complete.obs")^2)

pcoa.out <- pcoa(snp.cor.dist)
biplot(pcoa.out,col=rep(gps$color[index],2),cex=.1)

pdf(file = "calicorhclust.pdf",paper = "a4r",width = 11, height = 15)
plot(hclust(snp.cor.dist),cex = 0.3)
dev.off()


pdf(file = "bd6.NJcor.pdf",paper = "a4",width = 11, height = 15)
plot(nj(snp.cor.dist),cex = 0.3)
dev.off()

#_____________ old stuff

install.packages("phangorn")
library(phangorn)

snp.euc.dist <- dist(t(g.minor)) # problems with missing data
pel.upgma <- upgma(snp.euc.dist)

## operational
boot.euc.dist <- list()
upgmaboot <- list()
bootstraps <- matrix(nrow = ncol(g.minor),ncol=100)
for (i in 1:100){
boot.euc.dist[[i]] <- dist(t(g.minor[sample(1:nrow(g.minor),replace=T),]))
#snp.cor.dist <- as.dist(1-cor(g.nohet,use = "pairwise.complete.obs"))
upgmaboot[[i]] <- upgma(boot.euc.dist[[i]])
cat(i)
}

bstree <- plotBS(pel.upgma, upgmaboot)
pdf("pel4.pdf")
plot(pel.upgma,cex=0.4)
nodelabels(bstree$node.label,frame="none",cex=0.4)
dev.off()
