## ------------------------------------------------------------------------
## Read in the SNP data
snpData = read.table('GENE360snpData.txt',sep='\t',header=T)

## ------------------------------------------------------------------------
dim(snpData)

## ------------------------------------------------------------------------
names(snpData)

## ---- eval=FALSE---------------------------------------------------------
## ## Can view the data file within RStudio via
## View(snpData)

## ------------------------------------------------------------------------
table(snpData$Population)

## ------------------------------------------------------------------------
## Make a genotype frequency table for the first SNP
table(snpData$rs3826656)

## ------------------------------------------------------------------------
## Calculate proportions
prop.table(table(snpData$rs3826656))

## ------------------------------------------------------------------------
## Create contingeny table - genotypes across populations
table(snpData$Population, snpData$rs3826656)

## And calculate proportions (rounded)
round(prop.table(table(snpData$Population, snpData$rs3826656),1),2)

## ------------------------------------------------------------------------
snpFreqs = t(prop.table(table(snpData$Population, snpData$rs3826656), 1))
barplot(snpFreqs, beside=TRUE, legend.text=TRUE, ylim=c(0,1))

## ---- cache=TRUE---------------------------------------------------------
## Load ancestry data
snpAns = read.table('GENE360snpAncestry.txt', header=T, sep='\t')
## Check dimensionality of data
dim(snpAns)

## ---- eval=FALSE---------------------------------------------------------
## ## Can view the data file within RStudio via
## View(snpAns)

## ------------------------------------------------------------------------
table(snpAns$SubPopulation, snpAns$Population)

## ------------------------------------------------------------------------
## Create object containing only the SNP data - remove the first three columns
snpAnsDat = snpAns[,-c(1,2,3)]

## ---- cache=TRUE---------------------------------------------------------
## Load a custom function to convert to allele counts
source('alleleCounts.R')

## Apply the function to the genotype data, one column (SNP) at a time
## This will take 30 seconds or so...
snpAnsCount = apply(snpAnsDat, 2, alleleCounts)

## ---- eval=F-------------------------------------------------------------
## View(snpAnsCount)

## ---- cache=TRUE---------------------------------------------------------
## Load custom function to perform principal components analysis
source("pcaGenotypes.R")

## Apply this function to the allele count data
## This will take a couple of minutes to run
pca = pcaGenotypes(snpAnsCount)

## ------------------------------------------------------------------------
## Create an object relating to the population data
ansPop = snpAns$Population

## Generate colours to associate with each population
pCols = c("brown","red","purple","blue","green")
popCol = pCols[as.numeric(as.factor(ansPop))]
names(popCol) = ansPop

## Check that they correspond to populations
table(popCol,ansPop)

## ---- fig.height=7-------------------------------------------------------
## Load a custom plotting function to generate the plots:
source('plotPCA.R')

## Plot the PCA data that we have generated, along with the population information
plotPCA(pca, popCol)

## ------------------------------------------------------------------------
library(scatterplot3d)
scatterplot3d(pca[,1], pca[,2], pca[,3], color=popCol, pch=16,
              cex.symbols=0.5, xlab="PC1", ylab="PC2", zlab="PC3")

## ---- cache=TRUE---------------------------------------------------------
## Load custom function for calculating major homozygote frequency:
source('calcMajorFreq.R')

## Apply this function to the SNP count data
## Takes a few seconds
popFreqs = calcMajorFreq(snpAnsCount, ansPop)

## ------------------------------------------------------------------------
## Make a cluster tree based on these frequencies
plot(hclust(dist(popFreqs)),hang=-1)

## ---- eval=FALSE---------------------------------------------------------
## plot(pca[,1], pca[,2], col=popCol, pch=16, cex=0.5, xlim=1.5*range(pca[,1]))
## legend("topright", legend=names(table(snpAns$Population)), fill=pCols, cex=0.5)
## 
## ## Notes on parameters:
## ##  - "col" specifies the colors
## ##  - "pch" specifies the shape of the plotting symbol
## ##  - "cex" specifies the size of the points (or text, for the legend)
## ##  - "xlim" specifies the range of the x axis
## ##  - "fill" specifies the colours to be used in the legend

## ---- eval=FALSE---------------------------------------------------------
## load('subPopColours.RData')

