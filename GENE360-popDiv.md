GENE360 tutorial
================
Associate Professor Mik Black
6 May 2020

<!-- The following will produce markdown output that will be viewble on GitHub: -->

<!-- rmarkdown::render('GENE360-popDiv.Rmd', output_format="github_document") -->

## Overview

The objective of this module is to give you some familiarity with the
types of data and tools that are commonly used for genetic studies of
population diversity.

HMTL rendered version (slightly nicer looking) of this document is
available at:

<https://github.com/mikblack/GENE360-PopDiv/blob/master/GENE360-popDiv.md>

The data we will be using are publicly available high-throughput
sequencing data from the 1000 Genomes Project. The main web page for the
1000 Genomes Project is:

<http://www.internationalgenome.org/>

The following page provides background information about the 1000
Genomes Project (and is worth reading, particularly the information
about the populations being used in the study):

<http://www.internationalgenome.org/about/>

We will be using the R and RStudio software applications to analyze the
data.

R: <https://cran.rstudio.com/> <BR> RStudio:
<https://www.rstudio.com/product/rstudio/download/>

### Getting started - loading the data

To get started, we need to load some data. I have placed a zip file
containing the data files on the GENE360 Blackboard page. Download and
unzip this file - then we need to tell RStudio where to find the data.
This involves “setting the working directory” for R - choose “Set
Working Directory” and then “Choose Directory” from the “Session” menu
in RStudio, and then find the folder containing the module data and
click “Open”.

The first data set we will look at contains genetic data on just a few
loci (specific points in the genome) across the various populations
being studied. Looking at this small data set will give you a feel for
what we mean by “genotype data” for a particular locus.

Open RStudio and read in the first data set using the following command:

``` r
## Read in the SNP data
snpData = read.table('GENE360snpData.txt',sep='\t',header=T)
```

The `##` characters indicate a comment - anything typed after them is
ignored by R. The next line tells R to read in data from the file
“GENE360snpData.txt”, which it calls “snpData”.

In this handout, the “code” used to run commands in R (and any output
produced) is indicated by `courier font`. The “output” from each command
appears below the code, prefixed by the `##` symbols.

We can ask R about the characteristics of the snpData object. The
command:

``` r
dim(snpData)
```

    ## [1] 2504   11

will tell us how many rows (individuals) and columns the data set
contains. The output shows us that there are 2504 rows and 11 columns.

The column names can be found using the following command:

``` r
names(snpData)
```

    ##  [1] "SubjectID"   "Population"  "rs3826656"   "rs13387042"  "rs4779584"   "rs2398162"   "rs7299940"   "rs1344706"  
    ##  [9] "rs7659604"   "rs734553"    "rs113010081"

To look at the full data set, you can use the “View” command:

``` r
## Can view the data file within RStudio via 
View(snpData)
```

The second column is called “Population”. We can make a table of this
information to see how many individuals are present in each population
(the “$” sign tells R to use the “Population” column from the “snpData”"
object):

``` r
table(snpData$Population)
```

    ## 
    ## AFR AMR EAS EUR SAS 
    ## 661 347 504 503 489

The super-population codes are:

AFR, African; AMR, Ad-Mixed American; EAS, East Asian; EUR, European;
SAS, South Asian.

Additional information about the composition of these populations can be
found at:

<http://www.internationalgenome.org/faq/which-populations-are-part-your-study>

## Looking at SNP frequencies

The other nine columns in the data set relate to specific single
nucleotide polymorphisms (SNPs) in the genome - we have the genotype
data for each SNP for every individual in the data set.

We can use the “table” command again to summarize the genotype
information for each SNP:

``` r
## Make a genotype frequency table for the first SNP
table(snpData$rs3826656)
```

    ## 
    ##   AA   AG   GG 
    ## 1112  952  440

We can also calculate the proportions associated with each genotype:

``` r
## Calculate proportions
prop.table(table(snpData$rs3826656))
```

    ## 
    ##        AA        AG        GG 
    ## 0.4440895 0.3801917 0.1757188

and examine differences in genotype frequencies across populations:

``` r
## Create contingeny table - genotypes across populations
table(snpData$Population, snpData$rs3826656)
```

    ##      
    ##        AA  AG  GG
    ##   AFR 445 199  17
    ##   AMR 225 109  13
    ##   EAS  47 224 233
    ##   EUR 312 161  30
    ##   SAS  83 259 147

``` r
## And calculate proportions (rounded)
round(prop.table(table(snpData$Population, snpData$rs3826656),1),2)
```

    ##      
    ##         AA   AG   GG
    ##   AFR 0.67 0.30 0.03
    ##   AMR 0.65 0.31 0.04
    ##   EAS 0.09 0.44 0.46
    ##   EUR 0.62 0.32 0.06
    ##   SAS 0.17 0.53 0.30

These results can be plotted as a bar plot:

``` r
snpFreqs = t(prop.table(table(snpData$Population, snpData$rs3826656), 1))
barplot(snpFreqs, beside=TRUE, legend.text=TRUE, ylim=c(0,1))
```

![](GENE360-popDiv_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

These analyses can be repeated for different SNPs by changing the SNP ID
(e.g., rs3826656) in the above commands.

### Why look at these SNPs?

The nine SNPs contained in this data set were not just randomly chosen -
they are SNPs that have been found to be associated with altered disease
risk. That is, an individual’s genotype at a particular position in the
genome affects their risk of developing a particular disease. This is
not absolute though, it really just raises or lowers the probability of
disease - it doesn’t guarantee complete protection or susceptibility.

| SNP                                                         | Alleles | Nearby Gene | Disease                            |
| ----------------------------------------------------------- | ------- | ----------- | ---------------------------------- |
| [rs3826656](http://www.snpedia.com/index.php/Rs3826656)     | A/G     | CD33        | Alzheimer’s Disease                |
| [rs13387042](http://www.snpedia.com/index.php/Rs13387042)   | A/G     | DIRC3       | Breast Cancer                      |
| [rs4779584](http://www.snpedia.com/index.php/Rs4779584)     | T/C     | GREM1       | Colorectal Cancer                  |
| [rs2398162](http://www.snpedia.com/index.php/Rs2398162)     | A/G     | NR2F2       | Hypertension                       |
| [rs7299940](http://www.snpedia.com/index.php/Rs7299940)     | C/G     | \-          | Panic Disorder                     |
| [rs1344706](http://www.snpedia.com/index.php/Rs1344706)     | A/C     | ZNF408A     | Schizophrenia and Bipolar Disorder |
| [rs7659604](http://www.snpedia.com/index.php/Rs7659604)     | T/C     | TMEM155     | Type 2 Diabetes                    |
| [rs734553](http://www.snpedia.com/index.php/Rs734553)       | T/G     | SLC2A9      | Gout                               |
| [rs113010081](http://www.snpedia.com/index.php/Rs113010081) | C/T     | CCR5        | HIV/AIDs “resistance”              |

For each SNP, increased risk of disease (or for rs113010081, decreased
risk) is associated with the minor allele. Variation in genotype
frequencies across populations can help to explain some of the
population-specific differences in rates of different diseases.

Clicking on SNP IDs in the table above will link through to aditional
information about each variant.

## Population diversity

The second data set that we will look at can be used to investigate
population (genetic) diversity. These SNPs were chosen because they are
“ancestry informative” - variation in their genotype frequencies tends
to be associated with differences between population groups.

``` r
## Load ancestry data
snpAns = read.table('GENE360snpAncestry.txt', header=T, sep='\t')
## Check dimensionality of data
dim(snpAns)
```

    ## [1] 2504 2305

Now we have a lot more data: 2504 individuals, and 2302 SNP genotypes
for each. Note that there are 2305 columns in the data set, but the
first three relate to non-SNP data.

Have a look using the View command:

``` r
## Can view the data file within RStudio via 
View(snpAns)
```

Note that we also have data about the sub-populations now:

``` r
table(snpAns$SubPopulation, snpAns$Population)
```

    ##          
    ##           AFR AMR EAS EUR SAS
    ##   AFR_ACB  96   0   0   0   0
    ##   AFR_ASW  61   0   0   0   0
    ##   AFR_ESN  99   0   0   0   0
    ##   AFR_GWD 113   0   0   0   0
    ##   AFR_LWK  99   0   0   0   0
    ##   AFR_MSL  85   0   0   0   0
    ##   AFR_YRI 108   0   0   0   0
    ##   AMR_CLM   0  94   0   0   0
    ##   AMR_MXL   0  64   0   0   0
    ##   AMR_PEL   0  85   0   0   0
    ##   AMR_PUR   0 104   0   0   0
    ##   EAS_CDX   0   0  93   0   0
    ##   EAS_CHB   0   0 103   0   0
    ##   EAS_CHS   0   0 105   0   0
    ##   EAS_JPT   0   0 104   0   0
    ##   EAS_KHV   0   0  99   0   0
    ##   EUR_CEU   0   0   0  99   0
    ##   EUR_FIN   0   0   0  99   0
    ##   EUR_GBR   0   0   0  91   0
    ##   EUR_IBS   0   0   0 107   0
    ##   EUR_TSI   0   0   0 107   0
    ##   SAS_BEB   0   0   0   0  86
    ##   SAS_GIH   0   0   0   0 103
    ##   SAS_ITU   0   0   0   0 102
    ##   SAS_PJL   0   0   0   0  96
    ##   SAS_STU   0   0   0   0 102

To examine population diversity, we need to do two things:

1.  create a data object of ONLY the SNP genotype data (i.e., remove the
    first three columns).
2.  convert the genotypes to allele counts (e.g., TT, AT, AA to 0, 1, 2
    - count the number of A’s).

Step 1:

``` r
## Create object containing only the SNP data - remove the first three columns
snpAnsDat = snpAns[,-c(1,2,3)]
```

Step 2:

``` r
## Load a custom function to convert to allele counts
source('alleleCounts.R')

## Apply the function to the genotype data, one column (SNP) at a time
## This will take 30 seconds or so...
snpAnsCount = apply(snpAnsDat, 2, alleleCounts)
```

Now we have a data set of just the SNP data, with genotypes converted to
allele counts.

Have a look at the new data set using the View command:

``` r
View(snpAnsCount)
```

So, what are we going to do with this new data set?

### Dimension reduction (principal components analysis)

In terms of examining population diversity, we have 2302 dimensions of
data available - one dimension for each SNP.

Rather than trying to comprehend this huge amount of data in
2302-dimensional space, genetics researchers often use a statistical
tool called Principal Components Analysis (PCA) to reduce the
dimensionality of the data.

The idea is to find the most important variation in the data, and
examine the samples in terms of that variation, ignoring the rest. In
practice, this works fairly well, because genetic differences between
populations provide a strong (and relatively consistent) source of
variation across genomic loci (i.e., SNPs). Rather than looking at 2302
dimensions of data, we end up looking at variation across just 2 or 3
dimensions - each dimension is defined by a combination of SNPs which
vary in a similar way across the individuals in the study.

``` r
## Load custom function to perform principal components analysis
source("pcaGenotypes.R")

## Apply this function to the allele count data
## This will take a couple of minutes to run
pca = pcaGenotypes(snpAnsCount)
```

To visualise the populations across the principal components, we need to
define colours for each populations. My arbitrary choice (slightly
influenced by the Y-chomosome paper from last week) was: AFR (brown),
AMR (red), EAS (purple), EUR (blue), SAS (green).

``` r
## Create an object relating to the population data
ansPop = snpAns$Population

## Generate colours to associate with each population
pCols = c("brown","red","purple","blue","green")
popCol = pCols[as.numeric(as.factor(ansPop))]
names(popCol) = ansPop

## Check that they correspond to populations
table(popCol,ansPop)
```

    ##         ansPop
    ## popCol   AFR AMR EAS EUR SAS
    ##   blue     0   0   0 503   0
    ##   brown  661   0   0   0   0
    ##   green    0   0   0   0 489
    ##   purple   0   0 504   0   0
    ##   red      0 347   0   0   0

Now we can plot the principal components and colour the points based on
the population each sample belongs to.

``` r
## Load a custom plotting function to generate the plots:
source('plotPCA.R')

## Plot the PCA data that we have generated, along with the population information
plotPCA(pca, popCol)
```

![](GENE360-popDiv_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

From the plots we can see that samples from the same population tend to
cluster together, and that the first three principal components do a
reasonable job of capturing the genetic diversity between the
populations.

With the `scatterplot3d` package, you can plot the first three principal
components at once (i.e., combining the information from the three
scatterplots above). This shows that the European (EUR), East Asian
(EAS) and South Asian (SAS) super-populations are relatively
homogeneous, while the Ad-Mixed American (AMR) and African (AFR)
super-populations exhibit greater variation, suggesting admixture within
these groups.

``` r
library(scatterplot3d)
scatterplot3d(pca[,1], pca[,2], pca[,3], color=popCol, pch=16,
              cex.symbols=0.5, xlab="PC1", ylab="PC2", zlab="PC3")
```

![](GENE360-popDiv_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

## Phylogenetic trees

Another way to visualise genetic similarity is via “phylogenetic trees”,
also known more generally as “dendrograms”.

These tree diagrams group items together based on similarity scores. In
this setting our “items” are the populations, and the scores are
calculated based on the genetic similarities between each pair of
populations.

To calculate similarities, we need to create a set of “average”
genotypes for each population. One way to do this is to calculate the
frequency of the major homozygote for each SNP in each population.

``` r
## Load custom function for calculating major homozygote frequency:
source('calcMajorFreq.R')

## Apply this function to the SNP count data
## Takes a few seconds
popFreqs = calcMajorFreq(snpAnsCount, ansPop)
```

``` r
## Make a cluster tree based on these frequencies
plot(hclust(dist(popFreqs)),hang=-1)
```

![](GENE360-popDiv_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

The way in which the “leaves” of the tree cluster, reflects the
similarity between the items analysed. In this case the populations are
grouped based on their genetic similarity, as measured by these
particular loci. Here it appears that the AMR (Ad-Mixed American), and
SAS (South Asian) super-populations are most genetically similar, then
EUR (European), followed by EAS (East Asian). The AFR (African)
super-population is the most genetically dissimilar relative to the
others. It is likely that these groupings reflect (to some degree)
migration patterns and shared ancestry experienced by these populations.

## Admixture

We can also use the PCA data to explore admixture in the 1000 Genomes
populations.

Looking at he PCA plots we made above, it is fairly clear the the first
principal component (PC1) relates to genetic variation that distinguishs
the African super-population for the other super-populations.

A very simple approach to looking at admixture involving African
ancestry would be to rescale PC1 so that it provides an estimate of
African ancestry for each individual.

<b>Note:</b> this is a very crude and simplified method for exploring
admixture - in a research setting we would use more sophisiticated
approaches to investigate individuals’ genetic ancestries.

``` r
## Get the data for PC1 (first column of the PCA matrix)
pc1 = pca[,1]

## Rescale the PC1 data so that it is on the 0-1 scale
## (subtract the minimum, and divide through by the new maximum)
pc1Prop = ( pc1 - min(pc1) ) / max( pc1 - min(pc1) )

## Check the results
range(pc1Prop)
```

    ## [1] 0 1

``` r
## And plot
hist(pc1Prop, 50, xlab="African ancestry proportion", 
     main="Proportion of African ancestry across 1000 Genomes individuals")
```

![](GENE360-popDiv_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

We can also look at this per subpopulation. Here is a summary of the
admixture proportions based on PC1 for the AFR\_ESN subpopulation:

``` r
## Create new variable for subpopulations (just to make coding easier)
subPop <- snpAns$SubPopulation 

## Have a look at African ancestry proportion for AFR_ESN subpopulation
summary( pc1Prop[subPop=="AFR_ESN"] )
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.9613  0.9740  0.9781  0.9784  0.9829  0.9950

It looks like this subpopulation has very low levels of admixture (and
uniformly high levels of African ancestry).

How about the AFR\_ASW subpopulation?

``` r
## Have a look at African ancestry proportion for AFR_ASW subpopulation
summary( pc1Prop[subPop=="AFR_ASW"] )
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.1249  0.7066  0.8033  0.7558  0.8438  0.9318

Here it look like there is more evidence of admixture.

Lets try to reproduce the style of the admixture bar plots we saw last
week.

First we need to create a matrix of data where the rows sum to 1 for
each individual.

``` r
## Create a matrix where the rows sum to one per individual (i.e., column):
admix <- rbind(pc1Prop, 1-pc1Prop)

## Add row names
rownames(admix) <- c("African ancestry", "Non-African ancestry")

## Admixture data for first ten individuals 
admix[, 1:10]
```

    ##                           [,1]      [,2]       [,3]      [,4]      [,5]      [,6]       [,7]      [,8]      [,9]
    ## African ancestry     0.8207639 0.2205823 0.01520012 0.1343279 0.1182994 0.2106127 0.03185302 0.1336951 0.2145277
    ## Non-African ancestry 0.1792361 0.7794177 0.98479988 0.8656721 0.8817006 0.7893873 0.96814698 0.8663049 0.7854723
    ##                           [,10]
    ## African ancestry     0.01200789
    ## Non-African ancestry 0.98799211

Extract data for AFR\_ASW population, and sort it based on level of
African genetic ancestry\(\dots\)

``` r
## Extract data for ASW population
asw_admix <- admix[, subPop=="AFR_ASW"]

## Sort inidividuals based on level of African genetic ancestry
asw_admix <- asw_admix[, order(asw_admix[1,])]
```

\(\dots\)and now we can make the bar plot:

``` r
## Barplot of admixture proportions for ASW subpopulation
barplot(asw_admix, legend=TRUE, args.legend=list(x="topleft"),
        main="Admixture proportions for individuals in AFR_ASW subpopulation")
```

![](GENE360-popDiv_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

How about the AFR\_ESN subpopulation?

``` r
## Extract data and sort individuals
esn_admix <- admix[, subPop=="AFR_ESN"]
esn_admix <- esn_admix[, order(esn_admix[1,])]
```

``` r
## Barplot for AFR_ESN admixture
barplot(esn_admix, legend=TRUE, args.legend=list(x="bottomleft"),
        main="Admixture proportions for individuals in AFR_ESN population")
```

![](GENE360-popDiv_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

## Assessment tasks

Use the data and code from above to answer the following questions. Hand
in your code and the output generated, along with your answers to the
questions below.

1.  Select one of the nine SNPs from the first data set and investigate
    its variation (in terms of genotype frequencies) across the
    populations of the 1000 Genomes Project. Briefly overview the
    support for this SNP being a modifier of disease risk, and try to
    find evidence (in the literature) for variation in
    population-specific risk as a result of differences in genotype
    frequency.

2.  Create a plot of PC1 *versus* PC2, colouring points by
    sub-population rather than super-population. Rather than use my
    `plotPCA` function, you will need to use R’s generic plotting
    commands. The following code generates a PCA plot from scratch,
    using the super-population colours:

<!-- end list -->

``` r
plot(pca[,1], pca[,2], col=popCol, pch=16, cex=0.5, xlim=1.5*range(pca[,1]))
legend("topright", legend=names(table(snpAns$Population)), fill=pCols, cex=0.5)

## Notes on parameters:
##  - "col" specifies the colors
##  - "pch" specifies the shape of the plotting symbol
##  - "cex" specifies the size of the points (or text, for the legend)
##  - "xlim" specifies the range of the x axis
##  - "fill" specifies the colours to be used in the legend
```

Loading the `subPopColours.RData` file via:

``` r
load('subPopColours.RData')
```

will provide two colour objects: `sbCols` and `subPopCol`, which you
should be able to use to generate a PCA plot with the points coloured by
sub-population.

Comment on the distribution of the sub-populations with respect to
genetic variability. Which of the sub-populations exhibit evidence of
admixture, and which populations (sub- or super-) are the most likely
contributors? Is there support for these findings in the literature?
Provide a brief summary of any relevant publications you find.

Generate similar plots of PC1 *versus* PC3 and PC2 *versus* PC3. Based
on these plots, describe the main population features that are captured
by each of the first three principal components.

3.  Generate a phylogentic tree for the sub-population data. Comment on
    how well the structure of the tree agrees with the relationships
    between the super- and sub-populations of the 1000 Genomes Project
    data. Are there any super-populations that are split across
    different parts of the tree? If so, what is the likely cause?
    Discuss the shape of the tree (and the placement of the super- and
    sub-populations) in the context of the “out of Africa” theory of
    human migration. How does the shape of the tree support (or
    question) the proposed migration path? Also comment of the
    relationship between the shape of the tree, and the distribution of
    the populations in the PCA plot you generated in (2). What
    advantages and disadvantages do you feel the two plotting methods
    (PCA plot and phylogenetic tree) have when visualising genetic
    diversity data?

4.  Use PC2 to investigate admixture involving European and East Asian
    ancestries. Choose two sub-populations whose members exhibit
    substantial variability in their admixture proportions for these
    ancestries, and provide summary information for each. Make admixture
    bar plots for these two sub-populations. Comment on the range of
    variation that is present, and discuss these results in the context
    of historical human migration patterns.
