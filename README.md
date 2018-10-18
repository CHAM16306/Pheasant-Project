# Pheasant-Project


# Step1: Look at fastq data (in terminal)

Last login: Thu Oct 11 12:38:35 on ttys000

#So now when you open terminal, navigate to that folder using this command: cd /User/projects/PheasantProject

```eduroam-140-210:~ aminachabach$ cd /Users/aminachabach/Project/Fastq```

#where /User/projects/PheasantProject is whatever your file path is. 
You can check where your terminal is with the following command: pwd

````
eduroam-140-210:Fastq aminachabach$ pwd

/Users/aminachabach/Project/Fastq
`````

#You can then see if your PHE079.fastq.gz files in there with: ls

```
eduroam-140-210:Fastq aminachabach$ ls

Pheasant1.fastq	Pheasant2.fastq
````

#The following commands assume that you're unzipped the files. 
#This command shows you the first few lines of a document: head filename 

```
eduroam-140-210:Fastq aminachabach$ head Pheasant1.fastq
@1_1101_14596_1577

ATTGCAAAAAATTATGGTATGAAAATCACCTTCTACTCTCTCCCCAGAGTATTGGCTGGTTTCTTCATTTTTTTTAATTCCTCTTCATTTTCATTTTATAATCCCCAAGTTCCCATTAAGCATTTTGTAGAAGATCGAACTACAAAAG
+
BGGFEGGGGEEGCGFHF5DGHFHHHHHHFHHHHHHHHGHHHHHHGHEGHAGFHHHHHGGAGHHHHHHHGHHHGGGHBGHGGHGEFEGFHHHHHHHHHHHFHHFEGGHGAHHHFFDHHFGHHFHHHHB??BDGHGFHHFGFHGB2FDGH
@1_1101_15010_1581
GTGTGTAGAAGTCACAGTCAGAACACATAGGCTCTTCTGTTAATGCACCAGCACTGCATGAAGTTAAAGACACTACTCAGATGTTTCCAGCACATTACATGTACTTCTGTGTCAAGTAGTTAGCAGGTCAAATAGCCTCCAGGGATCT
+
BGFGGGGGGGGHHHHHHHHHHHHHHHHHGGHHHHHHHHHFHHHHHHHHHHHHHGHEHHHHHFFGHFHFHHHGHGGHHHGHHHHGHHHHGGGHHHHHHGHHGGHHHHHHHHHHHHHFGHHHHHHHHHHHHHHHHHHHHHHHHEE/FHEH
@1_1101_17524_1628
CTTCCTACTGCACGACACATCAGCTCCCAGGGCTGCTCCTGCAGGAGTTCCATGGGCTGCAGCCTCCTTCAAGCCACATCCACTGCTGCACCACGTGCTCCTCCACAGCTGCACATGGAGATGTGTTCCAGGTGGTGCCCATGGGCT
````

#This opens a page at a time. You can scroll through using spacebar, and quit this view with q: less filename

``` 
eduroam-140-210:Fastq aminachabach$ less Pheasant1.fastq
````

#shows the last few lines of a file: tail filename

```
eduroam-140-210:Fastq aminachabach$ tail Pheasant1.fastq
+
@FGGGFDGGBEDFHHGB4AG2AEGFFH5BCFEGEAADGGHGD5GGDGGGCGFCGEHHHHH5BCGFDGHHFHHH3CEHHFA1FAFHGHHHHH4DGGFGHH3@FGFHBGHHHHHHHBFFGGFHFHHGHHFHHHHBGHHHHFFFHHHD1F?
@1_2114_15195_28873
AATACTGCGGAGATAGTGTATGTGAAAATATGCAGCTACAAGCAGAGATGCCCAAACTCTTTGGGTTAGTTTCGATCAAGATTGGCAGAGAGGCAGAAATCTCAGAACTGTGAAGCACCCACAAGTCTGAAGGAAATGGGCAGCTG
+
BGGGFGGGGG?EEHFHHGGHHHHHFHHHHHHHHGHHEGHHHHHHHHGGHHHHHHEHHHHHHGHHGFGGHHHFHGGHHGHHFHHHHHHGGHGGGFGGGHHHHHHFHFHHHHHHHFHHHHHGGGGHFFHFFHHHHGHGHHFHHHGGHH
@1_2114_12421_28917
CCAGTTCAGGAATAACTGGAGGAGAGTCAAGTGCTGGATACCTGGAATCACTACAGAAACATGAATTCAGAAGCAGGATCACTGCTATAACCTGGTGTCTCAGAACTTTCTCTCCCTGTGGCTAATGGGGCTTCCCTCAAATCCCAGC
+
BCEGFGGFFCGGHHHHFHH4EGFF2EHFDGHBGFH5GC5F3GHHHFHHGEGGHHHHGHHGHBGHHHHHHHHGHBGGEFFHHGHHGB@GFHGHFFEGHCGHD4FGHBEFHHHGBH33BFBGFEAFEDHFFDEEC2EGHFFG23GFGHHG
eduroam-140-210:Fastq aminachabach$ 
````

# Step2: Filter the data in linux

```
eduroam-140-210:Fastq aminachabach$ cd /Users/aminachabach/Project/WPAallMerged_outfiles

eduroam-140-210:WPAallMerged_outfiles aminachabach$ pwd
/Users/aminachabach/Project/WPAallMerged_outfiles

eduroam-140-210:WPAallMerged_outfiles aminachabach$ ls
WPAallMerged.alleles.loci	WPAallMerged.nex		WPAallMerged.u.geno
WPAallMerged.geno		WPAallMerged.phy		WPAallMerged.u.snps.phy
WPAallMerged.gphocs		WPAallMerged.snps.map		WPAallMerged.ustr
WPAallMerged.hdf5		WPAallMerged.snps.phy		WPAallMerged.vcf
WPAallMerged.loci		WPAallMerged.str		WPAallMerged_stats.txt

eduroam-140-210:WPAallMerged_outfiles aminachabach$ less WPAallMerged.loci
````

# vcf Installation

```
eduroam-140-210:vcftools aminachabach$ cd /Users/aminachabach/Project/vcftools_0.1.13
eduroam-140-210:vcftools_0.1.13 aminachabach$ pwd
/Users/aminachabach/Project/vcftools_0.1.13
eduroam-140-210:vcftools_0.1.13 aminachabach$ ls
Makefile	README.md	README.txt	cpp		examples	perl

eduroam-140-210:vcftools_0.1.13 aminachabach$ less Makefile

eduroam-140-210:vcftools_0.1.13 aminachabach$ make
...
g++ -O2 -D_FILE_OFFSET_BITS=64  vcftools.o bcf_file.o vcf_file.o variant_file.o header.o bcf_entry.o vcf_entry.o entry.o entry_getters.o entry_setters.o vcf_entry_setters.o bcf_entry_setters.o entry_filters.o variant_file_filters.o variant_file_output.o parameters.o variant_file_format_convert.o variant_file_diff.o output_log.o bgzf.o gamma.o -o vcftools -lz 
cp /Users/aminachabach/Project/vcftools_0.1.13/cpp/vcftools /Users/aminachabach/Project/vcftools_0.1.13//bin/vcftools

````

# To find path if using the same
```
eduroam-140-210:vcftools_0.1.13 aminachabach$ echo $PATH
/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin

#eduroam-140-210:vcftools_0.1.13 aminachabach$ ls ~/
Applications	Assignmnet	Desktop		Downloads	Movies		Pictures	Project
Assignment	BirthdayStick	Documents	Library		Music		Practical 2	Public
````

## Put WPAallMerged.vcf file into vcftools folder
#We can get rid of a lot of the data that we know are not useful using maf etc
#How much data is left after this? Note that the input file here is the output file from before. 

```
eduroam-140-210:WPAallMerged_outfiles aminachabach$ cd /Users/aminachabach/Project/vcftools_0.1.13
````

## ALLELE FILTERING: Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value and less than or equal to the "--max-maf" value. 
#Although I've specified the name as filename.maf, vcftools adds the ## recode.vcf extension ## Remember to add the full name in your command.

```
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools  --vcf WPAallMerged.vcf --maf 0.05 --recode --recode-INFO-all --out filter.maf

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf WPAallMerged.vcf
	--recode-INFO-all
	--maf 0.05
	--out filter1
	--recode

After filtering, kept 69 out of 69 Individuals
Outputting VCF file...
After filtering, kept 125671 out of a possible 184264 Sites
Run Time = 23.00 seconds
```

## Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).
## Remember filename.recode.vcf extension ## 

```
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools --vcf filter.maf.recode.vcf --max-missing 0.5--recode --recode-INFO-all --out filter.missing

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf filter1.recode.vcf
	--recode-INFO-all
	--max-missing 0.5
	--out filter.missing

After filtering, kept 69 out of 69 Individuals
After filtering, kept 41503 out of a possible 125671 Sites
Run Time = 3.00 seconds
```

## First we'll use vcftools to tell us how much data is missing per individual, and then we'll plot it in the command line. This creates an output called out.imiss.

```
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools --vcf filter.missing.recode.vcf --missing-indv

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf filter.missing.recode.vcf
	--missing-indv

After filtering, kept 69 out of 69 Individuals
Outputting Individual Missingness
After filtering, kept 41503 out of a possible 41503 Sites
Run Time = 1.00 seconds
````

# Use the following to look at the data (quit with q)

````
cat out.imiss
less out.imiss
head out.imiss
tail out.imiss
````

# What we're interested in is the F_MISS column. This tells us the proportion of missing data for that individual.You can just copy and paste this next bit in your command line. This is a bit of linux code to draw of histogram of the F_MISS column: 

## Install gnuplot on Mac OSX
````
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null

...Password...
...Downloading and Installing....

brew install gnuplot
````

## Plot histogram 
````
eduroam-140-210:vcftools_0.1.13 aminachabach$ awk '!/IN/' out.imiss | cut -f5 > totalmissing
eduroam-140-210:vcftools_0.1.13 aminachabach$ gnuplot << \EOF 
> set terminal dumb size 120, 30
> set autoscale 
> unset label
> set title "Histogram of % missing data per individual"
> set ylabel "Number of Occurrences"
> set xlabel "% of missing data"
> #set yr [0:100000]
> binwidth=0.01
> bin(x,width)=width*floor(x/width) + binwidth/2.0
> plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
> pause -1
> EOF
````
##  For now we can identify individuals with a lot of missing data and remove them. Decide on a cutoff based on your results from above. You can now use vcftools to identify these individuals by name:

#This is a bit of linux code that first finds all the lines where F_MISS > 0.5. It then copies the names of these lines (i.e. the first column) to an output file called lowDP.indv

```
eduroam-140-210:vcftools_0.1.13 aminachabach$ awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
````

#we can look at this like before

````
eduroam-140-210:vcftools_0.1.13 aminachabach$ cat lowDP.indv
INDV
PHE113.control3
PHE133
PHE134
PHE135
PHE137
````
#and we can count the number of indivs we'll lose

````
eduroam-140-210:vcftools_0.1.13 aminachabach$ cat lowDP.indv |wc -l
       6
````

## Once you're happy with this cutoff, we can remove the individuals from the vcf file

````
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools --vcf filter.missing.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out filename.final

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf filter.missing.recode.vcf
	--exclude lowDP.indv
	--recode-INFO-all
	--out filename.final
	--recode

Excluding individuals in 'exclude' list
After filtering, kept 64 out of 69 Individuals
Outputting VCF file...
After filtering, kept 41503 out of a possible 41503 Sites
Run Time = 5.00 seconds
````
#and see what your data look like now
````
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools --vcf filename.final.recode.vcf

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf filename.final.recode.vcf

After filtering, kept 64 out of 64 Individuals
After filtering, kept 41503 out of a possible 41503 Sites
Run Time = 1.00 seconds
````

## Remove individual PHE126 (nachträglich gemacht) before filtering with Hardy and heterzygocity exces
#We did this bc seems to be outlier (ie. Malaysian that clustered with Mountain sp.)

````
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools --remove-indv PHE126 --vcf filename.final.recode.vcf --recode --recode-INFO-all --out filename.removed.vcf

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf filename.final.recode.vcf
	--out filename.removed.vcf
	--remove-indv PHE126

Excluding individuals in 'exclude' list
After filtering, kept 63 out of 64 Individuals
After filtering, kept 41503 out of a possible 41503 Sites
Run Time = 1.00 seconds
````


##   Summary stats in R
#https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html
#To run the package, you need to install the package and load it using the following command lines:
````
install.packages("pcadapt")
library(pcadapt)
````
## A. Reading genotype data
```
install.packages("pcadapt")
library(pcadapt)
````

#For example, assume your genotype file is called “foo.lfmm” and is located in the directory “path_to_directory”, use the following command lines:
#data into r (change type to vcf bc it's a vcf file)
```
path_to_file <- "/Users/aminachabach/Project/vcftools_0.1.13/filename.final.recode.vcf"
filename <- read.pcadapt(path_to_file, type = "vcf")
````

#filename : shows whta the data is :
````
filename
[1] "/private/var/folders/s4/wvv7tp457svc9s2p1tgv6rzh0000gn/T/RtmpcqN6VJ/fileedb23ec34b55.pcadapt.bed"
attr(,"n")
[1] 64
attr(,"p")
[1] 41503
attr(,"class")
[1] "pcadapt_bed"
````

## B. Choosing the number K of Principal Components
#NB: by default, data are assumed to be diploid. To specify the ploidy, use the argument ploidy (ploidy=2 for diploid species and  ploidy = 1 for haploid species) in the pcadapt function.
#Look at what the terms mean
````
x <- pcadapt(input = filename, K = 20, ploidy=2)
summary(x)
                Length Class  Mode   

scores            1280 -none- numeric	- is a (n,K) matrix corresponding to the projections of the individuals onto each PC.

singular.values     20 -none- numeric	-  is a vector containing the K ordered squared root of the proportion of variance explained by each PC.

loadings        830060 -none- numeric	-  is a (L,K) matrix containing the correlations between each genetic marker and each PC.

zscores         830060 -none- numeric	- is a (L,K) matrix containing the z-scores (A z-score (aka, a standard score) indicates how many standard deviations an element is from the mean. A z-score can be calculated from the following formula. where z is the z-score, X is the value of the element, μ is the population mean, and σ is the standard deviation.)

af               41503 -none- numeric	- is a vector of size L containing allele frequencies of derived alleles where genotypes of 0 are supposed to code for homozygous for the reference allele.

maf              41503 -none- numeric	- missing allele frequency

chi2.stat        41503 -none- numeric	- is a vector of size L containing minor allele frequencies.

stat             41503 -none- numeric	- is a vector of size L containing squared Mahalanobis distances by default.

gif                  1 -none- numeric	- is a numerical value corresponding to the genomic inflation factor estimated from   

pvalues          41503 -none- numeric	- is a vector containing L p-values.

pass             40992 -none- numeric	- A list of SNPs indices that are kept after exclusion based on the minor allele frequency threshold.

> 
`````


## B.1. Scree plot
#The ‘scree plot’ displays in decreasing order the percentage of variance explained by each PC. Up to a constant, it corresponds to the eigenvalues in decreasing order. The ideal pattern in a scree plot is a steep curve followed by a bend and a straight line. The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. We recommend to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule). In the provided example, K = 2 is the optimal choice for K. 
#Another option to choose the number of PCs is based on the ‘score plot’ that displays population structure. The score plot displays the projections of the individuals onto the specified principal components. Using the score plot, the choice of K can be limited to the values of K that correspond to a relevant level of population structure.
#The ‘scree plot’ displays in decreasing order the percentage of variance explained by each PC. Up to a constant, it corresponds to the eigenvalues in decreasing order. The ideal pattern in a scree plot is a steep curve followed by a bend and a straight line. The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. We recommend to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule). In the provided example, K = 2 is the optimal choice for K. The plot function displays a scree plot:
````
plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 10)
`````

## B.2. Score plot

#What we tried:
````

poplist.names <- pheasantpops[,7]
poplist.names <- data.frame(poplist.names)

poplist.colour <- pheasantpops[,8]
poplist.colour <- data.frame(poplist.colour)

str(poplist.names)
str(poplist.colour)

poplist.names2 <- c(poplist.colour, poplist.names)
str(poplist.names2)

plot(x, option = "scores", pop = poplist.names2)

....


poplist.names <- c(rep("POP1", 20),rep("POP2", 20),rep("POP3", 24))
plot(x, option = "scores", pop = poplist.names)
str(poplist.names)

str(poplist.colour)
as.character()
poplist.colour1 <- as.character(poplist.colour)
str(poplist.colour1)
`````

#What worked (iport csv file, ecel files seems to have been the problem):

````
poplist.csv <- read.csv("/Users/aminachabach/Project/pheasantpops.csv", header=TRUE)
head(poplist.csv)

poplist.names <- poplist.csv[,7]
str(poplist.names)

plot(x, option = "scores", pop = poplist.names)
````

### An introduction to adegenet 2.0.0
```
R.version.string
````

### 2 Getting started
## 2.1 Installing the package - stable version

#Install adegenet with dependencies using:
````
install.packages("adegenet", dep=TRUE)
````

#We can now load the package alongside other useful packages:
```
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
````

#If at some point you are unsure about the version of the package, you can check it using:
````
packageDescription("adegenet", fields = "Version")
````

#adegenet version should read 2.0.0. --> It doesn't, 2.1.1

## 2.2 Installing the package - devel version

#The development of adegenet is hosted on github:https://github.com/thibautjombart/adegenet.
#You can install this version using the package devtools and the following commands:
```
library("devtools")
install_github("thibautjombart/adegenet")
library("adegenet")
````

## 3.2 genpop objects

#These objects store genetic data at a population level, plus various meta-data. Their struture
#is nearly identical to genind objects, only simpler as they no longer store group membership
#for individuals. genpop objects are created from genind objects using genind2genpop:
```
file <- read.structure("/Users/aminachabach/Project/finalname.final.recode.str")
````
#type number in the console bit: Output = (btw toook absolutely ages to load???)
````
#/// GENIND OBJECT /////////

#// 64 individuals; 41,503 loci; 83,436 alleles; size: 42.4 Mb

#// Basic content
#@tab:  64 x 83436 matrix of allele counts
#@loc.n.all: number of alleles per locus (range: 2-4)
#@loc.fac: locus factor for the 83436 columns of @tab
#@all.names: list of allele names for each locus
#@ploidy: ploidy of each individual  (range: 2-2)
#@type:  codom
#@call: read.structure(file = "/Users/aminachabach/Project/finalname.final.recode.str")

#// Optional content
#@pop: population of each individual (group size range: 2-22)
````

````
data(file)
````

#still doesn't work btw
````
phepop <- genind2genpop(file)
phepop

/// GENPOP OBJECT /////////

 // 8 populations; 41,503 loci; 83,436 alleles; size: 24.6 Mb

 // Basic content
   @tab:  8 x 83436 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-4)
   @loc.fac: locus factor for the 83436 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: genind2genpop(x = file)

 // Optional content
   - empty -
````

#As in genind objects, data are stored as numbers of alleles, but this time for populations
#(here, phe colonies): 5 rows & 10 colums
````
phepop$tab[1:5,1:10]
````

## 3.3 Using accessors
#One advantage of formal (S4) classes is that they allow for interacting simply with possibly complex objects. This is made possible by using accessors, i.e. functions that extract information from an object, rather than accessing the slots directly. Another advantage of this approach is that as long as accessors remain identical on the user’s side, the internal structure of an object may change from one release to another without generating errors in old scripts. Although genind and genpop objects are fairly simple, we recommend using accessors whenever possible to access or modify their content.

#nancycats=file in my code

#indNames†: returns/sets labels for individuals; only for genind.
```
head(indNames(file),10)
````


#nInd: returns the number of individuals in the object; only for genind.
#DONT'T RUN OUR NAMES ARE FINE AS THEY ARE
````
indNames(nancycats) <- paste("cat", 1:nInd(nancycats),sep=".")
head(indNames(nancycats),10)
````

#locNames: returns/sets labels for l
````
locNames(file)
`````

#returns the names of the alleles in the form ’loci.allele’.
````
temp <- locNames(file, withAlleles=TRUE)
head(temp, 10)
````

#The slot ’pop’ can be retrieved and set using pop:
````
obj <- file[sample(1:50,10)]
pop(obj)
# [1] 2 1 5 1 2 2 3 2 1 1
#Levels: 1 2 3 5

pop(obj) <- rep("newPop",10)
pop(obj)
````

#Accessors make things easier. For instance, when setting new names for loci, the columns
#of @tab are renamed automatically:
````
head(colnames(tab(obj)),20)

locNames(obj)

locNames(obj)[1] <- "newLocusName"
locNames(obj)

head(colnames(tab(obj)),20)

obj@pop <- rep("newPop",10)
````

## 5.4 Measuring and testing population structure (a.k.a F statistics)
```
library(devtools)
install_github("jgx65/hierfstat")
````
#Haven't actually installed this, do we even need it bc the net line runs?

#Can we find any population structure in the pheasent colonies? The basic F statistics are provided by:
```
library("hierfstat")
fstat(file)
````

#This table provides the three F statistics F st (pop/total), F it (Ind/total), and F is
#(ind/pop). These are overall measures which take into account all genotypes and all loci.
#For more detail, pegas provides estimates by locus:

#This table provides the three F statistics F st (pop/total), F it (Ind/total), and F is
#(ind/pop). These are overall measures which take into account all genotypes and all loci.

````
	pop       Ind
Total 0.7958457 0.8203680
pop   0.0000000 0.1201165
````
#SO: FST = 0.7958457; FIT = 0.8203680; FIS = 0.1201165


#The measures FIS, FST, and FIT are related to the amounts of heterozygosity at various levels of population structure. Together, they are called F-statistics, and are derived from F, the inbreeding coefficient.

#FIT is the inbreeding coefficient of an individual (I) relative to the total (T) population, as above; FIS is the inbreeding coefficient of an individual (I) relative to the subpopulation (S)*, using the above for subpopulations and averaging them; and FST is the effect of subpopulations (S) compared to the total population (T)


#The fixation index (FST) is a measure of population differentiation due to genetic structure. It is frequently estimated from genetic polymorphism data, such as single-nucleotide polymorphisms (SNP) or microsatellites. Developed as a special case of Wright's F-statistics, it is one of the most commonly used statistics in population genetics.

#For more detail, pegas provides estimates by locus:
```
library(pegas)
Fst(as.loci(file))


                    Fit           Fst           Fis
SNP_1      1.0000000000  1.0000000000  1.000000e+00
SNP_2      1.0000000000  1.0000000000  1.000000e+00
SNP_3      1.0000000000  1.0000000000  1.000000e+00
SNP_4      0.6421931609 -0.1550400496  6.902213e-01
SNP_5      1.0000000000  1.0000000000  1.000000e+00
SNP_6      1.0000000000  1.0000000000  1.000000e+00
SNP_7      1.0000000000  1.0000000000  1.000000e+00
SNP_8      1.0000000000  1.0000000000  1.000000e+00
SNP_9      1.0000000000  1.0000000000  1.000000e+00
SNP_10     1.0000000000  1.0000000000  1.000000e+00
SNP_11     0.8867695280  0.7440658303  5.575797e-01
SNP_12     1.0000000000  1.0000000000           NaN
SNP_13     1.0000000000  1.0000000000  1.000000e+00
SNP_14     0.8640776699  0.7281553398  5.000000e-01
SNP_15     1.0000000000  1.0000000000  1.000000e+00
SNP_16     1.0000000000  1.0000000000  1.000000e+00
SNP_17     1.0000000000  1.0000000000  1.000000e+00
SNP_18     1.0000000000  1.0000000000           NaN
SNP_19     1.0000000000  1.0000000000           NaN
SNP_20     0.6812609457  0.3625218914  5.000000e-01
SNP_21     0.6812609457  0.3625218914  5.000000e-01
SNP_22     1.0000000000  1.0000000000           NaN
SNP_23     0.7669654289  0.5339308579  5.000000e-01
SNP_24     0.5766515243  0.4701918341  2.009401e-01
SNP_25     1.0000000000  1.0000000000           NaN
SNP_26     1.0000000000  1.0000000000  1.000000e+00
SNP_27     1.0000000000  1.0000000000  1.000000e+00
SNP_28     1.0000000000  1.0000000000  1.000000e+00
SNP_29     1.0000000000  1.0000000000  1.000000e+00
SNP_30     0.8319085778  0.8474725984 -1.020408e-01
SNP_31     0.7443786398  0.7879185418 -2.052980e-01
SNP_32     0.1671661829  0.3979426864 -3.833132e-01
SNP_33     0.3916322500  0.4686658287 -1.449814e-01
SNP_34     1.0000000000  1.0000000000  1.000000e+00
SNP_35     0.5543013092  0.0323774610  5.393879e-01
SNP_36     0.9098902554  0.9133560148 -4.000000e-02
SNP_37    -0.0327237103  0.1431797788 -2.052980e-01
SNP_38     0.5549385165  0.0357001190  5.384615e-01
SNP_39     0.5549385165  0.0357001190  5.384615e-01
SNP_40     0.5692819604  0.6462869185 -2.177046e-01
SNP_41     1.0000000000  1.0000000000  1.000000e+00
SNP_42     1.0000000000  1.0000000000  1.000000e+00
SNP_43     1.0000000000  0.9495401014  1.000000e+00
SNP_44     0.4428679899  0.2814425614  2.246521e-01
SNP_45     0.9663496219  0.9733601173 -2.631579e-01
SNP_46     0.5550970567  0.4315129058  2.173913e-01
SNP_47     1.0000000000  1.0000000000  1.000000e+00
SNP_48     0.2027156896  0.1702030567  3.918143e-02
SNP_49     0.3421278654  0.3236099831  2.737752e-02
SNP_50     0.6829286917  0.6740037216  2.737752e-02
SNP_51     0.5521696731  0.5395640787  2.737752e-02
SNP_52     0.1152743998  0.1348992222 -2.268501e-02
SNP_53     0.8374061609  0.7370249173  3.817139e-01
SNP_54     0.3756832536  0.1373077686  2.763158e-01
SNP_55     0.3365322601  0.1667879163  2.037229e-01
SNP_56     0.9477620429  0.9539356196 -1.340206e-01
SNP_57     1.0000000000  1.0000000000  1.000000e+00
SNP_58     1.0000000000  1.0000000000           NaN
SNP_59     1.0000000000  1.0000000000           NaN
SNP_60     0.8232040254  0.2511103844  7.639225e-01
SNP_61     1.0000000000  1.0000000000  1.000000e+00
SNP_62     0.8232040254  0.2511103844  7.639225e-01
SNP_63     1.0000000000  1.0000000000  1.000000e+00
SNP_64     1.0000000000  1.0000000000  1.000000e+00
SNP_65     0.9673958993  0.9738749193 -2.480000e-01
SNP_66     1.0000000000  1.0000000000  1.000000e+00
SNP_67     1.0000000000  1.0000000000           NaN
SNP_68     1.0000000000  1.0000000000  1.000000e+00
SNP_69     0.1359605833  0.1034725627  3.623762e-02
SNP_70     1.0000000000  1.0000000000  1.000000e+00
SNP_71     0.3623817980 -0.0464883058  3.907068e-01
SNP_72     1.0000000000  0.6574982481  1.000000e+00
SNP_73     1.0000000000 -0.0643274854  1.000000e+00
SNP_74     1.0000000000 -0.0643274854  1.000000e+00
SNP_75     1.0000000000  1.0000000000  1.000000e+00
SNP_76     1.0000000000  0.8867746129  1.000000e+00
SNP_77     1.0000000000  1.0000000000  1.000000e+00
SNP_78     1.0000000000  1.0000000000  1.000000e+00
SNP_79     0.7096078130  0.4923210784  4.280003e-01
SNP_80     0.7331418632  0.4860509958  4.807692e-01
SNP_81     0.5491891966  0.0644284438  5.181440e-01
SNP_82     1.0000000000  1.0000000000  1.000000e+00
SNP_83     1.0000000000  1.0000000000  1.000000e+00
SNP_84     1.0000000000  1.0000000000  1.000000e+00
SNP_85     0.7474480426  0.7757599187 -1.262570e-01
SNP_86     0.9714326164  0.9153181130  6.626506e-01
SNP_87     1.0000000000  1.0000000000  1.000000e+00
SNP_88     1.0000000000  0.8469529086  1.000000e+00
SNP_89     0.7474480426  0.7757599187 -1.262570e-01
SNP_90     1.0000000000  1.0000000000  1.000000e+00
SNP_91     1.0000000000  1.0000000000  1.000000e+00
SNP_92     1.0000000000  0.9593148666  1.000000e+00
SNP_93     1.0000000000  1.0000000000  1.000000e+00
SNP_94     1.0000000000  1.0000000000  1.000000e+00
SNP_95     0.5173975586  0.6338988094 -3.182214e-01
SNP_96     0.8730219048  0.3995414148  7.885315e-01
SNP_97     1.0000000000  1.0000000000  1.000000e+00
SNP_98     1.0000000000  1.0000000000  1.000000e+00
SNP_99     1.0000000000  1.0000000000  1.000000e+00
SNP_100    0.1849155893  0.2214188027 -4.688427e-02
SNP_101   -0.0312006319  0.0792851501 -1.200000e-01
SNP_102    1.0000000000  1.0000000000  1.000000e+00
SNP_103    1.0000000000  1.0000000000  1.000000e+00
SNP_104    1.0000000000  1.0000000000  1.000000e+00
SNP_105    1.0000000000  1.0000000000           NaN
SNP_106    1.0000000000  1.0000000000  1.000000e+00
SNP_107    1.0000000000  1.0000000000  1.000000e+00
SNP_108    1.0000000000  1.0000000000  1.000000e+00
SNP_109    0.8600128080  0.8694350229 -7.216495e-02
SNP_110    1.0000000000  1.0000000000           NaN
SNP_111    0.6866524557  0.6322311164  1.479770e-01
SNP_112    0.1656627733  0.4175948474 -4.325718e-01
SNP_113   -0.1008970008  0.1171237778 -2.469438e-01
SNP_114    0.4570336526  0.2015200773  3.200000e-01
SNP_115    0.4570336526  0.2015200773  3.200000e-01
SNP_116    0.3187442228  0.1906232803  1.582958e-01
SNP_117    0.3702866039 -0.0241417215  3.851306e-01
SNP_118   -0.0169412496  0.2843746762 -4.210526e-01
SNP_119    0.9512884086  0.9570452331 -1.340206e-01
SNP_120    1.0000000000  1.0000000000  1.000000e+00
SNP_121    0.2339639307  0.3078769150 -1.067917e-01
SNP_122    0.9757139877  0.9785841528 -1.340206e-01
SNP_123    0.5846152782  0.1380767024  5.180723e-01
SNP_124    1.0000000000  1.0000000000  1.000000e+00
SNP_125    1.0000000000  1.0000000000  1.000000e+00
SNP_126    1.0000000000  1.0000000000  1.000000e+00
SNP_127    1.0000000000 -0.0910297533  1.000000e+00
SNP_128    1.0000000000  1.0000000000  1.000000e+00
SNP_129    0.7819852977  0.5879731756  4.708726e-01
SNP_130    1.0000000000  1.0000000000           NaN
SNP_131    1.0000000000  0.6003734827  1.000000e+00
SNP_132    1.0000000000  1.0000000000  1.000000e+00
SNP_133    1.0000000000  1.0000000000  1.000000e+00
SNP_134    1.0000000000  1.0000000000           NaN
SNP_135    1.0000000000  1.0000000000           NaN
SNP_136    0.8772153212  0.7157062436  5.681063e-01
SNP_137    0.7178088371  0.7415056591 -9.167250e-02
SNP_138   -0.0644116020  0.0817987829 -1.592357e-01
SNP_139    0.5907398873  0.5500762222  9.037901e-02
SNP_140    0.9512522645  0.9601581008 -2.235294e-01
SNP_141    1.0000000000  1.0000000000  1.000000e+00
SNP_142    0.9505703422  0.9125475285  4.347826e-01
SNP_143    1.0000000000  1.0000000000  1.000000e+00
SNP_144    0.9381432671  0.9690716335 -1.000000e+00
SNP_145    0.9381432671  0.9690716335 -1.000000e+00
SNP_146    1.0000000000  1.0000000000  1.000000e+00
SNP_147    1.0000000000  1.0000000000  1.000000e+00
SNP_148    1.0000000000  1.0000000000  1.000000e+00
SNP_149    1.0000000000  1.0000000000  1.000000e+00
SNP_150    1.0000000000  1.0000000000  1.000000e+00
SNP_151    0.1368577940  0.2759615103 -1.921220e-01
SNP_152    0.6556609437 -0.0649004149  6.766467e-01
SNP_153    0.8059244472  0.7457250859  2.367491e-01
SNP_154    0.3045957625  0.4091985347 -1.770523e-01
SNP_155    0.9452004349  0.5178368931  8.863464e-01
SNP_156    0.3680460769  0.0127345348  3.598946e-01
SNP_157    0.3634001973  0.3668932003 -5.517241e-03
SNP_158    0.9676386666  0.9709274200 -1.131222e-01
SNP_159    0.0186250198  0.2356445276 -2.839248e-01
SNP_160    0.8272786718  0.4857415535  6.641352e-01
SNP_161    1.0000000000  1.0000000000  1.000000e+00
SNP_162    1.0000000000  1.0000000000  1.000000e+00
SNP_163    1.0000000000  1.0000000000  1.000000e+00
SNP_164    0.9680741504  0.9713186473 -1.131222e-01
SNP_165    1.0000000000  1.0000000000  1.000000e+00
SNP_166    1.0000000000  1.0000000000  1.000000e+00
SNP_167    1.0000000000  1.0000000000  1.000000e+00
SNP_168    1.0000000000  1.0000000000  1.000000e+00
SNP_169    0.7531451323  0.2650172270  6.641352e-01
SNP_170    1.0000000000  1.0000000000  1.000000e+00
SNP_171    0.7808234359  0.6219064488  4.203113e-01
SNP_172    0.7246729317  0.6046124195  3.036527e-01
SNP_173    0.9710965006  0.9697877325  4.331910e-02
SNP_174    0.7711888549  0.5780896533  4.576783e-01
SNP_175    1.0000000000  1.0000000000           NaN
SNP_176    1.0000000000  1.0000000000  1.000000e+00
SNP_177    1.0000000000  1.0000000000           NaN
SNP_178    1.0000000000  1.0000000000  1.000000e+00
SNP_179    1.0000000000  1.0000000000           NaN
SNP_180    1.0000000000  0.8446128858  1.000000e+00
SNP_181    1.0000000000  1.0000000000  1.000000e+00
SNP_182    0.4793349584  0.3641410995  1.811626e-01
SNP_183    1.0000000000  1.0000000000           NaN
SNP_184    1.0000000000  1.0000000000  1.000000e+00
SNP_185    0.7252946957  0.4796235105  4.721028e-01
SNP_186    0.4862924282  0.7431462141 -1.000000e+00
SNP_187    0.4664159458  0.4304749039  6.310704e-02
SNP_188    0.7864817584  0.7938295521 -3.563941e-02
SNP_189    0.6038121784  0.4182294620  3.189964e-01
SNP_190    0.0506345187  0.3874094316 -5.497553e-01
SNP_191    0.7895119248  0.6693547645  3.634021e-01
SNP_192    0.8503087365  0.8484906239  1.200000e-02
SNP_193    0.7122844756  0.7623726033 -2.107843e-01
SNP_194   -0.0948363241 -0.0726306491 -2.070207e-02
SNP_195    1.0000000000  1.0000000000  1.000000e+00
SNP_196    0.0753805738  0.3225384589 -3.648294e-01
SNP_197    0.6306834516  0.7378308452 -4.086957e-01
SNP_198    1.0000000000  1.0000000000  1.000000e+00
SNP_199    1.0000000000  1.0000000000  1.000000e+00
SNP_200    0.0064350064  0.2136642765 -2.635379e-01
SNP_201    0.6140107219  0.6247326463 -2.857143e-02
SNP_202    0.5309390839  0.3715142129  2.536650e-01
SNP_203    0.5309390839  0.3715142129  2.536650e-01
SNP_204    0.5309390839  0.3715142129  2.536650e-01
SNP_205    0.6692212897  0.5618786394  2.450067e-01
SNP_206    0.8723333885  0.8741571972 -1.449275e-02
SNP_207    0.2323322219  0.3219637619 -1.321928e-01
SNP_208    0.9237794608  0.9423460024 -3.220339e-01
SNP_209    1.0000000000  1.0000000000           NaN
SNP_210    0.2626895105  0.0741081737  2.036753e-01
SNP_211    0.8997858043  0.6608134915  7.045455e-01
SNP_212    1.0000000000  1.0000000000           NaN
SNP_213    1.0000000000  1.0000000000  1.000000e+00
SNP_214    0.8959009211  0.8995584563 -3.641457e-02
SNP_215    1.0000000000  1.0000000000  1.000000e+00
SNP_216    0.5420501365  0.6114031206 -1.784703e-01
SNP_217    0.5961842365  0.5790857312  4.062230e-02
SNP_218    0.7350980239  0.4090648226  5.517241e-01
SNP_219    1.0000000000  1.0000000000  1.000000e+00
SNP_220    1.0000000000  1.0000000000           NaN
SNP_221    1.0000000000  1.0000000000  1.000000e+00
SNP_222    1.0000000000  1.0000000000           NaN
SNP_223    1.0000000000  1.0000000000  1.000000e+00
SNP_224    0.9117890670  0.9006668189  1.119691e-01
SNP_225    1.0000000000  1.0000000000           NaN
SNP_226    0.5121853294  0.4977315614  2.877698e-02
SNP_227    0.5858307609  0.3026120901  4.061135e-01
SNP_228    1.0000000000  1.0000000000  1.000000e+00
SNP_229    0.2446235419  0.4184800283 -2.989691e-01
SNP_230   -0.0174488362  0.2441808646 -3.461538e-01
SNP_231    0.7501220765  0.3895839297  5.906433e-01
SNP_232   -0.0303353984  0.0036317026 -3.409091e-02
SNP_233   -0.0156731230 -0.0011635069 -1.449275e-02
SNP_234    0.9625157335  0.9659071671 -9.947644e-02
SNP_235    0.8359612814  0.8508028797 -9.947644e-02
SNP_236    1.0000000000  1.0000000000           NaN
SNP_237    1.0000000000  0.7333944114  1.000000e+00
SNP_238    1.0000000000  1.0000000000  1.000000e+00
SNP_239    1.0000000000  1.0000000000  1.000000e+00
SNP_240    1.0000000000  1.0000000000           NaN
SNP_241    1.0000000000  1.0000000000  1.000000e+00
SNP_242    0.9613221483  0.9648215730 -9.947644e-02
SNP_243    0.9613221483  0.9648215730 -9.947644e-02
SNP_244    1.0000000000  1.0000000000           NaN
SNP_245    1.0000000000  1.0000000000  1.000000e+00
SNP_246    0.8815063871  0.8922272378 -9.947644e-02
SNP_247    1.0000000000  1.0000000000  1.000000e+00
SNP_248   -0.0297464739  0.1527747265 -2.154341e-01
SNP_249    1.0000000000  1.0000000000  1.000000e+00
SNP_250    0.1994245310  0.3559011391 -2.429388e-01
SNP_251    1.0000000000  1.0000000000  1.000000e+00
SNP_252    0.3977209944  0.4945056209 -1.914653e-01
SNP_253    0.0704632913  0.2723664398 -2.774792e-01
SNP_254    0.7953948804  0.7197655629  2.698787e-01
SNP_255    1.0000000000  1.0000000000  1.000000e+00
SNP_256    0.6744338771  0.5816021900  2.218742e-01
SNP_257    1.0000000000  1.0000000000  1.000000e+00
SNP_258    0.1074634697  0.3396488987 -3.516091e-01
SNP_259    1.0000000000  0.8982315649  1.000000e+00
SNP_260    1.0000000000  1.0000000000  1.000000e+00
SNP_261    1.0000000000  1.0000000000  1.000000e+00
SNP_262    0.8354858549  0.8519372694 -1.111111e-01
SNP_263    1.0000000000  1.0000000000           NaN
SNP_264    1.0000000000  1.0000000000           NaN
SNP_265    1.0000000000  1.0000000000  1.000000e+00
SNP_266    1.0000000000  1.0000000000           NaN
SNP_267    1.0000000000  1.0000000000           NaN
SNP_268    1.0000000000  1.0000000000  1.000000e+00
SNP_269    1.0000000000  1.0000000000           NaN
SNP_270    0.8714405863  0.8971524691 -2.500000e-01
SNP_271    1.0000000000  1.0000000000  1.000000e+00
SNP_272    0.8714405863  0.8971524691 -2.500000e-01
SNP_273    1.0000000000  1.0000000000  1.000000e+00
SNP_274    0.8324045793  0.8463708643 -9.090909e-02
SNP_275    0.3164573728  0.2523752515  8.571429e-02
SNP_276    1.0000000000  1.0000000000  1.000000e+00
SNP_277    1.0000000000  1.0000000000           NaN
SNP_278    0.5666845004  0.6704996722 -3.150685e-01
SNP_279    1.0000000000  1.0000000000  1.000000e+00
SNP_280    1.0000000000  0.8041843740  1.000000e+00
SNP_281    0.2172203435  0.4047613028 -3.150685e-01
SNP_282    0.6054487737  0.6109286518 -1.408451e-02
SNP_283    0.5666845004  0.6704996722 -3.150685e-01
SNP_284    1.0000000000  1.0000000000  1.000000e+00
SNP_285    0.6054487737  0.6109286518 -1.408451e-02
SNP_286    0.9535867669  0.9502024687  6.796117e-02
SNP_287    1.0000000000  1.0000000000  1.000000e+00
SNP_288    0.8324045793  0.8463708643 -9.090909e-02
SNP_289    1.0000000000  1.0000000000  1.000000e+00
SNP_290    1.0000000000  1.0000000000  1.000000e+00
SNP_291    0.9535867669  0.9502024687  6.796117e-02
SNP_292    0.9466613934  0.9631709621 -4.482759e-01
SNP_293    0.7065095968  0.6184624759  2.307692e-01
SNP_294    0.9466613934  0.9631709621 -4.482759e-01
SNP_295    1.0000000000  1.0000000000  1.000000e+00
SNP_296    1.0000000000  1.0000000000  1.000000e+00
SNP_297    0.7065095968  0.6184624759  2.307692e-01
SNP_298    1.0000000000  1.0000000000  1.000000e+00
SNP_299    1.0000000000  1.0000000000           NaN
SNP_300    0.5637642278  0.4043703879  2.676056e-01
SNP_301    0.8275841410  0.7808881792  2.131148e-01
SNP_302    0.8462330221  0.7102083878  4.693878e-01
SNP_303    0.5963547700  0.2343814160  4.727855e-01
SNP_304    0.5289014276  0.5224302934  1.355014e-02
SNP_305    0.7102278526  0.6360517287  2.038095e-01
SNP_306    0.9055093624  0.9091436177 -4.000000e-02
SNP_307    0.2302952889  0.3043053573 -1.063830e-01
SNP_308    0.2636739882  0.0862999084  1.941272e-01
SNP_309    1.0000000000  1.0000000000  1.000000e+00
SNP_310    1.0000000000  1.0000000000  1.000000e+00
SNP_311    0.4514868748  0.4957218043 -8.771930e-02
SNP_312    0.7156065668  0.7504676973 -1.397059e-01
SNP_313    0.7156065668  0.7504676973 -1.397059e-01
SNP_314    1.0000000000  1.0000000000  1.000000e+00
SNP_315    0.6067415730  0.7247191011 -4.285714e-01
SNP_316    1.0000000000  1.0000000000  1.000000e+00
SNP_317    1.0000000000  1.0000000000  1.000000e+00
SNP_318    1.0000000000  1.0000000000  1.000000e+00
SNP_319    0.4152046784  0.5438596491 -2.820513e-01
SNP_320    1.0000000000  1.0000000000  1.000000e+00
SNP_321    1.0000000000  1.0000000000  1.000000e+00
SNP_322    1.0000000000  1.0000000000  1.000000e+00
SNP_323    0.7014050191  0.3221798230  5.594776e-01
SNP_324    1.0000000000  1.0000000000  1.000000e+00
SNP_325    0.9687287997  0.9666991785  6.094808e-02
SNP_326    1.0000000000  1.0000000000  1.000000e+00
SNP_327    0.6551974285  0.3938879156  4.311241e-01
SNP_328    1.0000000000  1.0000000000           NaN
SNP_329    0.6776974646  0.5184223453  3.307361e-01
SNP_330    1.0000000000  1.0000000000  1.000000e+00
SNP_331    1.0000000000  1.0000000000  1.000000e+00
SNP_332    0.0765299951  0.3393231736 -3.977636e-01
SNP_333    1.0000000000  1.0000000000           NaN
 [ reached getOption("max.print") -- omitted 41170 rows ]
 ````
 
#Are these values significant? This question can be addressed using the G-statistic test
#[3]; it is implemented for genind objects and produces a randtest object (package ade4).
#nsim is used to simulate a single or "n" prospective clinical trial(s) using the PK data and then link them to toxicity under a specified dose-toxicity configuration. The objective is to determine the maximum tolerated dose (MTD).
````
Gtest <- gstat.randtest(file,nsim=99)
Gtest
plot(Gtest)
````

````
> Gtest <- gstat.randtest(file,nsim=99)
There were 50 or more warnings (use warnings() to see the first 50)
> Gtest
Monte-Carlo test
Call: gstat.randtest(x = file, nsim = 99)

Observation: 0 

Based on 99 replicates
Simulated p-value: 1 
Alternative hypothesis: greater 

    Std.Obs Expectation    Variance 
        NaN           0           0 
````


##How to subset your data (in terminal!!)

#Create a plink file from the vcf file.

````
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools --vcf filename.final.recode.vcf --plink --out CHN.229.Neutral.plink

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf filename.final.recode.vcf
	--out CHN.229.Neutral.plink
	--plink

After filtering, kept 64 out of 64 Individuals
Writing PLINK PED and MAP files ... 

Unrecognized values used for CHROM: locus_7 - Replacing with 0.
....

````

#Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names (in r)

````
CHN.Neutral.loci.names.map <- read.table("/Users/aminachabach/Project/vcftools_0.1.13/CHN.229.Neutral.plink.map", header=F)
CHN.Neutral.loci.names <-  CHN.Neutral.loci.names.map$V2
CHN.Neutral.loci.names <- as.data.frame(CHN.Neutral.loci.names)
CHN.1000.loci.names <- CHN.Neutral.loci.names[sample(nrow(CHN.Neutral.loci.names),1000),]
CHN.1000.loci.names <- as.data.frame(CHN.1000.loci.names)
summary(CHN.1000.loci.names)
write.table(CHN.1000.loci.names, "CHN.1000.loci.names", quote=F, row.names=F, col.names=F)

>summary(CHN.1000.loci.names)
      CHN.1000.loci.names
 locus_10016:62 :  1     
 locus_10043:246:  1     
 locus_10066:2  :  1     
 locus_10127:16 :  1     
 locus_10127:225:  1     
 locus_10204:37 :  1     
 (Other)        :994 


````

#Subset the vcf file to get 1000 loci (in terminal)

````
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools --vcf filename.final.recode.vcf --snps CHN.1000.loci.names --recode --recode-INFO-all --out CHN.1000Neutral

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf filename.final.recode.vcf
	--recode-INFO-all
	--out CHN.1000Neutral
	--recode
	--snps CHN.1000.loci.names

After filtering, kept 64 out of 64 Individuals
Outputting VCF file...
After filtering, kept 0 out of a possible 41503 Sites
No data left for analysis!
Run Time = 1.00 seconds

````


## Hardy-Weinberg filter

# Installing plink
```
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/plink-1.07-mac-intel/plink --noweb

@----------------------------------------------------------@
|        PLINK!       |     v1.07      |   10/Aug/2009     |
|----------------------------------------------------------|
|  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
|----------------------------------------------------------|
|  For documentation, citation & bug-report instructions:  |
|        http://pngu.mgh.harvard.edu/purcell/plink/        |
@----------------------------------------------------------@

Skipping web check... [ --noweb ] 
Writing this text to log file [ plink.log ]
Analysis started: Thu Oct 18 14:32:58 2018

Options in effect:
	--noweb

Before frequency and genotyping pruning, there are 0 SNPs
0 founders and 0 non-founders found
0 SNPs failed missingness test ( GENO > 1 )
0 SNPs failed frequency test ( MAF < 0 )
After frequency and genotyping pruning, there are 0 SNPs

ERROR: Stopping as there are no SNPs left for analysis
```
#Dw about the error (--noweb not necessary either!)

#Create plink file (with the filename.removed so only 63 individuals)

````
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools --vcf filename.removed.recode.vcf --plink --out CHN.229.Neutral.plink

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf filename.removed.recode.vcf
	--out CHN.229.Neutral.plink
	--plink

After filtering, kept 63 out of 63 Individuals
Writing PLINK PED and MAP files ... 
After filtering, kept 41503 out of a possible 41503 Sites
Run Time = 3.00 seconds
````

# To exclude markers that failure the Hardy-Weinberg test at a specified significance threshold, use the option:
#The following output will appear in the console window and in plink.log, detailing how many SNPs failed the Hardy-Weinberg test, for the sample as a whole, and (when PLINK has detected a disease phenotype) for cases and controls separately:
#--hardy: Reports a p-value for each site from a Hardy-Weinberg Equilibrium test (as defined by Wigginton, Cutler and Abecasis (2005)). The resulting file (with suffix ".hwe") also contains the Observed numbers of Homozygotes and Heterozygotes and the corresponding Expected numbers under HWE.

````
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/plink-1.07-mac-intel/plink --file CHN.229.Neutral.plink --hwe 0.001

@----------------------------------------------------------@
|        PLINK!       |     v1.07      |   10/Aug/2009     |
|----------------------------------------------------------|
|  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
|----------------------------------------------------------|
|  For documentation, citation & bug-report instructions:  |
|        http://pngu.mgh.harvard.edu/purcell/plink/        |
@----------------------------------------------------------@

Web-based version check ( --noweb to skip )
Recent cached web-check found...Problem connecting to web

Writing this text to log file [ plink.log ]
Analysis started: Thu Oct 18 15:35:07 2018

Options in effect:
	--file CHN.229.Neutral.plink
	--hwe 0.001

41081 (of 41081) markers to be included from [ CHN.229.Neutral.plink.map ]
Warning, found 63 individuals with ambiguous sex codes
These individuals will be set to missing ( or use --allow-no-sex )
Writing list of these individuals to [ plink.nosex ]
63 individuals read from [ CHN.229.Neutral.plink.ped ] 
0 individuals with nonmissing phenotypes
Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
Missing phenotype value is also -9
0 cases, 0 controls and 63 missing
0 males, 0 females, and 63 of unspecified sex
Before frequency and genotyping pruning, there are 41081 SNPs
63 founders and 0 non-founders found
27447 markers to be excluded based on HWE test ( p <= 0.001 )
	0 markers failed HWE test in cases
	0 markers failed HWE test in controls
Total genotyping rate in remaining individuals is 0.753766
0 SNPs failed missingness test ( GENO > 1 )
0 SNPs failed frequency test ( MAF < 0 )
After frequency and genotyping pruning, there are 13634 SNPs
After filtering, 0 cases, 0 controls and 63 missing
After filtering, 0 males, 0 females, and 63 of unspecified sex

Analysis finished: Thu Oct 18 15:35:12 2018

````
## Heterozygosity excess filter
#--het: Calculates a measure of heterozygosity on a per-individual basis. Specfically, the inbreeding coefficient, F, is estimated for each individual using a method of moments. The resulting file has the suffix ".het".


````
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/plink-1.07-mac-intel/plink --file CHN.229.Neutral.plink --het

@----------------------------------------------------------@
|        PLINK!       |     v1.07      |   10/Aug/2009     |
|----------------------------------------------------------|
|  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
|----------------------------------------------------------|
|  For documentation, citation & bug-report instructions:  |
|        http://pngu.mgh.harvard.edu/purcell/plink/        |
@----------------------------------------------------------@

Web-based version check ( --noweb to skip )
Recent cached web-check found...Problem connecting to web

Writing this text to log file [ plink.log ]
Analysis started: Thu Oct 18 15:48:17 2018

Options in effect:
	--file CHN.229.Neutral.plink
	--het

41081 (of 41081) markers to be included from [ CHN.229.Neutral.plink.map ]
Warning, found 63 individuals with ambiguous sex codes
These individuals will be set to missing ( or use --allow-no-sex )
Writing list of these individuals to [ plink.nosex ]
63 individuals read from [ CHN.229.Neutral.plink.ped ] 
0 individuals with nonmissing phenotypes
Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
Missing phenotype value is also -9
0 cases, 0 controls and 63 missing
0 males, 0 females, and 63 of unspecified sex
Before frequency and genotyping pruning, there are 41081 SNPs
63 founders and 0 non-founders found
Total genotyping rate in remaining individuals is 0.753766
0 SNPs failed missingness test ( GENO > 1 )
0 SNPs failed frequency test ( MAF < 0 )
After frequency and genotyping pruning, there are 41081 SNPs
After filtering, 0 cases, 0 controls and 63 missing
After filtering, 0 males, 0 females, and 63 of unspecified sex
Converting data to Individual-major format
Writing individual heterozygosity information to [ plink.het ] 

Analysis finished: Thu Oct 18 15:48:22 2018
````




### (Tips and Tricks)

#sed = stream editor ...
#rm = remove


## Project ideas

#how do filters affect the analysis?
#what analysis will we be using?
#data that is not sensitive to analysis
#What data is available online? --> project question: given demographic history of population what inbreeding do we expect?
#simulated data set (can have pedigree)
#genomic data: look for paper on captive bred populations
#paper might have link to genomic data

#taxonomically different? compare if same demographic history













