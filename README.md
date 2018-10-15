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



## Summary stats in R
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


