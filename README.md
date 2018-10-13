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

# To find path if using the same (did this work?????)
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

````
???
`````






