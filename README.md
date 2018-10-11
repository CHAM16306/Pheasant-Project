# Pheasant-Project


# Step1: Look at fastq data (in terminal)

Last login: Thu Oct 11 12:38:35 on ttys000

#So now when you open terminal, navigate to that folder using this command: cd /User/projects/PheasantProject
```eduroam-140-210:~ aminachabach$ cd /Users/aminachabach/Project/Fastq```

#where /User/projects/PheasantProject is whatever your file path is. 
You can check where your terminal is with the following command: pwd
eduroam-140-210:Fastq aminachabach$ pwd
/Users/aminachabach/Project/Fastq

#You can then see if your PHE079.fastq.gz files in there with: ls
eduroam-140-210:Fastq aminachabach$ ls
Pheasant1.fastq	Pheasant2.fastq

#The following commands assume that you're unzipped the files. 
#This command shows you the first few lines of a document
head filename 
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

#This opens a page at a time. You can scroll through using spacebar, and quit this view with q: less filename
eduroam-140-210:Fastq aminachabach$ less Pheasant1.fastq


#shows the last few lines of a file: tail filename
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




#Step2: Filter the data in linux

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




#vcf Installation

eduroam-140-210:vcftools_0.1.13 aminachabach$ cd/Users/aminachabach/Project/vcftools_0.1.13
-bash: cd/Users/aminachabach/Project/vcftools_0.1.13: No such file or directory

eduroam-140-210:vcftools_0.1.13 aminachabach$ pwd
/Users/aminachabach/Project/vcftools_0.1.13

eduroam-140-210:vcftools_0.1.13 aminachabach$ ls
Makefile	README.md	README.txt	cpp		examples	perl

eduroam-140-210:vcftools_0.1.13 aminachabach$ less Makefile

#Installation
eduroam-140-210:vcftools_0.1.13 aminachabach$ make
g++ -c -O2 -D_FILE_OFFSET_BITS=64  vcftools.cpp -o vcftools.o


#To find path if using the same
eduroam-140-210:vcftools_0.1.13 aminachabach$ echo $PATH
/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin

#eduroam-140-210:vcftools_0.1.13 aminachabach$ ls ~/
Applications	Assignmnet	Desktop		Downloads	Movies		Pictures	Project
Assignment	BirthdayStick	Documents	Library		Music		Practical 2	Public


#Bin vcf tools (= path/bin/vcftools)
eduroam-140-210:vcftools_0.1.13 aminachabach$ /Users/aminachabach/Project/vcftools_0.1.13/bin/vcftools 
VCFtools (v0.1.13)
© Adam Auton and Anthony Marcketta 2009

Process Variant Call Format files

For a list of options, please go to:
	https://vcftools.github.io/examples.html

Questions, comments, and suggestions should be emailed to:
	vcftools-help@lists.sourceforge.net

eduroam-140-210:vcftools_0.1.13 aminachabach$ 

