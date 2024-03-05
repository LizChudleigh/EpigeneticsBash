# Epigenetics Bash Assignment
## Elizabeth Chudleigh
## 4. EN‐TEx ATAC‐seq data: downstream analyses
### 4.1)Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.

Starting Docker
```bash
cd epigenomics_uvic
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
mkdir ATAC-seq
cd ATAC-seq
```

Making the directories within the ATAC-seq folder
```bash
mkdir analyses
mkdir data
mkdir annotation
mkdir atac-seq
mkdir analyses/peaks.analysis
mkdri data/bed.files
mkdir data/bed.files
mkdir data/bigBed.files
```
### 4.2)Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections.  Make sure your md5sum values coincide with the ones provided by ENCODE.

Calling the code download.metadata.sh from inside the bin folder the metadata for ATAC-Seq for Sigmoid Colon and Stomach for ENCDO451RUA
```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=ATAC-seq&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&type=Experiment"
```

Getting the data on bigBed_narrowPeak and pseudoreplicated_peak data for humans in the GRCH38 assembly from the ATAC-Seq for Sigmoid Colon and Stomach for ENCDO451RUA
```bash
grep -F "bigBed_narrowPeak" metadata.tsv |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```
Making sure the md5sum values coincide with the ones provided by ENCODE.
In root@5ce471a4d155:/home/emchudleigh/epigenomics_uvic/ATAC-seq# 
```bash
head analyses/bigBed.peaks.ids.txt
```
Results
ENCFF287UHP     sigmoid_colon
ENCFF762IFP     stomach
Can manually check against the annotation file sections in https://www.encodeproject.org/experiments/ENCSR851SBY/ and https://www.encodeproject.org/experiments/ENCSR086OGH/

Here we expect no results if the values coincide 
```bash
for file_type in bigBed; do
  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
```
No results so the values of the MD5 from ENCODE and our md5sum are the same

### 4.3)For each tissue, run an intersection analysis using BEDTools: report 1) the number of peaks that intersect promoter regions, 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions).

Convert bigBed files with peaks for ATAC to regular bed files. This is our list of peaks for the two tissue types ENCFF762IFP for stomach and ENCFF287UHP for sigmoid colon
```bash
cut -f1 analyses/bigBed.peaks.ids.txt |\
 while read filename; do
   bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
 done
```

Retrieve genes with peaks of ATAC at the promoter region – This is our list of TSS site: gencode.v24.protein.coding.non.redundant.TSS.bed

For the overlap with TSS and ATAC-seq peaks 
```bash
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -b /home/emchudleigh/epigenomics_uvic/ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -a /home/emchudleigh/epigenomics_uvic/ATAC-seq/data/bed.files/"$filename".bed -u > /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis/genes.with.peaks."$tissue".ATACTSS.bed
done
```

Retrieve genes with peaks of ATAC at the promoter region – This is our list of gene bodies to NOT overlap: gencode.v24.protein.coding.gene.body.bed
For ATAC-seq peaks that fall OUTSIDE of the genes (DISTAL)
Unfiltered

```bash
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -b /home/emchudleigh/epigenomics_uvic/ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed -a /home/emchudleigh/epigenomics_uvic/ATAC-seq/data/bed.files/"$filename".bed -v > /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis/genes.with.peaks."$tissue".ATACDistal.bed
done
```

Here we get our 4 numbers for task 4 
2 for each tissue type
2 files for each overlapp type, 1 for peaks that overlap with TSS (-u) (ATACTSS), 1  for peaks that do not #(-v) overlap with gene body (ATACDistal)
```bash
wc -l analyses/peaks.analysis/*
```
Gives us
   37035 analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.ATACDistal.txt
   47871 analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.ATACTSS.txt
   34537 analyses/peaks.analysis/genes.with.peaks.stomach.ATACDistal.txt
   44749 analyses/peaks.analysis/genes.with.peaks.stomach.ATACTSS.txt
  164192 total

1) the number of peaks that intersect promoter regions (Unfiltered)
  37035 Sigmoid Colon
  34537 Stomach
  Total 71572

2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions) (Unfiltered)
  47871 Sigmoid Colon
  44749 Stomach
  Total 92620

# 5. Distal regulatory activity
Study distal regulatory activity
From section 4., you should have obtained a set of ATAC-seq peaks in stomach and sigmoid_colon that lie outside gene coordinates We will use these peaks as a starting point to build a catalogue of distal regulatory regions.

### Task 1: Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results.
In root@bf2c267898b8:/home/emchudleigh/epigenomics_uvic#
```bash
mkdir regulatory_element
install.dependecies.txt  nextflow  regulatory_element  test
```
### Task 2: Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. You will get a list of candidate distal regulatory elements for each tissue. How many are they?

Here we use the metadata we already have from the ChIP-seq analysis to index and retrieve peaks for each marker.

For H3K27ac
In root@bf2c267898b8:/home/emchudleigh/epigenomics_uvic#

```bash
grep -F H3K27ac ChIP-seq/metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > regulatory_element/bigBed.H3K27ac.peaks.ids.txt

cut -f1 regulatory_element/bigBed.H3K27ac.peaks.ids.txt|\
while read filename; do
  wget -P regulatory_element/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

```bash
head regulatory_element/bigBed.H3K27ac.peaks.ids.txt

ENCFF872UHN     sigmoid_colon   H3K27ac-human
ENCFF977LBD     stomach H3K27ac-human
```
For H3K4me1

```bash
grep -F H3K4me1 ChIP-seq/metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > regulatory_element/bigBed.H3K4me1.peaks.ids.txt

cut -f1 regulatory_element/bigBed.H3K4me1.peaks.ids.txt |\
while read filename; do
  wget -P regulatory_element/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

```bash
head regulatory_element/bigBed.H3K4me1.peaks.ids.txt

ENCFF724ZOF     sigmoid_colon   H3K4me1-human
ENCFF844XRN     stomach H3K4me1-human
```
Here we expect no results if the values coincide 

```bash
for file_type in bigBed; do

  ../bin/selectRows.sh <(cut -f1 analyses/*"$file_type".peaks.ids.txt) ../ChIP-seq/metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  awk '2!=3' data/"$file_type".files/md5sum.txt

done
```
No results so the values of the MD5 from ENCODE and our md5sum are the same

Converting the bigBed to bed files
```bash
cd regulatory_element
mkdir bed.files

cut -f1 regulatory_element/bigBed.H3K4me1.peaks.ids.txt |\
while read filename; do
bigBedToBed regulatory_element/bigBed.files/"$filename".bigBed regulatory_element/bed.files/"$filename".bed
done

cut -f1 regulatory_element/bigBed.H3K27ac.peaks.ids.txt |\
while read filename; do
bigBedToBed regulatory_element/bigBed.files/"$filename".bigBed regulatory_element/bed.files/"$filename".bed
done
```
For peaks which intersect promoter regions, we use the command bedtools intersect with our bed files as first argument, the protein coding non redundant TSS bed file as second argument (promoter region), we use parameter -u to look for unique peaks overlapping at least one TSS.

ATAC-seq peaks that intersect BOTH of the H3K27ac AND H3Kme1 

First we get a file where both peaks are present

For Sigmoid Colon
```bash
#from within regulatroy_element/data/bed.files directory
bedtools intersect -a ENCFF724ZOF.bed -b ENCFF872UHN.bed -u > /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis/H3Kme1_H3K27ac_overlap_sigmoid_colon.bed 
```
For Stomach 
```bash
#from within regulatroy_element/data/bed.files directory
bedtools intersect -a ENCFF844XRN.bed -b ENCFF977LBD.bed -u > /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis/H3Kme1_H3K27ac_overlap_stomach.bed
```
```bash
# quick check on row count in root@f9044d112704:/home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis#
wc -l *
   53615 H3Kme1_H3K27ac_overlap_sigmoid_colon.bed
   40992 H3Kme1_H3K27ac_overlap_stomach.bed
   37035 genes.with.peaks.sigmoid_colon.ATACDistal.bed
   47871 genes.with.peaks.sigmoid_colon.ATACTSS.bed
   34537 genes.with.peaks.stomach.ATACDistal.bed
   44749 genes.with.peaks.stomach.ATACTSS.bed
  258799 total
# To compare to the orginal bed files in regulatroy_element/data/bed.files directory
   97950 ENCFF724ZOF.bed
   68664 ENCFF844XRN.bed
   56661 ENCFF872UHN.bed
   57121 ENCFF977LBD.bed
  280396 total
```

Now we intersect the distal genes from the ATAC-seq peak analysis for each tissue with its combined H3Kme1_H3K27ac_overlap file
Sigmoid Colon
```bash
# from within /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis
bedtools intersect -a genes.with.peaks.sigmoid_colon.ATACDistal.bed  -b H3Kme1_H3K27ac_overlap_sigmoid_colon.bed -u > combined.peaks.sigmoid_colon.bed
```
Stomach
```bash
# from within /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis
bedtools intersect -a genes.with.peaks.stomach.ATACDistal.bed  -b H3Kme1_H3K27ac_overlap_stomach.bed -u > combined.peaks.stomach.bed
```
Attain counts of overlap unfiltered
14333 combined.peaks.sigmoid_colon.bed
8043 combined.peaks.stomach.bed

Remove Repeated coordinates
```bash
for file in combined.peaks.sigmoid_colon.bed combined.peaks.stomach.bed; \
  do sort -k 2,2 -k 3,3 "$file"| awk 'BEGIN{OFS="\t"} !seen[$1,$2,$3]++ {print $1, $2, $3, $4}' | sort -u | wc -l; done
```
8749 Filtered for Sigmoid Colon 
5189 Filtered for Stomach

### Task 3: Focus on regulatory elements that are located on chromosome 1 (hint: to parse a file based on the value of a specific column, have a look at what we did here), and generate a file regulatory.elements.starts.tsv that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.

Single File
```bash
grep -w 'chr1' combined.peaks."$tissue".bed | sort -k1,1 -k2,2n | awk '{print $4 "\t" $2}' > regulatory.elements.starts.tsv 
```
992 regulatory.elements.starts.tsv

Creating the regulatory.elements.starts.tsv unfiltered per tissue
```bash
#still within /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis#
for tissue in sigmoid_colon stomach; do
  grep -w 'chr1' combined.peaks."$tissue".bed | sort -k1,1 -k2,2n | awk '{print $4 "\t" $2}' > regulatory.elements.starts_"$tissue".tsv 
done
```
1537 regulatory.elements.starts_sigmoid_colon.tsv
992 regulatory.elements.starts_stomach.tsv

Filtered (This command removes duplicate entries based on the combination of chromosome, start position, and end position. It ensures that only unique peaks are retained. It also sets the output field separator (OFS) to a tab character.)
```bash
for tissue in sigmoid_colon stomach; do
  grep -w 'chr1' combined.peaks."$tissue".bed | sort -k 2,2 -k 3,3 | awk 'BEGIN{OFS="\t"} !seen[$1,$2,$3]++ {print $1, $2, $3, $4}' | awk '{print $4 "\t" $2}' | sort -u > regulatory.elements.starts_"$tissue"_filtered.tsv
done
```
953 regulatory.elements.starts_sigmoid_colon_filtered.tsv
649 regulatory.elements.starts_stomach_filtered.tsv

### Task 4: Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here, prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3'). Use the command below as a starting poi
```bash
awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}'
```

```bash
# from within /home/emchudleigh/epigenomics_uvic/ATAC-seq/annotation
 grep -w 'chr1' /home/emchudleigh/epigenomics_uvic/ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed| awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}'> gene.starts.tsv 
```

### Task 5: Download or copy this python script inside the epigenomics_uvic/bin folder. Have a look at the help page of this script to understand how it works:

```bash
python ../bin/get.distance.py -h
```
```bash
#move to bin folder
wget https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/get.distance.py
```
Need to add the input and start 
```bash
python ../bin/get.distance.py -h
Usage: get.distance.py [options]

Options:
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
  -s START, --start=START
```

The code itself shows more with this Python script seeming to be designed to find the closest gene to a given enhancer start position within a certain maximum distance threshold (x).

### This script takes as input two distinct arguments: 1) --input corresponds to the file gene.starts.tsv (i.e. the file you generated in Task #4); 2) --start corresponds to the 5' coordinate of a regulatory element. Complete the python script so that for a given coordinate --start the script returns the closest gene, the start of the gene and the distance of the regulatory element.

Edit python script
```bash
# in epigenomics_uvic/bin
nano get.distance.py
```
New Python Script
```python
#!/usr/bin/env python


#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser


#*****************
# OPTION PARSING *
#*****************

parser = OptionParser()
parser.add_option("-i", "--input", dest="input")
parser.add_option("-s", "--start", dest="start")
options, args = parser.parse_args()

open_input = open(options.input)
enhancer_start = int(options.start)


#********
# BEGIN *
#********

x=1000000 # set maximum distance to 1 Mb
selectedGene="" # initialize the gene as empty
selectedGeneStart=0 # initialize the start coordinate of the gene as empty

for line in open_input.readlines(): # for each line in the input file
	gene, y = line.strip().split('\t') # split the line into two columns based on a tab 
	position = int(y) # define a variable called position that correspond to the integer of the start of the gene
	distance = abs(position - enhancer_start) # compute the absolute value of the difference between position and enhancer_start

	if distance <x: # if this absolute value is lower than x
		x = distance # this value will now be your current x
		selectedGene = gene # save gene as selectedGene
		selectedGeneStart = position # save position as selectedGeneStart
print "\t".join([selectedGene, str(selectedGeneStart), str(x)])
```
### To make sure your script is working fine, run the following command:

Test python script
```bash
root@f9044d112704:/home/emchudleigh/epigenomics_uvic/ATAC-seq/annotation# python /home/emchudleigh/epigenomics_uvic/bin/get.distance.py --input gene.starts.tsv --start 980000
ENSG00000187642.9       982093  2093
```
### Task 6. For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above. Use the command below as a starting point:

```bash
cat regulatory.elements.starts.tsv | while read element start; do 
   python ../bin/get.distance.py ... # to be completed by you; 
done > regulatoryElements.genes.distances.tsv
```
```bash
cd  /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis

for tissue in sigmoid_colon stomach; do
   cat regulatory.elements.starts_"$tissue".tsv | while read element start; do 
      python /home/emchudleigh/epigenomics_uvic/bin/get.distance.py -i /home/emchudleigh/epigenomics_uvic/ATAC-seq/annotation/gene.starts.tsv --start "$start"
   done >> regulatoryElements.genes.distances."$tissue".tsv
done
```



### Task 7: Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv

First we filter the file for duplicates
```bash
for tissue in sigmoid_colon stomach; do
    awk '!seen[$0]++' regulatoryElements.genes.distances."$tissue".tsv > filtered.regulatoryElements.genes.distances."$tissue".tsv
done
```

Using Rscript -e which according to Rscript --help "Expressions (one or more '-e <expr>') may be used *instead* of 'file' from within R"
This script first uses awk to remove duplicate lines from the file 
This Rscript gives the mean and median of column 3 which is the distance calculated by the get.distance.py script
```bash
for tissue in sigmoid_colon stomach; do
  Rscript -e "x <- read.table('regulatoryElements.genes.distances."$tissue".tsv', header = FALSE, sep='\t'); cat('${tissue} Mean:', mean(x[,3]), '\n'); cat('${tissue} Median:', median(x[,3]), '\n')";
done

#Output
sigmoid_colon Mean: 72064.72
sigmoid_colon Median: 35827
stomach Mean: 45503.29
stomach Median: 27735
```




