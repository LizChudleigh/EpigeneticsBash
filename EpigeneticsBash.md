# Epigenetics Bash Assignment
## Elizabeth Chudleigh
## 4. EN‐TEx ATAC‐seq data: downstream analyses
#4.1)Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.

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
#4.2)Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections.  Make sure your md5sum values coincide with the ones provided by ENCODE.

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

#4.3)For each tissue, run an intersection analysis using BEDTools: report 1) the number of peaks that intersect promoter regions, 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions).

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
  bedtools intersect -b /home/emchudleigh/epigenomics_uvic/ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -a /home/emchudleigh/epigenomics_uvic/ATAC-seq/data/bed.files/"$filename".bed -u > /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis/genes.with.peaks."$tissue".ATACTSS.txt
done
```

Retrieve genes with peaks of ATAC at the promoter region – This is our list of gene bodies to NOT overlap: gencode.v24.protein.coding.gene.body.bed
For ATAC-seq peaks that fall OUTSIDE of the genes (DISTAL)
Unfiltered

```bash
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -b /home/emchudleigh/epigenomics_uvic/ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed -a /home/emchudleigh/epigenomics_uvic/ATAC-seq/data/bed.files/"$filename".bed -v > /home/emchudleigh/epigenomics_uvic/ATAC-seq/analyses/peaks.analysis/genes.with.peaks."$tissue".ATACDistal.txt
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

##Task 1: Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results.
In root@bf2c267898b8:/home/emchudleigh/epigenomics_uvic#
```bash
mkdir regulatory_element
install.dependecies.txt  nextflow  regulatory_element  test
```
##Task 2: Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. You will get a list of candidate distal regulatory elements for each tissue. How many are they?

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
List of candidate distal regulatory elements for each tissue. How many are they? Answer below:
```bash
wc -l regulatory_element/bigBed.files/*.bigBed
   15883 regulatory_element/bigBed.files/ENCFF724ZOF.bigBed
   12362 regulatory_element/bigBed.files/ENCFF844XRN.bigBed
    8668 regulatory_element/bigBed.files/ENCFF872UHN.bigBed
    9162 regulatory_element/bigBed.files/ENCFF977LBD.bigBed
   46075 total
```   

ATAC-seq peaks that intersect BOTH of the H3K27ac AND H3Kme1 

First we convert bigBed to bed files
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

##Task 3: Focus on regulatory elements that are located on chromosome 1 (hint: to parse a file based on the value of a specific column, have a look at what we did here), and generate a file regulatory.elements.starts.tsv that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.





