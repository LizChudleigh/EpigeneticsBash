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

Here we expect no results if the values conincide
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


