# Epigenetics Bash Assignment
## Elizabeth Chudleigh
## 4. EN‐TEx ATAC‐seq data: downstream analyses
#a)Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.

To download the metadata file for ATAC-Seq for Sigmoid and Colon for ENCDO451RUA
```bash
cd epigenomics_uvic
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
mkdir ATAC-seq
```
#b)Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. Hint: have a look at what we did here. Make sure your md5sum values coincide with the ones provided by ENCODE.

Calling the code download.metadata.sh from inside the bin folder
```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=ATAC-seq&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&type=Experiment"
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

Getting the data on bigBed_narrowPeak and pseudoreplicated_peak data for humans in the 38 assembly
```bash
grep -F "bigBed_narrowPeak" metadata.tsv |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt
```
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```


