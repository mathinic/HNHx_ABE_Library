#!/bin/bash
datum=$(date +"%Y%m%d")
touch $datum"_cutadaptlog.txt"
touch $datum"_flashlog.txt"
for filename in *R1_001.fastq.gz; do
	shortname="${filename:0:-21}"
	flash -o $shortname $filename $shortname"_L001_R2_001.fastq.gz" &>> $datum"_flashlog.txt"
	echo " " >> $datum"_flashlog.txt"
done
for filename in *.extendedFrags.fastq; do
	shortname="${filename:0:-20}"
	cutadapt -g AAAGGACGAAACACCG -o $shortname"_5trim.fastq.gz" $filename >>$datum"_cutadaptlog.txt"
	echo " " >>$datum"_cutadaptlog.txt"
	cutadapt -a gttttagagctagaaa -m 20 -M 20 -o $shortname"_Protospacer.fastq.gz" $shortname"_5trim.fastq.gz" >>$datum"_cutadaptlog.txt"
	echo " " >>$datum"_cutadaptlog.txt" 
	cutadapt -g gagtcggtgcTTTTTT -o $shortname"_5trim_temp_barcode1.fastq.gz" $filename >>$datum"_cutadaptlog.txt"
	echo " " >>$datum"_cutadaptlog.txt"
	cutadapt -g agaaatagcaagtta -o $shortname"_target_5trim.fastq.gz" $filename >>$datum"_cutadaptlog.txt"
	echo " " >>$datum"_cutadaptlog.txt"
	cutadapt -a ctgaaagcgagcgtga -o $shortname"_target.fastq.gz" $shortname"_target_5trim.fastq.gz" >>$datum"_cutadaptlog.txt"
	echo " " >>$datum"_cutadaptlog.txt"
	cutadapt -l 6 -o $shortname"_barcode1.fastq.gz" $shortname"_5trim_temp_barcode1.fastq.gz" >>$datum"_cutadaptlog.txt"
	echo " " >>$datum"_cutadaptlog.txt"
done


