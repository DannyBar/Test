#!/bin/bash
LibName=$1

#/home/db/Desktop/NIH/Software/BISCUIT/biscuit align -t 10 /home/db/Desktop/NIH/Software/BISCUIT/hg19.fa $FName $FName2 > $LibName.sam
#samtools sort -T . -O bam -o $LibName.bam $LibName.sam
#samtools index $LibName.bam

cat $LibName.sam | awk -v var="$mycol_new" -F $'\t' 'BEGIN {OFS = FS} {print $3, $4, $4+150}' > $LibName.Pos.bed
bedtools getfasta -fi hg19.fa -bed $LibName.Pos.bed -fo $LibName.MappedReads.seq
cat $LibName.Pos.bed | grep "chr">$LibName.AllReadPos.bed #A cleaner bed file to show read distribution in IGV.

perl GpCH2BED.pl $LibName > $LibName.bed #Locations of all GpCH in reads passing QC.
	#sort $LibName.bed | uniq > $LibName.U.bed
sortBed -i $LibName.bed > $LibName.Sorted.bed #Sorted. | awk 'seen[$0]++' 
perl MovingAvUP2.pl $LibName.Sorted.bed 1 1 > $LibName.Sorted.MovingAv.1.1.bed #Peak-center moving average of unique GpCH sites.
perl MovingAvUP2.pl $LibName.Sorted.bed 200 3 > $LibName.Sorted.MovingAv.200.3.bed #Peak-center moving average of unique GpCH sites.
perl MovingAvUP2.pl $LibName.Sorted.bed 200 5 > $LibName.Sorted.MovingAv.200.5.bed #Peak-center moving average of unique GpCH sites.

bedtools intersect -a $LibName.Sorted.MovingAv.1.1.bed -b Hela-DS11167.peaks.fdr0.01.hg19.expand7500.bed.gz -wo -bed | awk '{r=($2+$3)/2; p=($6+$7)/2; d=r-p; print d,$4}' > dist.$LibName.Sorted.MovingAv.Dir.1.1.bed.CTCFDS11167ex7500 #Average distance of moving average positions from CTCF binding sites.
Rscript MakePlot.r dist.$LibName.Sorted.MovingAv.Dir.1.1.bed.CTCFDS11167ex7500 #Histogram of said distance.

bedtools intersect -a $LibName.Sorted.MovingAv.200.3.bed -b Hela-DS11167.peaks.fdr0.01.hg19.expand7500.bed.gz -wo -bed | awk '{r=($2+$3)/2; p=($6+$7)/2; d=r-p; print d,$4}' > dist.$LibName.Sorted.MovingAv.Dir.200.3.bed.CTCFDS11167ex7500 #Average distance of moving average positions from CTCF binding sites.
Rscript MakePlot.r dist.$LibName.Sorted.MovingAv.Dir.200.3.bed.CTCFDS11167ex7500 #Histogram of said distance.

bedtools intersect -a $LibName.Sorted.MovingAv.200.5.bed -b Hela-DS11167.peaks.fdr0.01.hg19.expand7500.bed.gz -wo -bed | awk '{r=($2+$3)/2; p=($6+$7)/2; d=r-p; print d,$4}' > dist.$LibName.Sorted.MovingAv.Dir.200.5.bed.CTCFDS11167ex7500 #Average distance of moving average positions from CTCF binding sites.
Rscript MakePlot.r dist.$LibName.Sorted.MovingAv.Dir.200.5.bed.CTCFDS11167ex7500 #Histogram of said distance.


#/home/db/Desktop/NIH/Software/BISCUIT/biscuit align -t 10 /home/db/Desktop/NIH/Software/BISCUIT/hg19.fa D_S4_L001_R1_001.fastq.gz D_S4_L001_R2_001.fastq.gz | samtools sort -T . -O bam -o D.bam
#samtools index D.bam
#/home/db/Desktop/NIH/Software/BISCUIT/biscuit pileup -r /home/db/Desktop/NIH/Software/BISCUIT/hg19.fa -i D.bam -o D.vcf -q 20
#bgzip D.vcf
#tabix -p vcf D.vcf.gz
#/home/db/Desktop/NIH/Software/BISCUIT/biscuit vcf2bed -k 10 -t gch D.vcf.gz
#bedtools intersect -a Test.bed -b /home/db/Desktop/NIH/ChAMP/DS_Data/Commons/Hela-DS11167.peaks.fdr0.01.hg19.expand2500.bed.gz -wo -bed | awk '{r=($2+$3)/2; p=($6+$7)/2; d=r-p;if(d<0){d=-d} print d,$4}' > dist.Test.Sorted.MovingAvU.bed.CTCFDS11167ex2500
#Rscript /home/db/Desktop/NIH/ChAMP/DS_Data/Commons/MakePlot.r dist.Test.Sorted.MovingAvU.bed.CTCFDS11167ex2500 #Histogram of said distance.
