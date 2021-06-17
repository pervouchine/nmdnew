.PHONY : all
qth=0.5 
filter=awk '$$6>0.05 || $$6<-0.05' 

all :: data/upf1xrn1/upf1xrn1vscontrol.tsv data/upf1xrn1/smg6xrn1vscontrol.tsv

##### data download ######

data/shRNA_table.tsv data/upf1xrn1/M07/SRR1275413Aligned.sortedByCoord.out.M07.tsv :
	wget --no-check-certificate http://cb.skoltech.ru/dp/papers/nmd/data.tar.gz
	tar -xf data.tar.gz
	rm -f data.tar.gz

##### annotation ####

# genes.bed is a bed file containing gene names
data/genes.bed : data/gencode.v19.annotation.gtf
	awk '$$3=="gene"' data/gencode.v19.annotation.gtf | grep 'gene_type "protein_coding";' | perl Perl/print_gff_attributes.pl GEN gene_name | sort -u > data/genes.bed 

data/gene_names.bed : data/genes.bed
	cat data/upf1xrn1/M07/* | cut -f1 | sort -u | awk -v OFS="\t" '{split($$1,a,"_");print a[1],a[2],a[2]+1,$$1,1,a[3]}' | intersectBed -a stdin -b data/genes.bed -s -wa -wb | cut -f4,13 > data/gene_names.bed

data/annotated_gencode.tsv : ~/db/annotation/gencode.v19.annotation.gtf
	awk '$$3=="exon"' ~/db/annotation/gencode.v19.annotation.gtf | awk '{if($$7=="+"){print $$1"_"$$4"_"$$7"_A\n"$$1"_"$$5"_"$$7"_D"}else{print $$1"_"$$4"_"$$7"_D\n"$$1"_"$$5"_"$$7"_A"}}' | sort -u > data/annotated_gencode.tsv

data/annotated_refgene.tsv : data/refGene.txt.gz
	zcat data/refGene.txt.gz | awk '{split($$10,a,",");split($$11,b,",");for(i=1;i<length(a);i++){if($$4=="+"){print $$3"_"a[i]"_+_A\n"$$3"_"b[i]"_+_D"}else{print $$3"_"a[i]"_-_D\n"$$3"_"b[i]"_-_A"}}}' | sort -u > data/annotated_refgene.tsv

data/annotated_all.tsv : data/annotated_refgene.tsv data/annotated_gencode.tsv
	cat data/annotated_refgene.tsv data/annotated_gencode.tsv | sort -u > data/annotated_all.tsv

clean ::
	rm -f data/genes.bed data/gene_names.bed data/annotated_sites.tsv

#### eCLIP ####

# exon_peaks.bed contains eCLIP peaks in 5000-nt window around exons
# cols 1-6 are exon windows (#4 is exon id); #7 is gene name
# cols 7-17 are eCLIP peaks
data/eCLIP/selected_peaks.bed : data/eCLIP/dump.tsv data/gene_names.bed
	awk -v OFS="\t" -v w=5000 '{split($$1,a,"_");if(a[2]<w){a[2]=w};print a[1],a[2]-w,a[2]+w,$$1,1,a[3],$$2}' data/gene_names.bed | intersectBed -a stdin -b data/eCLIP/dump.tsv -wa -wb -s > data/eCLIP/selected_peaks.bed

# This file tabulates the distance to cognate eCLIP peaks for all exons in RBP
#data/eCLIP/exon_peaks_dist.tsv : data/eCLIP/exon_peaks.bed
#	awk -v OFS="\t" '{split($$11,a,"_");if(a[1]==$$7){d=(($$9+$$10)-($$2+$$3))/2;if(d<0){d=-d};print $$4,a[1],d,$$14}}' data/eCLIP/exon_peaks.bed > data/eCLIP/exon_peaks_dist.tsv

# self_peaks.tsv contains three columns: exon_id, peak_id, gene name
#data/eCLIP/self_peaks.bed : data/eCLIP/dump.tsv data/genes.bed
#	intersectBed -a data/eCLIP/dump.tsv -b data/genes.bed -wa -wb -s | awk '{split($$4,a,"_");if(a[1]==$$17){print}}' | sort -k1,1 -k2,2n > data/eCLIP/self_peaks.bed

# this is a bed file for track hub with all self peaks
#hub/eCLIP_peaks.bed :  data/eCLIP/self_peaks.bed
#	cut -f1-6 data/eCLIP/self_peaks.bed | awk -v OFS="\t" '{split($$4,a,"_");$$4=a[1];print}' | bedtools merge -s -c 4 -o distinct -i stdin | awk -v OFS="\t" '{print $$1,$$2,$$3,$$5,1000,$$4}' | sort -k1,1 -k2,2n > hub/eCLIP_peaks.bed

# A bigBed version of it for track hub
#hub/hg19/eCLIP.bb : hub/eCLIP_peaks.bed hub/hg19/hg19.chrom.sizes
#	./bedToBigBed hub/eCLIP_peaks.bed hub/hg19/hg19.chrom.sizes  hub/hg19/eCLIP.bb
#	git add hub/hg19/eCLIP.bb

#all :: hub/hg19/eCLIP.bb data/eCLIP/exon_peaks_dist.tsv

#### shRNA-KD ####

# This step creates a make file to compute deltaPSI for all shRNA-KD
# In the make file, a script called deltaPSIc.r is used to compute deltaPSI, correct for the expression changes, and estimate p-values and q-values
shRNA.mk : data/shRNA_table.tsv
	awk -v dir=data/shRNA/ '{out=dir"N07/"$$6"_"$$5".tsv "; print  out ": "dir"M07/"$$1".M07.tsv "dir"M07/"$$2".M07.tsv "dir"M07/"$$3".M07.tsv "dir"M07/"$$4".M07.tsv\n\tmkdir -p "dir"N07/\n\tRscript deltaPSIc.r "dir,$$1,$$2,$$3,$$4,$$5,$$6"\n"; all=all out}END{print "all :" all}' data/shRNA_table.tsv > shRNA.mk

# At this step, the makefile is executed, and all deltaPSI are pooled for exons in RBPs
data/shRNA/deltaPSI.tsv : shRNA.mk
	make -f shRNA.mk all
	awk '$$1==$$4' data/shRNA/N07/*.tsv > data/shRNA/deltaPSI.tsv

# this is a bed file for track hub
hub/shRNA-KD.bed : data/shRNA/deltaPSI.tsv
	mkdir -p hub/
	${filter} data/shRNA/deltaPSI.tsv | awk -f trackB.awk | sort -k1,1 -k2,2n > hub/shRNA-KD.bed 

# and its bigBed version
hub/hg19/shRNA.bb : hub/shRNA-KD.bed hub/hg19/hg19.chrom.sizes
	./bedToBigBed hub/shRNA-KD.bed hub/hg19/hg19.chrom.sizes hub/hg19/shRNA.bb
	git add hub/hg19/shRNA.bb

all :: hub/hg19/shRNA.bb
clean :: 
	rm -f hub/hg19/shRNA.bb hub/shRNA-KD.bed data/shRNA/deltaPSI.tsv shRNA.mk

####### UPF1/XRN1 KD ######

# this step generates a pdf, png, and tsv table for UPF1/XRN1 KD vs control
data/upf1xrn1/upf1xrn1vscontrol.tsv : data/upf1xrn1/M07/SRR1275413Aligned.sortedByCoord.out.M07.tsv data/upf1xrn1/M07/SRR1275416Aligned.sortedByCoord.out.M07.tsv plot.r data/gene_names.bed
	Rscript plot.r data/upf1xrn1/M07/ SRR1275413Aligned.sortedByCoord.out SRR1275416Aligned.sortedByCoord.out ${qth} UPF1XRN1 data/upf1xrn1/upf1xrn1vscontrol 

# same for SMG6
data/upf1xrn1/smg6xrn1vscontrol.tsv : data/upf1xrn1/M07/SRR1275413Aligned.sortedByCoord.out.M07.tsv data/upf1xrn1/M07/SRR1275415Aligned.sortedByCoord.out.M07.tsv plot.r data/gene_names.bed
	Rscript plot.r data/upf1xrn1/M07/ SRR1275413Aligned.sortedByCoord.out SRR1275415Aligned.sortedByCoord.out ${qth} SMG6XRN1 data/upf1xrn1/smg6xrn1vscontrol

# delta PSI are pooled for UPF1/XRN1 and SMG6/XRN1 expts
data/upf1xrn1/deltaPSI.tsv : upf1andsmg6.r data/upf1xrn1/upf1xrn1vscontrol.tsv data/upf1xrn1/smg6xrn1vscontrol.tsv
	Rscript upf1andsmg6.r data/upf1xrn1/upf1xrn1vscontrol.tsv data/upf1xrn1/smg6xrn1vscontrol.tsv data/upf1xrn1/deltaPSI.tsv

# This script plots relative positions of reactive exons
#data/upf1xrn1/relpos.pdf : data/upf1xrn1/upf1xrn1vscontrol.tsv data/genes.bed
#	cat data/upf1xrn1/upf1xrn1vscontrol.tsv | grep chr | awk -v OFS="\t" '{split($$1,a,"_");print a[1],a[2],a[2]+1,$$1,1,a[3],$$4}' | sort -k1,1 -k2,2n | intersectBed -a stdin -b data/genes.bed -s -wa -wb -f 1 | awk -v OFS="\t" '{p = ($$2-$$9)/(($$10-$$9)-($$3-$$2));if($$6=="-"){p=1-p};print p,$$7}' | Rscript -e 'df = read.delim("stdin", header=F);df$$Group = factor(df$$V2>0, levels=c(T,F), labels=c("dPSI>0","dPSI<0"));library(ggplot2);pdf("data/upf1xrn1/relpos.pdf");ggplot(df,aes(x=100*V1,fill=Group)) + geom_histogram(position="identity",alpha=0.5,bins=20,aes(y=..density..)) + theme_bw() + xlab("Relative position, %")'

# this is a bed file for track hub
hub/nmd.bed: data/upf1xrn1/upf1xrn1vscontrol.tsv data/upf1xrn1/smg6xrn1vscontrol.tsv
	tail -n+2 data/upf1xrn1/upf1xrn1vscontrol.tsv | ${filter} | awk -f trackA.awk | sort -k1,1 -k2,2n > hub/upf1xrn1.bed
	tail -n+2 data/upf1xrn1/smg6xrn1vscontrol.tsv | ${filter} | awk -f trackA.awk | sort -k1,1 -k2,2n > hub/smg6xrn1.bed 
	cat hub/upf1xrn1.bed hub/smg6xrn1.bed | sort -k1,1 -k2,2n > hub/nmd.bed

# and its bigBed version
hub/hg19/nmd.bb : hub/nmd.bed
	./bedToBigBed hub/nmd.bed hub/hg19/hg19.chrom.sizes hub/hg19/nmd.bb
	git add hub/hg19/nmd.bb

all :: hub/hg19/nmd.bb data/upf1xrn1/deltaPSI.tsv
clean ::
	rm -f hub/hg19/nmd.bb hub/nmd.bed data/upf1xrn1/deltaPSI.tsv data/upf1xrn1/smg6xrn1vscontrol.tsv data/upf1xrn1/upf1xrn1vscontrol.tsv

####################

# Here all the pieces are combined together in one table
data/combined.pdf: data/upf1xrn1/deltaPSI.tsv data/shRNA/deltaPSI.tsv combined.r data/eCLIP/exon_peaks_dist.tsv
	Rscript combined.r data/upf1xrn1/deltaPSI.tsv data/shRNA/deltaPSI.tsv data/eCLIP/exon_peaks_dist.tsv data/combined

# Script for Figure 2C
#data/poison.pdf : poison.r data/upf1xrn1/upf1xrn1vscontrol.tsv data/stop.tsv
	#Rscript poison.r data/upf1xrn1/upf1xrn1vscontrol.tsv data/poison.pdf

# Script for Figure 2D
#data/essential.pdf : essential.r data/upf1xrn1/upf1xrn1vscontrol.tsv data/stop.tsv
#	Rscript essential.r data/upf1xrn1/upf1xrn1vscontrol.tsv data/essential.pdf

all :: data/combined.pdf 
	# data/essential.pdf data/poison.pdf

# End of pipeline
# Should you have any questions please don't hesitate to contact Dmitri Pervouchine pervouchine@gmail.com
