args <- commandArgs(trailingOnly = TRUE)
print(args)
# Usage: Rscript plot.r target1 target2 control1 control2 pvalue.cutoff title output

require(plyr)
library(qvalue)
require(dplyr)
require(ggplot2)
require(ggrepel)


compute_psi <- function(df) {
    df$sj_count = df$inc + df$exc + df$ret
    df$psi   = df$inc / df$sj_count
    subset(df, sj_count>=20)
}

pool_biorep <- function(df1,df2) {
    merge(df1,df2,by="id", all=F) -> df
    df[is.na(df)] <- 0
    with(df, data.frame(id=id, inc=inc.x+inc.y, exc=exc.x+exc.y, ret=ret.x + ret.y))
}

dir = args[1]

target1 = read.delim(paste(dir,"M07/",args[2],".M07.tsv", sep=""), col.names=c('id','inc','exc','ret'))
target2 = read.delim(paste(dir,"M07/",args[3],".M07.tsv", sep=""), col.names=c('id','inc','exc','ret'))
target = compute_psi(pool_biorep(target1,target2))

control1 = read.delim(paste(dir,"M07/",args[4],".M07.tsv", sep=""), col.names=c('id','inc','exc','ret'))
control2 = read.delim(paste(dir,"M07/",args[5],".M07.tsv", sep=""), col.names=c('id','inc','exc','ret'))
control = compute_psi(pool_biorep(control1, control2))

gene.names = read.delim(pipe("grep -v orf data/gene_names.bed"), header=F)
colnames(gene.names) = c('id','gene.name')
merge(merge(control,target,by='id'), gene.names, by='id') -> df

df$cell = args[6]
df$KD   = args[7]
df$q.cutoff = args[8]

output = paste(dir,"N07/",args[7],"_",args[6],".tsv",sep="")
print(output)

df$deltaPSI = round(with(df, psi.y-psi.x), digits=2)
df$log10sjcount = round(with(df, log10(sj_count.x + sj_count.y)), digits=4)
df$log10FC = round(with(df, log10(sj_count.x/sj_count.y)), digits=4)
model = lm(deltaPSI~log10FC, df)
print(summary(model))
df$deltaPSIc = round(model$residuals, digits=2)

df1 = subset(df, deltaPSI!=0)
df1$bin = cut_number(df1$log10sjcount, 10)

ddply(df1,.(bin), summarise, mean=mean(deltaPSIc), sd=sd(deltaPSIc)) -> sd_bins
merge(df1, sd_bins)-> df2
sd_bins$center=unlist(lapply(strsplit(as.character(sd_bins$bin),"[\\[,\\]\\(\\)]",perl=T),function(x){(as.numeric(x[3])+as.numeric(x[2]))/2}))
m = lm(log10(sd)~log10(center),sd_bins)
df2$SD = 10^m$coefficients[1]*df2$log10sjcount^m$coefficients[2]
df2$z = with(df2, round((deltaPSIc-mean)/SD, digits=2))
df2$p = pnorm(-abs(df2$z))
df2$q = qvalue_truncp(df2$p)$qvalues
df2$p = round(-log10(df2$p),2)
df2$q = round(-log10(df2$q),2)

fileds = c('KD', 'cell', 'id', 'gene.name', 'deltaPSI', 'deltaPSIc', 'z', 'p', 'q')
write.table(df2[,fileds], file=output, col.names=T, row.names=F, quote=F, sep="\t")

