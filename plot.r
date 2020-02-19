args <- commandArgs(trailingOnly = TRUE)
print(args)
# Usage: Rscript plot.r target1 target2 control1 control2 pvalue.cutoff title output

require(plyr, quietly =T)
library(qvalue, quietly =T)
require(dplyr, quietly = T)
require(ggplot2, quietly =T)
require(ggrepel, quietly =T)


compute_psi <- function(df) {
    df$expr = df$inc + df$exc + df$ret
    df$psi   = df$inc / df$expr
    subset(df, expr>=20)
}

dir = args[1]

control = compute_psi(read.delim(paste(dir,args[2],".M07.tsv", sep=""), col.names=c('id','inc','exc','ret')))

target  = compute_psi(read.delim(paste(dir,args[3],".M07.tsv", sep=""), col.names=c('id','inc','exc','ret')))

gene.names = read.delim(pipe("grep -v orf data/gene_names.bed"), header=F)
colnames(gene.names) = c('id','gene.name')
merge(merge(control,target,by='id'), gene.names, by='id') -> df

q.cutoff = args[4]
title    = args[5]
output   = args[6]

df$KD = title;
df$cell = 'HEK293';

df$deltaPSI = round(with(df, psi.y-psi.x), digits=2)
df$log10expt = round(with(df, log10(expr.x + expr.y)), digits=4)
df$log10FC = round(with(df, log10(expr.x/expr.y)), digits=4)
model = lm(deltaPSI~log10FC, df)
print(summary(model))
df$deltaPSIc = round(model$residuals, digits=2)

df1 = subset(df, deltaPSI!=0)
df1$bin = cut_number(df1$log10expt, 40)

ddply(df1,.(bin), summarise, mean=mean(deltaPSIc), sd=sd(deltaPSIc)) -> sd_bins
merge(df1, sd_bins)-> df2
sd_bins$mid=unlist(lapply(strsplit(as.character(sd_bins$bin),"[\\[,\\]\\(\\)]",perl=T),function(x){(as.numeric(x[3])+as.numeric(x[2]))/2}))
m = lm(log10(sd)~log10(mid),sd_bins)
print(summary(m))
df2$SD = 10^(m$coefficients[1] + log10(df2$log10expt)*m$coefficients[2])
df2$z = with(df2, round((deltaPSIc-mean)/SD, digits=2))
df2$p = pnorm(-abs(df2$z))
df2$q = qvalue_truncp(df2$p)$qvalues
df2$p = round(-log10(df2$p),2)
df2$q = round(-log10(df2$q),2)

fileds = c('KD', 'cell', 'id', 'gene.name', 'deltaPSI', 'deltaPSIc', 'z', 'p', 'q')
write.table(df2[,fileds], file=paste(output, ".tsv", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

df3 = subset(df2, q>q.cutoff) %>% group_by(gene.name) %>% slice(which.max(abs(deltaPSIc)))
df3 = head(df3[order(-df3$q),], n=70)

p = ggplot(df2, aes(x=log10expt, y=deltaPSIc)) + geom_point(size=0.5, alpha=0.5,aes(color=(q>q.cutoff))) 
p = p + geom_abline(slope=0,lty="dashed") + geom_label_repel(size=3, data=df3, aes(label = gene.name, color=(q>q.cutoff)), min.segment.length = 0)
p = p + xlab("log(expression)") + ylab(expression(paste(Delta,Psi,"'(KD-Control)")))
md = data.frame(x=seq(0.5,5.5,0.01))
md$y = 10^(m$coefficients[1]+log10(md$x)*m$coefficients[2])

p = p + theme_classic() + theme(legend.position="none") + scale_color_brewer(palette="Set2") + ggtitle(title) + geom_line(data=md,aes(x=x,y=y),color='blue',linetype="dashed") + geom_line(data=md,aes(x=x,y=-y),color='blue',linetype="dashed")

pdf(width=5,height=6, paste(output, ".pdf", sep=""))
print(p)
ggsave(width=5,height=6,paste(output, ".png", sep=""))
dev.off()

