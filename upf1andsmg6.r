args <- commandArgs(trailingOnly = TRUE)

upf1 = args[1]
smg6 = args[2]
out = args[3]

fields = c('id','gene.name','deltaPSI','deltaPSIc','z','p','q')
df1 = read.delim(upf1)[,fields]
df2 = read.delim(smg6)[,fields]
df3 = merge(df1,df2,by=c('id','gene.name'),suffixes=c('.UPF1','.SMG6'),all=T)
df3$deltaPSIc = apply(df3[,c('deltaPSIc.UPF1','deltaPSIc.SMG6')],1,mean,na.rm=T)
df3$p = apply(df3[,c('p.UPF1','p.SMG6')],1,sum,na.rm=T)
df3$q = apply(df3[,c('q.UPF1','q.SMG6')],1,sum,na.rm=T)
write.table(df3, file=args[3], sep="\t", col.names=T, row.names=F, quote=F)
