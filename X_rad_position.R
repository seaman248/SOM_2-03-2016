ggplot(FISHdata, aes(x=X_efvMax*100, y=..density..))+
  #geom_histogram(aes(fill=Sp), alpha=0.4)+
  #geom_density(aes(col=Sp))+
  stat_density()+
  
  theme_bw()+
  facet_wrap(~Sp_Tissue, scales = 'free')+
  scale_fill_manual(values=c('orange', 'lightblue'))+
  theme(legend.position='none', legend.title=element_blank())+
  
  xlab('Radial position [%]')+
  ylab('Density [%]')

aggregate(X_efvMax~Sp_Tissue, data=FISHdata, mean)
aggregate(X_efvMax~Sp_Tissue, data=FISHdata, sd)

maxRadPosComparison <- dunnTest(X_efvMax~as.factor(Sp_Tissue), FISHdata, method='bonf')
maxRadPosComparison$res$check <- maxRadPosComparison$res$P.adj < 0.01
subset(maxRadPosComparison$res, select=c(Comparison, P.adj, check))
