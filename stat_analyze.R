source('~/prog/r/spatial_organization_measurement_2-03-2016/fview.R')

# Libraries
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library("MASS", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library("FSA", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library("dunn.test", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")

# Data without Germ
FISHdata <- subset(FISHdata, Tissue != 'Germ c.')

FISHdata <- subset(FISHdata, XN>0.03)

plevel <- 0.05

# Devide into different namespase
atr_NC <- subset(FISHdata, Sp == 'An. atroparvus' & Tissue == 'Nurse c.')
lab_NC <- subset(FISHdata, Sp == 'An. labranchiae' & Tissue == 'Nurse c.')
atr_FE <- subset(FISHdata, Sp == 'An. atroparvus' & Tissue == 'Fol ep.')
lab_FE <- subset(FISHdata, Sp == 'An. labranchiae' & Tissue == 'Fol ep.')

# Radial position of X-chromosome

aggregate(X_rad_pos ~ Sp+Tissue, FISHdata, mean)
radDenPlot <- ggplot(FISHdata)+
  geom_density(aes(x=X_rad_pos*100, y=..density..*100, fill=Sp), alpha=0.8)+theme_light() +
  xlab('Radial position of X CT [%]') +
  ylab('Relative DNA content [%]') +
  geom_vline(xintercept = mean(FISHdata$X_rad_pos)*100, col='red', linetype='dashed', size=1.2)+
  facet_wrap(~Sp_Tissue)+
  theme(legend.position='bottom', legend.title=element_blank(), plot.title=element_text())+
  scale_fill_manual(values=c('orange', 'lightblue'))+
  ggtitle('Radial position of X CT')+
  scale_x_continuous(breaks=seq(0, 100, 10))+
  scale_y_continuous(breaks=seq(0, 6, 1))

ggsave('Radial position Density plot.png', plot=radDenPlot, dpi=300, width = 9, height = 5)

radPosPlot <- ggplot(FISHdata)+
  geom_boxplot(
    aes(x = Sp, y = X_rad_pos * 100, fill = Sp)
  )+
  geom_hline(yintercept = mean(FISHdata$X_rad_pos*100), linetype='dashed', size=1.2, col='red') +
  facet_wrap(~Tissue)+
  xlab('')+
  ylab('Radial position of X CT [%]') +
  theme_bw()+
  theme(legend.position='bottom', legend.title=element_blank())+
  scale_fill_manual(values=c('orange', 'lightblue'))+
  ggtitle('Radial position of X CT')+
  scale_y_continuous(breaks=seq(0, 100, 10))

ggsave('Radial position plot.png', plot=radPosPlot, dpi=300, width = 9, height = 5)
radPos_comparison <- dunnTest(X_rad_pos~as.factor(Sp_Tissue), data=FISHdata, method='bonf')
radPos_comparison$res$check <- radPos_comparison$res$P.adj < plevel
subset(radPos_comparison$res, select=c(Comparison, P.adj, check))

# X volumes
xVolumesPlot <- ggplot(FISHdata, aes(x=Sp, y=X_volume/Nuc_volume*100, fill=Sp))+
  geom_boxplot()+
  facet_wrap(~Tissue)+
  theme_bw()+
  ylab('Relative volume of X [%]')+
  xlab('')+
  geom_hline(yintercept = mean(FISHdata$X_volume/FISHdata$Nuc_volume*100), linetype='dashed', size=1.2, col='red')+
  scale_fill_manual(values=c('orange', 'lightblue'))+
  ggtitle('Relative volume of X CT')+
  theme(legend.position='bottom', legend.title=element_blank(), plot.title=element_text())+
  scale_y_continuous(breaks=seq(0, 40, 5))

ggsave('Volumes of X.png', plot=xVolumesPlot, dpi=300, width = 9, height = 5)

tapply(FISHdata$X_volume/FISHdata$Nuc_volume, FISHdata$Sp_Tissue, mean)
tapply(FISHdata$X_volume/FISHdata$Nuc_volume, FISHdata$Sp_Tissue, sd)
tapply(FISHdata$X_surface/FISHdata$Nuc_volume, FISHdata$Sp_Tissue, mean)

vol_comparisons <- dunnTest(XN~as.factor(Sp_Tissue), data=FISHdata)
vol_comparisons$res$check <- vol_comparisons$res$P.adj < plevel
subset(vol_comparisons$res, select=c(Comparison, P.adj, check))
surf_comparisons <- dunnTest(X_surface/Nuc_volume~as.factor(Sp_Tissue), FISHdata)
surf_comparisons$res$check <- surf_comparisons$res$P.adj < plevel
subset(surf_comparisons$res, select=c(Comparison, P.adj, check))


# Compacity / Elongation

comPlot <- ggplot(FISHdata, aes(x=Sp, y=X_compacity*100))+
  geom_boxplot(aes(fill=Sp))+
  facet_wrap(~Tissue)+
  xlab('')+
  ylab('Compacity [%]')+
  scale_fill_manual(values=c('orange', 'lightblue'))+
  theme_bw()


com_comparison <- dunnTest(X_compacity~as.factor(Sp_Tissue), data=FISHdata, method='bonf')
com_comparison$res$check <- com_comparison$res$P.adj < plevel
subset(com_comparison$res, select = c(Comparison, P.adj, check))

elongPlot <- ggplot(FISHdata, aes(x=Sp, fill=Sp, y=X_rel_ellong))+
  geom_boxplot()+
  facet_wrap(~Tissue)+
  xlab('')+
  ylab('Relative elongation [%]') +
  scale_fill_manual(values=c('orange', 'lightblue'))+
  theme_bw()

elong_comparison <- dunnTest(X_rel_ellong~as.factor(Sp_Tissue), FISHdata, method='bonf')
elong_comparison$res$check <- elong_comparison$res$P.adj < plevel
subset(elong_comparison$res, select = c(Comparison, P.adj, check))

mComEl <- grid.arrange(comPlot, elongPlot)
ggsave('Compacity and elongation plot.png', plot=mComEl, dpi=300, width=9, height=9)

# DC avg

dcPlot <- ggplot(FISHdata, aes(x=X_DC_avg/X_volume))+
  
  geom_density(aes(y=..density.., fill=Sp))+
  facet_wrap(~Sp_Tissue)+
  theme_bw()+
  xlab('Relative average distance from center to the surface of X CT [μm^-1]')+
  ylab('Density [%]') +
  geom_vline(xintercept = mean(FISHdata$X_DC_avg/FISHdata$X_volume), size=1.3, col='red', linetype='dashed')+
  scale_fill_manual(values=c('orange', 'lightblue'))+
  theme(legend.position='bottom', legend.title=element_blank())+
  xlim(c(0, 0.4))


dc_comparison <- dunnTest(X_DC_avg/X_volume~as.factor(Sp_Tissue), FISHdata, method='bonf')
dc_comparison$res$check <- dc_comparison$res$P.adj < plevel
subset(dc_comparison$res, select = c(Comparison, check))
  
  
ggsave('DC plot.png', plot=dcPlot, dpi=300, width = 9, height = 5)

