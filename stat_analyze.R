source('~/prog/r/spatial_organization_measurement_2-03-2016/fview.R')
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library("MASS", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")

FISHdata <- subset(FISHdata, Tissue != 'Germ c.')

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

# wilcox.test(atr_NC$X_rad_pos, lab_NC$X_rad_pos, alternative = 'less', conf.level=0.95)
# # p-value = 0.102
# wilcox.test(atr_FE$X_rad_pos, lab_FE$X_rad_pos)
# # p-value = 0.03693
# wilcox.test(atr_NC$X_rad_pos, atr_FE$X_rad_pos)
# # p-value = 0.5912
# wilcox.test(lab_NC$X_rad_pos, lab_FE$X_rad_pos)
# # p-value = 0.4378
# shapiro.test(atr_NC$X_rad_pos)
# # p-value = 0.3819
# shapiro.test(lab_NC$X_rad_pos)
# # p-value = 0.7698
# t.test(atr_NC$X_rad_pos, lab_NC$X_rad_pos)
# # p-value = 0.259
pairwise.wilcox.test(FISHdata$X_rad_pos, FISHdata$Sp_Tissue, p.adj='bonf')
ks.test(atr_FE$X_rad_pos, lab_FE$X_rad_pos)
ks.test(atr_NC$X_rad_pos, lab_NC$X_rad_pos)

radPosPlot <- ggplot(FISHdata)+
  geom_boxplot(aes(x=Sp, y=X_rad_pos, fill=Sp))+
  geom_hline(yintercept = mean(FISHdata$X_rad_pos), linetype='dashed', size=1.2, col='red') +
  facet_wrap(~Tissue)+
  xlab('')+
  ylab('Radial position of X CT [%]') +
  theme_light()+
  theme(legend.position='bottom', legend.title=element_blank())+
  scale_fill_manual(values=c('orange', 'lightblue'))+
  ggtitle('Radial position of X CT')+
  scale_y_continuous(breaks=seq(0, 1, 0.1))

ggsave('Radial position plot.png', plot=radPosPlot, dpi=300, width = 9, height = 5)

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

# Compacity / Elongation

comElongPlot <- ggplot(FISHdata, aes(x=X_compacity*100, y=X_rel_ellong*100, col=Sp))+
  geom_point(size=3)+
  stat_smooth(aes(col=Sp), method=glm, se=FALSE)+
  facet_wrap(~Tissue)+
  xlab('Compacity of X CT[%]')+
  ylab('Relative elongation of X CT [%]')+
  scale_colour_manual(values=c('orange', 'lightblue'))+
  ggtitle('Compacity ~ Elongation ratio')+
  theme(legend.position='bottom', legend.title=element_blank())+
  theme_bw()

ggsave('Compacity-Elongation ratio plot.png', plot = comElongPlot, dpi=300, width = 9, height = 5)

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
  
ggsave('DC plot', plot=dcPlot, dpi=300, width = 9, height = 5)

