ggplot(FISHdata, aes(x=X_efvMax, y=..density..))+
  facet_wrap(~Sp_Tissue)+
  geom_histogram(aes(fill=Sp))+
  geom_density(aes(col=Sp))

table(FISHdata$X_efvMax, FISHdata$Sp_Tissue)
