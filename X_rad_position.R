ggplot(FISHdata, aes(x=X_efvMax, y=..density.., fill=Sp))+
  facet_wrap(~Sp_Tissue)+
  geom_histogram()+
  geom_density()

table(FISHdata$X_efvMax, FISHdata$Sp_Tissue)
