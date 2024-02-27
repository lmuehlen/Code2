save_vis<-function(p,name,height=15,width=30,units="cm",path="C:/Users/leona/Dropbox/Apps/Overleaf/Fiscal Rules Budget Composition/Graphics"){
  
  ggsave(paste(name,".pdf"), p, device = "pdf", height = 15, width = 30, units = "cm", path = "Output/")
}