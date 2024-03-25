save_vis<-function(p,name,height=12,width=30,units="cm",path="C:/Users/leona/Dropbox/Apps/Overleaf/Fiscal Rules Budget Composition/Graphics"){
  
  ggsave(paste0(name,".pdf"), p, device = "pdf", height = height, width = width, units = "cm", path =path)
}
