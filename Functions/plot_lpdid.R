plot_lpdid<-function (reg, conf = 0.9, segments = TRUE, add = FALSE, xlab = NULL, 
          ylab = NULL, main = "", x.shift = 0, pch = 19, cex = 1, 
          col = "black", opacity = 1) 
{
  if (nrow(reg$coeftable) != length(reg$window)) 
    stop("coeftable and window are not the same length.  It is likely that pooled=TRUE in the lpdid function.  An event study cannot be plotted when pooled=TRUE.")
  coeftable <- reg$coeftable
  coeftable$t <- reg$window
  conf_z <- abs(qnorm((1 - conf)/2))
  coeftable$uCI <- coeftable$Estimate + conf_z * coeftable$`Std. Error`
  coeftable$lCI <- coeftable$Estimate - conf_z * coeftable$`Std. Error`
 
coeftable%>%
  ggplot(aes(t))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = -1,linetype=2)+
  geom_segment(aes(y=uCI,yend=lCI,xend=t))+
  geom_point(aes(y=Estimate))+
  theme_tufte(base_size = 19)+
  theme(panel.grid.major.y = element_line(linetype = "dotted",color = "grey",linewidth=0.5),axis.ticks = element_blank(),legend.position = "bottom")+
  labs(x="Years since fiscal rule change",y="Effect")

}
