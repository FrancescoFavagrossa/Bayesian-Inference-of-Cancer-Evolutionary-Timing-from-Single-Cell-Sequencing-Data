setwd("/Users/francescofavagrossa/Documents/Projects/Thesis")
suppressMessages({library(ggplot2);library(dplyr)})
source("robustness_checks/_plot_style.R")
df<-readRDS("/tmp/hp_sweep.rds")
df$clone<-factor(df$clone, levels=c("Clone 4","Clone 6"), labels=c("Clone 4 (chain)","Clone 6 (branching)"))
pal<-c("Clone 4 (chain)"="#E45756","Clone 6 (branching)"="#4C78C8")
base_thm<-theme_thesis(15)+theme(legend.position="top",legend.title=element_blank())
d1<-df%>%filter(sweep=="sigma_rate")
p1<-ggplot(d1,aes(value,tbias,colour=clone))+
  geom_vline(xintercept=1,linetype="dashed",colour="grey55")+
  annotate("text",x=1.05,y=max(d1$tbias),label="chosen~lambda==1",parse=TRUE,hjust=0,size=4.3,colour="grey40")+
  geom_line(linewidth=1)+geom_point(size=3)+scale_colour_manual(values=pal)+
  labs(title="Event-time error vs random-walk step prior",
       x=expression(sigma[log~w]~"~ Exponential("*lambda*")   ( larger "*lambda*" = tighter )"),
       y="Mean |bias|  (pseudotime)")+base_thm
ggsave("Thesis_final/immagini/robust_sigma_rate.png",p1,width=8,height=5.4,dpi=150,bg="white")
d2<-df%>%filter(sweep=="omega_mean")
p2<-ggplot(d2,aes(value,tbias,colour=clone))+
  geom_vline(xintercept=1.25,linetype="dashed",colour="grey55")+
  annotate("text",x=1.29,y=max(d2$tbias),label="chosen~omega[0]==1.25",parse=TRUE,hjust=0,size=4.3,colour="grey40")+
  geom_line(linewidth=1)+geom_point(size=3)+scale_colour_manual(values=pal)+
  labs(title="Event-time error vs clock-factor prior mean",
       x=expression(omega~"~ LogNormal(log "*omega[0]*", 0.5)"),
       y="Mean |bias|  (pseudotime)")+base_thm
ggsave("Thesis_final/immagini/robust_omega_mean.png",p2,width=8,height=5.4,dpi=150,bg="white")
cat("FIX_DONE\n")
