setwd("/Users/francescofavagrossa/Documents/Projects/Thesis")
suppressMessages({library(cmdstanr);library(posterior);library(stringr)})
m<-cmdstan_model("stan_models/base_rw.stan")
sm<-function(f,v){s<-f$summary(v);s[order(as.integer(str_extract(s$variable,"[0-9]+"))),]}
go<-function(sd){m$sample(data=sd,chains=4,parallel_chains=4,iter_warmup=2000,iter_sampling=2500,adapt_delta=0.99,max_treedepth=14,seed=42,refresh=0,show_messages=FALSE)}
truth<-list(
 "Clone 6"=list(file="robustness_checks/templates/c6d_sd.rds",te=c(0,0,0.291,0.602,0.602,0.715,0.715,0.742,0.760),tg=c(0.0493,0.092,0.0946,0.0969,0.1879,0.1165),ne=9),
 "Clone 4"=list(file="robustness_checks/templates/c4d_sd.rds",te=c(0,0,0.383,0.530,0.530,0.730),tg=c(0.0493,0.0983,0.135,0.1883),ne=6))
rows<-list(); k<-1
for(cl in names(truth)){ T<-truth[[cl]]
  base<-readRDS(T$file);base$sigma_nu_rate<-NULL
  ## sweep A: omega prior mean (sd 0.5, rate 1)
  for(om in c(1.0,1.25,1.5,2.0)){ sd<-base;sd$sigma_log_w_rate<-1;sd$omega_prior_mean<-om;sd$omega_prior_sd<-0.5
    f<-go(sd);t<-sm(f,"tevents");w<-sm(f,"w_node_out");ae<-abs(t$median-T$te)[3:T$ne]
    ct<-mean(T$te[3:T$ne]>=t$q5[3:T$ne]&T$te[3:T$ne]<=t$q95[3:T$ne])
    rows[[k]]<-data.frame(clone=cl,sweep="omega_mean",value=om,tbias=mean(ae),tmax=max(ae),tcov=ct,gbias=mean(abs(w$median-T$tg)));k<-k+1 }
  ## sweep B: sigma_log_w rate (omega LN 1.25,0.5)
  for(rt in c(0.5,1,2,3)){ sd<-base;sd$sigma_log_w_rate<-rt;sd$omega_prior_mean<-1.25;sd$omega_prior_sd<-0.5
    f<-go(sd);t<-sm(f,"tevents");w<-sm(f,"w_node_out");ae<-abs(t$median-T$te)[3:T$ne]
    ct<-mean(T$te[3:T$ne]>=t$q5[3:T$ne]&T$te[3:T$ne]<=t$q95[3:T$ne])
    rows[[k]]<-data.frame(clone=cl,sweep="sigma_rate",value=rt,tbias=mean(ae),tmax=max(ae),tcov=ct,gbias=mean(abs(w$median-T$tg)));k<-k+1 }
  cat("done",cl,"\n")
}
df<-do.call(rbind,rows); saveRDS(df,"/tmp/hp_sweep.rds"); print(df); cat("HP_DONE\n")
