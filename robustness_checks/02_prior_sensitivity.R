setwd("/Users/francescofavagrossa/Documents/Projects/Thesis")
suppressMessages({library(cmdstanr); library(posterior); library(stringr)})
tmpl <- readRDS("robustness_checks/templates/c6d_sd.rds"); tmpl$sigma_nu_rate <- NULL
HY <- list(w_bar_prior_mean=0.06,w_bar_prior_sd=0.35,sigma_log_w_rate=1,omega_prior_mean=1.25,omega_prior_sd=0.5,phi=100)
sim_dgp <- function(tmpl,HY,Ncells_target=3000,seed=NA){
  if(!is.na(seed)) set.seed(seed)
  Nn<-tmpl$Nnodes;Ne<-tmpl$Nevents;Neq<-tmpl$Nequations;par<-tmpl$parents;noe<-tmpl$number_of_events
  ife<-tmpl$idx_first_event;ifm<-tmpl$idx_first_mutation;len<-tmpl$lengths;cf<-tmpl$coeff;mu<-tmpl$mu;dt<-tmpl$Tf_real-tmpl$t0_real
  wbar<-rlnorm(1,log(HY$w_bar_prior_mean),HY$w_bar_prior_sd);sig<-min(max(rexp(1,HY$sigma_log_w_rate),0.12),0.6)
  omega<-rlnorm(1,log(HY$omega_prior_mean),HY$omega_prior_sd);u<-runif(Ne);te<-numeric(Ne)
  for(node in 1:Nn){i0<-ife[node];p<-par[node];for(ev in 1:noe[node]){idx<-i0+ev-1
    if(p==0){te[idx]<-if(ev==1)u[idx] else u[idx]*te[i0]}else{tp<-te[ife[p]];te[idx]<-if(ev==1)tp+u[idx]*(1-tp) else tp+u[idx]*(te[i0]-tp)}}}
  lw<-numeric(Nn);lw[1]<-log(wbar);for(node in 2:Nn)lw[node]<-lw[par[node]]+rnorm(1,0,sig);w<-exp(lw);nu<-omega*w
  mut<-integer(Neq)
  for(node in 1:Nn){i0<-ife[node];p<-par[node];tp<-if(p==0)0 else te[ife[p]];for(ev in 1:noe[node]){idx<-i0+ev-1
    if(ev==1){im<-ifm[i0];r<-nu[node]*mu*dt*len[i0]*cf[im]*(te[i0]-tp);mut[im]<-rpois(1,max(r,1e-9))}else{
      ia<-ifm[idx];ib<-ia+1;tc<-te[idx];ra<-nu[node]*mu*dt*len[idx]*cf[ia]*(tc-tp);rb<-nu[node]*mu*dt*len[idx]*cf[ib]*(te[i0]-tc)
      mut[ia]<-rpois(1,max(ra,1e-9));mut[ib]<-rpois(1,max(rb,1e-9))}}}
  M<-numeric(Nn);for(node in 1:Nn){tn<-te[ife[node]];M[node]<-exp(min(w[node]*dt*(1-tn),200))}
  for(ch in 2:Nn){p<-par[ch];tc<-te[ife[ch]];M[p]<-M[p]-exp(min(w[p]*dt*(1-tc),200))};M[M<1e-12]<-1e-12;pp<-M/sum(M)
  Nc<-rnbinom(Nn,size=HY$phi,mu=Ncells_target*pp)
  sd<-tmpl;sd$mutations<-as.integer(mut);sd$Ncells_per_clone<-as.integer(Nc);sd$Ncells<-as.integer(sum(Nc))
  sd$w_bar_prior_mean<-HY$w_bar_prior_mean;sd$w_bar_prior_sd<-HY$w_bar_prior_sd;sd$sigma_log_w_rate<-HY$sigma_log_w_rate
  sd$omega_prior_mean<-HY$omega_prior_mean;sd$omega_prior_sd<-HY$omega_prior_sd;sd$phi<-HY$phi
  list(stan_data=sd,true_te=te,true_w=w)
}
m<-cmdstan_model("stan_models/base_rw.stan")
sm<-function(f,v){s<-f$summary(v);s[order(as.integer(str_extract(s$variable,"[0-9]+"))),]}
NR<-12
settings<-list(
  base    =function(sd){sd},
  om_tight=function(sd){sd$omega_prior_sd<-0.2;sd},
  om_wide =function(sd){sd$omega_prior_sd<-1.0;sd},
  sig_rate2=function(sd){sd$sigma_log_w_rate<-2;sd})
acc<-setNames(lapply(names(settings),function(x)list(tb=c(),t90=c(),wb=c(),w90=c(),rh=c())),names(settings))
for(r in 1:NR){
  S<-sim_dgp(tmpl,HY,seed=200+r)
  for(nm in names(settings)){
    sd<-settings[[nm]](S$stan_data)
    f<-m$sample(data=sd,chains=2,parallel_chains=2,iter_warmup=1200,iter_sampling=1500,adapt_delta=0.95,max_treedepth=12,seed=1,refresh=0,show_messages=FALSE)
    te<-sm(f,"tevents");w<-sm(f,"w_node_out")
    acc[[nm]]$tb<-c(acc[[nm]]$tb,mean(abs(te$median-S$true_te)))
    acc[[nm]]$t90<-c(acc[[nm]]$t90,(S$true_te>=te$q5 & S$true_te<=te$q95))
    acc[[nm]]$wb<-c(acc[[nm]]$wb,mean(abs(w$median-S$true_w)))
    acc[[nm]]$w90<-c(acc[[nm]]$w90,(S$true_w>=w$q5 & S$true_w<=w$q95))
    acc[[nm]]$rh<-c(acc[[nm]]$rh,max(c(te$rhat,w$rhat),na.rm=TRUE))
  }
  cat(sprintf("rep %2d/%d done\n",r,NR))
}
saveRDS(acc,"/tmp/dgp_part2_res.rds")
cat(sprintf("\n=== PART 2 (prior sensitivity, %d repliche fisse dal DGP baseline) ===\n",NR))
cat(sprintf("%-10s | TEMPI |bias|=%s cov90=%s | GROWTH |bias|=%s cov90=%s | maxrhat=%s\n","setting","","","","",""))
for(nm in names(settings)){a<-acc[[nm]]
  cat(sprintf("%-10s | TEMPI |bias|=%.3f cov90=%.0f%% | GROWTH |bias|=%.3f cov90=%.0f%% | maxrhat=%.2f\n",
    nm, mean(a$tb), 100*mean(a$t90), mean(a$wb), 100*mean(a$w90), max(a$rh)))}
cat("PART2_DONE\n")
