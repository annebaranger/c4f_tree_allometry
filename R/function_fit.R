#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_fit.R  
#' @description R script with all functions for fitting allometric models
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - data fit ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param tree.merge file where trees were merges
#' @param plot.merge Idem for plots

data.fit<-function(tree,
                   plot,
                   frac=1){
  data_tree=tree %>%
    left_join(plot) %>%
    filter(!(H==1&dbh>5)) %>% # filter weird tree
    filter(!(H<15&dbh>100)) %>%
    filter(dbh<200) %>% 
    filter(is.na(ba_tot)==FALSE) %>%
    filter(is.na(bio01)==FALSE) %>%
    mutate(system=ordered(system,levels=c("forest","secondary_forest","plantation","agroforestry")),
           sys=case_when(system=="forest"~1,
                         system=="secondary_forest"~2,
                         system=="plantation"~3,
                         system=="agroforestry"~4),
           ori=as.numeric(case_when(origin=="remnant"~1,
                                    origin=="recruited"~2,
                                    origin=="planted"~3)),
           id_plot=as.factor(as.character(id_plot)),
           num_plot=as.numeric(id_plot),
           g_s=as.factor(paste0(genusCorr,"_",speciesCorr)),
           num_species=as.numeric(g_s),
           systori=as.numeric(as.factor(paste0(system,"_",origin))),
           v_compet=as.numeric(scale(v_compet,center=FALSE)),
           lba=log(1+ba_tot),
           ba_tot=as.numeric(scale(ba_tot,center=FALSE)),
           bio01=bio01/100,
           bio05=bio05/100,
           bio10=bio10/100,
           bio12=bio12/1000,
           bio17=bio17/100)%>%
    mutate_at(c("database","id_plot","id_tree","system","origin","genus","species"),as.factor)
  if(frac!=1){
    data_tree<-data_tree |> 
      sample_frac(size=frac)
  }
  return(data_tree)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - model shape selection ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param tree.merge file where trees were merges
#' @param plot.merge Idem for plots

shape.select <- function(folder="mod_shape_select",
                         data.mod){
  if(!dir.exists(folder)){
    dir.create(path=folder)
  }
  data_HD = list(
    N = dim(data.mod)[1],
    H = data.mod$H ,
    dbh=data.mod$dbh
  )
  
  ## LOG-LIN
  
  ### Fit
  if(!file.exists(file.path(folder,"HD_loglin.rdata"))){
    mod_HD_loglin=stan(file="stan/model_HD_loglin.stan", # stan program
                       data = data_HD,         # dataset
                       warmup = 1000,          # number of warmup iterations per chain
                       iter = 2000)   
    save(mod_HD_loglin,file=file.path(folder,"HD_loglin.rdata"))
  }else{
    load(file.path(folder,"HD_loglin.rdata"))
  }
  
  ### Prediction 
  model <- function(x,alpha,beta) {
    alpha + beta * log(x)}
  par_mod<- rstan::extract(mod_HD_loglin)
  loglin<-summary(mod_HD_loglin)[[1]]
  data.lin<-data.mod |> 
    rowwise() |> 
    mutate(pred=median(model(dbh,
                             par_mod$alpha,
                             par_mod$beta)),
           q05=quantile(model(dbh,
                              par_mod$alpha,
                              par_mod$beta),
                        probs = 0.05),
           q95=quantile(model(dbh,
                              par_mod$alpha,
                              par_mod$beta),
                        probs = 0.95)) |>
    mutate(mod="loglin") |> 
    ungroup()
  
  ## LOG - LOG
  
  ### Fit
  if(!file.exists(file.path(folder,"HD_loglog.rdata"))){
    mod_HD_loglog=stan(file="stan/model_HD_loglog.stan", # stan program
                       data = data_HD,         # dataset
                       warmup = 1000,          # number of warmup iterations per chain
                       iter = 2000)   
    save(mod_HD_loglog,file=file.path(folder,"HD_loglog.rdata"))
  }else{
    load(file.path(folder,"HD_loglog.rdata"))
  }
  
  ### Prediction 
  model <- function(x,alpha,beta) {
    exp(alpha + beta * log(x))}
  par_mod<- rstan::extract(mod_HD_loglog)
  loglog<-summary(mod_HD_loglog)[[1]]
  data.log<-data.mod |> 
    rowwise() |> 
    mutate(pred=median(model(dbh,
                             par_mod$alpha,
                             par_mod$beta)),
           q05=quantile(model(dbh,
                              par_mod$alpha,
                              par_mod$beta),
                        probs = 0.05),
           q95=quantile(model(dbh,
                              par_mod$alpha,
                              par_mod$beta),
                        probs = 0.95)) |>
    mutate(mod="loglog") |> 
    ungroup()
  
  ## WEIBULL
  
  ### Fit
  if(!file.exists(file.path(folder,"HD_weibull.rdata"))){
    mod_HD_weibull=stan(file="stan/model_HD_weibull.stan", # stan program
                       data = data_HD,         # dataset
                       warmup = 1000,          # number of warmup iterations per chain
                       iter = 2000)   
    save(mod_HD_weibull,file=file.path(folder,"HD_weibull.rdata"))
  }else{
    load(file.path(folder,"HD_weibull.rdata"))
  }
  
  ### Prediction
  model <- function(x,alpha,beta) {
    alpha *(1-exp(-x/beta))}
  par_mod<- rstan::extract(mod_HD_weibull)
  weibull<-summary(mod_HD_weibull)[[1]]
  data.weibull<-data.mod |> 
    rowwise() |> 
    mutate(pred=median(model(dbh,
                             par_mod$alpha,
                             par_mod$beta)),
           q05=quantile(model(dbh,
                              par_mod$alpha,
                              par_mod$beta),
                        probs = 0.05),
           q95=quantile(model(dbh,
                              par_mod$alpha,
                              par_mod$beta),
                        probs = 0.95)) |>
    mutate(mod="weibull") |> 
    ungroup()
  
  ## MICHAELIS - MENTEN
  
  ### Fit
  if(!file.exists(file.path(folder,"HD_micmen.rdata"))){
    mod_HD_micmen=stan(file="stan/model_HD_micmen.stan", # stan program
                        data = data_HD,         # dataset
                        warmup = 1000,          # number of warmup iterations per chain
                        iter = 2000)   
    save(mod_HD_micmen,file=file.path(folder,"HD_micmen.rdata"))
  }else{
    load(file.path(folder,"HD_micmen.rdata"))
  }
  
  ### Predictions
  model <- function(x,alpha,beta) {
    (alpha *x)/((1/beta)+x)}
  par_mod<- rstan::extract(mod_HD_micmen)
  micmen<-summary(mod_HD_micmen)[[1]]
  data.micmen<-data.mod |> 
    rowwise() |> 
    mutate(pred=median(model(dbh,
                             par_mod$alpha,
                             par_mod$beta)),
           q05=quantile(model(dbh,
                              par_mod$alpha,
                              par_mod$beta),
                        probs = 0.05),
           q95=quantile(model(dbh,
                              par_mod$alpha,
                              par_mod$beta),
                        probs = 0.95)) |>
    mutate(mod="micmen") |> 
    ungroup()
  
  data_shape=rbind(data.lin,
                   data.log,
                   data.weibull,
                   data.micmen)
  mod_shape=list(loglin=loglin,
                 loglog=loglog,
                 weibull=weibull,
                 micmen=micmen)
  
  return(list(data_shape=data_shape,
              mod_shape=mod_shape))
}