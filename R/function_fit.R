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
                   frac=1,
                   ba_require=TRUE){
  if(ba_require){
    tree<-tree |> 
      filter(!is.na(ba_tot))
  }
  if(frac!=1){
    tree<-tree |> 
      sample_frac(size=frac)
  }
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
  data.mod<-data.mod |> 
    filter(!is.na(H),!is.na(dbh))
  
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




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - model covariable selection ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' select relevant cofactor 
#' @note cofactor effect is set on both alpha and beta
#' @param tree.merge file where trees were merges
#' @param plot.merge Idem for plots
#' @param folder folder path

cofactor.select<-function(mod.data,
                          folder="mod_cof_select"){
  if(!dir.exists(folder)){
    dir.create(path=folder)
  }
  
  mod.bic=data.frame(model=c("nul","sys","sys","ori","ori","systori","systori"),
                     effect=c("nul",rep(c("a","b"),3)),
                     npar=NA,
                     lp=NA)
  mod.data<-mod.data |> 
    filter(!is.na(H),!is.na(dbh),!is.na(sys),!is.na(ori),!is.na(systori))
  # Nul model
  data_nul = list(
    N = dim(mod.data)[1],
    p =nlevels(mod.data$id_plot),
    sp=nlevels(mod.data$g_s),
    plot=mod.data$num_plot,
    species=mod.data$num_species,
    H = mod.data$H ,
    dbh=mod.data$dbh
  )
  
  if(!file.exists(file.path(folder,"cof_nul.rdata"))){
    HD_cof_nul=stan(file="stan/model_cof_nul.stan", # stan program
                    data = data_nul,         # dataset
                    warmup = 1000,          # number of warmup iterations per chain
                    iter = 2000,
                    cores = 2)   
    save(HD_cof_nul,file=file.path(folder,"cof_nul.rdata"))
  }else{
    load(file.path(folder,"cof_nul.rdata"))
  }
  mod.bic[mod.bic$model=="nul",c("npar","lp")]=
    c(data_nul$p+data_nul$sp+5,
      summary(HD_cof_nul)$summary["lp__",1])
  
  cofactors=mod.data %>% 
    select(sys,ori,systori)
  
  for (i in 1:dim(cofactors)[2]){
    print(i)
    cof=colnames(cofactors)[i]
    data_cof = list(
      N = dim(mod.data)[1],
      p =nlevels(mod.data$id_plot),
      sp=nlevels(mod.data$g_s),
      cof=max(cofactors[,i]),
      plot=mod.data$num_plot,
      species=mod.data$num_species,
      cofactor=pull(cofactors,colnames(cofactors)[i]),
      H = mod.data$H ,
      dbh=mod.data$dbh
    )
    if(!file.exists(file.path(folder,paste0("cof_a_",cof,".rdata")))){
      HD_cofa=stan(file="stan/model_cof_a.stan",
                   data=data_cof,
                   warmup = 1000,
                   iter=2000,
                   core=4)
      save(HD_cofa,file=file.path(folder,paste0("cof_a_",cof,".rdata")))
    }else{
      load(file.path(folder,paste0("cof_a_",cof,".rdata")))
    }
    mod.bic[mod.bic$model==cof&mod.bic$effect=="a",c("npar","lp")]=
      c(data_cof$p+data_cof$sp+data_cof$cof+5,
        summary(HD_cofa)$summary["lp__",1])
    
    if(!file.exists(file.path(folder,paste0("cof_b_",cof,".rdata")))){
      HD_cofb=stan(file="stan/model_cof_b.stan",
                   data=data_cof,
                   warmup = 1000,
                   iter=2000,
                   core=4)
      save(HD_cofb,file=file.path(folder,paste0("cof_b_",cof,".rdata")))
    }else{
      load(file.path(folder,paste0("cof_b_",cof,".rdata")))
    }
    mod.bic[mod.bic$model==cof&mod.bic$effect=="b",c("npar","lp")]=
      c(data_cof$p+data_cof$sp+data_cof$cof+5,
        summary(HD_cofb)$summary["lp__",1])

  }
  mod.bic<-mod.bic |> 
    mutate(bic=log(dim(mod.data)[1])*npar-2*lp)

}


#' select relevant covariables 
#' @note cofactor effect is set on both alpha and beta
#' @param sub.mod.data input dataset, without plots for which BA is not available
#' @param folder folder path

covariable.select<-function(sub.mod.data,
                            folder="mod_cov_select"){
  if(!dir.exists(folder)){
    dir.create(path=folder)
  }
  
  # mod.bic=data.frame(model=c("nul","sys","sys","ori","ori","systori","systori"),
  #                    effect=c("nul",rep(c("a","b"),3)),
  #                    npar=NA,
  #                    lp=NA)
  mod.data<-sub.mod.data |> 
    filter(!is.na(H),!is.na(dbh),!is.na(sys),!is.na(ori),!is.na(systori),
           !is.na(bio01),!is.na(bio05),!is.na(bio12),!is.na(bio17),
           !is.na(ba_tot),!is.na(n_tree10))
  
  # Nul model
  data_nul = list(
    N = dim(mod.data)[1],
    p =nlevels(mod.data$id_plot),
    sp=nlevels(mod.data$g_s),
    plot=mod.data$num_plot,
    species=mod.data$num_species,
    H = mod.data$H ,
    dbh=mod.data$dbh
  )
  
  if(!file.exists(file.path(folder,"nul.rdata"))){
    HD_nul=stan(file="stan/model_cof_nul.stan", # stan program
                    data = data_nul,         # dataset
                    warmup = 1000,          # number of warmup iterations per chain
                    iter = 2000,
                    cores = 4)   
    save(HD_nul,file=file.path(folder,"nul.rdata"))
  }else{
    load(file.path(folder,"nul.rdata"))
  }

  cofactors=c("sys","ori","systori")
  
  for (cof in cofactors){
    print(cof)
    data_cof = list(
      N = dim(mod.data)[1],
      p =nlevels(mod.data$id_plot),
      sp=nlevels(mod.data$g_s),
      ncof=max(mod.data[,cof]),
      plot=mod.data$num_plot,
      species=mod.data$num_species,
      cof=pull(mod.data,cof),
      H = mod.data$H ,
      dbh=mod.data$dbh
    )
    if(!file.exists(file.path(folder,paste0("nul_",cof,".rdata")))){
      HD_nul_cof=stan(file="stan/model_cov_nul.stan",
                   data=data_cof,
                   warmup = 1000,
                   iter=2000,
                   core=4)
      save(HD_nul_cof,file=file.path(folder,paste0("nul_",cof,".rdata")))
    }else{
      load(file.path(folder,paste0("nul_",cof,".rdata")))
    }
  }
  
  
  covariables=c("ba_tot","n_tree10","bio01","bio05","bio12","bio17")
  
  for (cov in covariables){
    print(cov)
    data_systori_cov=list(
      N = dim(mod.data)[1],
      p =nlevels(mod.data$id_plot),
      sp=nlevels(mod.data$g_s),
      so=nlevels(as.factor(mod.data$systori)),
      plot=mod.data$num_plot,
      species=mod.data$num_species,
      systori=mod.data$systori,
      H = mod.data$H,
      dbh=mod.data$dbh,
      cov=pull(mod.data,cov)
    )
    if(!file.exists(file.path(folder,paste0("so_",cov,".rdata")))){
      HD_so_cov=stan(file="stan/model_cov_nul.stan",
                   data=data_systori_cov,
                   warmup = 1000,
                   iter=2000,
                   core=4)
      save(HD_so_cov,file=file.path(folder,paste0("so_",cov,".rdata")))
    }else{
      load(file.path(folder,paste0("so_",cov,".rdata")))
    }
  }
 
  
}

















# data_nul= list(
#   N = dim(mod.data)[1],
#   p =nlevels(mod.data$id_plot),
#   sp=nlevels(mod.data$g_s),
#   so=nlevels(as.factor(mod.data$systori)),
#   plot=mod.data$num_plot,
#   species=mod.data$num_species,
#   systori=as.numeric(mod.data$systori),
#   H = mod.data$H,
#   dbh=mod.data$dbh
# )
# 
# if(!file.exists(file.path(folder,"cov_nul.rdata"))){
#   HD_cov_nul=stan(file="stan/model_cof_nul.stan", # stan program
#                   data = data_nul,         # dataset
#                   warmup = 1000,          # number of warmup iterations per chain
#                   iter = 2000)   
#   save(HD_cov_nul,file=file.path(folder,"cov_nul.rdata"))
# }else{
#   load(file.path(folder,"cov_nul.rdata"))
# }
# 

