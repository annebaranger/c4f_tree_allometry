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

#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
#' @author Julien Barrere
create_dir_if_needed <- function(file.in){
  
  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

model_nul <- function(x,alpha,beta,gamma_sp,gamma_plot) {
  gamma_sp*gamma_plot*(alpha * x)/ 
    (beta + x)
}

model_origin <- function(x,alpha_ori,beta_ori,gamma_sp,gamma_plot){
  gamma_sp*gamma_plot*(alpha_ori * x)/ 
    (beta_ori+x)
}

model_system <- function(x,alpha_sys,beta_sys,gamma_sp,gamma_plot){
  gamma_sp*gamma_plot*(alpha_sys * x)/ 
    (beta_sys+x)
}

model_systori <- function(x,alpha_systori,beta_systori,gamma_sp,gamma_plot){
  gamma_sp*gamma_plot*(alpha_systori * x)/ 
    (beta_systori + x)
}

model_complete <- function(x,ba,precmin,alpha_sys,beta_sys,beta_ba,beta_precmin){
   alpha_sys * x /
    (beta_sys*(ba^beta_ba) * (precmin^beta_precmin) + x)}

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
    filter(!is.na(origin)) |> 
    filter(!is.na(system)) |> 
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
#' @note fit models to be compared
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
           !is.na(ba_tot),!is.na(v_compet))
  
  # Nul model
  print("nul model")
  data_nul = list(
    N = dim(mod.data)[1],
    p =nlevels(mod.data$id_plot),
    sp=nlevels(mod.data$g_s),
    plot=mod.data$num_plot,
    species=mod.data$num_species,
    H = mod.data$H ,
    dbh=mod.data$dbh
  )
  
  if(!file.exists(file.path(folder,"mod_nul.rdata"))){
    HD_nul=stan(file="stan/model_cof_nul.stan", # stan program
                    data = data_nul,         # dataset
                    warmup = 1000,          # number of warmup iterations per chain
                    iter = 2000,
                    include = FALSE,
                    pars=c("gamma_plot","gamma_sp"),
                    core=4)   
    save(HD_nul,file=file.path(folder,"mod_nul.rdata"))
  }else{
    print("model fitted")
    # load(file.path(folder,"mod_nul.rdata"))
  }
  
  # cofactors test

  cofactors=c("sys","ori","systori")
  
  for (cof in cofactors){
    print(cof)
    data_cof = list(
      N = dim(mod.data)[1],
      p =nlevels(mod.data$id_plot),
      sp=nlevels(mod.data$g_s),
      ncof=nlevels(as.factor(mod.data[[cof]])),
      plot=mod.data$num_plot,
      species=mod.data$num_species,
      cof=as.numeric(as.factor(mod.data[[cof]])),
      H = mod.data$H ,
      dbh=mod.data$dbh
    )
    if(!file.exists(file.path(folder,paste0("nul_",cof,".rdata")))){
      HD_nul_cof=stan(file="stan/model_cov_nul.stan",
                   data=data_cof,
                   warmup = 1000,
                   iter=2000,
                   include = FALSE,
                   pars=c("gamma_plot","gamma_sp"),
                   core=4)
      save(HD_nul_cof,file=file.path(folder,paste0("nul_",cof,".rdata")))
    }else{
      print("model fitted")
      # load(file.path(folder,paste0("nul_",cof,".rdata")))
    }
  }
  
  # covariables tests
  covariables=c("ba_tot","v_compet","bio01","bio05","bio12","bio17")
  
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
      HD_so_cov=stan(file="stan/model_systori_cov.stan",
                   data=data_systori_cov,
                   warmup = 1000,
                   iter=2000,
                   include = FALSE,
                   pars=c("gamma_plot","gamma_sp"),
                   core=4)
      save(HD_so_cov,file=file.path(folder,paste0("so_",cov,".rdata")))
    }else{
      print("model fitted")
      # load(file.path(folder,paste0("so_",cov,".rdata")))
    }
  }
 
  # complete model
  data_complete = list(
    N = dim(mod.data)[1],
    p =nlevels(mod.data$id_plot),
    sp=nlevels(mod.data$g_s),
    so=nlevels(as.factor(mod.data$systori)),
    plot=mod.data$num_plot,
    species=mod.data$num_species,
    systori=mod.data$systori,
    H = mod.data$H,
    dbh=mod.data$dbh,
    ba=mod.data$ba_tot,
    precmin=mod.data$bio17
  )
  
  if(!file.exists(file.path(folder,"complete.rdata"))){
    HD_complete=stan(file="stan/model_total_ba_prec.stan", # stan program
                data = data_complete,         # dataset
                warmup = 1000,          # number of warmup iterations per chain
                iter = 2000,
                include = FALSE,
                pars=c("gamma_plot","gamma_sp"),
                core=4)   
    save(HD_complete,file=file.path(folder,"complete.rdata"))
  }else{
    print("model fitted")
    # load(file.path(folder,"complete.rdata"))
  }
  
}



#' compare models
#' @note fit models to be compared
#' @param files.list files of the models to compare
#' @param list.names names of the model, same order as files

select.mod <- function(files.list,
                       list.names) {
  files <- files.list
  loo.list <- list()
  
  for (i in seq_along(files)) {
    new_name <- paste0("mod.",list.names[i])
    load.rename(files[i], new_name)
    loo.list[[new_name]] <- loo(get(new_name))
    rm(list = new_name, envir = .GlobalEnv)
  }
  
  names(loo.list) <- paste0("loo.",list.names)
  comp<-loo_compare(loo.list)
  return(comp)
}


#' Function to load and rename objects
#' @param dir.object directory of the object to load
#' @param new_name names for the object

load.rename<-function(dir.object,new_name){
  # Load the RData file and capture the name of the loaded object
  loaded_names <- load(dir.object)
  
  # Check if only one object was loaded and rename it
  if (length(loaded_names) == 1) {
    original_name <- loaded_names[1]
    # Assign the object to the new name and remove the original object
    assign(new_name, get(original_name),envir = parent.env(as.environment(-1)))
    rm(list = original_name)
  } else {
    cat("More than one object found in the RData file. Please specify which object to rename.")
  }
  return(invisible(NULL))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 4 - spatial cross validation ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' select relevant cofactor 
#' @note cofactor effect is set on both alpha and beta
#' @param tree.merge file where trees were merges
#' @param plot.merge Idem for plots
#' @param folder folder path

make_subdataset<-function(mod.data,
                          folder="rds/"){
  create_dir_if_needed(folder)
  plot.data<-mod.data |> 
    dplyr::select(database,id_plot,system,lat,long) |> 
    unique() 
  km=kmeans(as.matrix(plot.data[c("long","lat")]),centers=7)
  km.center=as.data.frame(km$centers) |> tibble::rownames_to_column(var="cls") |> 
    mutate(cluster_near=NA,
           cls=as.numeric(cls))
  for (i in 1:dim(km.center)[1]){
    print(i)
    if(km.center$lat[i]>quantile(mod.data$lat,probs = 0.99)[[1]]){
      ind=which.min(st_distance(st_as_sf(km.center[i,],
                                         coords=c("long","lat")),
                                st_as_sf(km.center[-i,],
                                         coords=c("long","lat"))
      ))
      km.center$cluster_near[i]=as.numeric(km.center[-i,][ind,"cls"])
    }
  }
  plot.class=cbind(plot.data,
                   cluster=km$cluster) |> 
    left_join(km.center[,c("cls","cluster_near")],by=c("cluster"="cls")) |> 
    mutate(cluster=case_when(!is.na(cluster_near)~cluster_near,
                             TRUE~cluster),
           cluster=as.numeric(as.factor(as.character(cluster)))) |> 
    group_by(cluster) |> 
    mutate(n_clust=n()) |> 
    ungroup()
  unique(plot.class$cluster)
  list.subdataset <- vector("list", length(unique(plot.class$cluster)))
  for (i in unique(plot.class$cluster)){
    print(i)
    plot.select=plot.class|>
      filter(cluster==i)|>
      pull(id_plot)
    subdata=mod.data|>
      filter(id_plot%in%plot.select) |> 
      mutate(g_s=as.factor(as.character(g_s)),
             id_plot=as.factor(as.character(id_plot)))
    saveRDS(subdata,
            file=paste0(folder,"subdataset_",i,".rds"))
    list.subdataset[[i]]=paste0(folder,"subdataset_",i,".rds")
  }
  return(list.subdataset)
}


#' Create blocks for spatial cross-validation
#' @note use blockCV package
#' @param mod.data
#' @param folder folder path

make_blockCV<-function(mod.data,
                       nfold=5){
  point_data <- sf::st_as_sf(mod.data, coords = c("long", "lat"), crs = 4326)
  scv <- cv_cluster(x = point_data,
                    column = "sys", # optional: counting number of train/test records
                    k = nfold)
  
  return(cbind(mod.data,fold=scv$folds_ids))
}

#' Fit model on a subselection of the dataset
#' @param sim.plan table with all model to run on which dataset
#' @param id.sim id of the model x dataset to run
#' @param mod.data.block data with block from blockCV
#' @param folder where to save models

make_model<-function(sim.plan,
                     id.sim,
                     mod.data.block,
                     folder="mod_spatial_cross_valid/"){
  print(id.sim)
  file_data=paste0(folder,"subdata_",sim.plan$fold[id.sim])
  file_mod=paste0(file_data,"/",sim.plan$model[id.sim],".RData")
  if(!dir.exists(folder))dir.create(folder)
  if(!dir.exists(file_data))dir.create(file_data)
  mod.fold=mod.data.block[mod.data.block$fold!=sim.plan$fold[id.sim],]
  
  ## model nul
  if(sim.plan$model[id.sim]=="nul"){
    data_nul = list(
      N = dim(mod.fold)[1],
      p =nlevels(mod.fold$id_plot),
      sp=nlevels(mod.fold$g_s),
      plot=as.numeric(mod.fold$id_plot),
      species=as.numeric(mod.fold$g_s),
      H = mod.fold$H ,
      dbh=mod.fold$dbh
    )
    HD_nul=stan(file="stan/model_cof_nul.stan", # stan program
                    data = data_nul,         # dataset
                    warmup = 1000,          # number of warmup iterations per chain
                    iter = 2000,
                    include = FALSE,
                    pars=c("gamma_plot","gamma_sp","log_lik"),
                    cores = 4)
    save(HD_nul,file=file_mod)
  }
  
  ## model origin
  if(sim.plan$model[id.sim]=="origin"){
    cof="origin"
    data_origin = list(
      N = dim(mod.fold)[1],
      p =nlevels(mod.fold$id_plot),
      sp=nlevels(mod.fold$g_s),
      ncof=nlevels(as.factor(mod.fold[[cof]])),
      plot=as.numeric(mod.fold$id_plot),
      species=as.numeric(mod.fold$g_s),
      cof=as.numeric(as.factor(mod.fold[[cof]])),
      H = mod.fold$H ,
      dbh=mod.fold$dbh
    )
    HD_origin=stan(file="stan/model_cov_nul.stan",
                      data=data_origin,
                      warmup = 1000,
                      iter=2000,
                      include = FALSE,
                      pars=c("gamma_plot","gamma_sp","log_lik"),
                      core=4)
    save(HD_origin,file=file_mod)
  }
  ## model system
  if(sim.plan$model[id.sim]=="system"){
    cof="system"
    data_system = list(
      N = dim(mod.fold)[1],
      p =nlevels(mod.fold$id_plot),
      sp=nlevels(mod.fold$g_s),
      ncof=nlevels(as.factor(mod.fold[[cof]])),
      plot=as.numeric(mod.fold$id_plot),
      species=as.numeric(mod.fold$g_s),
      cof=as.numeric(as.factor(mod.fold[[cof]])),
      H = mod.fold$H ,
      dbh=mod.fold$dbh
    )
    HD_system=stan(file="stan/model_cov_nul.stan",
                   data=data_system,
                   warmup = 1000,
                   iter=2000,
                   include = FALSE,
                   pars=c("gamma_plot","gamma_sp","log_lik"),
                   core=4)
    save(HD_system,file=file_mod)
  }
  ## model system origin
  if(sim.plan$model[id.sim]=="systori"){
    cof="systori"
    data_systori = list(
      N = dim(mod.fold)[1],
      p =nlevels(mod.fold$id_plot),
      sp=nlevels(mod.fold$g_s),
      ncof=nlevels(as.factor(mod.fold[[cof]])),
      plot=as.numeric(mod.fold$id_plot),
      species=as.numeric(mod.fold$g_s),
      cof=as.numeric(as.factor(mod.fold[[cof]])),
      H = mod.fold$H ,
      dbh=mod.fold$dbh
    )
    HD_systori=stan(file="stan/model_cov_nul.stan",
                   data=data_systori,
                   warmup = 1000,
                   iter=2000,
                   include = FALSE,
                   pars=c("gamma_plot","gamma_sp","log_lik"),
                   core=4)
    save(HD_systori,file=file_mod)
  }
  ## model complete
  if(sim.plan$model[id.sim]=="complete"){
    data_complete = list(
      N = dim(mod.fold)[1],
      p =nlevels(mod.fold$id_plot),
      sp=nlevels(mod.fold$g_s),
      so=nlevels(as.factor(mod.fold$systori)),
      plot=as.numeric(mod.fold$id_plot),
      species=as.numeric(mod.fold$g_s),
      systori=as.numeric(as.factor(mod.fold[["systori"]])),
      H = mod.fold$H,
      dbh=mod.fold$dbh,
      ba=mod.fold$ba_tot,
      precmin=mod.fold$bio17
    )
    HD_complete=stan(file="stan/model_total_ba_prec.stan", # stan program
                     data = data_complete,         # dataset
                     warmup = 1000,          # number of warmup iterations per chain
                     iter = 2000,
                     include = FALSE,
                     pars=c("gamma_plot","gamma_sp","log_lik"),
                     core=4)   
    save(HD_complete,file=file_mod)
  }
  return(file_mod)
}

#' Fit model on a subselection of the dataset
#' @param sim.plan table with all model to run on which dataset
#' @param id.sim id of the model x dataset to run
#' @param mod.data.block data with block from blockCV
#' @param folder where to save models

get_prediction<-function(spatial.cross.val,
                         sim.plan,
                         id.sim,
                         mod.data.block){
  print(id.sim)
  file_mod=spatial.cross.val[id.sim]
  mod=sim.plan$model[id.sim]
  mod.train=mod.data.block[mod.data.block$fold!=sim.plan$fold[id.sim],]
  mod.test=mod.data.block[mod.data.block$fold==sim.plan$fold[id.sim],]
  mod_fit=loadRData(file_mod)
  par_mod<- as.data.frame(mod_fit) |> 
    dplyr::select(!matches(c("log_lik","lp__"))) 
  
  ## model nul
  if(mod=="nul"){
    
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
    data_nul = list(
      N = dim(mod.fold)[1],
      p =nlevels(mod.fold$id_plot),
      sp=nlevels(mod.fold$g_s),
      plot=as.numeric(mod.fold$id_plot),
      species=as.numeric(mod.fold$g_s),
      H = mod.fold$H ,
      dbh=mod.fold$dbh
    )
    HD_nul=stan(file="stan/model_cof_nul.stan", # stan program
                data = data_nul,         # dataset
                warmup = 1000,          # number of warmup iterations per chain
                iter = 2000,
                cores = 4)
    save(HD_nul,file=file_mod)
  }
  
  ## model origin
  if(sim.plan$model[id.sim]=="origin"){
    cof="origin"
    data_origin = list(
      N = dim(mod.fold)[1],
      p =nlevels(mod.fold$id_plot),
      sp=nlevels(mod.fold$g_s),
      ncof=nlevels(as.factor(mod.fold[[cof]])),
      plot=as.numeric(mod.fold$id_plot),
      species=as.numeric(mod.fold$g_s),
      cof=as.numeric(as.factor(mod.fold[[cof]])),
      H = mod.fold$H ,
      dbh=mod.fold$dbh
    )
    HD_origin=stan(file="stan/model_cov_nul.stan",
                   data=data_origin,
                   warmup = 1000,
                   iter=2000,
                   include = FALSE,
                   pars=c("gamma_plot","gamma_sp"),
                   core=4)
    save(HD_origin,file=file_mod)
  }
  ## model system
  if(sim.plan$model[id.sim]=="system"){
    cof="system"
    data_system = list(
      N = dim(mod.fold)[1],
      p =nlevels(mod.fold$id_plot),
      sp=nlevels(mod.fold$g_s),
      ncof=nlevels(as.factor(mod.fold[[cof]])),
      plot=as.numeric(mod.fold$id_plot),
      species=as.numeric(mod.fold$g_s),
      cof=as.numeric(as.factor(mod.fold[[cof]])),
      H = mod.fold$H ,
      dbh=mod.fold$dbh
    )
    HD_system=stan(file="stan/model_cov_nul.stan",
                   data=data_system,
                   warmup = 1000,
                   iter=2000,
                   include = FALSE,
                   pars=c("gamma_plot","gamma_sp"),
                   core=4)
    save(HD_system,file=file_mod)
  }
  ## model system origin
  if(sim.plan$model[id.sim]=="systori"){
    cof="systori"
    data_systori = list(
      N = dim(mod.fold)[1],
      p =nlevels(mod.fold$id_plot),
      sp=nlevels(mod.fold$g_s),
      ncof=nlevels(as.factor(mod.fold[[cof]])),
      plot=as.numeric(mod.fold$id_plot),
      species=as.numeric(mod.fold$g_s),
      cof=as.numeric(as.factor(mod.fold[[cof]])),
      H = mod.fold$H ,
      dbh=mod.fold$dbh
    )
    HD_systori=stan(file="stan/model_cov_nul.stan",
                    data=data_systori,
                    warmup = 1000,
                    iter=2000,
                    include = FALSE,
                    pars=c("gamma_plot","gamma_sp"),
                    core=4)
    save(HD_systori,file=file_mod)
  }
  ## model complete
  if(mod=="complete"){
    n=dim(mod.test)[1]
    for (i in 1:n){
      print(i)
      cof=mod.test$systori[i]
      dbh=mod.test$dbh[i]
      ba=mod.test$ba_tot[i]
      prec=mod.test$bio17[i]
      H=model_complete(dbh,ba, prec, 
                       alpha_sys=par_mod[[paste0("alpha[",cof,"]")]],
                       beta_sys=par_mod[[paste0("beta[",cof,"]")]],
                       beta_ba=par_mod[["beta_ba"]],
                       beta_precmin=par_mod[["beta_precmin"]] )
      mod.train$H_med[i]=median(H)
      mod.train$H_q05[i]=quantile(H,probs=0.05)
      mod.train$H_q95[i]=quantile(H,probs=0.95)
    }
  }
}
