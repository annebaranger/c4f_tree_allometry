library(targets)
library(dplyr)
library(terra)
library(stringr)
library(ggplot2)
tar_load(c(data_bhkamani,
           data_nguessan,
           data_tene,
           data_esanial,
           data_lataha,
           data_mopri.sangoue,
           data_oibt,
           data_akouassi))



data_env=bind_rows(data_bhkamani$plot,
                   data_nguessan$plot,
                   data_tene$plot) |> 
  mutate(plot=as.factor(plot)#,
         # prec=as.character(prec)
  ) |> 
  bind_rows(data_esanial$plot,
            data_lataha$plot,
            data_mopri.sangoue$plot,
            data_oibt$plot |> select(-prec),
            data_akouassi$plot) |> 
  select(database,id_plot,X,Y,lat,long,area_plot.ha,system,age) |> 
  mutate(system=case_when(system=="secondary_forest"&age<=20~"young_sf",
                          system=="secondary_forest"&age>20~"old_sf",
                          TRUE~system))
# summary(data_env)
# str(data_env)

# Load bioclim data 
files=paste0("data/chelsa/",list.files("data/chelsa/",pattern=".tif"))
names(files)=paste0("bio",str_sub(files,start=-6,end=-5))
clim=rast(lapply(names(files),
                 function(z){
                   r=rast(files[[z]])
                   names(r)=z
                   return(r)
                 }
)
)
# extract data for each point
data_env=cbind(data_env,
               terra::extract(clim,
                              y=data.frame(x=data_env$long,
                                           y=data_env$lat))[,-1])
plot=data_env
tree<-tar_read(tree) |> 
  select(-system) |> 
  left_join(data_env[,c("id_plot","system")],by="id_plot")

tree |> filter(system%in%c("young_sf","old_sf","forest")) |> 
  ggplot(aes(dbh,H,color=system))+
  geom_point()+
  geom_smooth(method="gam")+
  facet_grid(~origin)


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
  # filter(is.na(ba_tot)==FALSE) %>%
  filter(is.na(bio01)==FALSE) %>%
  filter(!is.na(origin)) |> 
  filter(!is.na(system)) |> 
  mutate(system=ordered(system,levels=c("forest","old_sf","young_sf","plantation","agroforestry")),
         sys=as.numeric(system),
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


data_tree |> filter(system%in%c("young_sf","old_sf","forest")) |> 
  ggplot(aes(dbh,H,color=system))+
  geom_point()+
  geom_smooth(method="gam")+
  facet_grid(~origin)

sub.mod.data.3=data_tree


#%%%%%%% test of models with stratified sampling

tar_load(sub.mod.data.3)
sub.mod.data.3=sub.mod.data.3 |> 
  mutate(dbh_cat=cut(dbh,
                     breaks=c(-Inf,10,20,30,40,50,60,70,Inf))) 
sub.mod.data.3|> 
  ggplot(aes(dbh_cat,fill=system)) +
  geom_histogram(stat="count")



n_samp<-as.data.frame(table(sub.mod.data.3[,c("system","dbh_cat")])[,8]) |> 
  tibble::rownames_to_column()
colnames(n_samp)=c("system","n_samp")
data_strat=sub.mod.data.3 |> 
  dplyr::select(id_tree,id_plot,g_s,system,origin,systori,dbh,H,dbh_cat,lat,long) |> 
  left_join(n_samp) |> 
  group_by(system,dbh_cat) |> 
  sample_n(size=unique(n_samp),replace=TRUE)

table(data_strat[,c("system","dbh_cat")])

data_strat |> 
  ggplot(aes(dbh,H,color=origin))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~system)

sub.mod.data.1 |> ggplot(aes(dbh,H,color=origin))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~system)

cof="systori"
data_systori = list(
  N = dim(data_strat)[1],
  p =nlevels(data_strat$id_plot),
  sp=nlevels(data_strat$g_s),
  ncof=nlevels(as.factor(data_strat[[cof]])),
  plot=as.numeric(data_strat$id_plot),
  species=as.numeric(data_strat$g_s),
  cof=as.numeric(as.factor(data_strat[[cof]])),
  H = data_strat$H ,
  dbh=data_strat$dbh
)
library(rstan)
HD_systori=stan(file="stan/model_cov_nul.stan",
                data=data_systori,
                warmup = 500,
                iter=1000,
                include = FALSE,
                pars=c("gamma_plot","gamma_sp","log_lik"),
                core=3)  
library(shinystan)

data_strat_2=sub.mod.data.3 |> 
  dplyr::select(id_tree,id_plot,g_s,system,origin,systori,dbh,H,dbh_cat,lat,long) |> 
  sample_n(size=dim(data_strat)[1],replace=TRUE)
data_systori_2 = list(
  N = dim(data_strat_2)[1],
  p =nlevels(data_strat_2$id_plot),
  sp=nlevels(data_strat_2$g_s),
  ncof=nlevels(as.factor(data_strat_2[[cof]])),
  plot=as.numeric(data_strat_2$id_plot),
  species=as.numeric(data_strat_2$g_s),
  cof=as.numeric(as.factor(data_strat_2[[cof]])),
  H = data_strat_2$H ,
  dbh=data_strat_2$dbh
)
HD_systori_2=stan(file="stan/model_cov_nul.stan",
                data=data_systori_2,
                warmup = 500,
                iter=1000,
                include = FALSE,
                pars=c("gamma_plot","gamma_sp","log_lik"),
                core=3)  
summary(HD_systori)$summary[,1]
summary(HD_systori_2)$summary[,1]
data_strat |>ungroup() |>  select(systori,system,origin) |>  unique()
