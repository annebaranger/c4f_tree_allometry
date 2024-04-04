# libraries
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

# load models

load("mod_cov_select/mod_nul.rdata")
post_nul<-rstan::extract(HD_nul)
rm(HD_nul)
gc()
load("mod_cov_select/complete.rdata")
post_complete<-rstan::extract(HD_complete)
rm(HD_complete)
gc()
load("mod_cov_select/nul_ori.rdata")
post_ori<-rstan::extract(HD_nul_cof)
rm(HD_nul_cof)
gc()
load("mod_cov_select/nul_sys.rdata")
post_sys<-rstan::extract(HD_nul_cof)
rm(HD_nul_cof)
gc()
load("mod_cov_select/nul_systori.rdata")
post_systori<-rstan::extract(HD_nul_cof)
rm(HD_nul_cof)
gc()


tar_load(mod.data)
library(BIOMASS)
dataWD <- getWoodDensity(
  genus = mod.data$genusCorr,
  species = mod.data$speciesCorr,
  stand = mod.data$id_plot)
dataHchave <- retrieveH(
  D = mod.data$dbh,
  coord = mod.data[, c("long", "lat")]
)

## PREDICTIONS ##
#################
model_system <- function(x,alpha_sys,beta_sys,gamma_plot,gamma_species) {
  gamma_species *gamma_plot * alpha_sys * x /
    (beta_sys + x)}
model_cov <- function(x,ba,precmin,alpha_sys,beta_sys,beta_ba,beta_precmin,gamma_plot,gamma_species) {
  gamma_species *gamma_plot * alpha_sys * x /
    (beta_sys*(ba^beta_ba) * (precmin^beta_precmin) + x)}
mod.data<-mod.data |> 
  mutate(Hnul=NA_real_,
         Hsystem=NA_real_,
         Hsystori=NA_real_,
         Hori=NA_real_,
         Hcov=NA_real_)
n=dim(mod.data)[1]
for (i in 1:n){
  print(i)
  sys=mod.data$sys[i]
  systori=mod.data$systori[i]
  ori=mod.data$ori[i]
  # plot=mod.data$num_plot[i]
  # species=mod.data$num_species[i]
  mod.data$Hnul[i]=median(model_system(mod.data$dbh[i],post_nul$alpha_0,post_nul$beta_0,1,1))#post_system$gamma_plot[,plot],post_system$gamma_sp[,species]))
  mod.data$Hsystem[i]=median(model_system(mod.data$dbh[i],post_sys$alpha[,sys],post_sys$beta[,sys],1,1))#post_system$gamma_plot[,plot],post_system$gamma_sp[,species]))
  mod.data$Hsystori[i]=median(model_system(mod.data$dbh[i],post_systori$alpha[,systori],post_systori$beta[,systori],1,1))#post_systori$gamma_plot[,plot],post_systori$gamma_sp[,species]))
  mod.data$Hori[i]=median(model_system(mod.data$dbh[i],post_ori$alpha[,ori],post_ori$beta[,ori],1,1))
  mod.data$Hcov[i]=median(model_cov(mod.data$dbh[i],mod.data$ba_tot[i],mod.data$bio17[i],post_complete$alpha[,systori],post_complete$beta[,systori],post_complete$beta_ba,post_complete$beta_precmin,1,1))#post_systori_cov$gamma_plot[,plot],post_systori_cov$gamma_sp[,species]))
}

mod.data=cbind(mod.data,Hchave=dataHchave$H)  

## Comparison of height prediction ##
#####################################
mod.data%>% 
  ggplot(aes(x=log(H), y=log(Hcov),color=system)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope=1,color="darkred") +
  labs(x="H observed (log)",
       y="H predicted (log)")

mod.data %>% 
  #filter(dbh<60) %>% 
  mutate(id=row_number(),
         Esystem=(H-Hsystem)/H,
         Esystori=(H-Hsystori)/H,
         Ecov=(H-Hcov)/H,
         Echave=(H-Hchave)/H) %>% 
  ggplot()+
  geom_density(aes(Esystori),color="green")+
  geom_density(aes(Ecov),color="red")+
  geom_density(aes(Esystem),color="blue")+
  geom_density(aes(Echave),color="black")+
  xlim(-5,5)
# geom_point(aes(id,Esystori),color="green",size=0.2)+
# geom_point(aes(id,Ecov),color="red",size=0.2)+
# geom_point(aes(id,Esystem),color="blue",size=0.2)

data_chave

AGBobs<- summaryByPlot(
  computeAGB(D = mod.data$dbh,
             WD = dataWD$meanWD,
             H = mod.data$H),
  mod.data$id_plot)
AGBchave<- summaryByPlot(
  computeAGB(D = mod.data$dbh,
             WD = dataWD$meanWD,
             H = mod.data$Hchave),
  mod.data$id_plot)
AGBnul<-summaryByPlot(
  computeAGB(D = mod.data$dbh,
             WD = dataWD$meanWD,
             H = mod.data$Hnul),
  mod.data$id_plot)
AGBsys<- summaryByPlot(
  computeAGB(D = mod.data$dbh,
             WD = dataWD$meanWD,
             H = mod.data$Hsystem),
  mod.data$id_plot)
AGBori<- summaryByPlot(
  computeAGB(D = mod.data$dbh,
             WD = dataWD$meanWD,
             H = mod.data$Hori),
  mod.data$id_plot)
AGBsystori<- summaryByPlot(
  computeAGB(D = mod.data$dbh,
             WD = dataWD$meanWD,
             H = mod.data$Hsystori),
  mod.data$id_plot)
AGBcov<- summaryByPlot(
  computeAGB(D = mod.data$dbh,
             WD = dataWD$meanWD,
             H = mod.data$Hcov),
  mod.data$id_plot)

# create df with agb estimates
AGB=data.frame(plot=AGBobs$plot,
               AGBobs=AGBobs$AGB,
               AGBchave=AGBchave$AGB,
               AGBnul=AGBnul$AGB,
               AGBsys=AGBsys$AGB,
               AGBori=AGBori$AGB,
               AGBsystori=AGBsystori$AGB,
               AGBcov=AGBcov$AGB
               ) %>% 
  left_join(mod.data[,c("id_plot","area_plot.ha","system")],by=c("plot"="id_plot")) %>% 
  mutate_at(vars(AGBobs,AGBchave,AGBnul,AGBori,AGBsys,AGBsystori,AGBcov),~./area_plot.ha) %>% 
  mutate(
         Echave=100*(AGBchave-AGBobs)/AGBobs,
         Enul=100*(AGBnul-AGBobs)/AGBobs,
         Eori=100*(AGBori-AGBobs)/AGBobs,
         Esys=100*(AGBsys-AGBobs)/AGBobs,
         Esystori=100*(AGBsystori-AGBobs)/AGBobs,
         Ecov=100*(AGBcov-AGBobs)/AGBobs
         # Echave=AGBchave-AGBobs,
          # Enul=AGBnul-AGBobs,
          # Eori=AGBori-AGBobs,
          # Esys=AGBsys-AGBobs,
          # Esystori=AGBsystori-AGBobs,
          # Ecov=AGBcov-AGBobs
         )

# distribution of errors by model
AGB %>% 
  select(plot,area_plot.ha,Echave,Enul, Eori,Esys,Esystori) %>% #,Ecov
  reshape2::melt(id.vars=c("plot","area_plot.ha")) %>% 
  ggplot(aes(value,color=variable))+geom_density()+
  labs(color="Model")

# distribution of errors by model and system
AGB %>% 
  select(plot,area_plot.ha,system,Echave,Enul, Eori,Esys,Esystori,Ecov) %>% 
  mutate(system=recode_factor(system,
                              `forest`="Forest",
                              `secondary_forest`="Secondary forest",
                              `plantation`="Plantation",
                              `agroforestry`="Agroforestry system")) %>% 
  rename("Chave"=`Echave`,
         "Nul"=`Enul`,
         "System only"=`Esys`,
         "Origin only"=`Eori`,
         "System*origin"=`Esystori`) %>% #"Complete model"=`Ecov`
  reshape2::melt(id.vars=c("plot","area_plot.ha","system")) %>% 
  ggplot(aes(value,color=system))+geom_density()+facet_wrap(~variable)+
  scale_color_manual(values=c("forestgreen","royalblue2","red3","goldenrod1"))+
  geom_vline(xintercept=0,color="black")+
  labs(x="Biomass estimation error (%)",
       color="System")

# box plots by system and models
library(multcomp)
data_for_anova <- AGB %>%
  dplyr::select(plot, area_plot.ha, system, Echave, Enul, Esys, Ecov) %>%
  mutate(system = recode_factor(system,
                                `forest` = "Forest",
                                `secondary_forest` = "Secondary forest",
                                `plantation` = "Plantation",
                                `agroforestry` = "Agroforestry system")) %>%
  rename("Chave" = `Echave`, "Nul" = `Enul`, "System*origin" = `Esys`, "Complete" = Ecov) %>%
  pivot_longer(cols = -c("plot", "area_plot.ha", "system"), names_to = "Model", values_to = "Error") %>%
  mutate(Model = factor(Model, levels = c("Chave", "Nul", "System*origin", "Complete")))

anova <-  aov(Error ~ Model*system, data = data_for_anova)


# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

# table with factors and 3rd quantile
dt <- group_by(data_for_anova, Model,system) %>%
  summarise(q95=quantile(Error,probs=0.99)[[1]]) %>%
  arrange(desc(q95)) |> ungroup() 
cld <- as.data.frame.list(cld$`Model:system`)
dt$cld <- cld$Letters
print(dt)

AGB %>% 
  dplyr::select(plot,area_plot.ha,system,Echave,Enul, Esys,Ecov) %>% #,Eori,Esystori
  mutate(system=recode_factor(system,
                              `forest`="Forest",
                              `secondary_forest`="Secondary forest",
                              `plantation`="Plantation",
                              `agroforestry`="Agroforestry system")) %>% 
  rename("Chave"=`Echave`,
         "Nul"=`Enul`,
         "System*origin"=`Esys`,
         # "Origin only"=`Eori`,
         # "System*origin"=`Esystori`,
         "Complete"=Ecov) |> 
  pivot_longer(cols = -c("plot","area_plot.ha","system")) |> 
  mutate(name=factor(name,levels=c("Chave","Nul","System*origin","Complete"))) |> 
  ggplot()+
  geom_boxplot(aes(name,value),outlier.colour = NA)+
  geom_text(data=dt,aes(x=Model,y=q95+20,label=cld))+
  geom_hline(yintercept = 0)+
  facet_wrap(~system)+
  labs(y="Biomass prediction error  ",x="Model")

# agb obs by system
AGB %>% 
  select(plot,system,AGBobs) %>% 
  ggplot(aes(AGBobs,color=system))+geom_density()+
  ylim(0,0.025)+xlim(0,1000)

# fig tuned for m2 report
AGB %>% 
  select(plot,area_plot.ha,system,Echave,Esys,Esystori,Ecov) %>% 
  mutate(system=recode_factor(system,
                              `forest`="Forêt",
                              `secondary_forest`="Forêt secondaire",
                              `plantation`="Plantation",
                              `agroforestry`="Agroforesterie")) %>% 
  rename("Modèle de Chave (2014)"=`Echave`,
         "Modèle système"=`Esys`,
         "Modèle système*origine"=`Esystori`,
         "Modèle complet"=`Ecov`) %>% 
  reshape2::melt(id.vars=c("plot","area_plot.ha","system")) %>% 
  ggplot(aes(value,color=system))+geom_density()+facet_wrap(~variable)+
  scale_color_manual(values=c("forestgreen","royalblue2","red3","goldenrod1"))+
  geom_vline(xintercept=0,color="black")+
  theme(#axis.title.y = element_blank(),
    legend.position = "bottom")+
  labs(x="Erreur sur la prédiction de biomasse (%)",
       y="Densité",
       color="Système")+
  xlim(-200,200)#+
  # ggsave("A08-figures/biomasspredictions-prez.png",width=7,height=4.3,units = "in")


## plot je sais plus prq

plotBA=data_tree %>%
  left_join(data_env %>% select(id_plot,long,lat,area_plot.ha),by=c("id_plot")) %>%
  group_by(id_plot) %>%
  summarise(dbh_mean=mean(dbh))

AGB %>% 
  left_join(plotBA, by=c("plot"="id_plot"))%>% 
  select(plot,system,dbh_mean,Echave,Esys,Esystori,Ecov) %>% 
  reshape2::melt(id.vars=c("plot","system","dbh_mean")) %>% 
  ggplot(aes(dbh_mean,value,color=system))+geom_point()+facet_wrap(~variable)



## trajectories

trajectory=mod.data %>% 
  dplyr::select(system,origin,systori) %>% 
  unique() %>% 
  mutate(systori=as.numeric(systori)) %>% 
  tidyr::crossing(dbh = seq(0, 200, 1))
n=dim(trajectory)[1]
for (i in 1:n){
  systori=trajectory$systori[i]
  trajectory$pred[i]=median(median(model_system(trajectory$dbh[i],post_systori$alpha[,systori],post_systori$beta[,systori],1,1)))
  trajectory$mu5[i]=quantile(model_system(trajectory$dbh[i],post_systori$alpha[,systori],post_systori$beta[,systori],1,1),prob=0.05)
  trajectory$mu95[i]=quantile(model_system(trajectory$dbh[i],post_systori$alpha[,systori],post_systori$beta[,systori],1,1),prob=0.95)
}


trajectory %>% 
  mutate(system=recode_factor(system,
                              `forest`="Forest",
                              `secondary_forest`="Secondary forest",
                              `plantation`="Plantation",
                              `agroforestry`="Agroforestry system")) %>% 
  ggplot()+
  geom_point(data=mod.data %>% 
               mutate(system=recode_factor(system,
                                           `forest`="Forest",
                                           `secondary_forest`="Secondary forest",
                                           `plantation`="Plantation",
                                           `agroforestry`="Agroforestry system")),
             aes(x=dbh,y=H,color=system),size=0.5,alpha=.40)+
  geom_ribbon(aes(x=dbh,ymin=mu5,ymax=mu95,fill=system), alpha = .25)+
  geom_line(aes(x=dbh,y=pred,col=system))+ 
  scale_y_log10()+
  scale_fill_manual(values=c("forestgreen","royalblue2","red3","goldenrod1"))+
  scale_color_manual(values=c("forestgreen","royalblue2","red3","goldenrod1"))+
  #theme(legend.position=c(0.6,0.15))+
  facet_wrap(~origin)+
  labs(color="System",fill="System")

tar_read(comp.cofactor) |> 
  as.data.frame() |> 
  tibble::rownames_to_column(var="model") |>
  mutate(model=case_match(model,
                          "loo.nul"~"Nul",
                          "loo.ori"~"Origin only",
                          "loo.sys"~"System only",
                          "loo.systori"~"Sytem origin"),
         model=factor(model,levels=c("Sytem origin","System only","Origin only","Nul")))|> 
  ggplot(aes(elpd_diff,model))+geom_point()+labs(y="",x="Expected pointwise \n log likelihood (in relative)")


tar_read(comp.covar) |> 
  as.data.frame() |> 
  tibble::rownames_to_column(var="model") |>
  mutate(model=case_match(model,
                          "loo.mpdq"~"MPDQ",
                          "loo.complete"~"MPDQ + density",
                          "loo.map"~"MAP",
                          "loo.mat"~"MAT",
                          "loo.mtwm"~"MTWM",
                          "loo.vcomp"~"Density",
                          "loo.ba"~"BA"),
         model=factor(model,levels=c("MPDQ","MPDQ + density","MAP","MAT","MTWM","Density","BA")))|> 
  ggplot(aes(elpd_diff,model))+geom_point()+labs(y="",x="Expected pointwise \n log likelihood (in relative)")
