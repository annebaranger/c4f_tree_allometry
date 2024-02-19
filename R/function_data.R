#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script with all functions formatting databases
#' relative to trees and plots
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - OIBT data ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param dir.data directory of oibt files
oibt <- function(dir.data="data/oibt"){
  OIBTtree<-read.csv2(file.path(dir.data,"Ivoire_tdc_cor.csv"))|> 
    filter(!is.na(DBH)) |> # erase blank rows
    ## Rename variables for consistency with other database
    rename_with(tolower,.cols=c("DBH","Date","Pointeur")) |> 
    rename(site=Foret,
           name=NomCom,
           name_sc=NomSc,
           n_ech=numEch,
           visee_b=ViseeB,
           visee_h=ViseeH,
           H=Ltot,
           H_fut=Lfut) |>  #formating names of variables
    unique() |>  
    ## Compute new variables database, id_plot, origin and system
    mutate(database=as.factor("OIBT"),
           site=tools::toTitleCase(tolower(gsub(pattern=" |-","_",site))),
           name=tools::toTitleCase(tolower(name)),
           id_plot=as.factor(paste0("OIBT_",site)), #definition of id_plot
           origin=as.factor(case_when(type=="Naturelle"~"remnant",
                                      type=="Plantation"~"planted")), # definition of origin
           system=as.factor(case_when(type=="Naturelle"~"forest",
                                      type=="Plantation"~"plantation"))) |> 
    ## Compute new variable id_tree
    group_by(id_plot) |> 
    mutate(id_tree=as.factor(paste0(id_plot,"_",row_number())),
           density=NA) |> 
    ungroup() |> 
    ## Tidy
    relocate(database,id_plot,id_tree,.before=site)
  
  # Transitory dataset  
  datatree_OIBT=OIBTtree |> 
    ## Correct species
    # mutate(name_sc=recode_factor(as.factor(name_sc),
    #                              `Chlorophora_regia-C_excelsa`="Chlorophora_regia",
    #                              `Lovoa_trichilio\xefdes`="Lovoa_trichiliodes",
    #                              `Lovoa_trichilio\x95des`="Lovoa_trichiliodes")) |> 
    tidyr::separate(name_sc,c("genus","species"),sep="_",remove=TRUE) |> #creating genus and species
    ## Filter inconsistent data
    filter(pointeur!="KOUASSI") |>  # weird data for these samples
    filter(H>0 & H<100)
  
  # original dataset
  OIBTenv<-read.csv2(file.path(dir.data,"env.csv")) |> 
    ## Rename for consistency
    rename(long=Long,
           lat=Lat) |> 
    ## Compute new variable id_plot database and area
    mutate(site=tools::toTitleCase(tolower(gsub(pattern=" |-","_",site))),
           id_plot=as.factor(paste0("OIBT_",site)),
           database=as.factor("OIBT"),
           area_plot.ha=as.numeric(NA)) |> 
    ## get system
    left_join(OIBTtree |> 
                select(id_plot,system) |> 
                distinct(),
              by=c("id_plot")) |> 
    ## Tidy
    relocate(database,id_plot,.before=site) 

  
  ## coordinates conversion into UTM
  crs="+proj=longlat +datum=WGS84"#"+proj=utm +zone=30"
  OIBT_longlat30 <- data.frame(id_plot=OIBTenv$id_plot,long=OIBTenv$long,lat=OIBTenv$lat) |> 
    filter(long>-6) |>
    vect(geom = c("long", "lat"), crs = crs) |> 
    project("+proj=utm +zone=30 ellps=WGS84") |> 
    as.data.frame(geom="XY") |> 
    rename_with(.cols=c("x","y"),
                toupper)
  OIBT_longlat29 <- data.frame(id_plot=OIBTenv$id_plot,long=OIBTenv$long,lat=OIBTenv$lat) |> 
    filter(long<(-6)) |>
    vect(geom = c("long", "lat"), crs = crs) |> 
    project("+proj=utm +zone=29 ellps=WGS84") |> 
    as.data.frame(geom="XY") |> 
    rename_with(.cols=c("x","y"),
                toupper)
  
  ### Gathering the 2 utm zones
  OIBT_longlat=bind_rows(OIBT_longlat29,OIBT_longlat30)
  
  # Transitory
  dataenv_OIBT=OIBTenv |> 
    # select(-long,-lat) |> 
    ##Gathering dataenv with newly calculated coordinates
    left_join(OIBT_longlat,by=c("id_plot")) 
  
  return(list(tree=datatree_OIBT,
              plot=dataenv_OIBT))
}
                       
                       

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - Tene data ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param dir.data directory of tene files
tene <- function(dir.data="data/tene"){
  #original datasets, all trees
  TeneTotal <- read.csv2(file.path(dir.data,"Global_Tene_census.csv")) |> 
    select(-X.1) |> # remove row number
    unique() |>  # remove duplicated rows
    # filter(year>1987) |> # because below this date too many mistakes (duplicata)
    filter(year==2018) |> # because anyhow height measures were performed in 2020
    # select randomly duplicated tree in 2018 year
    group_by(Parcelle_ID,Subplot,Id_Arbre,year) |> 
    sample_n(1) |> 
    ungroup() |>  
    # rename variables
    rename(plot=Parcelle_ID,
           subplot=Subplot,
           n_ech_old=Former_Id_tree,
           genus=Genus,
           sp=Species,
           n_ech=Id_Arbre) |> 
    mutate(species=paste0(substr(tolower(genus),1,8),substr(str_to_title(sp),1,8))) |> 
    # biotic covariates (need all tree on plots)
    group_by(plot) |>  # plot because surface of subplot is not available
    mutate(v_compet=n()/0.2,
           ba_tot=sum(pi*((dbh/200)^2),na.rm=TRUE)/0.2,
           n_tree10=n()) |> 
    ungroup() |>  
    mutate(dbh_cat=cut(dbh,
                       breaks=seq(0,200,by=10))) |> 
    group_by(plot,dbh_cat) |> 
    mutate(h_compet=n()/0.2) |> 
    ungroup() |> 
    filter(!duplicated(n_ech)) 
  
  
  # height measures
  TeneTree <- read.csv2(file.path(dir.data,"teneHD.csv")) |> 
    select(-X) |> 
    # rename variables
    rename(n_ech=IDtree,
           H=H_total,
           H_fut=H_stem,
           subplot=plot)  |> 
    # Compute new variables database, id_plot, origin and system
    mutate(plot=str_sub(subplot, end=-3),
           database=as.factor("Tene"),
           id_plot=as.factor(paste0("Tene_",plot)), # definition of id_plot
           origin=as.factor("remnant"), # definition of origin
           system=as.factor("forest"),
           dbh=100*dbh)|> 
    # Compute new variable id_tree
    group_by(id_plot) |> 
    mutate(id_tree=as.factor(paste0(id_plot,"_",row_number()))) |> 
    ungroup() |> 
    # add species
    left_join(TeneTotal |> 
                select(n_ech,genus,sp,species,dbh,dbh_cat,v_compet,h_compet,ba_tot,Treatment,X,Y),
              by=c("n_ech")) |> 
    # Tidy
    relocate(database,id_plot,id_tree,.before=n_ech) |> 
    relocate(plot,.before=subplot)
  
  
  
  # Tree final
  datatree_Tene=TeneTree |> 
    filter(!is.na(species)) |> 
    filter(abs(dbh.y-dbh.x)<10) |> 
    select(-dbh.y) |> 
    rename(dbh=dbh.x)
  
  # plot n
  plotn_tene<-TeneTotal |> 
    select(plot,n_tree10) |> 
    unique() |> 
    mutate(id_plot=paste0("Tene_",plot)) |> 
    select(-plot)
  
  dataenv_Tene=read.csv2(file.path(dir.data,"tene_data_env.csv")) |> 
    ## Rename for consistency
    rename(long=Longitude,
           lat=Latitude,
           alt=altitude,
           # area_plot.ha="area_plot-ha",
           id_plot=id) |> 
    ## Compute new variable id_plot database and area
    mutate(database=as.factor("Tene"),
           system=as.factor("forest"),
           area_plot.ha=0.2,
           lat=as.numeric(lat),
           long=as.numeric(long)) |> 
    ## Tidy
    relocate(database,id_plot,.before=site) |> 
    left_join(plotn_tene)
  
  return(list(tree=datatree_Tene,
              plot=dataenv_Tene))
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - Amani data ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param dir.data directory of amani files

bhkamani<-function(dir.data="data/bhkamani"){
  
  ## Load all subdatasets
  dataenv_badenou=read.csv2(file.path(dir.data,"badenou_data_env.csv"))
  dataenv_foumbou=read.csv2(file.path(dir.data,"foumbou_data_env.csv"))
  dataenv_HtSassandra=read.csv2(file.path(dir.data,"HtSassandra_data_env.csv"))
  dataenv_Irobo=read.csv2(file.path(dir.data,"irobo_data_env.csv"))
  dataenv_Niegre=read.csv2(file.path(dir.data,"niegre_data_env.csv"))
  dataenv_AmTene=read.csv2(file.path(dir.data,"tene_data_env.csv"))
  dataenv_Yaya=read.csv2(file.path(dir.data,"yaya_data_env.csv"))
  
  ## gather all datasets
  dataenv_Amani=bind_rows(dataenv_badenou,dataenv_foumbou,dataenv_HtSassandra,dataenv_Irobo,dataenv_Niegre,dataenv_AmTene,dataenv_Yaya) |> 
    ## Rename data for consistency with other database
    rename(lat="Latitude",
           long="Longitude",
           id_plot=id,
           alt="altitude") |> 
    ## Compute new variables: database, system, area
    mutate(database=as.factor("Amani"),
           id_plot=paste0("Am_",id_plot),
           system=as.factor(if_else(categorie=="foretanc",
                                    "forest",
                                    "secondary_forest")),
           lat=as.numeric(lat),
           long=as.numeric(long),
           id_plot=as.factor(id_plot)) |> 
    relocate(database,.before=id_plot)
  
  # Clean env
  # rm(dataenv_badenou,dataenv_foumbou,dataenv_HtSassandra,dataenv_Irobo,dataenv_Niegre,dataenv_AmTene,dataenv_Yaya) 
  
  
  
  # Load tree data
  AmaniTree_foumbou <- read.csv2(file.path(dir.data,"foumbou_data_tree.csv")) |> 
    mutate(dbh=circ/pi)
  AmaniTree_badenou <- read.csv2(file.path(dir.data,"badenou_data_tree.csv"))  |>  
    mutate(dbh=circ/pi)
  AmaniTree_oth <- read.csv2(file.path(dir.data,"total.tree.csv")) |> 
    filter(!site%in%c("Foumbou","Badenou")) # for some reason two sites were measured indeptly
  
  AmaniTree<-bind_rows(AmaniTree_badenou,
                       AmaniTree_foumbou,
                       AmaniTree_oth) |> 
    # Rename data for consistency with other database
    rename(id_tree=id) |> 
    # recompute plot number to match data_env
    group_by(site) |> 
    mutate(plot_new=plot-min(plot)+1) |> # because weird plot assignment
    ungroup() |> 
    # Compute new variables: id_plot, database, origin, system
    mutate(id_plot=paste0("Am_",site,"_",plot_new),
           database=as.factor("Amani"),
           origin=case_when(rem==0~"recruited",
                            rem==1~"remnant")) |> 
    # Compute system type
    left_join(dataenv_Amani |> 
                select(id_plot,system),
              by="id_plot") |> 
    # use directly corrected values because i have the flemme to correct les accents
    mutate(species=speciesCorr,
           genus=genusCorr) |> 
    # correct remnant in forest system
    mutate(## put secondary_forest where data is missing
      system=as.character(system),
      system=as.factor(if_else(is.na(system)==TRUE,
                               "secondary_forest",
                               system)),
      ## set origin
      origin=if_else(system=="forest",
                     "remnant",
                     origin),
      origin=as.factor(origin)) |> 
    # dbh cat and horizontal competition
    mutate(dbh_cat=cut(dbh,
                       breaks=seq(0,200,by=10))) |> 
    group_by(plot,dbh_cat) |> 
    mutate(h_compet=n()/0.2) |> 
    ungroup() |> 
    ## Tidy 
    relocate(database,id_plot,.before="id_tree") 
  
  # competition
  AmaniDensity <- AmaniTree |> 
    filter(dbh>10) |> 
    group_by(id_plot) |> 
    mutate(v_compet=n()/0.2,
           ba_tot=sum(pi*((dbh/200)^2),na.rm=TRUE)/0.2,
           n_tree10=n()) |> 
    ungroup() |> 
    select(id_plot,v_compet,ba_tot,n_tree10) |> 
    unique()
  
  # Final tree
  datatree_Amani=AmaniTree |> 
    left_join(AmaniDensity,by=c("id_plot")) |> 
    filter(dbh<200) |> 
    filter(H<40) |> 
    filter(!is.na(species))
  
  # Final env
  dataenv_Amani<-dataenv_Amani |> 
    left_join(AmaniDensity)
  
  return(list(tree=datatree_Amani,
              plot=dataenv_Amani))
  
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 4 - Nguessan data ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param dir.data directory of amani files

nguessan<-function(dir.data="data/nguessan"){
  AnnyEnv=read.csv2(file.path(dir.data,"data_env.csv")) |> 
    left_join(read_excel(file.path(dir.data,"coord.xlsx")),by=c("plot")) |> 
    select(-alt,-X) |> 
    # Rename data for consistency with other database
    rename(prox=Proximiteforesti.re.,
           density=Densiteforesti.re.,
           remnant=nbr_remanent,
           alt=Altitude,
           Y=lat,
           X=long,
           year_crop=Annedeculture,
           type_crop=Precedentcultural,
           type_soil=Typedesol,
           topo=Topographie) |> 
    ## Compute new variables: database,id_plot,system, area
    mutate(id_plot=as.factor(paste0("Nguessan_",plot)),
           database=as.factor("Nguessan"),
           system=as.factor(if_else(categorie%in%c("foretanci ","foretexpl "),
                                    "forest",
                                    "secondary_forest")),
           area_plot.ha=0.2) |> 
    ## Tidy
    relocate(X,Y,alt,.before="age") |> 
    relocate(database,id_plot,.before="plot") 
  
  
  ## coordinates conversion into latlong
  ### Compute for zone 30
  crs="+proj=utm +zone=30"
  Anny_UTM<-data.frame(id_plot=AnnyEnv$id_plot,X=AnnyEnv$X,Y=AnnyEnv$Y) |> 
    vect(geom = c("X", "Y"), crs = crs) |> 
    project("+proj=longlat") |> 
    as.data.frame(geom="XY") |> 
    rename(long=x,
           lat=y)
  
  ### Gathering dataenv with newly calculated coordinates
  dataenv_Anny=AnnyEnv |> 
    select(-X,-Y) |> 
    left_join(Anny_UTM,by=c("id_plot")) 
  
  # original dataset
  Annytree=read.csv2(file.path(dir.data,"data_anny_anne.csv")) |> 
    # Rename data for consistency with other database
    rename(dbh=D) |> 
    # Compute new variables: database,id_plot, system, origin
    mutate(database=as.factor("Nguessan"),
           id_plot=paste0(database,"_",plot),
           origin=case_when(rem==0~"recruited",
                            rem==1~"remnant")) |> 
    # Compute id_tree
    group_by(plot) |> 
    mutate(id_tree=paste0(id_plot,"_",row_number())) |> 
    ungroup() |> 
    # Compute system type
    left_join(dataenv_Anny |> 
                select(id_plot,system),
              by="id_plot") |> 
    # correct remnant in forest system
    mutate(origin=if_else(system=="forest",
                          "remnant",
                          origin),
           origin=as.factor(origin)) |>
    # dbh cat and horizontal competition
    mutate(dbh_cat=cut(dbh,
                       breaks=seq(0,200,by=10))) |> 
    group_by(plot,dbh_cat) |> 
    mutate(h_compet=n()/0.2) |> 
    ungroup() |> 
    ## tidy
    select(-X) |> 
    relocate(database,id_plot,id_tree,.before="plot")
  
  AnnyDensity <- Annytree |> 
    filter(dbh>10) |> 
    group_by(id_plot) |> 
    mutate(v_compet=n()/0.2,
           ba_tot=sum(pi*((dbh/200)^2),na.rm=TRUE)/0.2,
           n_tree10=n()) |> 
    ungroup() |> 
    select(id_plot,v_compet,ba_tot,n_tree10) |> 
    unique()
  
  
  # Final tree
  datatree_Anny=Annytree |> 
    left_join(AnnyDensity,by=c("id_plot")) |> 
    filter(dbh<200)
  
  
  # Final env
  dataenv_Anny<-dataenv_Anny |> 
    left_join(AnnyDensity)
  
  return(list(tree=datatree_Anny,
              plot=dataenv_Anny))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 5 - Lataha data ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param dir.data directory of lataha files

lataha<-function(dir.data="data/lataha"){
  #Total dataset
  dataLataha_tot=read.csv2(file.path(dir.data,"ForestInnov_Lataha.csv")) |> 
    # Useless columns
    select(-Plot,-X,x,y,-Color,-Number,-POM_cm,-Splot_name,-idtree) |> 
    # Rename data for consistency with other database
    rename(dbh=DBH,
           genus=Genus,
           species=Species,
           X=X_UTM,
           Y=Y_UTM,
           name_sc=Name,
           n_ech=ID,
           planting_density=Planting_Density,
           planting_year=Planting_Year,
           dir_deg=Dir_deg,
           H=Height_Total_cm) |> 
    # Compute new variables: database,id_plot, & correct H and n_ech
    mutate(database=as.factor("ForestInnov_Lataha"),
           origin=as.factor("planted"),
           system=as.factor("plantation"),
           id_plot=as.factor(paste0(database,"_",plot)),
           H=H/100,
           n_ech=as.character(n_ech),
    ) |> 
    # Compute new variable id_tree
    group_by(id_plot) |> 
    mutate(id_tree=paste0(id_plot,"_",row_number())) |> 
    ungroup() |> 
    ## tidy
    relocate(database,id_plot,id_tree,.before="plot")
  
  # project coordinates
  crs="+proj=utm +zone=30"
  Lataha_UTM<-data.frame(id_plot=dataLataha_tot$id_plot,X=dataLataha_tot$X,Y=dataLataha_tot$Y) |> 
    unique() |> 
    na.omit() |> 
    vect(geom = c("X", "Y"), crs = crs) |> 
    project("+proj=longlat") |> 
    as.data.frame(geom="XY") |> 
    rename(long=x,
           lat=y)
  
  #environmental data
  dataenv_Lataha=dataLataha_tot |> 
    group_by(id_plot) |> 
    mutate(area_plot.ha=(max(x)*max(y))/10000)|> 
    select(database,id_plot,plot,X,Y,planting_density,planting_year,area_plot.ha,system)|> 
    left_join(Lataha_UTM,by="id_plot") |> 
    mutate(plot=as.factor(plot),
           # fill missing coordinates with mean coordinates of plantation
           X=case_when(is.na(X)~220244,TRUE~X),
           Y=case_when(is.na(Y)~1058451,TRUE~Y),
           long=case_when(is.na(long)~-5.548512,TRUE~long),
           lat=case_when(is.na(lat)~9.565892,TRUE~lat)) |> 
    unique() |> 
    ungroup()
  
  # competition
  LatahaDensity <- dataLataha_tot |> 
    group_by(id_plot) |> 
    mutate(area_plot.ha=(max(x)*max(y))/10000)|> 
    filter(dbh>10) |> 
    filter(alive==1) |> 
    mutate(v_compet=n()/area_plot.ha,
           ba_tot=sum(pi*((dbh/200)^2),na.rm=TRUE)/area_plot.ha,
           n_tree10=n()) |> 
    ungroup() |> 
    select(id_plot,v_compet,ba_tot,n_tree10) |> 
    unique()
  
  LatahaHcompet <- dataLataha_tot |> 
    group_by(id_plot) |> 
    mutate(area_plot.ha=(max(x)*max(y))/10000)|> 
    filter(alive==1) |>
    ungroup() |> 
    ## dbh cat and horizontal competition
    mutate(dbh_cat=cut(dbh,
                       breaks=seq(0,200,by=10))) |> 
    group_by(plot,dbh_cat) |> 
    mutate(h_compet=n()/area_plot.ha) |> 
    ungroup() |> 
    select(id_tree,dbh_cat,h_compet)
  
  
  datatree_Lataha=dataLataha_tot |> 
    left_join(LatahaDensity,by=c("id_plot")) |> 
    left_join(LatahaHcompet,by=c("id_tree")) |> 
    select(database,id_tree,id_plot,plot,n_ech,dbh,genus,species,name_sc,alive,growth,H,origin,system,dbh_cat,h_compet,v_compet,ba_tot,n_tree10) |>
    ## Filter useless individuals : dead trees and trees with no height measures
    filter(alive==1) |> 
    filter(!is.na(H)) 
  
  return(list(tree=datatree_Lataha,
              plot=dataenv_Lataha))
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 5 - Mopri-Sangoue data ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param dir.data directory of lataha files

mopri.sangoue<-function(dir.data="data/mopri_sangoue"){
  # Species correspondance
  species_MopriSangoue=data.frame(
    vernacular=c("aboudikro","acajou","adjoaba","assamela","badi","bete","cordia",
                 "frake","framire","fromager","iroko","irvingia","koto","makore",
                 "movinguie","pouo","samba","sipo","sougue","spondias","tiama"),
    genus=c("Entandrophragma","Khaya","Dacryodes","Pericopsis","Nauclea","Mansonia",
            "Cordia","Terminalia","Terminalia","Ceiba","Milicia","Irvingia","Pterygota",
            "Tieghemella","Distemonanthus","Funtumia","Triplochiton","Entandrophragma",
            "Parinari","Spondias","Entandrophragma"),
    species=c("cylindricum","anthoteca","klaineana","elata","diderrichii","altissima",
              "platythyrsa","superba","ivorensis","pentandra","excelsa","gabonensis",
              "macrocarpa","heckelii","benthamianus","africana","scleroxylon","utile",
              "excelsa","mombin","angolense")
  )
  
  ## All dataset for competition index computation
  ################################################
  dataMopriSangoue_tot=read.csv2(file.path(dir.data,"dataAnne_Compet.csv"))  |> 
    # Useless columns
    select(-X,color,-pom) |>    
    # Rename data for consistency with other database
    rename(plot=plot_name,        
           location=splot) |> 
    # Compute new variables: database, system, origin, id_plot, id_tree, dbh in cm
    mutate(database=as.factor(paste0("ForestInnov_",location)), 
           system=as.factor("plantation"),
           origin=as.factor("planted"),
           id_plot=as.factor(paste0(database,"_",tolower(plot))),
           n_ech=as.character(number),
           id_tree=paste0(id_plot,"_",tolower(vernacular),tolower(color),n_ech),
           dbh=dbh/10,
           long=case_when(location=="mopri"~-4.883175,
                          location=="sangoue"~-5.515998),
           lat=case_when(location=="mopri"~5.839006,
                         location=="sangoue"~6.307813),
           X=case_when(location=="mopri"~291492.748,
                       location=="sangoue"~221628.194),
           Y=case_when(location=="mopri"~645756.684,
                       location=="sangoue"~697902.930)) |> 
    # compute plot area
    group_by(id_plot) |> 
    mutate(area_plot.ha=(max(x)*max(y))/10000) |> 
    ungroup() |> 
    # dead trees
    filter(is.na(dbh)==FALSE) |> 
    filter(alive!=0) |> 
    # tidy
    relocate(database,id_plot,id_tree,.before="plot") 
  
  # competition
  MopriSangoue_hcompet=dataMopriSangoue_tot |> 
    # h_compet
    mutate(dbh_cat=cut(dbh,
                       breaks=seq(0,200,by=10))) |> 
    group_by(plot,dbh_cat) |> 
    mutate(h_compet=n()/area_plot.ha) |> 
    ungroup() |> 
    select(id_tree,dbh_cat,h_compet)
  
  MopriSangoue_vcompet=dataMopriSangoue_tot |> 
    group_by(id_plot) |> 
    #v_compet
    filter(dbh>10) |> 
    summarise(v_compet=n()/area_plot.ha,
              ba_tot=sum(pi*((dbh/200)^2),na.rm=TRUE)/area_plot.ha,
              n_tree10=n()) |> 
    ungroup() |>
    unique()
  
  
  ## Dataset with only H/DBH measures
  ###################################
  datatree_MopriSangoue=read.csv2(file.path(dir.data,"dataAnne_Height.csv")) |> 
    select(-X,-species_number,-pom,-dbh_pred,-coef,-vol,-growth_alpha,-growth_beta,-dbh_pred,-Height_Bole_cm,-growth) |>
    rename(plot=plot_name,        
           location=splot,
           H=Height_Total_cm) |>
    # Compute new variables: database, system, origin, id_plot, id_tree, dbh in cm
    mutate(database=as.factor(paste0("ForestInnov_",location)), 
           system=as.factor("plantation"),
           origin=as.factor("planted"),
           id_plot=as.factor(paste0(database,"_",tolower(plot))),
           n_ech=as.character(number),
           id_tree=paste0(id_plot,"_",tolower(vernacular),tolower(color),n_ech),
           dbh=dbh/10,
           H=H/100,
           long=case_when(location=="mopri"~-4.883175,
                          location=="sangoue"~-5.515998),
           lat=case_when(location=="mopri"~5.839006,
                         location=="sangoue"~6.307813),
           X=case_when(location=="mopri"~291492.748,
                       location=="sangoue"~221628.194),
           Y=case_when(location=="mopri"~645756.684,
                       location=="sangoue"~697902.930)) |> 
    relocate(database,id_plot,id_tree,.before="plot") |> 
    left_join(species_MopriSangoue,by=c("vernacular")) |> 
    left_join(MopriSangoue_hcompet,by=c("id_tree")) |> 
    left_join(MopriSangoue_vcompet,by=c("id_plot")) |> 
    filter(is.na(dbh_cat)==FALSE) |> 
    filter(H<60) |> 
    filter(dbh<200)
  
  dataenv_MopriSangoue=dataMopriSangoue_tot |>
    left_join(MopriSangoue_vcompet) |> 
    select(database,id_plot,plot,location,inventory,system,long,lat,X,Y,area_plot.ha,n_tree10) |> 
    unique()
  
  return(list(tree=datatree_MopriSangoue,
              plot=dataenv_MopriSangoue))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 6 - Elsa Sanial data ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param dir.data directory of ESanial files

esanial<-function(dir.data="data/esanial"){
  ### Environmental data ###
  ##########################
  # original dataset
  ElsaEnv=read.csv2(file.path(dir.data,"plot.csv")) |> 
    rename(plot=code,
           lat=latitude,
           long=longitude,
           area_plot.ha=sup_totale,
           age=age_champ,
           density=densit) |> 
    mutate(database=as.factor("Sanial"),
           id_plot=as.factor(paste0(database,"_",plot)),
           system=as.factor("agroforestry"))
  
  
  
  ## coordinates conversion into UTM
  ### Compute for zone 30
  crs="+proj=longlat +datum=WGS84"#"+proj=utm +zone=30"
  Elsa_longlat30<-data.frame(id_plot=ElsaEnv$id_plot,long=ElsaEnv$long,lat=ElsaEnv$lat) |> 
    filter(long>-6) |> 
    vect(geom = c("long", "lat"), crs = crs) |> 
    project("+proj=utm +zone=30 ellps=WGS84") |> 
    as.data.frame(geom="XY") |> 
    rename_with(.cols=c("x","y"),
                toupper)
  
  
  ### Compute for zone 29
  crs="+proj=longlat +datum=WGS84"#"+proj=utm +zone=30"
  Elsa_longlat29<-data.frame(id_plot=ElsaEnv$id_plot,long=ElsaEnv$long,lat=ElsaEnv$lat) |> 
    filter(long<(-6)) |> 
    vect(geom = c("long", "lat"), crs = crs) |> 
    project("+proj=utm +zone=29 ellps=WGS84") |> 
    as.data.frame(geom="XY") |> 
    rename_with(.cols=c("x","y"),
                toupper)
  
  ### Gathering the 2 utm zones
  Elsa_longlat=bind_rows(Elsa_longlat29,Elsa_longlat30)
  
  # Transitory
  dataenv_Elsa=ElsaEnv |> 
    # select(long,lat) |> 
    ##Gathering dataenv with newly calculated coordinates
    left_join(Elsa_longlat,by=c("id_plot")) 
  
  ### Tree data ###
  #################
  
  # Original dataset
  ElsaTree=read.csv2(file.path(dir.data,"treedata.csv")) |> 
    ## Rename data for consistency with other database
    rename(n_ech=X,
           name=nom_arb,
           H=taille,
           strata=strate,
           circ=peri,
           plot=code,
           age=age_arb,
           family=Famille,
           genus=Genre,
           species=Espece,
           genusCorr=GenreCorr,
           speciesCorr=EspeceCorr) |> 
    # Compute new variables database, origine, id_plot, system and correct n_ech
    mutate(database=as.factor("Sanial"),
           origin=as.factor(case_when(meth_intro==1~"planted",
                                      meth_intro==2~"recruited",
                                      meth_intro==3~"remnant")),
           id_plot=as.factor(paste0(database,"_",plot)),
           system=as.factor("agroforestry"),
           n_ech=as.character(n_ech)) |> 
    # Compute id_tree
    group_by(id_plot) |> 
    mutate(id_tree=paste0(id_plot,"_",row_number())) |> 
    ungroup() |> 
    
    # Tidy
    relocate(database,id_plot,id_tree,.before="n_ech")
  
  ElsaDensity=ElsaTree |> 
    left_join(dataenv_Elsa[,c("id_plot","area_plot.ha")],by=c("id_plot")) |> 
    filter(dbh>10) |> 
    group_by(id_plot) |> 
    mutate(v_compet=n()/area_plot.ha,
           ba_tot=sum(pi*((dbh/200)^2),na.rm=TRUE)/area_plot.ha,
           n_tree10=n()) |> 
    ungroup() |> 
    select(id_plot,v_compet,ba_tot,n_tree10) |> 
    unique()
  
  ElsaHcompet=ElsaTree |> 
    left_join(dataenv_Elsa |> select(id_plot,area_plot.ha),by=c("id_plot")) |> 
    mutate(dbh_cat=cut(dbh,
                       breaks=seq(0,200,by=10))) |> 
    group_by(id_plot,dbh_cat) |> 
    mutate(h_compet=n()/area_plot.ha) |> 
    ungroup() |> 
    select(id_tree,dbh_cat,h_compet)
  
  
  # Transitory
  datatree_Elsa=ElsaTree |> 
    left_join(ElsaDensity,by=c("id_plot")) |> 
    left_join(ElsaHcompet,by=c("id_tree")) |> 
    select(database,id_plot,id_tree,n_ech,name,H,strata,circ,age,family,genus,
           species,dbh,dbh_cat,genusCorr,speciesCorr,familyAPG,origin,system,
           v_compet,h_compet,ba_tot,n_tree10)
  
  
  dataenv_Elsa=dataenv_Elsa |> 
    left_join(ElsaDensity)
  
  return(list(tree=datatree_Elsa,
              plot=dataenv_Elsa))
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 7 - Elsa Sanial data ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param dir.data directory of ESanial files
akouassi<-function(dir.data="data/akouassi"){
  load(file.path(dir.data,"tree_hd.Rdata"))
  tree_hd |> 
    rename(H=htot,
           n_ech=stemID,
           site=clus) |>
    mutate(id_plot=paste0("ak_",site,"_",type),
           id_tree=paste0("ak_",n_ech),
           database="AimeK",
           system="agroforestry",
           origin=case_when(origin=="spontaneous"~"recruited",
                            TRUE~origin))
  
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 8 - Merge data ####
#' @author Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get tree data
#' @description
#' Merge all tree data of different inventories and check taxonomy
#' @param dataset datasets to merge
#' @param output.file output to be saved

get.tree <- function(data_bhkamani,
                     data_nguessan,
                     data_esanial,
                     data_lataha,
                     data_tene,
                     data_mopri.sangoue,
                     data_oibt,
                     output.file){
  #DATATREE
  data_tree<-bind_rows(data_bhkamani$tree,
                       data_nguessan$tree,
                       data_esanial$tree,
                       data_oibt$tree) |> #datatree_OIBT,
    mutate(plot=as.character(plot),
           subplot=as.character(subplot)) |> 
    bind_rows(data_lataha$tree,
              data_tene$tree,
              data_mopri.sangoue$tree) |> 
    mutate(origin=as.character(origin),
           origin_2=as.factor(if_else(system=="forest",
                                      "forest",
                                      origin)),
           origin=as.factor(origin)) |> 
    select(database,id_plot,id_tree,system,origin,genus,species,dbh,H,dbh_cat,
           v_compet,h_compet,ba_tot,n_tree10) # selecting variables of interest
  
  ## Check taxo
  Taxo=data_tree |> 
    select(genus,species) |> 
    unique() |> 
    arrange(genus,species) 
  
  Taxonomy <- correctTaxo(genus = Taxo$genus, species = Taxo$species)
  Taxo$genusCorr <- Taxonomy$genusCorrected
  Taxo$speciesCorr <- Taxonomy$speciesCorrected
  Taxo$modif <- Taxonomy$nameModified
  
  data_tree=data_tree |> 
    left_join(Taxo,by=c("genus","species")) 
  
  
  write.csv2(data_tree,file.path("data",output.file))
  return(data_tree)
}

#' Get environnemental data
#' @description
#' Merge all tree data of different inventories and check taxonomy
#' @param dataset datasets to merge
#' @param output.file output to be saved

get.env <- function(data_bhkamani,
                    data_nguessan,
                    data_esanial,
                    data_lataha,
                    data_tene,
                    data_mopri.sangoue,
                    data_oibt,
                    output.file){
  data_env=bind_rows(data_bhkamani$plot,
                     data_nguessan$plot,
                     data_tene$plot) |> 
    mutate(plot=as.factor(plot)#,
           # prec=as.character(prec)
           ) |> 
    bind_rows(data_esanial$plot,
              data_lataha$plot,
              data_mopri.sangoue$plot,
              data_oibt$plot |> select(-prec)) |> 
    select(database,id_plot,X,Y,lat,long,area_plot.ha,system)
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

  
  write.csv2(data_env,file.path("data",output.file))
  
  return(data_env)
}
