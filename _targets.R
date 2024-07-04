library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("dplyr","tidyr","stringr",#"readxl","BIOMASS",
               "sf","terra",
               "rstan","loo","blockCV",
               "future"), # packages that your targets need to run
  format = "rds", # default storage format
  memory="transient"
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")
future::plan(future::multisession, workers = 6)
# tar_make_future() configuration (okay to leave alone):

# tar_source()
source("R/function_data.R")
source("R/function_fit.R")
# source("other_functions.R") # Source other scripts as needed. # nolint

list(
  ## Data formating
  tar_target(
    data_oibt,
    oibt(dir.data="data/oibt")
  ),
  tar_target(
    data_tene,
    tene(dir.data="data/tene")
  ),
  tar_target(
    data_bhkamani,
    bhkamani(dir.data="data/bhkamani")
  ),
  tar_target(
    data_nguessan,
    nguessan(dir.data="data/nguessan")
  ),
  tar_target(
    data_lataha,
    lataha(dir.data="data/lataha")
  ),
  tar_target(
    data_mopri.sangoue,
    mopri.sangoue(dir.data="data/mopri_sangoue")
  ),
  tar_target(
    data_esanial,
    esanial(dir.data="data/esanial/")
  ),
  tar_target(
    data_akouassi,
    akouassi(dir.data="data/akouassi/")
  ),
  tar_target(
    tree,
    get.tree(data_bhkamani,
             data_nguessan,
             data_esanial,
             data_lataha,
             data_tene,
             data_mopri.sangoue,
             data_oibt,
             data_akouassi,
             output.file="all_tree_corr.csv")
  ),
  tar_target(
    plot,
    get.env(data_bhkamani,
            data_nguessan,
            data_esanial,
            data_lataha,
            data_tene,
            data_mopri.sangoue,
            data_oibt,
            data_akouassi,
            output.file="all_plot_corr.csv")
  ),
  
  # Fit models
  tar_target(
    mod.data,
    data.fit(tree,
             plot,
             frac=1)
  ),
  tar_target(
    sub.mod.data.1,
    data.fit(tree,
             plot,
             frac=0.1,
             ba_require=FALSE)
  ),
  tar_target(
    sub.mod.data.3,
    data.fit(tree,
             plot,
             frac=0.3,
             ba_require=FALSE)
  ),
  tar_target(
    sub.mod.data.3.ba,
    data.fit(tree,
             plot,
             frac=0.3,
             ba_require=TRUE)
  ),
  tar_target(
    shape.selection,
    shape.select(folder="mod_shape_select",
                 sub.mod.data.1)
  ),
  tar_target(
    cofactor.selection,
    cofactor.select(folder="mod_cof_select",
                    sub.mod.data.3)
  ),
  tar_target(
    covariable.models,
    covariable.select(sub.mod.data.3.ba,
                      folder="mod_cov_select")
  ),
  tar_target(
    comp.cofactor,
    select.mod(files.list=c("mod_cov_select/mod_nul.rdata",
                            "mod_cov_select/nul_ori.rdata",
                            "mod_cov_select/nul_sys.rdata",
                            "mod_cov_select/nul_systori.rdata"),
               list.names=c("nul","ori","sys","systori"))
  ),
  tar_target(
    comp.covar,
    select.mod(files.list=c("mod_cov_select/so_bio01.rdata",
                            "mod_cov_select/so_bio05.rdata",
                            "mod_cov_select/so_bio12.rdata",
                            "mod_cov_select/so_bio17.rdata",
                            "mod_cov_select/so_v_compet.rdata",
                            "mod_cov_select/so_ba_tot.rdata",
                            "mod_cov_select/complete.rdata"),
               list.names=c("mat","mtwm","map","mpdq","vcomp","ba","complete"))
  ),
  tar_target(
    list.subdataset,
    make_subdataset(mod.data,
                    folder="rds/")
  ),
  tar_target(
    mod.data.block,
    make_blockCV(mod.data,
                 nfold=4)
  ),
  tar_target(
    sim.plan,
    data.frame(model=c("nul","system","origin","systori","complete")) |> 
      crossing(fold=1:4)
  ),
  tar_target(
    id.sim,
    1:dim(sim.plan)[1]
  ),
  tar_target(
    spatial.cross.val,
    make_model(sim.plan,
               id.sim,
               mod.data.block,
               folder="mod_spatial_cross_valid/"),
    pattern=map(id.sim),
    iteration="vector",
    format="file"
  ),
  tar_target(
    predict_spatial_crossval,
    get_prediction(spatial.cross.val,
                   sim.plan,
                   id.sim,
                   mod.data.block),
    pattern=map(id.sim),
    iteration="vector"
  ),
  NULL
)
