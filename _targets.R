library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("dplyr","stringr","readxl","sf","terra","BIOMASS"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(
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
    tree,
    get.tree(data_bhkamani,
             data_nguessan,
             data_esanial,
             data_lataha,
             data_tene,
             data_mopri.sangoue,
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
            output.file="all_plot_corr.csv")
  ),
  NULL
)