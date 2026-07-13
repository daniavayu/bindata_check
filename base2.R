library(fastverse)
library(data.table)
library(haven)
library(fs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("functions.R")
source("duplicatehousehols.R")

op <- options(joyn.reportvar = "report")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Paths
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_dir <- "C:/Users/wb661551/OneDrive - WBG/Desktop/Internship/Bottom Censoring/01-input"

output_dir <- "C:/Users/wb661551/OneDrive - WBG/Desktop/Internship/Bottom Censoring/02-output"


fs::dir_create(output_dir, recurse = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find .dta files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

files <- fs::dir_ls(
  data_dir,
  recurse = TRUE,
  regexp = "\\.dta$"
)

length(files)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nq <- 1000


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Process files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dlt <- lapply(cli::cli_progress_along(files), \(i) {
  
  x <- files[i]
  
  message("Processing: ", x)
  
  
  #----------------------------------------
  # Read data
  #----------------------------------------
  
  dt <- haven::read_dta(x) |>
    data.table::as.data.table()
  
  
  #----------------------------------------
  # Create Lorenz table
  #----------------------------------------
  
  lt <- lorenz_table(
    dt,
    nq = nq
  )
  
  
  #----------------------------------------
  # Censoring
  #----------------------------------------
  
  lt[
    bin >= nq,
    quantile := NA_real_
  ]
  
  
  #----------------------------------------
  # Create ID
  #----------------------------------------
  
  # Extract filename
  fname <- basename(x)
  
  lt[
    ,
    id := fname
  ]
  
  
  lt
  
})


# Combine all files

dlt_final <- rbindlist(
  dlt,
  fill = TRUE
)


# Order columns

cols <- c(
  "id",
  "reporting_level",
  "bin"
)

setorderv(
  dlt_final,
  cols
)

setcolorder(
  dlt_final,
  cols
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save outputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

haven::write_dta(
  dlt_final,
  fs::path(
    output_dir,
    "1kbins.dta"
  )
)


fst::write_fst(
  dlt_final,
  fs::path(
    output_dir,
    "1kbins.fst"
  )
)


dlt_final