#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   Subfunctions   ---------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Grouped Data ----
# 1. Grouped Data ----
gd_pop_wlf <- function(pl) {
  dt <- pipload::pip_load_cache(
    pl$country_code,
    pl$surveyid_year,
    verbose = FALSE,
    version = pl$version
  )
  
  levels <- dt[, unique(reporting_level)]
  
  welfare <- vector("list", length(levels))
  population <- vector("list", length(levels))
  
  for (i in seq_along(levels)) {
    nn <- levels[[i]]
    welfare[[i]] <- dt[reporting_level == nn, mean]
    population[[i]] <- dt[reporting_level == nn, weight]
  }
  
  names(welfare) <- names(population) <- levels
  
  return(list(
    welfare = welfare,
    population = population
  ))
}


get_gd_calcs <- function(level, vctr, mean, id) {
  welfare <- vctr$welfare[[level]]
  population <- vctr$population[[level]]
  mean <- mean[[level]]
  
  params <- wbpip:::get_gd_quantiles(
    welfare,
    population,
    complete = TRUE,
    mean = mean,
    popshare = popshare,
    lorenz = lorenz
  )
  
  povlines <- params$dist_stats$quantiles
  lorenz <- params$selected_lorenz$for_dist
  
  wlf_share <-
    get_gd_wlf_share_by_qtl(params = params, lorenz = lorenz, n = nq) |>
    {
      \(.) .$dist_stats$welfare_share
    }()
  
  pop_share <- c(popshare[1], diff(popshare))
  
  avg_wlf_qtl <- (wlf_share * mean) / pop_share
  
  dt <- data.table(
    quantile = povlines,
    welfare_share = wlf_share,
    pop_share = pop_share,
    avg_welfare = avg_wlf_qtl,
    reporting_level = level,
    id = id
  )
  
  dt[, bin := .I][bin == max(bin), quantile := NA_real_]
  
  return(dt)
}

poss_get_gd_calcs <- purrr::possibly(.f = get_gd_calcs, otherwise = NULL)


fmt_sve <- function(dt, id) {
  dt <- copy(dt)
  nvars <- c("country_code", "year", "welfare_type")
  id <- unique(dt[, id])
  
  dt[,
     (nvars) := tstrsplit(id, split = "_")
  ][, id := NULL]
  
  id <- glue("{id}_{nq}bin")
  
  # save
  # haven::write_dta(dt, fs::path("data/singles", id, ext = "dta"))
  fs::path(singles_dir, id, ext = ext) |>
    qs::qsave(x = dt, file = _)
}


# 2. Micro Data ----
get_micro_dist <- function(pl) {
  
  df <- pipload::pip_load_cache(
    pl$country_code,
    pl$surveyid_year,
    verbose = FALSE,
    version = pl$version,
    welfare_type = pl$wt_call
  )
  
  
  setorder(
    df,
    imputation_id,
    reporting_level,
    welfare_type,
    mean,
    weight
  )
  
  
  popct <- pop[
    country_code == pl$country_code &
      reporting_year == pl$reporting_year
  ]
  
  
  n_ids <- df[, uniqueN(imputation_id)]
  
  wt <- funique(df$welfare_type)
  
  
  # Renombrar mean -> welfare para usar lorenz_table()
  df <- df |>
    fselect(
      reporting_level,
      welfare = mean,
      weight
    )
  
  
  df[, weight := weight / n_ids]
  
  
  lt <- lorenz_table(df, nq = nq)
  
  
  lt <- lt[
    popct,
    on = "reporting_level",
    pop := i.pop
  ]
  
  
  lt <- lt[bin >= nq, quantile := NA_real_]
  
  
  lt[
    ,
    id := paste(
      pl$country_code,
      pl$reporting_year,
      wt,
      sep = "_"
    )
  ]
  
  
  return(lt)
}


# 2. Rur/Urb to National ----
# Functions for surveys with urban/rural reporting level only
create_synth_bins <- function(vctr, mean, pop, nbins, id) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load synthetic distributions ----
  dt <- get_synth_level_rbind(
    vctr = vctr,
    mean = mean,
    pop = pop
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Bins and totals ----
  dt[,
     bin := wbpip:::md_compute_bins(
       welfare,
       weight,
       nbins = nbins,
       output = "simple"
     )
  ]
  dt[, `:=`(
    tot_pop = sum(weight, na.rm = TRUE),
    tot_wlf = sum(welfare * weight, na.rm = TRUE)
  )]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Main measures ----
  dt[,
     c("avg_welfare", "pop_share", "welfare_share", "quantile") := {
       avg_welfare <- weighted.mean(welfare, weight, na.rm = TRUE)
       pop_share <- sum(weight, na.rm = TRUE) / tot_pop
       welfare_share <- sum(welfare * weight, na.rm = TRUE) / tot_wlf
       quantile <- max(welfare, na.rm = TRUE)
       .(avg_welfare, pop_share, welfare_share, quantile)
     },
     by = .(bin)
  ]
  
  # Getting means ----
  dt <- dt[,
           lapply(.SD, mean, na.rm = TRUE),
           by = .(bin),
           .SDcols = c("avg_welfare", "pop_share", "welfare_share", "quantile")
  ]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Add reporting level variable ----
  dt[,
     reporting_level := 'national'
  ]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Censoring --------
  
  dt <- dt[bin == 100, quantile := NA_real_]
  
  return(dt)
}

poss_create_synth_bins <- purrr::possibly(
  .f = create_synth_bins,
  otherwise = NULL
)

get_synth_level_rbind <- function(vctr, mean, pop) {
  # Extract levels
  levels <- names(mean)
  
  # Create synth bins for each level
  synth_level <- purrr::map(
    levels,
    ~ {
      level <- .x
      # Extract data for the current level
      welfare <- vctr$welfare[[level]]
      population <- vctr$population[[level]]
      mean_level <- mean[[level]]
      pop_level <- pop[[level]]
      
      # Create synthetic vector for this level
      synth_vector_level <- wbpip:::sd_create_synth_vector(
        welfare = welfare,
        population = population,
        mean = mean_level,
        pop = pop_level
      )
      return(synth_vector_level)
    }
  )
  
  # Combine all the results into a single data.table
  combined_results <- purrr::reduce(synth_level, rbind) |>
    setorder(welfare)
  
  return(combined_results)
}


# Variation on gd_pop_wlf which gets only urban/rural: maybe we can adapt the
# original later to take an additional argument.
gd_pop_wlf_rur_urb <- function(pl) {
  
  dt <- pipload::pip_load_cache(
    pl$country_code,
    pl$surveyid_year,
    verbose = FALSE,
    version = pl$version
  )
  
  
  if ("national" %in% dt$reporting_level) {
    return(NULL)
  }
  
  
  levels <- dt[, unique(reporting_level)]
  
  welfare <- vector("list", length(levels))
  population <- vector("list", length(levels))
  
  
  for (i in seq_along(levels)) {
    
    nn <- levels[[i]]
    
    welfare[[i]] <-
      dt[reporting_level == nn, mean]
    
    population[[i]] <-
      dt[reporting_level == nn, weight]
  }
  
  
  names(welfare) <- names(population) <- levels
  
  
  return(
    list(
      welfare = welfare,
      population = population
    )
  )
}

get_micro_dist_rur_urb <- function(pl) {
  
  
  dt <- pipload::pip_load_cache(
    pl$country_code,
    pl$surveyid_year,
    verbose = FALSE,
    version = pl$version,
    welfare_type = pl$wt_call
  )
  
  
  setorder(
    dt,
    imputation_id,
    welfare_type,
    mean,
    weight
  )
  
  
  unique_levels <- unique(dt$reporting_level)
  
  if ("national" %in% unique_levels) {
    return(NULL)
  }
  
  
  
  dt[,
     bin := wbpip:::md_compute_bins(
       mean,
       weight,
       nbins = nq,
       output = "simple"
     ),
     by = c(
       "imputation_id",
       "welfare_type"
     )
  ][,
    `:=`(
      tot_pop = sum(weight),
      tot_wlf = sum(mean * weight)
    ),
    by = c(
      "imputation_id",
      "welfare_type"
    )
  ]
  
  
  
  dt[,
     c(
       "avg_welfare",
       "pop_share",
       "welfare_share",
       "quantile"
     ) :=
       {
         avg_welfare <-
           weighted.mean(mean, weight)
         
         pop_share <-
           sum(weight) / tot_pop
         
         welfare_share <-
           sum(mean * weight) / tot_wlf
         
         quantile <-
           max(mean)
         
         list(
           avg_welfare,
           pop_share,
           welfare_share,
           quantile
         )
       },
     by = .(
       imputation_id,
       welfare_type,
       bin
     )
  ]
  
  
  
  dt <-
    dt[,
       lapply(
         .SD,
         mean
       ),
       by = .(
         imputation_id,
         welfare_type,
         bin
       ),
       .SDcols =
         c(
           "avg_welfare",
           "pop_share",
           "welfare_share",
           "quantile"
         )
    ][,
      lapply(
        .SD,
        mean
      ),
      by = .(
        welfare_type,
        bin
      ),
      .SDcols =
        c(
          "avg_welfare",
          "pop_share",
          "welfare_share",
          "quantile"
        )
    ]
  
  
  
  dt[
    bin == nq,
    quantile := NA_real_
  ]
  
  
  dt[
    ,
    id := paste(
      pl$country_code,
      pl$reporting_year,
      welfare_type,
      sep="_"
    )
  ]
  
  
  return(dt)
}

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Bins ant totals --------
  
  dt[,
     # get bins and total pop and welfare
     bin := wbpip:::md_compute_bins(
       welfare_ppp,
       weight,
       nbins = nq,
       output = "simple"
     ),
     by = c("imputation_id", "welfare_type") # remove reporting level
  ][,
    `:=`(
      tot_pop = sum(weight),
      tot_wlf = sum(welfare_ppp * weight)
    ),
    by = c("imputation_id", "welfare_type") # remove reporting level
  ]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Main measures --------
  
  dt[,
     # get avg wlf, pop and wlf shared
     c("avg_welfare", "pop_share", "welfare_share", "quantile") := {
       avg_welfare <- weighted.mean(welfare_ppp, weight)
       pop_share <- sum(weight) / tot_pop
       welfare_share <- sum(welfare_ppp * weight) / tot_wlf
       quantile <- max(welfare_ppp)
       list(avg_welfare, pop_share, welfare_share, quantile)
     },
     by = .(imputation_id, welfare_type, bin)
  ]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Getting means --------
  # Note: same as micro dist but no grouping by reporting level.
  
  dt <-
    dt[,
       # mean by imputation and bin
       lapply(.SD, mean),
       by = .(imputation_id, welfare_type, bin),
       .SDcols = c("avg_welfare", "pop_share", "welfare_share", "quantile")
    ][,
      # mean by bin
      lapply(.SD, mean),
      by = .(welfare_type, bin),
      .SDcols = c("avg_welfare", "pop_share", "welfare_share", "quantile")
    ]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Censoring --------
  
  dt <- dt[bin == 100, quantile := NA_real_]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Creating id --------
  
  dt[,
     id := paste(pl$country_code, pl$reporting_year, welfare_type, sep = "_")
  ]
  
  return(dt)



fmt_sve_synth <- function(dt, id) {
  dt <- copy(dt)
  nvars <- c("country_code", "year", "welfare_type")
  id <- unique(dt[, id])
  
  dt[,
     (nvars) := tstrsplit(id, split = "_")
  ][, id := NULL]
  
  id <- glue("{id}_synth_{nq}bin") # added synth to name to avoid overwriting
  
  fs::path(singles_dir, id, ext = ext) |>
    qs::qsave(x = dt, file = _)
}


# Fix welfare share

fix_welfare_share <- \(x, nq) {
  ori_name <- copy(names(x))
  
  nbybins <- x[, .N, by = reporting_level]
  
  sv <- vector("list", nrow(nbybins))
  
  for (i in seq_len(nrow(nbybins))) {
    nbins <- nbybins[i, N]
    sv[[i]] <-
      seq(log(1e-7), log(1e-6), length.out = nbins) |>
      exp()
  }
  sv <- unlist(sv)
  
  y <- x |>
    add_vars(small = sv) |>
    #  1. create equal sized bins
    ftransform(
      tot_pop = fsum(pop, reporting_level, TRA = "fill"),
      bin_total = fnobs(pop, reporting_level, TRA = "fill")
    ) |>
    ftransform(pop_adj = tot_pop / bin_total) |>
    # 2. create total welfare of each bin and
    # aggregate welfare in the population
    ftransform(
      tot_welf_orign = avg_welfare * pop,
      tot_welf_equal = avg_welfare * pop_adj
    ) |>
    ftransform(
      agg_welf_orign_sum = fsum(tot_welf_orign, reporting_level, TRA = "fill"),
      agg_welf_equal_sum = fsum(tot_welf_equal, reporting_level, TRA = "fill")
    ) |>
    #  3. generate the discrepancy (residual) of aggregate welfare
    # in population using original or equalized population  distribute
    # the discrepancy using the original distribution of average welfare
    ftransform(
      residual = agg_welf_orign_sum - agg_welf_equal_sum,
      tot_avg_welf = fsum(avg_welfare, reporting_level, TRA = "fill")
    ) |>
    ftransform(
      tot_welf_adj = tot_welf_equal +
        small +
        (residual * avg_welfare / tot_avg_welf)
    ) |>
    # first adjust to avg welfare. change name of var to avg_welf_adj
    # to compare
    #. 4. re-create a new average welfare
    ftransform(avg_welfare = tot_welf_adj / pop_adj) |>
    # Second adjust to welfare share. change name of var to welf_share_adj
    # to compare
    ftransform(welfare_share = tot_welf_adj / agg_welf_orign_sum) |>
    get_vars(ori_name)
  
  y
}
