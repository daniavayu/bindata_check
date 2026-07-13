# library(ggplot2)
library(data.table)
library(collapse)

# Step 1: Identify Outliers
find_outliers <- \(DT, weight = "weight", threshold = 2.5) {
  mean_w <- fmean(DT[[weight]])
  sd_w   <- fsd(DT[[weight]])
  
  DT[, is_outlier := .SD[[1]] > (mean_w + threshold * sd_w), .SDcols = weight]
  DT
}


optimize_ratio <- \(y, m) {
  # Candidate optimum (smallest x possible given that
  # 1. y/x < x, whichi x > y^(1/2)
  # 2. x > m, hwere m is the min or mean of the distribution
  
  opt_x <- pmax(m, sqrt(y))
  
  # Check the extra condition: we need x_opt <= y/2 to ensure y/x >= 2.
  to_one <- (opt_x > y/2)
  opt_x[to_one] <- 1
  
  round(y / opt_x)
  
}


# Step 2: Calculate Replications and Partitioning for Outliers
duplicate_obs <- \(DT, weight = "weight") {
  min_w <- fmin(DT[[weight]])
  
  # Get houdehold ID
  DT[, hhindex := .I]
  # Get replication count
  DT[, rep_count := optimize_ratio(.SD[[1]],  min_w),
     .SDcols = weight
  ][is_outlier == FALSE,
    rep_count := 1]
  
  Y <- DT[rep(1:.N,rep_count)]
  Y
}


split_weights <- \(x, rep) {
  base <- round(x / rep)
  rem <- x[1] - fsum(base)
  
  add_to_base <- floor(abs(rem)/base[1])
  
  base[length(base)] <- base[1] + rem
  base
}

add_new_weights <- \(DT, weight = "weight") {
  X <- DT[is_outlier == TRUE]
  Y <- DT[is_outlier == FALSE]
  # X[, x := split_weights(.SD[[1]], rep_count),
  #    by = hhindex,
  #    .SDcols = weight]
  X[, x := .SD[[1]]/ rep_count,
    .SDcols = weight]
  setnames(X, c(weight, "x"), c("x", weight)) # reverse names
  X[, x := NULL]
  
  rowbind(Y,X, fill = TRUE)
}

clean_new_weights <- \(DT, ori_names) {
  DT <- DT[, ..ori_names]
}

lorenz_table <- \(df, nq = 100, tolerance = 1e-6) {
  # number of data labels
  no_dl <- funique(df$reporting_level) |>
    length()
  
  
  if (no_dl > 1) {
    ## Bins at national level ----
    df2 <-  copy(df)
    df2[, reporting_level := "national"]
    df <- rowbind(df, df2)
  }
  
  ## Bins at reporting level -------
  # sort according to bins calculation method
  setorder(df, reporting_level, welfare)
  df[, id := .I]
  
  dt <-
    df[, new_bins(welfare = welfare,
                  weight = weight,
                  id = id,
                  nbins = nq,
                  tolerance = tolerance),
       by = reporting_level]  |>
    ftransform(wt_welfare = welfare*weight) |>
    fgroup_by(reporting_level) |>
    ## totals by reporting level
    fmutate(tot_pop = fsum(weight),
            tot_wlf = fsum(wt_welfare)) |>
    fungroup() |>
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Main measures --------
  ## shares at the observation level ---------
  ftransform(pop_share = weight/tot_pop,
             welfare_share = (wt_welfare)/tot_wlf) |>
    ## aggregate
    fgroup_by(reporting_level, bin) |>
    fsummarise(avg_welfare    = fmean(welfare, w = weight),
               pop_share      = fsum(pop_share),
               welfare_share  = fsum(welfare_share),
               quantile       = fmax(welfare),
               pop            = fsum(weight)) |>
    fungroup() |>
    setorder(reporting_level, bin)
  
  dt
  
}

#' Assign households to weighted bins (percentiles)
#'
#' This function assigns each household to a percentile bin using a weighted distribution,
#' ensuring that each bin has approximately the same total weight. Households with large
#' weights may be duplicated and split across bins to preserve monotonicity of welfare shares.
#'
#' @param welfare Numeric vector of household welfare (e.g., income).
#' @param weight  Numeric vector of sampling weights (same length as welfare).
#' @param nbins   Number of bins to divide the population into (default = 100).
#' @param id      Optional ID vector for each observation. If NULL, it will be generated.
#'
#' @return A `data.table` with columns: `id`, `welfare`, `weight`, and `bin`.
#'   The output may have more rows than the input due to household splitting.
#'
#' @export
new_bins1 <- function(welfare, weight, nbins = 100, tolerance = 1e-6, id = NULL) {
  stopifnot(length(welfare) == length(weight))
  
  # Handle NAs
  if (anyNA(welfare) || anyNA(weight)) {
    valid <- !is.na(welfare) & !is.na(weight)
    welfare <- welfare[valid]
    weight  <- weight[valid]
    if (!is.null(id)) id <- id[valid]
  } else {
    id <- if (is.null(id)) seq_along(welfare) else id
  }
  
  # Sort by welfare
  o <- order(welfare)
  welfare <- welfare[o]
  weight  <- weight[o]
  id      <- id[o]
  
  total_weight <- collapse::fsum(weight)
  bin_size     <- total_weight / nbins
  
  # Preallocate result vectors
  out_id      <- integer(2 * length(welfare))
  out_welfare <- numeric(2 * length(welfare))
  out_weight  <- numeric(2 * length(welfare))
  out_bin     <- integer(2 * length(welfare))
  
  cur_bin    <- 1
  cur_weight <- 0
  out_index  <- 1
  
  for (i in seq_along(welfare)) {
    w  <- welfare[i]
    wt <- weight[i]
    id_i <- id[i]
    
    while (wt > 0 && cur_bin <= nbins) {
      room <- bin_size - cur_weight
      
      if (nearly_equal(room, 0, tolerance)) {
        cur_bin <- cur_bin + 1
        cur_weight <- 0
        next
      }
      
      take <- min(wt, room)
      
      out_id[out_index]      <- id_i
      out_welfare[out_index] <- w
      out_weight[out_index]  <- take
      out_bin[out_index]     <- cur_bin
      out_index <- out_index + 1
      
      wt <- wt - take
      cur_weight <- cur_weight + take
      
      if (wt < room || nearly_equal(wt, room, tolerance)) {
        cur_bin <- cur_bin + 1
        cur_weight <- 0
      }
    }
  }
  
  data.table::data.table(
    id   = out_id[1:(out_index - 1)],
    bin  = out_bin[1:(out_index - 1)],
    # You can drop these two if you only want bin assignment:
    weight = out_weight[1:(out_index - 1)],
    welfare = out_welfare[1:(out_index - 1)]
  )
}


new_bins2 <- function(welfare, weight, nbins = 100, tolerance = 1e-6, id = NULL) {
  stopifnot(length(welfare) == length(weight))
  
  # Handle NAs
  if (anyNA(welfare) || anyNA(weight)) {
    valid <- !is.na(welfare) & !is.na(weight)
    welfare <- welfare[valid]
    weight  <- weight[valid]
    if (!is.null(id)) id <- id[valid]
  } else {
    id <- if (is.null(id)) seq_along(welfare) else id
  }
  
  # Sort by welfare
  o <- order(welfare)
  welfare <- welfare[o]
  weight  <- weight[o]
  id      <- id[o]
  
  total_weight <- fsum(weight)
  bin_size     <- total_weight / nbins
  
  # Result containers
  out_id      <- integer(2 * length(welfare))
  out_welfare <- numeric(2 * length(welfare))
  out_weight  <- numeric(2 * length(welfare))
  out_bin     <- integer(2 * length(welfare))
  
  cur_bin     <- 1
  cur_weight  <- 0
  out_index   <- 1
  bin_people  <- integer(nbins)  # Tracks number of people per bin
  
  for (i in seq_along(welfare)) {
    w  <- welfare[i]
    wt <- weight[i]
    id_i <- id[i]
    
    while (wt > 0 && cur_bin <= nbins) {
      room <- bin_size - cur_weight
      
      if (abs(room) < tolerance) {
        cur_bin <- cur_bin + 1
        cur_weight <- 0
        next
      }
      
      take <- fmin(wt, room)
      
      out_id[out_index]      <- id_i
      out_welfare[out_index] <- w
      out_weight[out_index]  <- take
      out_bin[out_index]     <- cur_bin
      
      out_index <- out_index + 1
      wt <- wt - take
      cur_weight <- cur_weight + take
      bin_people[cur_bin] <- bin_people[cur_bin] + 1
      
      # Optional person-count check: Ensure growth over bins
      if (cur_bin > 1 &&
          bin_people[cur_bin] <= bin_people[cur_bin - 1] &&
          cur_weight >= bin_size - tolerance) {
        # If this bin isn't growing in person count, force more splits
        next
      }
      
      if (wt < room || abs(wt - room) < tolerance) {
        cur_bin <- cur_bin + 1
        cur_weight <- 0
      }
    }
  }
  
  data.table::data.table(
    id      = out_id[1:(out_index - 1)],
    bin     = out_bin[1:(out_index - 1)],
    weight  = out_weight[1:(out_index - 1)],
    welfare = out_welfare[1:(out_index - 1)]
  )
}

new_bins3 <- function(welfare, weight, nbins = 100, tolerance = 1e-6, id = NULL) {
  stopifnot(length(welfare) == length(weight))
  
  
  if (anyNA(welfare) || anyNA(weight)) {
    valid <- !is.na(welfare) & !is.na(weight)
    welfare <- welfare[valid]
    weight  <- weight[valid]
    if (!is.null(id)) id <- id[valid]
  } else {
    id <- if (is.null(id)) seq_along(welfare)
  }
  
  # Order the data by welfare
  o <- order(welfare)
  welfare <- welfare[o]
  weight  <- weight[o]
  if (!is.null(id)) id <- id[o]
  
  # Define bin target size
  total_weight <- fsum(weight)
  bin_size     <- total_weight / nbins
  
  # Initialize result containers (oversized)
  out_id      <- integer(2 * length(welfare))
  out_welfare <- numeric(2 * length(welfare))
  out_weight  <- numeric(2 * length(welfare))
  out_bin     <- integer(2 * length(welfare))
  
  cur_bin     <- 1
  cur_weight  <- 0
  out_index   <- 1
  
  for (i in seq_along(welfare)) {
    w  <- welfare[i]
    wt <- weight[i]
    original_id <- if (!is.null(id)) id[i] else i
    
    while (wt > 0 && cur_bin <= nbins) {
      room <- bin_size - cur_weight
      
      # Check if within tolerance of bin end
      if (abs(room) < tolerance) {
        cur_bin <- cur_bin + 1
        cur_weight <- 0
        next
      }
      
      out_id[out_index]      <- original_id
      out_welfare[out_index] <- w
      out_bin[out_index]     <- cur_bin
      out_index <- out_index + 1
      
      if (wt <= room + tolerance) {
        # Fits in current bin
        out_weight[out_index]  <- wt
        cur_weight <- cur_weight + wt
        break
      } else {
        # Split into current bin and remainder
        out_weight[out_index]  <- room
        wt <- wt - room
        cur_bin <- cur_bin + 1
        cur_weight <- 0
      }
    }
  }
  
  # Build result
  res <- data.table::data.table(
    id = out_id[1:(out_index - 1)],
    bin = out_bin[1:(out_index - 1)]
  )
  return(res)
}

# this one works
new_bins4 <- function(welfare, weight, nbins = 100, tolerance = 1e-6, id = NULL) {
  stopifnot(length(welfare) == length(weight))
  
  # Handle NAs
  if (anyNA(welfare) || anyNA(weight)) {
    valid <- !is.na(welfare) & !is.na(weight)
    welfare <- welfare[valid]
    weight  <- weight[valid]
    if (!is.null(id)) id <- id[valid]
  } else {
    id <- if (is.null(id)) seq_along(welfare) else id
  }
  
  # Sort by welfare
  o <- order(welfare)
  welfare <- welfare[o]
  weight  <- weight[o]
  id      <- id[o]
  
  total_weight <- collapse::fsum(weight)
  bin_size     <- total_weight / nbins
  
  # Preallocate result vectors
  out_id      <- integer(2 * length(welfare))
  out_welfare <- numeric(2 * length(welfare))
  out_weight  <- numeric(2 * length(welfare))
  out_bin     <- integer(2 * length(welfare))
  
  cur_bin    <- 1
  cur_weight <- 0
  out_index  <- 1
  
  for (i in seq_along(welfare)) {
    w  <- welfare[i]
    wt <- weight[i]
    id_i <- id[i]
    
    while (wt > 0 && cur_bin <= nbins) {
      room <- bin_size - cur_weight
      
      if (abs(room) < tolerance) {
        cur_bin <- cur_bin + 1
        cur_weight <- 0
        next
      }
      
      take <- min(wt, room)
      
      out_id[out_index]      <- id_i
      out_welfare[out_index] <- w
      out_weight[out_index]  <- take
      out_bin[out_index]     <- cur_bin
      out_index <- out_index + 1
      
      wt <- wt - take
      cur_weight <- cur_weight + take
      
      if (cur_weight >= bin_size - tolerance) {
        cur_bin <- cur_bin + 1
        cur_weight <- 0
      }
    }
  }
  
  data.table::data.table(
    id   = out_id[1:(out_index - 1)],
    bin  = out_bin[1:(out_index - 1)],
    # You can drop these two if you only want bin assignment:
    weight = out_weight[1:(out_index - 1)],
    welfare = out_welfare[1:(out_index - 1)]
  )
}


new_bins <- function(welfare, weight, nbins = 100, tolerance = 1e-6, id = NULL) {
  stopifnot(length(welfare) == length(weight))
  
  # Handle NAs
  if (anyNA(welfare) || anyNA(weight)) {
    valid <- !is.na(welfare) & !is.na(weight)
    welfare <- welfare[valid]
    weight  <- weight[valid]
    if (!is.null(id)) id <- id[valid]
  } else {
    id <- if (is.null(id)) seq_along(welfare) else id
  }
  
  # Sort by welfare
  o <- order(welfare)
  welfare <- welfare[o]
  weight  <- weight[o]
  id      <- id[o]
  
  total_weight <- fsum(weight)
  bin_size     <- total_weight / nbins
  
  # Preallocate result vectors
  out_id      <- integer(2 * length(welfare))
  out_welfare <- numeric(2 * length(welfare))
  out_weight  <- numeric(2 * length(welfare))
  out_bin     <- integer(2 * length(welfare))
  
  cur_bin    <- 1
  cur_weight <- 0 # current among ot people in a bin
  out_index  <- 1 # index of the new vectors
  
  # loop over each row of the original data.
  for (i in seq_along(welfare)) {
    w  <- welfare[i]
    wt <- weight[i]
    id_i <- id[i]
    
    # as long as the weight to allocate is positive and
    # the value of the current bin is smaller or equal to
    # the number of bins proceed
    while (wt > 0 && cur_bin <= nbins) {
      # room in the bin available after subtracting
      # the people that is already in there.
      room <- bin_size - cur_weight
      
      # Move to next bin if very close to full
      if (abs(room) < tolerance) {
        cur_bin <- cur_bin + 1
        cur_weight <- 0
        next
      }
      
      #
      take <- min(wt, room)
      
      out_id[out_index]      <- id_i
      out_welfare[out_index] <- w
      out_weight[out_index]  <- take
      out_bin[out_index]     <- cur_bin
      out_index <- out_index + 1
      
      wt <- wt - take
      cur_weight <- cur_weight + take
      
      if (cur_weight >= bin_size - tolerance) {
        cur_bin <- cur_bin + 1
        cur_weight <- 0 # reset and go to the next bin
      }
    }
  }
  
  # Handle any leftover weight that couldn't be assigned due to tolerance
  if (cur_bin > nbins && wt > 0) {
    # Add remaining weight to the last bin
    out_id[out_index]      <- id_i
    out_welfare[out_index] <- w
    out_weight[out_index]  <- wt
    out_bin[out_index]     <- nbins
    out_index <- out_index + 1
  }
  
  # Return result
  data.table::data.table(
    id      = out_id[1:(out_index - 1)],
    bin     = out_bin[1:(out_index - 1)],
    weight  = out_weight[1:(out_index - 1)],
    welfare = out_welfare[1:(out_index - 1)]
  )
}



nearly_equal <- function(x, y, tol = 1e-12) abs(x - y) < tol








# Wrapper Function
duplicate_households <- function(DT,
                                 weight = "weight",
                                 threshold = 2.5,
                                 i = 0,
                                 li = 5,
                                 super_limit = 20) {
  R <- copy(DT)  # work on a copy to avoid modifying original DT
  
  # First iteration to make avoid goind through the whole
  # transformation of weights
  if (i == 0) {
    lt <- lorenz_table(R) # this is very inefficient, but that's what I have for now
    welfare_share_bad <-
      lt[, diff(lt$welfare_share),
         by = reporting_level
      ][, any(V1 < 0)]
    
    setattr(R, "welfare_share_OK", !welfare_share_bad)
    setattr(R, "threshold", threshold)
    setattr(R, "iterations", i)
    setattr(R, "lorenz", lt)
    if (!welfare_share_bad) {
      return(R)
    }
  }
  
  if (i >= super_limit)
    return(R)
  
  
  i = i + 1
  ori_names <- R |>
    names() |>
    copy()
  R <- find_outliers(R, weight, threshold)
  R <- duplicate_obs(R, weight)
  R <- add_new_weights(R, weight)
  R <- clean_new_weights(R, ori_names)
  lt <- lorenz_table(R) # this is very inefficient, but that's what I have for now
  
  welfare_share_bad <-
    lt[, diff(lt$welfare_share),
       by = reporting_level
    ][, any(V1 < 0)]
  
  setattr(R, "lorenz", lt)
  setattr(R, "welfare_share_OK", !welfare_share_bad)
  setattr(R, "threshold", threshold)
  setattr(R, "iterations", i)
  
  # if after li iteration still does nor, then lower the threshold and
  # increase the number of iteration
  if (welfare_share_bad && i >= li && threshold > 0) {
    threshold <- max(threshold - .5, 0)
    li <- li*2
  }
  
  if (welfare_share_bad) {
    R <- duplicate_households(R, weight = weight,
                              threshold = threshold,
                              i = i,
                              li = li,
                              super_limit = super_limit)
  }
  
  R
  
}




new_bins_old <- \(welfare, weight, nbins) {
  # deal with NAs -----
  if (anyNA(welfare)) {
    ina      <- !is.na(welfare)
    weight   <- weight[ina]
    welfare  <- as.numeric(welfare)[ina]
  }
  
  if (anyNA(weight)) {
    ina      <- !is.na(weight)
    weight   <- weight[ina]
    welfare  <- as.numeric(welfare)[ina]
  }
  
  # Sort data ------
  
  if (is.unsorted(welfare)) {
    o       <- order(welfare) # this is faster than collapse::radixorder
    welfare <- welfare[o]
    weight  <- weight[o]
  }
  
  p <- fcumsum(weight)/fsum(weight)
  bins <- ceiling(p * nbins) # Bins
  bins[bins > nbins] <- nbins
  bins
  # bins <-  cut(p, c(0, probs), labels = FALSE)
}




# civ <- pipload::pip_load_cache("CIV", 2002)
# civ2 <- replicate_households(civ)
#
# attributes(civ2)
# nrow(civ2)
# nrow(civ)
#
# wbpip::md_compute_gini(civ$welfare, civ$weight)
# wbpip::md_compute_gini(civ2$welfare, civ2$weight)
#
#
# ury <- pipload::pip_load_cache("PRY", 2018)
# ury2 <- replicate_households(ury)
#
#
#
# civ2  |>
#   ggplot(aes(x = weight)) +
#     geom_histogram(bins = 100) +
#     geom_vline(aes(xintercept=fsd(weight)*2.5+fmean(weight)),
#                color="blue", linetype="dashed", linewidth=1)
#
#
#
#
#
# civp1 <- lorenz_table(civ)
# civp2 <- lorenz_table(civ2)
#
# perr <- which(diff(civp1$welfare_share) < 0)
# perr2 <- which(diff(civp2$welfare_share) < 0)
#
# civp1[perr]
# civp1[perr2]
# civp2[c(perr2, perr2+1) |> sort()]
#
# dim(civ2)