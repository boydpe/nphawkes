# Peter Boyd
# Mass Shooting Functions

# library(Rcpp)
# library(ggplot2)
# library(gridExtra)
# nah....library(grid)
# library(lubridate)
# library(tidyverse)




# Nonparametric Hawkes ----------------------------------------------------

# maybe rename misd?
# only uses lubridate
#' Model Independent Stochastic Declustering
#'
#' This function uses nonparametric procedures to analyze a Hawkes
#' process in temporal or spatio-temporal domain, with or without marks.
#'
#'
#' @param dates a vector of dates, as LIST DATE FORMATS ALLOWED
#' @param lat a vector of latitudes, omit if not using spatial data
#' @param lon a vector of longitudes, omit if not using spatial data
#' @param marks a vecotr of marks, or magnitudes, omit if not using marked process
#' @param g_bins a vector of cutoff values for temporal bins of time differences
#' @param h_bins a vector of cutoff values for spatial bins of distance differences
#' @param k_bins a vector of cutoff values for magnitude bins
#' @param g_quantile FALSE by default to use defined temporal bins, TRUE to establish unifrom bins
#' @param g_length a scalar to define the number of desired uniform temporal bins CHANGE THIS IN CODE
#' @param h_quantile FALSE by default to use defined spatial bins, TRUE to establish uniform bins
#' @param h_length a scalar to define the number of desired unifrom spatial bins CHANGE THIS
#' @param k_quantile FALSE by default to use defined magnitude bins, TRUE to establish unifrom bins
#' @param k_length a scalar to define the number of desired magnitude bins
#' @param ref_date a date to serve as time 0, defaults to earliest observation
#' @param time_unit character string that specifies the desired unit of time
#' @param dist_unit character string that specifies the desired unit of distance: meter, kilometer, or mile
#' @param stop_when scalar that serves as conversion criterion, 1e-3 as default
#'
#' @return Probability matrix \code{p0} containing the probabilities that event
#' \code{i} is an offspring of event \code{j}, \code{i > j}. Diagonal elements
#' represent the probability that event \code{i} is a background event.
#' @return \code{g} is a vector of the estimated values for each bin of the temporal triggering component
#' @return \code{h} is a vector of the estimated values for each bin of the spatial triggering component
#' @return \code{k} is a vector of the estimataed values for each bin of the magnitude triggering component
#' @return \code{br} is the estimated background rate of the process
#' @return \code{perc_diag} is the proportion of mass lying on the diagonal of matrix \code{p0}
#' @return \code{time_bins} is a matrix containing the temporal bin of each pair of events
#'  @return \code{dist_bins} is a matrix containing the spatial bin of each pair of events
#'  @return \code{mark_bins} is a vector containing the magnitude bin of each event
#'  @return \code{input} is a list of all inputs


#' @export
nph <- function(dates, ref_date = min(dates),
                lat = rep(0, length(dates)),
                lon = rep(0, length(dates)),
                marks = rep(0, length(dates)),
                g_bins = c(0,1), h_bins = c(0,1),
                k_bins = c(0,1), g_quantile = FALSE,
                k_quantile = FALSE, h_quantile = FALSE,
                k_length = 6, h_length = 6,
                g_length = 6, stopwhen = 1e-3,
                time_unit = "day", dist_unit = "mile"){

  usethis::use_package("Rcpp")
  usethis::use_package("lubridate")

  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/0_bins.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/0_dist_matrix1.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/0_time_matrix.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/0_dist_bins.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/0_time_bins.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/0_mark_matrix.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/1_p0.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/2_br.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/3_trig_marks.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/3_trig_time.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/3_trig_space.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/4_update1.cpp")
  Rcpp::sourceCpp("~/OSU/PhD/Packages/nphawkes/R/5_check.cpp")

  # create times from dates
  df = data.frame(dates, lat, lon, marks)
  # put dates in correct format
  if(class(dates) == "Date"){
    dates_clean = dates
  } else {
    dates_clean = lubridate::mdy(dates)
  }
  # put refernce date in correct format
  if(class(ref_date) == "Date") {
    ref_date = ref_date
  } else {
    ref_date = lubridate::mdy(ref_date)
  }

  times = lubridate::time_length(lubridate:interval(ref_date, dates_clean), time_unit)
  times = times + runif(length(times), 0, 1)
  df$times = times
  df = df[order(df$times),]
  df = df[which(df$times > 0),]

  # output ordered events
  lat = df$lat
  lon = df$lon
  marks = df$marks
  times = df$times
  dates = df$dates

  # establish time, space differences for each event
  time_mat = get_time(times)
  dist_mat = get_dist(lat, lon)

  # convert distances to correct units
  if (dist_unit == "mile") {
    dist_mat = dist_mat*0.000621371
  } else if (dist_unit == "kilometer") {
    dist_mat = dist_mat*0.001
  } else {
    dist_mat = dist_mat
  }

  # implement option of bin by quantile
  if(g_quantile == TRUE){
    time_mat[lower.tri(time_mat, diag = TRUE)] = NA
    g_bins = as.vector(quantile(time_mat, na.rm = TRUE,
                                probs = seq(0,1, length.out= g_length)))
    g_bins[length(g_bins)] = g_bins[length(g_bins)] + 1
    g_bins[1] = 0
  } else {g_bins = g_bins}

  if(k_quantile == TRUE){
    k_bins = unique(as.vector(quantile(marks,
                                       probs = seq(0, 1, length.out = k_length))))
    k_bins[length(k_bins)] = k_bins[length(k_bins)] + 1
    k_bins[1] = 0 # accounts for simulating points below lowest bin
  } else {k_bins = k_bins}

  if(h_quantile == TRUE){
    dist_mat[lower.tri(dist_mat, diag = TRUE)] = NA
    h_bins = as.vector(quantile(dist_mat, na.rm = TRUE,
                                probs = seq(0,1, length.out= h_length)))
    h_bins[length(h_bins)] = h_bins[length(h_bins)] + 1
    h_bins[1] = 0
  } else {h_bins = h_bins}

  # place time, dist, and marks in bins
  time_bins = get_time_bins(time_mat, g_bins)
  mark_mat = get_mark(marks, k_bins)
  dist_bins = get_dist_bins(dist_mat, h_bins)

  # calculate br and trig components
  # update matrix, check if converge
  p0 = init_p0(times)
  max_diff = 1

  while( max_diff > stopwhen){
    br = calc_br(p0, times)
    g = get_g(p0, g_bins, time_mat)
    h = get_h(p0, h_bins, dist_mat)
    k = get_k(p0, marks, k_bins)
    p = update_p(p0, time_mat, dist_mat, mark_mat,
                 g, h, k,
                 h_bins, g_bins, k_bins,
                 br, time_bins, dist_bins)
    max_diff = check_p(p0, p)
    p0 = p
  }
  # is the max value per row background or trig?
  max_event = c()
  for ( i in 1:nrow(p0)){
    if (p0[i,i] == max(p0[i,])){
      max_event[i] = 1
    } else{
      max_event[i] = 0
    }
  }

  area = diff(k_bins)*k
  k_std = area / sum(area)

  perc_br = sum(max_event) / length(max_event)
  perc_diag = sum(diag(p0)) / nrow(p0)

  out = list(p0 = p0, g= g, h = h, k = k, br = br, k_std = k_std,
             time_bins = time_bins, mark_bins = mark_mat,
             dist_bins = dist_bins, perc_br = perc_br, perc_diag = perc_diag,
             g_bins = g_bins, k_bins = k_bins, h_bins = h_bins, data = df,
             ref_date = ref_date,
             input =   mget(names(formals()),sys.frame(sys.nframe())))
  # out = p0
  return(out)
}



# Conditional Intensity ---------------------------------------------------

#' Conditional Intensity
#'
#' This function estimates the conditional intensity function of the observed process.
#'
#' @param model the output from \code{nph()}
#'
#' @return a data frame containing the time, location, marks, and estimated conditional intensity
#' @export
cond_int = function(model) {

  times = model$data$times
  marks = model$data$marks
  lat = model$data$lat
  lon = model$data$lon

  # binning function
  g_bins = model$g_bins
  h_bins = model$h_bins
  k_bins = model$k_bins

  bin_f <- function(u, v) {
    x <- rep(0, length(u))
    for (j in 1:length(u)) {
      for (i in 1:(length(v) - 1)) {
        if (v[i] < u[j] & u[j] <= v[i + 1]) {
          x[j] <- i
        }
      }
    }
    return(x)
  }

  cond_int = c()
  # calculate br for obs 1 based on all data
  cond_int[1] = sum(diag(model$p0)) / (max(times) - min(times))

  n = length(times)
  g_values = matrix(nrow = n,
                    ncol = n,
                    data = rep(0, n * n))
  k_values = matrix(nrow = n,
                    ncol = n,
                    data = rep(0, n * n))
  h_values = matrix(nrow = n,
                    ncol = n,
                    data = rep(0, n * n))

  # create matrix of bin values for each obs in
  # 3 components. add 1 to convert from c++
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      g_values[i, j] = model$g[model$time_bins[j, i] + 1]
    }
  }

  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      h_values[i, j] = model$h[model$dist_bins[j, i] + 1]
    }
  }

  for (i in 2:n) {
      for (j in 1:(i - 1)) {
        k_values[i, j] = model$k[model$mark_bins[j] + 1]
      }
    }

  # get overall trig component per event
  trig = g_values * k_values * h_values
  br = model$br
  # get br and conditional intensity
  for (i in 2:n) {
    cond_int[i] = br + sum(trig[i, ])
  }

  ci = data.frame(times, lat, lon, marks, cond_int)
  ci$Date = as.Date(model$ref_date + ci$times)
  return(ci)
}

# Super-thinning ----------------------------------------------------------

#' Super-thinning
#'
#' This function performs super-thinning to assess model fit. Thinning only can be executed by
#' simply considering events of type "thin" or "retain", while superpositioning only can be executed
#' by simply treating "thinned" points as "retained".
#'
#' @param K a constant or character string that governs the amount of thinning and
#' superposing that is implemented. Can be a constant value, the median, mean, minimum, or maximum
#'  conditional intensity: "median_ci", "mean_ci", "min_ci", or "max_ci", respectively.
#'
#' @param model the output from \code{nph()}
#'
#' @return a data frame which includes the time, location, and estimted conditional intensity of events.
#' The type of event, either observed or simulated, is noted along with the probability that the event was kept
#' and wheter or not the point was in fact retained.
#' @export
super_thin = function(K = "median_ci",
                      model){

  # FIRST: CI FOR DATA
  times = model$data$times
  marks = model$data$marks
  lat = model$data$lat
  lon = model$data$lon
  # binning function
  g_bins = model$g_bins
  h_bins = model$h_bins
  k_bins = model$k_bins

  bin_f <- function(u,v){
    x <- rep(0,length(u))
    for (j in 1:length(u)){
      for (i in 1:(length(v)- 1)){
        if(v[i] < u[j] & u[j] <= v[i + 1]){
          x[j] <- i
        }
      }
    }
    return(x)
  }

  ci_data = cond_int(model)

  # SECOND: THIN DATA
  keep = c()
  p = c()
  type = c()

  if( K == "median_ci"){
    K = median(ci_data$cond_int)
  } else if(K == "mean_ci"){
    K = mean(ci_data$cond_int)
  } else if(K == "min_ci"){
    K = min(ci_data$cond_int)
  } else if(K == "max_ci"){
    K = max(ci_data$cond_int)
  } else {
    K = K
  }


  for(i in 1:nrow(ci_data)){
    p[i] = min(K/ci_data$cond_int[i],1)
    keep[i] = rbinom(1,1,p[i])
  }

  ci_data$p = p
  ci_data$keep = keep

  for(i in 1:nrow(ci_data)){
    type[i] = ifelse(keep[i] == 1, "retain", "thin")
  }
  ci_data$type = type
  ci_data$lat = lat
  ci_data$lon = lon

  # THIRD: SIMULATE POINTS

  # simulate, assign based on k.
  M1 = min(times)
  M2 = max(times)
  N = rpois(1, K*M2)
  N = ifelse(N ==0, 1, N) # not sure about adding this feature
  sim_time = sort(runif(N, M1, M2))
  sim_pp = data.frame(Time = sim_time, Obs = rep(0, N), marks = rep(NA, N))
  n2 = nrow(sim_pp)
  rate = (round(mean(marks)) - min(marks))^(-1)
  sim_pp$marks = round(rexp(n2, rate)) + min(marks)
  # not actually using simulated marks - so no parametric assumptions
  sim_lat = runif(n2, 0, 1)
  sim_lon = runif(n2, 0, 1)
  if(sum(model$lat) == 0){
    sim_pp = cbind(sim_pp, lat = rep(0,n2), lon = rep(0,n2))
  } else{
    sim_pp = cbind(sim_pp, lat = sim_lat, lon = sim_lon)
  }

  data_pp = data.frame(Time = times, marks = marks,
                       lat = lat, lon = lon)
  data_pp$Obs = rep(1, nrow(ci_data))
  # combine sim and real data
  all_pp = rbind(data_pp, sim_pp)
  all_pp = all_pp[order(all_pp$Time),]
  all_pp$row = c(1:(nrow(ci_data) + n2))
  r = 0
  for (i in 1:(nrow(ci_data) + n2)){
    r = ifelse(all_pp$Obs[i] == 1, r + 1, r)
    all_pp$row1[i] = r
  }
  sim_pp = all_pp[which(all_pp$Obs == 0),]


  cond_int_sim = c()
  br = model$br

  if(sum(sim_pp$lat) == 0){
    sim_pp$lat = rep(1,n2)
  } else {
    sim_pp$lat = sim_pp$lat
  }

  for (i in 1:n2) {
    sim_trig = 0
    # get triggering component for each
    for (j in 1:(sim_pp$row1[i])){# cut -1 from end
      td = sim_pp$Time[i] - data_pp$Time[j] # last bit was sim_data$row1[j]
      # change second to last obs time, need just obs data
      gb = bin_f(td, g_bins)
      gg = model$g[gb]

      kb = bin_f(data_pp$marks[j], k_bins) # was all_data, changed to data_pp
      kk = model$k[kb]

      #gk = gg*kk
      #sim_trig = sim_trig + gk

      hd = (sim_pp$lat[i] - data_pp$lat[j])^2 +
        (sim_pp$lon[i] - data_pp$lon[j])^2
      hb = bin_f(hd, h_bins)
      hh = model$h[hb]

      ghk = gg*kk*hh
      sim_trig = sim_trig + ghk
    }

    # get ci for each then remove the simulated point and correct row numbers
    cond_int_sim[i] = br + sim_trig
  }

  # FOURTH: THIN SIMULATED DATA

  ci_sim = data.frame(times = sim_pp$Time,
                      lat = sim_pp$lat, lon = sim_pp$lon,
                   cond_int = cond_int_sim, p = NA,
                   keep = NA, type = NA)
  # simulated data
  ci_sim$type = "sim"
  ci_sim$Date = as.Date(model$ref_date + ci_sim$times)
  for(i in 1:nrow(ci_sim)){
    ci_sim$p[i] = max((K - cond_int_sim[i])/K, 0)
    ci_sim$keep[i] = rbinom(1,1,ci_sim$p[i])
  }

  ci_sim = ci_sim[which(ci_sim$keep == 1),]
  ci_st = rbind(ci_data, ci_sim)
  ci_st = ci_st[order(ci_st$times),]
  #ci_st$rv = runif(nrow(ci_st), 0, 1)

  return(ci_st)
}


# Standard Error Bars -----------------------------------------------------

#' Standard Error Estimation
#'
#' This function estimates the variance of each bin for all utilized triggering components.
#'
#' @param model the output from \code{nph()}
#'
#' @return a list of variance estimates for all utilized triggering component


#' @export
se_bars = function(model){

  g_bins = model$g_bins
  h_bins = model$h_bins
  k_bins = model$k_bins

  nt = sum(model$p0) - sum(diag(model$p0))
  thetag = rep(0, length(model$g))
  thetak = rep(0, length(model$k))
  thetah = rep(0, length(model$h))

  # temporal standard error
  for(i in 2: nrow(model$p0)){
    for(j in 1:(i-1)){
      x_g = model$time_bins[j,i] + 1
      #should it  be + 1?
      thetag[x_g] = thetag[x_g] + model$p0[i,j]
    }
  }

  thetag = thetag / nt
  varg = thetag*(1 - thetag)/ (nt*(diff(g_bins))^2)

  # magnitude standard error
  for(i in 2: nrow(model$p0)){
    for(j in 1:(i-1)){
      x_k = model$mark_bins[j] + 1
      thetak[x_k] = thetak[x_k] + model$p0[i,j]
    }
  }

  thetak = thetak / nt
  vark = as.vector(nt * thetak*(1 - thetak)/ (table(model$mark_bins)^2))

  # spatial standard error
  for(i in 2: nrow(model$p0)){
    for(j in 1:(i-1)){
      x_h = model$dist_bins[j,i] + 1
      #should it  be + 1?
      thetag[x_h] = thetag[x_h] + model$p0[i,j]
    }
  }

  thetah = thetah / nt
  varh = thetah*(1 - thetah)/ (nt*(diff(h_bins))^2)

  out = list(varg = varg, vark = vark, varh = varh)
  return(out)
}


# Trig Plots ------------------------------------------------------------

# ggplot2 and gridExtra
#' Triggering Plots
#'
#' This function exports histogram estimators for all utilized triggering components. Plots show
#' estimated value of triggering functions for each bin along with \pm 2 standard errors bars
#'
#' @param model the output from \code{nph()}
#' @param g_max the maximum x axis value shown in the temporal plot
#' @param h_max the maximum x axis value shown in the spatial plot
#' @param k_min the minimum x axis value shown in the magnitude plot
#' @param k_max the maximum x axis value shown in the magnitude plot
#' @param mag_label character string representing what the magnitude measures
#'
#' @return histogram estimaotrs for all utilized triggering components
#' @export
trig_plots = function(model, g_max = max(model$g_bins),
                      k_max = max(model$k_bins),
                      k_min = min(model$k_bins),
                      h_max = max(model$h_bins),
                      mag_label = "magnitude"){

  g_bins = model$g_bins
  n1 = length(g_bins)
  k_bins = model$k_bins
  n2 = length(k_bins)
  h_bins = model$h_bins
  n3 = length(h_bins)

  ser = se_bars(model)
  se_g = sqrt(ser$varg)
  se_k = sqrt(ser$vark)
  se_h = sqrt(ser$varh)
  time_df = data.frame(g = c(model$g, model$g[length(model$g)]),
                       time = g_bins,
                       se = c(se_g, se_g[length(se_g)]))
  #time_df = time_df[-nrow(time_df),] may actually need this one

  mag_df = data.frame(k = c(model$k, model$k[length(model$k)]),
                      magnitude = k_bins,
                      se = c(se_k, se_k[length(se_k)]))

  dist_df = data.frame(h = c(model$h, model$h[length(model$h)]),
                        dist = h_bins,
                       se = c(se_h, se_h[length(se_h)]))

  # Time Plot
  trig_g = ggplot2::ggplot(time_df, aes(time, g)) +
    ggplot2::coord_cartesian(xlim = c(0, g_max)) +
    ggplot2::xlab(paste("t (time in ", model$input$time_unit, "s)", sep = "")) +
    ggplot2::ylab("g(t)")

  for (i in 1:(n1-1)){
    trig_g = trig_g +
      ggplot2::geom_polygon(aes(x = x, y = y),
                                   data = data.frame(
                                     x = c(time_df$time[i], time_df$time[i+1],
                                           time_df$time[i+1], time_df$time[i]),
                                     y = c(max(time_df$g[i] - 2*time_df$se[i], 0),
                                           max(time_df$g[i] - 2*time_df$se[i], 0),
                                           time_df$g[i] + 2*time_df$se[i],
                                           time_df$g[i] + 2*time_df$se[i]) ),
                                   fill = "grey")
  }

  trig_g = trig_g + ggplot2::geom_step(aes(group=1)) +
    ggplot2::scale_x_continuous(breaks = model$g_bins) +
    ggplot2::geom_point(aes(x = x, y = y, fill = type),
               shape = 21,
               data = data.frame(
                 x = sort(rep(time_df$time, 2))[-c(1,2,2*n1 -1 ,2*n1)],
                 y = c(rep(time_df$g[1:(n1-2)], each = 2)[-1], time_df$g[n1]),
                 type = c(rep(c("a", "b"), times = (n1-2)*2))
               )) +
    ggplot2::scale_fill_manual(values = c("black", "white")) +
    ggplot2::theme(legend.position = "none")

  # Magnitude Plot
  if (sum(model$mark_bins) != 0){
  trig_k = ggplot2::ggplot(mag_df, aes(magnitude, k)) +
    ggplot2::coord_cartesian(xlim = c(k_min, k_max)) +
    ggplot2::xlab(paste("m (", mag_label, ")", sep = "")) +
    ggplot2::ylab("k(m)")

  for (i in 1:(n2-1)){
    trig_k = trig_k + ggplot2::geom_polygon(aes(x = x, y = y),
                                   data = data.frame(
                                     x = c(mag_df$magnitude[i], mag_df$magnitude[i+1],
                                           mag_df$magnitude[i+1], mag_df$magnitude[i]),
                                     y = c(mag_df$k[i] - 2*mag_df$se[i],
                                           mag_df$k[i] - 2*mag_df$se[i],
                                           mag_df$k[i] + 2*mag_df$se[i],
                                           mag_df$k[i] + 2*mag_df$se[i])),
                                   fill = "grey")
  }

  trig_k = trig_k + ggplot2::geom_step(aes(group=1)) +
    ggpplot2::scale_x_continuous(breaks =
                         c(min(model$input$marks), model$k_bins[-1])) +
    ggplot2::geom_point(aes(x = x, y = y, fill = type),
               shape = 21,
               data = data.frame(
                 x = sort(rep(mag_df$magnitude, 2))[-c(1,2,2*n2 -1 ,2*n2)],
                 y = c(rep(mag_df$k[1:(n2-2)], each = 2)[-1], mag_df$k[n2]),
                 type = c(rep(c("a", "b"), times = (n2-2)*2))
               )) +
    ggplot2::scale_fill_manual(values = c("black", "white")) +
    ggplot2::theme(legend.position = "none")
  }

  #Space Plot
  if(sum(model$dist_bins) != 0) {

  trig_h = ggplot2::ggplot(dist_df, aes(dist, g)) +
    ggplot2::coord_cartesian(xlim = c(0, h_max)) +
    ggplot2::xlab(paste("s (distance in ", model$input$dist_unit, "s)", sep = "")) +
    ggplot2::ylab("h(s)")

  for (i in 1:(n1-1)){
    trig_h = trig_h + ggplot2::geom_polygon(aes(x = x, y = y),
                                   data = data.frame(
                                     x = c(dist_df$dist[i], dist_df$dist[i+1],
                                           dist_df$dist[i+1], dist_df$dist[i]),
                                     y = c(max(dist_df$h[i] - 2*dist_df$se[i], 0),
                                           max(dist_df$h[i] - 2*dist_df$se[i], 0),
                                           dist_df$h[i] + 2*dist_df$se[i],
                                           dist_df$h[i] + 2*dist_df$se[i]) ),
                                   fill = "grey")
  }

  trig_h = trig_h + ggplot2::geom_step(aes(group=1)) +
    ggplot2::scale_x_continuous(breaks = model$h_bins) +
    ggplot2::geom_point(aes(x = x, y = y, fill = type),
               shape = 21,
               data = data.frame(
                 x = sort(rep(dist_df$dist, 2))[-c(1,2,2*n3 -1 ,2*n3)],
                 y = c(rep(dist_df$h[1:(n3-2)], each = 2)[-1], dist_df$h[n3]),
                 type = c(rep(c("a", "b"), times = (n1-2)*2))
               )) +
    ggplot2::scale_fill_manual(values = c("black", "white")) +
    ggplot2::theme(legend.position = "none")
  }

  if (sum(model$mark_bins) == 0 & sum(model$dist_bins) == 0){
    out = gridExtra::grid.arrange(trig_g, ncol = 1)
  } else if (sum(model$mark_bins) != 0 & sum(model$dist_bins) == 0) {
    out = gridExtra::grid.arrange(trig_g, trig_k, ncol = 2)
  } else if (sum(model$mark_bins) == 0 & sum(model$dist_bins) != 0) {
    out = gridExtra::grid.arrange(trig_g, trig_h, ncol = 2)
  } else {
    out = gridExtra::grid.arrange(trig_g, trig_h, trig_k, ncol = 3)
  }
  return(out)
}



# Super-thinning Plot -----------------------------------------------------


# ggplot2
#' Super-thinning Plot
#'
#' This function exports a plot displaying the performance of the super-thinning procedure.
#'
#' @param superthin the output from \code{super_thin()}
#'
#' @return a tiered plot with superposed points on the top tier, points that were retained, not thinned,
#' on the middle tier, and points that were thinned on the bottom tier.
#' @export
st_plot = function(superthin){

  superthin$y = rep(0, nrow(superthin))
  for (i in 1:nrow(superthin)) {
    if (superthin$type[i] == "retain") {
      superthin$y[i] = 1
    } else if (superthin$type[i] == "thin") {
      superthin$y[i] = 0
    } else{
      superthin$y[i] = 2
    }
  }


  labs = c("thinned", "retained", "superimposed")
  labs = c("a", "b", "c")
  # plot1 = ggplot2::ggplot(data = superthin, aes(x = times, y = rv, col = type)) +
  #   ggplot2::geom_point() +
  #   ggplot2::theme(axis.title.y=element_blank(),
  #         axis.text.y=element_blank(),
  #         axis.ticks.y=element_blank())
  plot2 = ggplot2::ggplot(data = superthin, aes(x = Date, col = type)) +
    ggplot2::geom_segment(aes(x = Date, y = y, yend = y + 1, xend = Date),
                 show.legend = FALSE) +
    ggplot2::theme(axis.title.y=element_blank(),
          axis.ticks.y=element_blank()) +
          #axis.text.y=element_blank()) +
    ggplot2::scale_color_manual(breaks = c("sim", "retain", "thin"),
                       values=c("grey3", "grey40", "grey70")) +
    ggplot2::xlab("Year") +
    ggplot2::scale_y_continuous(breaks = c(0,1,2) + 0.5,
                       labels = c("thinned", "retained", "simulated"))
  # plot3 = ggplot2::ggplot(data = superthin, aes(x = Time, y = cond_int,
  #                                      col = type, shape = type)) +
  #   ggplot2::geom_point()

  #out = grid.arrange(plot2)
  out = plot2
  return(out)
}



# Condtitional Intensity Histogram -----------------------------------------

#' Histogram of super-thinned process
#'
#' This function exports a histogram of the super-thinned process.
#'
#' @param superthin the output from \code{super_thin()}
#' @param nbins scalar of the number of bins for the histogram
#' @param date_break character string of number followed by a time unit for the
#' desired time difference in between x axis labels
#' @param date_labels character string of % followed by first letter of time unit, i.e.
#' %Y for year, for desired label on x axis tick marks
#'
#' @return a histogram of the super-thinned process.
#' @export
ci_hist = function(superthin, nbins = 30,
                   date_break = "1 years", date_labels = "%Y"){
  st = superthin[which(superthin$type != "thin"),]
  out = ggplot2::ggplot(data = st) +
    ggplot2::geom_histogram(aes(x = Date), bins = nbins) +
    ggplot2::scale_x_date(date_breaks = date_break,
                          date_labels = date_label) +
    ggplot2::ylab("frequency")

  return(out)
}


# Conditional Intensity Plot ----------------------------------------------

#' Plot of Conditional Intensity Over Time
#'
#' In this function, the estimated conditional intensity is plotted against
#' the number of events. The estimated conditional intensity is found by taking the
#' median value in each month and multiplying it by the number of days in that month.
#'
#' @param model the output from \code{nph()}
#' @param superthin the output from \code{super_thin()}
#' @param min_date the minimum date to be shown on the plot
#' @param max_date the maximum date to be shown on the plot
#'
#' @return a plot that shows the estimated conditional intensity with the number of monthly events
#' @export
ci_plot = function(model, min_date = model$input$ref_date,
                   max_date = max(model$data$dates),
                   superthin){

  df = model$data
  df$Year = year(df$dates)
  df$Month = month(df$dates)

  min = as.Date(min_date)
  max = as.Date(max_date)

  df_counts  = df %>% dplyr::group_by(Year, Month) %>%
    dplyr::summarize(n = n(), Day = 1) %>%
    dplyr::mutate(date = lubridate::make_date(Year, Month, Day))

  ci_one = superthin %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::summarize(cond_int_sum = sum(cond_int),
            cond_int_mean = mean(cond_int),
            cond_int_med = median(cond_int),
            Day = 1) %>%
    dplyr::mutate(Date = lubridate::make_date(Year, Month, Day)) %>%
    dplyr::mutate(ndays = lubridate::days_in_month(Date)) %>%
    dplyr::mutate(n = ndays*cond_int_med) %>%
    dplyr::select(-c(cond_int_sum, cond_int_mean, cond_int_med, ndays)) %>%
    dplyr::mutate(type = "est")

  tmp = data.frame(Date = seq(min, max, by = "month"))
  tmp %>% dplyr::left_join(df_stanford_counts) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
    dplyr::mutate(type = "data") %>%
    dplyr::bind_rows(ci_one) %>%
    ggplot2::ggplot() +
      ggplot2::geom_line(aes(x = Date, y = n, linetype = type)) +
      ggplot2::theme(axis.title.y=element_blank(),
        legend.position = "right") +
        ggplot2::scale_linetype_manual(name = "Number of \nMonthly Events",
                        labels = c("Observed", "Estimated"),
                        values = c("solid", "dashed"))
}

