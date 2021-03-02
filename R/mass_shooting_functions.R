# Peter Boyd
# Mass Shooting Functions


# Rcpp::compileAttributes() #when changing c++ code
# devtools::document() #to update whole thing
# devtools::load_all() #to temp install
# library(nphawkes) #to use

# Nonparametric Hawkes ----------------------------------------------------

# maybe rename misd?
# only uses lubridate

#' Model Independent Stochastic Declustering
#'
#' This function uses nonparametric procedures to analyze a Hawkes
#' process in temporal or spatio-temporal domain, with or without marks
#' through the Model independent stochastic declustering algorithm.
#'
#'
#' @param dates a vector of dates as "yyyy-mm-dd".
#' @param lat a vector of latitudes, omit if not using spatial data
#' @param lon a vector of longitudes, omit if not using spatial data
#' @param marks a vecotr of marks, or magnitudes, omit if not using marked data
#' @param time_breaks a vector of cutoff values for temporal bins of time differences
#' @param space_breaks a vector of cutoff values for spatial bins of distance differences
#' @param mark_breaks a vector of cutoff values for magnitude bins
#' @param time_quantile FALSE by default to use defined temporal bins, TRUE to establish uniform bins
#' @param g_length a scalar to define the number of desired uniform temporal bins
#' @param space_quantile FALSE by default to use defined spatial bins, TRUE to establish uniform bins
#' @param h_length a scalar to define the number of desired unifrom spatial bins
#' @param mark_quantile FALSE by default to use defined magnitude bins, TRUE to establish uniform bins
#' @param k_length a scalar to define the number of desired magnitude bins
#' @param ref_date a date to serve as time 0, defaults to earliest observation
#' @param time_of_day character string that lists the time of day of events, as hour:minute:second
#' @param just_time TRUE or FALSE object. TRUE if \code{dates} object is a vector of only times,
#' indicating time elapsed since beginning of catalog.
#' @param time_unit character string that specifies the desired unit of time
#' @param dist_unit character string that specifies the desired unit of distance: meter, kilometer, or mile
#' @param stop_when scalar that serves as conversion criterion, 1e-3 as default
#'
#' This function can only be applied to data that contains a temporal feature. It can also be applied to data
#' consisting of time and space, time and marks, or time and space and marks.
#'
#' For each triggering component used (time, space, marks), a binning structure will be applied. The user may define
#' these right continuous bins as a vector ([c(1,5,10)] creates two bins for 1 < x \leq 5, and 5 < x \leq 10),
#' or may be generated automatically by specifying \code{time_breaks = TRUE}, with the same applying for marks and space.
#' This method will establish breaks such that the allocation of time differences within the data will be roughly
#' equal in each created bin. Further specifying \code{g_length} will control the number of bins established, with the
#' default producing 6 bins. Uniform binning methods may unattainable for discrete marks containing many replicates of values.
#'
#' If no time of day is provided, events will be randomly assigned a time during the event's date.
#'
#'
#' @return Probability matrix \code{p0} containing the probabilities that event
#' \code{i} is an offspring of event \code{j}, \code{i > j}. Diagonal elements
#' represent the probability that event \code{i} is a background event.
#' @return \code{g} is a vector of the estimated values for each bin of the temporal triggering component
#' @return \code{h} is a vector of the estimated values for each bin of the spatial triggering component
#' @return \code{k} is a vector of the estimataed values for each bin of the magnitude triggering component
#' @return \code{br} is the estimated background rate of the process
#' @return \code{perc_diag} is the proportion of mass lying on the diagonal of matrix \code{p0}
#' @return \code{perc_br} is the proportion of events in which the maximum probabilistic assignment
#' is as a background event
#' @return \code{time_bins} is a matrix containing the temporal bin of each pair of events
#' @return \code{dist_bins} is a matrix containing the spatial bin of each pair of events
#' @return \code{mark_bins} is a vector containing the magnitude bin of each event
#' @return \code{input} is a list of all inputs


#' @export
nph <- function(dates, ref_date = min(dates),
                lat = rep(0, length(dates)),
                lon = rep(0, length(dates)),
                marks = rep(0, length(dates)),
                time_breaks = c(0,1), space_breaks = c(0,1),
                mark_breaks = c(0,1), time_quantile = FALSE,
                mark_quantile = FALSE, space_quantile = FALSE,
                k_length = 6, h_length = 6,
                g_length = 6, stopwhen = 1e-3,
                time_of_day = NA,
                just_times = FALSE,
                time_unit = "day", dist_unit = "mile"){

  # create times from dates
  df = data.frame(dates, lat, lon, marks)
  # put dates in correct format
  if (just_times == FALSE) {
    if (class(dates) == "Date") {
      dates_clean = dates
    } else {
      dates_clean = lubridate::as_date(dates)
    }

    ref_date = lubridate::as_date(ref_date)

    times = lubridate::time_length(lubridate::interval(ref_date, dates_clean), time_unit)
    if (is.na(time_of_day[1]) == TRUE) {
      times = times + runif(length(times), 0, 1)
    } else {
      times_clean = lubridate::hms(time_of_day)
      h = lubridate::hour(times_clean) / 24
      m = lubridate::minute(times_clean) / 1440
      s = lubridate::second(times_clean) / 86400
      times = times + h + m + s
    }
  } else {
    times = dates
  }

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

  n = length(lat)
  dist_mat = matrix(data = NA, nrow = n, ncol = n)
  R = 6371e3

  for (i in 1:n) {
    for (j in 1:n) {
      phi1 = lat[i] * pi / 180
      phi2 = lat[j] * pi / 180
      latdiff = (lat[j] - lat[i]) * pi / 180
      londiff = (lon[j] - lon[i]) * pi / 180

      a = sin(latdiff / 2) * sin(latdiff / 2) +
        cos(phi1) * cos(phi2) * sin(londiff / 2) * sin(londiff / 2)
      b = 2 * atan2(sqrt(a), sqrt(1 - a))
      d = R * b

      dist_mat[i, j] = d
    }
  }


  # convert distances to correct units
  if (dist_unit == "mile") {
    dist_mat = dist_mat*0.000621371
  } else if (dist_unit == "kilometer") {
    dist_mat = dist_mat*0.001
  } else {
    dist_mat = dist_mat
  }

  # implement option of bin by quantile
  if(time_quantile == TRUE){
    time_mat1 = time_mat
    time_mat1[lower.tri(time_mat1, diag = TRUE)] = NA
    time_breaks = as.vector(quantile(time_mat1, na.rm = TRUE,
                                probs = seq(0,1, length.out= g_length + 1)))
    time_breaks[length(time_breaks)] = time_breaks[length(time_breaks)] + 1
    time_breaks[1] = 0
  } else {time_breaks = time_breaks}

  if(mark_quantile == TRUE){
    mark_breaks = unique(as.vector(quantile(marks,
                                       probs = seq(0, 1, length.out = k_length + 1))))
    mark_breaks[length(mark_breaks)] = mark_breaks[length(mark_breaks)] + 1
    mark_breaks[1] = 0 # accounts for simulating points below lowest bin
  } else {mark_breaks = mark_breaks}

  if(space_quantile == TRUE){
    dist_mat1 = dist_mat
    dist_mat1[lower.tri(dist_mat1, diag = TRUE)] = NA
    space_breaks = as.vector(quantile(dist_mat1, na.rm = TRUE,
                                probs = seq(0,1, length.out= h_length + 1)))
    space_breaks[length(space_breaks)] = space_breaks[length(space_breaks)] + 1
    space_breaks[1] = 0
  } else {space_breaks = space_breaks}

  # place time, dist, and marks in bins
  time_bins = get_time_bins(time_mat, time_breaks)
  mark_mat = get_mark(marks, mark_breaks)
  dist_bins = get_dist_bins(dist_mat, space_breaks)

  # calculate br and trig components
  # update matrix, check if converge
  p0 = init_p0(times)
  max_diff = 1

  while( max_diff > stopwhen){
    br = calc_br(p0, times)
    g = get_g(p0, time_breaks, time_mat)
    h = get_h(p0, space_breaks, dist_mat)
    k = get_k(p0, marks, mark_breaks)
    p = update_p(p0, time_mat, dist_mat, mark_mat,
                 g, h, k,
                 space_breaks, time_breaks, mark_breaks,
                 br, time_bins, dist_bins, lat)
    max_diff = check_p(p0, p)
    p0 = p
  }

  max_event = c()
  for ( i in 1:nrow(p0)){
    if (p0[i,i] == max(p0[i,])){
      max_event[i] = 1
    } else{
      max_event[i] = 0
    }
  }

  perc_br = sum(max_event) / length(max_event)
  perc_diag = sum(diag(p0)) / nrow(p0)

  out = list(p0 = p0, g= g, h = h, k = k, br = br,
             time_bins = time_bins, mark_bins = mark_mat,
             dist_bins = dist_bins, perc_br = perc_br, perc_diag = perc_diag,
             time_breaks = time_breaks, mark_breaks = mark_breaks, space_breaks = space_breaks, data = df,
             ref_date = ref_date,
             input =   mget(names(formals()),sys.frame(sys.nframe())))
  return(out)
}



# Conditional Intensity ---------------------------------------------------

#' Conditional Intensity
#'
#' This function estimates the conditional intensity function of the observed process.
#'
#' @param model the output from \code{nph()}
#'
#' This function is to be used in conjunction with the \code{nph()} function from the \code{nphawkes} library.
#' Using the output from the \code{nph()} function,
#'
#' @return a data frame containing the time, location, marks, and estimated conditional intensity
#' @export
cond_int = function(model) {

  times = model$data$times
  marks = model$data$marks
  lat = model$data$lat
  lon = model$data$lon

  # binning function
  time_breaks = model$time_breaks
  space_breaks = model$space_breaks
  mark_breaks = model$mark_breaks

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
#' @param model the output from \code{nph()}
#' @param method character string that defines residual analysis method as "superthin", "thin", or "superpose"
#' @param map name of map provided by the maps package, defaults to maps::world
#' @param region name of subregion to include, defaults to entirety of map
#' @param sim_grid TRUE if geographical coordinates do not necessarily pertain to
#' certain area, such as 0x1 by 0x1 grid, FALSE if using specific region or not using spatial information
#' @param lat_bounds vector containing minimum and maximum latitude bounds if sim_grid is TRUE
#' @param lon_bounds vector containing minimum and maximum longitude bounds if sim_gird is TRUE
#'
#' This function is to be used in conjunction with the \code{nph()} function from the \code{nphawkes} library.
#'
#' To simulate spatial data, the user may define the \code{map} and \code{region} to easily simulate points within
#' set political borders. Otherwise, simulated points may be established when \code{sim_grid = TRUE}.
#'
#' The parameter \code{K}
#'
#' @return a data frame which includes the time, location, and estimated conditional intensity of events.
#' The type of event, either observed or simulated, is noted along with the probability that the event was kept
#' and wheter or not the point was in fact retained.
#' @export
super_thin = function(K = "median_ci",
                      model, method = "superthin",
                      map = world, region = ".",
                      sim_grid = TRUE,
                      lat_bounds = c(min(model$data$lat), max(model$data$lat)),
                      lon_bounds = c(min(model$data$lon), max(model$data$lon))) {

  # FIRST: CI FOR DATA
  times = model$data$times
  marks = model$data$marks
  lat = model$data$lat
  lon = model$data$lon
  # binning function
  time_breaks = model$time_breaks
  space_breaks = model$space_breaks
  mark_breaks = model$mark_breaks

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
  rate = (mean(marks) - min(marks))^(-1)
  sim_pp$marks = round(rexp(n2, rate)) + min(marks)
  # not actually using simulated marks
  if (sum(model$lat) == 0) {
    sim_pp = cbind(sim_pp, lat = rep(0, n2), lon = rep(0, n2))
  } else{
    if (sim_grid == TRUE) {
      sim_lat = runif(n2, min = lat_bounds[1], max = lat_bounds[2])
      sim_lon = runif(n2, min = lon_bounds[1], max = lat_bounds[2])
      sim_pp = cbind(sim_pp, lat = sim_lat, lon = sim_lon)
    } else {
      region = ggplot2::map_data(maps::map(), region = region)
      region = region[, c('long', 'lat')] %>%
        sp::Polygon() %>%
        sp::spsample(n = n2, type = "random")
      sim_lat = region$y
      sim_lon = region$x
      sim_pp = as.data.frame(cbind(sim_pp, lat = sim_lat, lon = sim_lon))
    }
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
    for (j in 1:(sim_pp$row1[i])){
      td = sim_pp$Time[i] - data_pp$Time[j]
      gb = bin_f(td, time_breaks)
      gg = model$g[gb]

      kb = bin_f(data_pp$marks[j], mark_breaks)
      kk = model$k[kb]

      hd = (sim_pp$lat[i] - data_pp$lat[j])^2 +
        (sim_pp$lon[i] - data_pp$lon[j])^2
      hb = bin_f(hd, space_breaks)
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
                   keep = NA, type = NA, marks = sim_pp$marks)
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

  if (method == "superthin") {
    ci_st = ci_st
  } else if (method == "thin") {
    ci_st = ci_st %>% dplyr::filter(type != "sim")
  } else {
    ci_st = ci_st %>% dplyr::mutate(replace(type, type == "thin", "retain"))
  }
  return(ci_st)
}


# Standard Error Bars -----------------------------------------------------

#' Standard Error Estimation
#'
#' This function estimates the variance of each bin for all utilized triggering components.
#'
#' @param model the output from \code{nph()}
#'
#' This function is to be used in conjunction with the \code{nph()} function from the \code{nphawkes} library.
#' and is used within the function that produces triggering plots.
#'
#' @return a list of variance estimates for all utilized triggering component


#' @export
se_bars = function(model){

  time_breaks = model$time_breaks
  space_breaks = model$space_breaks
  mark_breaks = model$mark_breaks

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
  varg = thetag*(1 - thetag)/ (nt*(diff(time_breaks))^2)

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
      thetah[x_h] = thetah[x_h] + model$p0[i,j]
    }
  }

  thetah = thetah / nt
  varh = thetah*(1 - thetah)/ (nt*(diff(space_breaks))^2)

  out = list(varg = varg, vark = vark, varh = varh)
  return(out)
}


# Trig Plots ------------------------------------------------------------

#' Triggering Plots
#'
#' This function exports histogram estimators for all utilized triggering components. Plots show
#' estimated value of triggering functions for each bin along with \pm 2 standard errors bars
#'
#' @param model the output from \code{nph()}
#' @param g_xlim vector of minimum and maximum x-axis value shown in the temporal plot
#' @param h_xlim vector of minimum and maximum x-axis value shown in the spatial plot
#' @param k_xlim vector of minimum and maximum x-axis value shown in the magnitude plot
#' @param g_ylim vector of minimum and maximum y-axis value shown in the temporal plot
#' @param h_ylim vector of minimum and maximum y-axis value shown in the spatial plot
#' @param k_ylim vector of minimum and maximum y-axis value shown in the magnitude plot
#' @param mag_label character string representing what the magnitude measures
#' @param title character string for the title of the collection of triggering plots
#'
#' @return \code{all_plots} is a cowplot object of histogram estimators for all utilized triggering components
#' @return \code{time_plot} is a ggplot object of the histogram estimator for the temporal effect
#' @return \code{space_plot} is a ggplot object of the histogram estimator for the spatial effect,
#' NULL if data provided is not a spatial process
#' @return \code{mark_plot} is a ggplot object of the histogram estimator for the mark effect,
#' NULL if data provided is not a marked process

#' @export
trig_plots = function(model,
                      g_xlim = c(min(model$time_breaks), max(model$time_breaks)),
                      h_xlim = c(min(model$space_breaks), max(model$space_breaks)),
                      k_xlim = c(min(model$mark_breaks), max(model$mark_breaks)),
                      g_ylim = c(0,NA),
                      h_ylim = c(0,NA),
                      k_ylim = c(0,NA),
                      mag_label = "magnitude",
                      plot_title = NULL){

  time_breaks = model$time_breaks
  n1 = length(time_breaks)
  mark_breaks = model$mark_breaks
  mark_breaks[1] = min(model$data$marks)

  n2 = length(mark_breaks)
  space_breaks = model$space_breaks
  n3 = length(space_breaks)

  ser = se_bars(model)
  se_g = sqrt(ser$varg)
  se_k = sqrt(ser$vark)
  se_h = sqrt(ser$varh)
  time_df = data.frame(g = c(model$g, model$g[length(model$g)]),
                       time = time_breaks,
                       se = c(se_g, se_g[length(se_g)]))

  mag_df = data.frame(k = c(model$k, model$k[length(model$k)]),
                      magnitude = mark_breaks,
                      se = c(se_k, se_k[length(se_k)]))

  dist_df = data.frame(h = c(model$h, model$h[length(model$h)]),
                        dist = space_breaks,
                       se = c(se_h, se_h[length(se_h)]))

  g_ylim[2] = ifelse(is.na(g_ylim[2]) == TRUE,
                     max(time_df$g + time_df$se), g_ylim[2])
  h_ylim[2] = ifelse(is.na(h_ylim[2]) == TRUE,
                     max(dist_df$h + dist_df$se), h_ylim[2])
  k_ylim[2] = ifelse(is.na(k_ylim[2]) == TRUE,
                     max(mag_df$k + mag_df$se), k_ylim[2])


  # Time Plot
  trig_g = ggplot2::ggplot(time_df, ggplot2::aes(time, g)) +
    ggplot2::coord_cartesian(xlim = g_xlim, ylim = g_ylim) +
    ggplot2::xlab(paste("t (time in ", model$input$time_unit, "s)", sep = "")) +
    ggplot2::ylab("g(t)") +
    ggplot2::theme(axis.title.x = element_text(size = 20),
                   axis.title.y = element_text(size = 20),
                   axis.text = element_text(size = 20))

  for (i in 1:(n1-1)){
    trig_g = trig_g +
      ggplot2::geom_polygon(ggplot2::aes(x = x, y = y),
                                   data = data.frame(
                                     x = c(time_df$time[i], time_df$time[i+1],
                                           time_df$time[i+1], time_df$time[i]),
                                     y = c(max(time_df$g[i] - 2*time_df$se[i], 0),
                                           max(time_df$g[i] - 2*time_df$se[i], 0),
                                           time_df$g[i] + 2*time_df$se[i],
                                           time_df$g[i] + 2*time_df$se[i]) ),
                                   fill = "grey")
  }

  trig_g = trig_g + ggplot2::geom_step(ggplot2::aes(group=1)) +
    ggplot2::scale_x_continuous(breaks = round(model$time_breaks, 1),
                                minor_breaks = time_breaks) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y, fill = type),
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
  trig_k = ggplot2::ggplot(mag_df, ggplot2::aes(magnitude, k)) +
    ggplot2::coord_cartesian(xlim = k_xlim, ylim = k_ylim) +
    ggplot2::xlab(paste("m (", mag_label, ")", sep = "")) +
    ggplot2::ylab("k(m)") +
    # delete this
    ggplot2::theme(axis.title.x = element_text(size = 20),
                   axis.title.y = element_text(size = 20),
                   axis.text = element_text(size = 20))

  for (i in 1:(n2-1)){
    trig_k = trig_k + ggplot2::geom_polygon(ggplot2::aes(x = x, y = y),
                                   data = data.frame(
                                     x = c(mag_df$magnitude[i], mag_df$magnitude[i+1],
                                           mag_df$magnitude[i+1], mag_df$magnitude[i]),
                                     y = c(max(mag_df$k[i] - 2*mag_df$se[i], 0),
                                           max(mag_df$k[i] - 2*mag_df$se[i], 0),
                                           mag_df$k[i] + 2*mag_df$se[i],
                                           mag_df$k[i] + 2*mag_df$se[i])),
                                   fill = "grey")
  }

  trig_k = trig_k + ggplot2::geom_step(ggplot2::aes(group=1)) +
    ggplot2::scale_x_continuous(breaks =
                         c(round(min(model$input$marks),1),
                           round(model$mark_breaks[-1],1)),
                         minor_breaks = c(round(min(model$input$marks),1),
                                          round(model$mark_breaks[-1],1))) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y, fill = type),
               shape = 21,
               data = data.frame(
                 x = sort(rep(mag_df$magnitude, 2))[-c(2, 2*n2 -1 ,2*n2)],
                 y = c(rep(mag_df$k[1:(n2-2)], each = 2), mag_df$k[n2]),
                 type = c("a", rep(c("a", "b"), times = (n2-2)))
               )) +
    ggplot2::scale_fill_manual(values = c("black", "white")) +
    ggplot2::theme(legend.position = "none")
  } else {
    trig_k = NULL
  }

  #Space Plot
  if(sum(model$dist_bins) != 0) {

  trig_h = ggplot2::ggplot(dist_df, ggplot2::aes(dist, h)) +
    ggplot2::coord_cartesian(xlim = h_xlim, ylim = h_ylim) +
    ggplot2::xlab(paste("s (distance in ", model$input$dist_unit, "s)", sep = "")) +
    ggplot2::ylab("h(s)")

  for (i in 1:(n1-1)){
    trig_h = trig_h + ggplot2::geom_polygon(ggplot2::aes(x = x, y = y),
                                   data = data.frame(
                                     x = c(dist_df$dist[i], dist_df$dist[i+1],
                                           dist_df$dist[i+1], dist_df$dist[i]),
                                     y = c(max(dist_df$h[i] - 2*dist_df$se[i], 0),
                                           max(dist_df$h[i] - 2*dist_df$se[i], 0),
                                           dist_df$h[i] + 2*dist_df$se[i],
                                           dist_df$h[i] + 2*dist_df$se[i]) ),
                                   fill = "grey")
  }

  trig_h = trig_h + ggplot2::geom_step(ggplot2::aes(group=1)) +
    ggplot2::scale_x_continuous(breaks = round(model$space_breaks, 1)) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y, fill = type),
               shape = 21,
               data = data.frame(
                 x = sort(rep(dist_df$dist, 2))[-c(1,2,2*n3 -1 ,2*n3)],
                 y = c(rep(dist_df$h[1:(n3-2)], each = 2)[-1], dist_df$h[n3]),
                 type = c(rep(c("a", "b"), times = (n3-2)*2))
               )) +
    ggplot2::scale_fill_manual(values = c("black", "white")) +
    ggplot2::theme(legend.position = "none")

  } else {
    trig_h = NULL
  }

  # if (sum(model$mark_bins) == 0 & sum(model$dist_bins) == 0){
  #   out = gridExtra::grid.arrange(trig_g, ncol = 1)
  # } else if (sum(model$mark_bins) != 0 & sum(model$dist_bins) == 0) {
  #   out = gridExtra::grid.arrange(trig_g, trig_k, ncol = 2)
  # } else if (sum(model$mark_bins) == 0 & sum(model$dist_bins) != 0) {
  #   out = gridExtra::grid.arrange(trig_g, trig_h, ncol = 2)
  # } else {
  #   out = gridExtra::grid.arrange(trig_g, trig_h, trig_k, ncol = 3)
  # }
  if (sum(model$mark_bins) == 0 & sum(model$dist_bins) == 0){
    out = cowplot::plot_grid(trig_g, ncol = 1)
  } else if (sum(model$mark_bins) != 0 & sum(model$dist_bins) == 0) {
    out = cowplot::plot_grid(trig_g, trig_k, ncol = 2)
  } else if (sum(model$mark_bins) == 0 & sum(model$dist_bins) != 0) {
    out = cowplot::plot_grid(trig_g, trig_h, ncol = 2)
  } else {
    out = cowplot::plot_grid(trig_g, trig_h, trig_k, ncol = 3)
  }

  title = cowplot::ggdraw() +
    cowplot::draw_label(plot_title, x = 0, hjust = 0) +
    ggplot2::theme_set(cowplot::theme_cowplot(font_size=20)) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 7))

  out = cowplot::plot_grid(title, out, ncol = 1, rel_heights = c(0.15, 1))
  return(all_plots = out, time_plot = trig_g,
         space_plot = trig_h, mark_plot = trig_k)
}



# Super-thinning Plot -----------------------------------------------------

#' Super-thinning Plot
#'
#' This function exports a plot displaying the performance of the super-thinning procedure.
#'
#' @param superthin the output from \code{super_thin()}
#' @param method chacter string denoting if the residual analysis method utilized was
#' "superthin", "thin", or "superpose"
#' @param time_label character string that provides the frequency of tick marks on the x-axis
#' @param plot_title character string for the plot title
#'
#' @return a tiered plot with superposed points on the top tier, points that were retained, not thinned,
#' on the middle tier, and points that were thinned on the bottom tier.
#' @export
st_plot = function(superthin, method = "superthin",
                   time_label = "Year", plot_title = NULL){

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

  if (method == "superthin"){
    p = ggplot2::ggplot(data = superthin, ggplot2::aes(x = Date, col = type)) +
    ggplot2::geom_segment(ggplot2::aes(x = Date, y = y, yend = y + 1, xend = Date),
                 show.legend = FALSE) +
    ggplot2::theme(axis.title.y= ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank()) +
          #axis.text.y=element_blank()) +
    ggplot2::scale_color_manual(breaks = c("sim", "retain", "thin"),
                       values=c("grey3", "grey40", "grey70")) +
    ggplot2::xlab(time_label) +
    ggplot2::scale_y_continuous(breaks = c(0,1,2) + 0.5,
                       labels = c("Thinned", "Retained", "Simulated"))

  } else if (method == "thin"){
    p = superthin %>% dplyr::filter(type != "sim") %>%
    ggplot2::ggplot(ggplot2::aes(x = Date, col = type)) +
      ggplot2::geom_segment(ggplot2::aes(x = Date, y = y, yend = y + 1, xend = Date),
                            show.legend = FALSE) +
      ggplot2::theme(axis.title.y= ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank()) +
      #axis.text.y=element_blank()) +
      ggplot2::scale_color_manual(breaks = c("retain", "thin"),
                                  values=c("grey40", "grey70")) +
      ggplot2::xlab(time_label) +
      ggplot2::scale_y_continuous(breaks = c(0,1,2) + 0.5,
                                  labels = c("Thinned", "Retained"))

  } else {
    p = superthin %>% dplyr::mutate(type = replace(type, type == "thin", "retain")) %>%
      dplyr::mutate(replace(type, type == "retain", "observed")) %>%
      ggplot2::ggplot(data = superthin, ggplot2::aes(x = Date, col = type)) +
      ggplot2::geom_segment(ggplot2::aes(x = Date, y = y, yend = y + 1, xend = Date),
                            show.legend = FALSE) +
      ggplot2::theme(axis.title.y= ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank()) +
      ggplot2::scale_color_manual(breaks = c("sim", "observed"),
                                  values=c("grey3", "grey40")) +
      ggplot2::xlab(time_label) +
      ggplot2::scale_y_continuous(breaks = c(0,1,2) + 0.5,
                                  labels = c("Observed", "Simulated")) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12,
                                                        margin = ggplot2::margin(t = 8, b = -20)))
  }
  return(p)
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
#' @param plot_title character string for the title of the plot
#'
#' @return a histogram of the super-thinned process.
#' @export
ci_hist = function(superthin, nbins = 30,
                   date_break = "1 year", date_label = "%Y",
                   plot_title = NULL){
  st = superthin[which(superthin$type != "thin"),]
  out = ggplot2::ggplot(data = st) +
    ggplot2::geom_histogram(ggplot2::aes(x = Date), bins = nbins) +
    # ggplot2::scale_x_date(date_breaks = date_break,
    #                       date_labels = date_label) +
    ggplot2::ylab("frequency") +
    ggplot2::ggtitle(plot_title)
    # ggplot2::theme(plot.title = ggplot2::element_text(size = 12,
    #                                                   margin = ggplot2::margin(t = 8, b = -20)))

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
#' @param plot_title character string for the title of the plot
#'
#' @return a plot that shows the estimated conditional intensity as a black line with the number of monthly events as
#' gray vertical lines
#' @export
ci_plot = function(model, min_date = model$input$ref_date,
                   max_date = model$data$dates[nrow(model$data)],
                   superthin,
                   plot_title = NULL){

  df = model$data
  df$Year = lubridate::year(df$dates)
  df$Month = lubridate::month(df$dates)

  min = as.Date(min_date)
  max = as.Date(max_date)

  df_counts  = df %>% dplyr::group_by(Year, Month) %>%
    dplyr::summarize(n = dplyr::n(), Day = 1) %>%
    dplyr::mutate(Date = lubridate::make_date(Year, Month, Day))

  ci_one = superthin %>%
    dplyr::mutate(Year = lubridate::year(superthin$Date)) %>%
    dplyr::mutate(Month = lubridate::month(superthin$Date)) %>%
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

  tmp = data.frame(Date = seq.Date(min, max, by = "month"))
  tmp %>% dplyr::left_join(df_counts) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
    dplyr::mutate(type = "data") %>%
    dplyr::bind_rows(ci_one) %>%
    ggplot2::ggplot() +
      #ggplot2::geom_line(ggplot2::aes(x = Date, y = n, linetype = type)) +
    ggplot2::geom_segment(ggplot2::aes(x = Date, y = n,
                              xend = Date, yend = 0),
                          alpha = 0.5) +
    ggplot2::geom_line(data = ci_one,
                       ggplot2::aes(x = Date, y = n),
                       linetype = "solid") +
    ggplot2::theme(axis.title.y=ggplot2::element_blank(),
        legend.position = "right") +
        ggplot2::scale_linetype_manual(name = "Number of \nMonthly Events",
                        labels = c("Observed", "Estimated"),
                        values = c("solid", "dashed")) +
    #ggplot2::ylab("Number of Monthly Events") +
    ggplot2::ggtitle(plot_title)
    ggplot2::scale_x_date(breaks = seq(
      from = lubridate::year(model$ref_date),
      to = lubridate::year(max)))
    labels = lubridate::year(model$ref_date):lubridate::year(max)
}

