# Peter Boyd
# nphawkes functions


# Rcpp::compileAttributes() #when changing c++ code
# devtools::document() #to update whole thing
# devtools::load_all() #to temp install
# library(nphawkes) #to use

# Nonparametric Hawkes ----------------------------------------------------

#' Model Independent Stochastic Declustering
#'
#' This function uses nonparametric procedures to analyze a Hawkes
#' process in temporal or spatio-temporal domain, with or without marks
#' through the Model independent stochastic declustering algorithm.
#'
#' This function can only be applied to data that contains a temporal feature. It can also be applied to data
#' consisting of time and space, time and marks, or time and space and marks.
#'
#' For each triggering component used (time, space, marks), a binning structure will be applied. The user may define
#' these right continuous bins as a vector (\code{c(1,5,10)} creates two bins for \eqn{1 < x \le 5}, and \eqn{5 < x \le 10}),
#' or n time bins may be generated automatically by specifying \code{time_quantile = n}, with the same applying for marks and space.
#' This method will establish breaks such that the allocation of time differences within the data will be roughly
#' equal in each created bin. Uniform binning methods may unattainable for discrete marks containing many replicates of values.
#'
#' If no time of day is provided, events will be randomly assigned a time during the event's date.
#'
#' @param dates a vector of dates as "yyyy-mm-dd".
#' @param lat a vector of latitudes, omit if not using spatial data
#' @param lon a vector of longitudes, omit if not using spatial data
#' @param marks a vecotr of marks, or magnitudes, omit if not using marked data
#' @param time_breaks a vector of cutoff values for temporal bins of time differences
#' @param space_breaks a vector of cutoff values for spatial bins of distance differences
#' @param mark_breaks a vector of cutoff values for magnitude bins
#' @param time_quantile NA by default to use defined temporal bins, integer \code{n} value will
#' establish n time bins containing roughly equal number of pairwise time differences
#' @param space_quantile NA by default to use defined spatial bins, integer \code{n} value will
#' establish n spatial bins containing roughly equal number of pairwise spatial differences
#' @param mark_quantile NA by default to use defined magnitude bins, integer \code{n} value will
#' establish n mark bins containing roughly equal number of events
#' @param ref_date a date to serve as time 0, defaults to earliest observation
#' @param time_of_day character string that lists the time of day of events, as hour:minute:second
#' @param just_time TRUE or FALSE object. TRUE if \code{dates} object is a vector of only times,
#' indicating time elapsed since beginning of catalog.
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
#' @return \code{perc_br} is the proportion of events in which the maximum probabilistic assignment
#' is as a background event
#' @return \code{time_bins} is a matrix containing the temporal bin of each pair of events
#' @return \code{dist_bins} is a matrix containing the spatial bin of each pair of events
#' @return \code{mark_bins} is a vector containing the magnitude bin of each event
#' @return \code{n_iterations} is the number of iterations executed until convergence
#' @return \code{input} is a list of all inputs
#'
#' @examples
#' data("hm.csv")
#' out = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_breaks = c(0,0.5, 1, 10, 25, 100),
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = T)
#'
#'  out1 = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_quantile = 7,
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = T,
#'    nonstat_br = T)


#' @export
misd <- function(dates, ref_date = min(dates),
                lat = rep(0, length(dates)),
                lon = rep(0, length(dates)),
                marks = rep(0, length(dates)),
                time_breaks = c(0,1), space_breaks = c(0,1),
                mark_breaks = c(0,1), time_quantile = NA,
                mark_quantile = NA, space_quantile = NA,
                stopwhen = 1e-3, time_of_day = NA,
                just_times = FALSE, nonstat_br = FALSE,
                lon_lim = c(min(lon), max(lon), (max(lon) - min(lon))/10),
                lat_lim = c(min(lat), max(lat), (max(lat) - min(lat))/10),
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

  # sort distances. easier computation of number of
  # points for nonstationary background rate
  dist_mat2 = matrix(data = NA, nrow = n, ncol = n)
  for (i in 1:n){
    dist_mat2[i,] = sort(dist_mat[i,])
  }

  if(sum(lat) != 0) {
    dist_mat = ifelse(dist_mat == 0, 0.0001, dist_mat)
  }

  # THREE: nonstationary background rate
  diag(dist_mat) = 0

  # create spatial grid to contain events
  x_grid = seq(lon_lim[1], lon_lim[2], lon_lim[3])
  y_grid = seq(lat_lim[1], lat_lim[2], lat_lim[3])
  # grid is entire window, pix is midpoints of grid
  x_pix = seq(lon_lim[1] + lon_lim[3]*0.5,
              lon_lim[2] - lon_lim[3]*0.5,
              lon_lim[3])
  y_pix = seq(lat_lim[1] + lat_lim[3]*0.5,
              lat_lim[2] - lat_lim[3]*0.5,
              lat_lim[3])

  # calculate radii such that np number of events are within
  # distance d_i
  di = calc_d(dist_mat2, np = floor(0.05*n))# was doing 24

  # put all events into a pixel
  # select proper br in prob calcs
  pix = get_pix(x_grid, y_grid, lat, lon, x_pix, y_pix)
  pix = data.frame(pix)
  names(pix) = c("lon", "lat")

  # implement option of bin by quantile
  if(is.na(time_quantile) == FALSE){
    time_mat1 = time_mat
    time_mat1[lower.tri(time_mat1, diag = TRUE)] = NA
    time_breaks = as.vector(quantile(time_mat1, na.rm = TRUE,
                                probs = seq(0,1, length.out= time_quantile + 1)))
    time_breaks[length(time_breaks)] = time_breaks[length(time_breaks)] + 1
    time_breaks[1] = 0
  } else {time_breaks = time_breaks}

  if(is.na(mark_quantile) == FALSE){
    mark_breaks = unique(as.vector(quantile(marks,
                                       probs = seq(0, 1, length.out = mark_quantile + 1))))
    mark_breaks[length(mark_breaks)] = mark_breaks[length(mark_breaks)] + 1
    mark_breaks[1] = 0 # accounts for simulating points below lowest bin
  } else {mark_breaks = mark_breaks}

  if(is.na(space_quantile) == FALSE){
    dist_mat1 = dist_mat
    dist_mat1[lower.tri(dist_mat1, diag = TRUE)] = NA
    space_breaks = as.vector(quantile(dist_mat1, na.rm = TRUE,
                                probs = seq(0,1, length.out= space_quantile + 1)))
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
  n_iterations = 0

  while( max_diff > stopwhen){
    if(nonstat_br == TRUE) {
      tau = calc_tau(x_pix, y_pix, p0, di,
                     lon, lat, times)
      z = sum(tau[,1] * diff(x_pix)[1] * diff(y_pix)[1] * ceiling(max(times)))
      br = calc_br_nonstat(p0, times, z, tau)

      br_grid = data.frame(br)
      names(br_grid) = c("br", "lon", "lat")

      locs = dplyr::left_join(pix, br_grid, by = c("lon", "lat"))
      locs$br = ifelse(is.na(locs$br) == TRUE, 0, locs$br)
      br = as.vector(locs$br)
    } else {
      br = calc_br(p0, times)
      br = rep(br, n)
      br_grid = NA
    }
    g_vals = get_g(p0, time_breaks, time_mat)
    g = g_vals[,1] / g_vals[,2]
    h = get_h(p0, space_breaks, dist_mat)
    k = get_k(p0, marks, mark_breaks)
    p = update_p(p0, time_mat, dist_mat, mark_mat,
                 g, h, k,
                 space_breaks, time_breaks, mark_breaks,
                 br, time_bins, dist_bins, lat)
    # rounding difference may cause insignificant negatives
    # replace trace neagtives by 0
    p = ifelse(p<0, 0, p)
    max_diff = check_p(p0, p)
    p0 = p
    n_iterations = n_iterations + 1
  }

  max_diag = c()
  for ( i in 1:nrow(p0)){
    if (p0[i,i] == max(p0[i,])){
      max_diag[i] = 1
    } else{
      max_diag[i] = 0
    }
  }
  df$mainshock = max_diag

  max_event = c()
  for (i in 1:nrow(p0)) {
    max_event[i] = which.max(p0[i,])
  }
  df$parent = max_event

  df$dates = lubridate::ymd(ref_date) +
    lubridate::days(floor(df$times))

  perc_br = sum(max_diag) / length(max_diag)
  perc_diag = sum(diag(p0)) / nrow(p0)

  out = list(p0 = p0, g= g, h = h, k = k, br = br,  br_grid = br_grid,
             time_bins = time_bins, mark_bins = mark_mat,
             dist_bins = dist_bins, perc_br = perc_br, perc_diag = perc_diag,
             time_breaks = time_breaks, mark_breaks = mark_breaks, space_breaks = space_breaks, data = df,
             ref_date = ref_date, n_iterations = n_iterations,
             input = mget(names(formals()),sys.frame(sys.nframe())))
  return(out)
}



# Conditional Intensity ---------------------------------------------------

#' Conditional Intensity
#'
#' This function estimates the conditional intensity function of the observed process.
#'
#' This function is to be used in conjunction with the \code{misd()} function from the \code{nphawkes} library.
#' Using the output from the \code{misd()} function. The exported data frame will contain, for each event,
#' the time elapsed (in days) since the beginning of the observation window, the date of the event,
#' the conditional intensity at the event's location in time (or space-time), as well as the
#' coordinates and marks (if provided).
#'
#' @param model the output from \code{misd()}
#'
#' @return a data frame containing the time, location, marks, and estimated conditional intensity
#'
#' @examples
#' data("hm.csv")
#' out = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_breaks = c(0,0.5, 1, 10, 25, 100),
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = T)
#'
#' ci = cond_int(out)
#'
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
    cond_int[i] = br[i] + sum(trig[i, ])
  }

  ci = data.frame(times, lat, lon, marks, cond_int)
  ci$dates = lubridate::ymd(model$ref_date) +
    lubridate::days(floor(ci$times))
  return(ci)
}

# Super-thinning ----------------------------------------------------------

#' Super-thinning
#'
#' This function performs super-thinning to assess model fit. Thinning only can be executed by
#' simply considering events of type "thin" or "retain", while superpositioning only can be executed
#' by simply treating "thinned" points as "retained".
#'
#' This function is to be used in conjunction with the \code{misd()} function from the \code{nphawkes} library.
#'
#' To simulate spatial data, the user may define the \code{map} and \code{region} to easily simulate points within
#' set political borders. This feature is not quite implementable, but may be soon.
#' Otherwise, simulated points may be established when \code{sim_grid = TRUE}, based on minimum and maximum observed
#' latitude and longitude coordinates.
#'
#'
#' @param K a constant or character string that governs the amount of thinning and
#' superposing that is implemented. Can be a constant value, the median, mean, minimum, or maximum
#'  conditional intensity: "median_ci", "mean_ci", "min_ci", or "max_ci", respectively.
#' @param model the output from \code{misd()}
#' @param method character string that defines residual analysis method as "superthin", "thin", or "superpose"
#' @param map name of map provided by the maps package, defaults to maps::world
#' @param region name of subregion to include, defaults to entirety of map
#' @param sim_grid TRUE if geographical coordinates do not necessarily pertain to
#' certain area, such as 0x1 by 0x1 grid, FALSE if using specific region or not using spatial information
#' @param lat_lim vector containing minimum and maximum latitude bounds if sim_grid is TRUE
#' @param lon_lim vector containing minimum and maximum longitude bounds if sim_gird is TRUE
#'
#' @return a data frame which includes the time, location, and estimated conditional intensity of events.
#' The type of event, either observed or simulated, is noted along with the probability that the event was kept
#' and whether or not the point was in fact retained.
#'
#' @examples
#' data("hm.csv")
#' out = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_breaks = c(0,0.5, 1, 10, 25, 100),
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = T)
#' st = super_thin(K = "max_ci",
#'     model = out,
#'     method = "superthin",
#'     )
#'
#' @export

super_thin = function(K = "median_ci",
                      model, method = "superthin",
                      map = world, region = ".",
                      sim_grid = TRUE,
                      lat_lim = model$input$lat_lim,
                      lon_lim = model$input$lon_lim) {

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
  if (sum(model$data$lat) == 0) {
    sim_pp = cbind(sim_pp, lat = rep(0, n2), lon = rep(0, n2))
  } else{
    if (sim_grid == TRUE) {
      sim_lat = runif(n2, min = lat_lim[1], max = lat_lim[2])
      sim_lon = runif(n2, min = lon_lim[1], max = lon_lim[2])
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
  if(model$input$nonstat_br == TRUE) {
    x_grid = seq(model$input$lon_lim[1], model$input$lon_lim[2], model$input$lon_lim[3])
    y_grid = seq(model$input$lat_lim[1], model$input$lat_lim[2], model$input$lat_lim[3])
    # grid is entire window, pix is midpoints of grid
    x_pix = seq(model$input$lon_lim[1] + model$input$lon_lim[3]*0.5,
                model$input$lon_lim[2] - model$input$lon_lim[3]*0.5,
                model$input$lon_lim[3])
    y_pix = seq(model$input$lat_lim[1] + model$input$lat_lim[3]*0.5,
                model$input$lat_lim[2] - model$input$lat_lim[3]*0.5,
                model$input$lat_lim[3])

    pix = get_pix(x_grid, y_grid, sim_pp$lat, sim_pp$lon, x_pix, y_pix)
    pix = data.frame(pix)
    names(pix) = c("lon", "lat")

    locs = dplyr::left_join(pix, model$br_grid, by = c("lon", "lat"))
    locs$br = ifelse(is.na(locs$br) == TRUE, 0, locs$br)
    br = as.vector(locs$br)
  } else {
    br = rep(model$br[1], n2)
  }

  if(sum(sim_pp$lat) == 0){
    sim_pp$lat = rep(1,n2)
  } else {
    sim_pp$lat = sim_pp$lat
  }

  data_pp$marks = ifelse(data_pp$marks == mark_breaks[1],
                         data_pp$marks + 0.0001, data_pp$marks)

  for (i in 1:n2) {
    sim_trig = 0
    dist_vec = rep(NA, n2)
    R = 6371e3

    for (j in 1:(sim_pp$row1[i])) {
      phi1 = sim_pp$lat[i] * pi / 180
      phi2 = data_pp$lat[j] * pi / 180
      latdiff = (data_pp$lat[j] - sim_pp$lat[i]) * pi / 180
      londiff = (data_pp$lon[j] - sim_pp$lon[i]) * pi / 180

      a = sin(latdiff / 2) * sin(latdiff / 2) +
        cos(phi1) * cos(phi2) * sin(londiff / 2) * sin(londiff / 2)
      b = 2 * atan2(sqrt(a), sqrt(1 - a))
      d = R * b

      dist_vec[j] = d
    }

    # convert distances to correct units
    if (model$input$dist_unit == "mile") {
      dist_vec = dist_vec*0.000621371
    } else if (model$input$dist_unit == "kilometer") {
      dist_vec = dist_vec*0.001
    } else {
      dist_vec = dist_vec
    }

    # get triggering component for each
    for (j in 1:(sim_pp$row1[i])){
      td = sim_pp$Time[i] - data_pp$Time[j]
      td = ifelse(td == 0, 0.00001, td)
      gb = bin_f(td, time_breaks)
      gg = model$g[gb]

      kb = bin_f(data_pp$marks[j], mark_breaks)
      kk = model$k[kb]

      hd = dist_vec[j]
      hd = ifelse(hd == 0, 0.00001, hd)
      hb = bin_f(hd, space_breaks)
      hh = model$h[hb]

      ghk = gg*kk*hh
      sim_trig = sim_trig + ghk
    }

    # get ci for each then remove the simulated point and correct row numbers
    cond_int_sim[i] = br[i] + sim_trig
  }

  # FOURTH: THIN SIMULATED DATA

  ci_sim = data.frame(times = sim_pp$Time,
                      lat = sim_pp$lat, lon = sim_pp$lon,
                   cond_int = cond_int_sim, p = NA,
                   keep = NA, type = NA, marks = sim_pp$marks)
  # simulated data
  ci_sim$type = "sim"
  ci_sim$dates = lubridate::ymd(model$ref_date) +
    lubridate::days(floor(ci_sim$times))
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
#' This function is to be used in conjunction with the \code{misd()} function from the \code{nphawkes} library,
#' and is used within the function that produces triggering plots.
#'
#' @param model the output from \code{misd()}
#'
#' @return a list of variance estimates for each bin of all utilized triggering components
#'
#' @examples
#' data("hm.csv")
#' out = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_breaks = c(0,0.5, 1, 10, 25, 100),
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = T)
#' se = se_bars(out)
#'
#' @export
#'
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
#' estimated value of triggering functions for each bin along with \eqn{+/- 2} standard errors bars
#'
#' Depending on bin size, x-axis labels may overlap, impairing readability. Axis tick marks may be changed for plots
#' using \code{scale_x_continuous()} within the \code{ggplot2} library.
#'
#' @param model the output from \code{misd()}
#' @param time_xlim vector of minimum and maximum x-axis value shown in the temporal plot
#' @param space_xlim vector of minimum and maximum x-axis value shown in the spatial plot
#' @param mark_xlim vector of minimum and maximum x-axis value shown in the magnitude plot
#' @param time_ylim vector of minimum and maximum y-axis value shown in the temporal plot
#' @param space_ylim vector of minimum and maximum y-axis value shown in the spatial plot
#' @param mark_ylim vector of minimum and maximum y-axis value shown in the magnitude plot
#' @param mag_label character string representing what the magnitude measures
#' @param se_include TRUE to include standard error estimates for bins of triggering components, FALSE
#' to exclude them from the visualizations
#'
#' @return \code{all_plots} is a cowplot object of histogram estimators for all utilized triggering components
#' @return \code{time_plot} is a ggplot object of the histogram estimator for the temporal effect
#' @return \code{space_plot} is a ggplot object of the histogram estimator for the spatial effect,
#' NULL if data provided is not a spatial process
#' @return \code{mark_plot} is a ggplot object of the histogram estimator for the mark effect,
#' NULL if data provided is not a marked process
#'
#' @examples
#' data("hm.csv")
#' out = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_breaks = c(0,0.5, 1, 10, 25, 100),
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = T)
#' tp = trig_plots(out,
#'     time_xlim = c(0, 2),
#'     space_xlim = c(0, 10),
#'     mark_xlim = c(3.9, 5.1))
#'
#' @export
trig_plots = function(model,
                      time_xlim = c(min(model$time_breaks), max(model$time_breaks)),
                      space_xlim = c(min(model$space_breaks), max(model$space_breaks)),
                      mark_xlim = c(min(model$mark_breaks), max(model$mark_breaks)),
                      time_ylim = c(0,NA),
                      space_ylim = c(0,NA),
                      mark_ylim = c(0,NA),
                      mag_label = "magnitude",
                      se_include = TRUE){

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

  time_ylim[2] = ifelse(is.na(time_ylim[2]) == TRUE,
                     max(time_df$g + time_df$se), time_ylim[2])
  space_ylim[2] = ifelse(is.na(space_ylim[2]) == TRUE,
                     max(dist_df$h + dist_df$se), space_ylim[2])
  mark_ylim[2] = ifelse(is.na(mark_ylim[2]) == TRUE,
                     max(mag_df$k + mag_df$se), mark_ylim[2])


  # Time Plot
  trig_g = ggplot2::ggplot(time_df, ggplot2::aes(time, g)) +
    ggplot2::coord_cartesian(xlim = time_xlim, ylim = time_ylim) +
    ggplot2::xlab(paste("t (time in ", model$input$time_unit, "s)", sep = "")) +
    ggplot2::ylab("g(t)")

  if(se_include == TRUE) {
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
    ggplot2::coord_cartesian(xlim = mark_xlim, ylim = mark_ylim) +
    ggplot2::xlab(paste("m (", mag_label, ")", sep = "")) +
    ggplot2::ylab("k(m)")

  if(se_include == TRUE){
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
    ggplot2::coord_cartesian(xlim = space_xlim, ylim = space_ylim) +
    ggplot2::xlab(paste("s (distance in ", model$input$dist_unit, "s)", sep = "")) +
    ggplot2::ylab("h(s)")

  if(se_include == TRUE){
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

  if (sum(model$mark_bins) == 0 & sum(model$dist_bins) == 0){
    all_plots = cowplot::plot_grid(trig_g, ncol = 1)
  } else if (sum(model$mark_bins) != 0 & sum(model$dist_bins) == 0) {
    all_plots = cowplot::plot_grid(trig_g, trig_k, ncol = 2)
  } else if (sum(model$mark_bins) == 0 & sum(model$dist_bins) != 0) {
    all_plots = cowplot::plot_grid(trig_g, trig_h, ncol = 2)
  } else {
    all_plots = cowplot::plot_grid(trig_g, trig_h, trig_k, ncol = 3)
  }

  out = list(all_plots = all_plots, time_plot = trig_g,
         space_plot = trig_h, mark_plot = trig_k)
  return(out)
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
#'
#' @return a tiered plot with superposed points on the top tier, points that were retained, not thinned,
#' on the middle tier, and points that were thinned on the bottom tier.
#'
#' @examples
#' data("hm.csv")
#' out = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_breaks = c(0,0.5, 1, 10, 25, 100),
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = T)
#' st = super_thin(out)
#' sp = st_plot(st, "superthin")
#' @export
st_plot = function(superthin, method = "superthin",
                   time_label = "Year"){

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
    p = ggplot2::ggplot(data = superthin, ggplot2::aes(x = dates, col = type)) +
    ggplot2::geom_segment(ggplot2::aes(x = dates, y = y, yend = y + 1, xend = dates),
                 show.legend = FALSE) +
    ggplot2::theme(axis.title.y= ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank()) +
    ggplot2::scale_color_manual(breaks = c("sim", "retain", "thin"),
                       values=c("grey3", "grey40", "grey70")) +
    ggplot2::xlab(time_label) +
    ggplot2::scale_y_continuous(breaks = c(0,1,2) + 0.5,
                       labels = c("Thinned", "Retained", "Simulated"))

  } else if (method == "thin"){
    p = superthin %>% dplyr::filter(type != "sim") %>%
    ggplot2::ggplot(ggplot2::aes(x = dates, col = type)) +
      ggplot2::geom_segment(ggplot2::aes(x = dates, y = y, yend = y + 1, xend = dates),
                            show.legend = FALSE) +
      ggplot2::theme(axis.title.y= ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank()) +
      ggplot2::scale_color_manual(breaks = c("retain", "thin"),
                                  values=c("grey40", "grey70")) +
      ggplot2::xlab(time_label) +
      ggplot2::scale_y_continuous(breaks = c(0,1,2) + 0.5,
                                  labels = c("Thinned", "Retained"))

  } else {
    p = superthin %>% dplyr::mutate(type = replace(type, type == "thin", "retain")) %>%
      dplyr::mutate(replace(type, type == "retain", "observed")) %>%
      ggplot2::ggplot(data = superthin, ggplot2::aes(x = dates, col = type)) +
      ggplot2::geom_segment(ggplot2::aes(x = dates, y = y, yend = y + 1, xend = dates),
                            show.legend = FALSE) +
      ggplot2::theme(axis.title.y= ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank()) +
      ggplot2::scale_color_manual(breaks = c("sim", "observed"),
                                  values=c("grey3", "grey40")) +
      ggplot2::xlab(time_label) +
      ggplot2::scale_y_continuous(breaks = c(0,1,2) + 0.5,
                                  labels = c("Observed", "Simulated"))
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
#'
#' @return a histogram of the super-thinned process.
#'
#' @examples
#' data("hm.csv")
#' out = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_breaks = c(0,0.5, 1, 10, 25, 100),
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = T)
#' st = super_thin(out, K = "max_ci")
#' ci_h = ci_hist(st, nbins = 40)

#' @export
ci_hist = function(superthin, nbins = 30){
  st = superthin[which(superthin$type != "thin"),]
  out = ggplot2::ggplot(data = st) +
    ggplot2::geom_histogram(ggplot2::aes(x = dates), bins = nbins) +
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
#' @param model the output from \code{misd()}
#' @param superthin the output from \code{super_thin()}
#' @param min_date the minimum date to be shown on the plot
#' @param max_date the maximum date to be shown on the plot
#'
#' @return a plot that shows the estimated conditional intensity as a black line with the number of monthly events as
#' gray vertical lines
#'
#' @examples
#' data("hm.csv")
#' out = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_breaks = c(0,0.5, 1, 10, 25, 100),
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = T)
#' st = super_thin(out)
#' ci_p = ci_plot(out, st)
#'
#' @export
ci_plot = function(model, min_date = model$input$ref_date,
                   max_date = model$data$dates[nrow(model$data)],
                   superthin){

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
    ggplot2::ggtitle(plot_title)
    ggplot2::scale_x_date(breaks = seq(
      from = lubridate::year(model$ref_date),
      to = lubridate::year(max)))
    labels = lubridate::year(model$ref_date):lubridate::year(max)
}


# Background Rate Plot ----------------------------------------------------


#' Plot of nonstationary background rate over space
#'
#' This function exports a heat map of the background rate of a nonstationary background rate over space, and is paired
#' with the \code{misd()} function when a nonstationary background rate is specified.
#'
#' @param model the output from \code{misd()}
#' @param vals_include if TRUE, the background rate values will be printed on the map
#'
#' @return a ggplot of the background rates as a heat map
#'
#' @examples
#' data("hm.csv")
#' out = misd(dates = hm$t,
#'    ref_date = "1999-10-16",
#'    lat = hm$lat,
#'    lon = hm$lon,
#'    marks = hm$m,
#'    time_breaks = c(0,0.1, 0.5, 1,7,93,600),
#'    space_breaks = c(0,0.5, 1, 10, 25, 100),
#'    mark_breaks = c(3, 3.1,3.3, 4, 5, 8),
#'    just_times = TRUE,
#'    nonstat_br = TRUE)
#' br_p = br_plot(out, vals_include = TRUE)
#'
#' @export
br_plot = function(model, vals_include = FALSE) {
  br_grid = model$br_grid
  ggplot2::ggplot(data = br_grid, aes(x = lon, y = lat, fill = br)) +
    ggplot2::scale_fill_gradient(low = "white", high = "black") +
    ggplot2::geom_tile() +
    ggplot2::scale_x_continuous(minor_breaks = NULL) +
    ggplot2::scale_y_continuous(minor_breaks = NULL)
}
