# some data transforms -------------------------------------------------------------

circ_null <- circular(0)
base_radians <- circular((0:180) * 2 * pi/360)



# nlminb data: fitting radians -------------------------------------------------
prep_data <- function(data, replace_zero = circular(.Machine$double.xmin)) {
  # if (!missing(replace_zero)) {
  #    data$error_0 <- ifelse(data$error_0 == 0, replace_zero, data$error_0)
  # }
  
  dl <- split(circular(data$error_0), f = data$set_size)
  id <- unique(data$id)
  cvid <- unique(data$cvid)
  set_sizes <- sort(unique(data$set_size))
  exp <- unique(data$exp)
  leftout <- unique(data$leftout)
  stopifnot(names(dl) == as.character(set_sizes))
  return(list(
    datalist = dl,
    set_sizes = set_sizes,
    id = id,
    cvid = cvid,
    exp = exp,
    leftout = leftout
  ))
}

# GA data: fitting degrees -------------------------------------------------
prep_data_index <- function(data) {
  
  dl <- split((abs(data$deg_error_0) + 1), f = data$set_size)
  set_sizes <- sort(unique(data$set_size))
  id <- unique(data$id)
  cvid <- unique(data$cvid)
  exp <- unique(data$exp)
  leftout <- unique(data$leftout)
  stopifnot(names(dl) == as.character(set_sizes))
  return(list(
    datalist = dl,
    set_sizes = set_sizes,
    id = id,
    cvid = cvid,
    exp = exp,
    leftout = leftout
  ))
}


# Get start parameters ---------------------------------------------

get_start_vp <- function(model) {
  
  
  
  #if (is.null(res_ga)) {
  if (model == "MK_RNminus") {
    start <- c(
      mkappa1 = runif(1, 30, 60),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 10, 40)
    )
  } else if (model == "MK_RNplus") {
    start <- c(
      mkappa1 = runif(1, 100, 200),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 50, 130),
      kappa_r = runif(1, 10, 70)
    )
  } else if (model == "J_RNminus") {
    start <- c(
      J1bar = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 10, 40)
    )
  } else if (model == "J_RNplus") {
    start <- c(
      J1bar = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 10, 40),
      kappa_r = runif(1, 30, 60)
    )
    #  }
  }
  return(start)
}
