
ll_vp_full <- function(pars, model, error_list, set_sizes, ll_fun, type,...) {

  precision <- pars[1]/(set_sizes^pars[2])
  parscont <- c(pars[3])
  
  if (grepl("RNplus",model)){
    parscont <- c(parscont,pars[4])
  }
  
  out <- vector("numeric", length(set_sizes))
  for (i in seq_along(error_list)) {
   i <- 1
    out[[i]] <- ll_fun(pars = c(precision[i],parscont), 
                       errors = error_list[[i]],
                       model = model,
                       type = type)
  }
  
  return(sum(out))
  
}

ll_vp3 <- function(pars, errors, model, type) {
  
  if (type == "r") {
    coreFunction <- paste0("int_fun_",model)
  } else {
    coreFunction <- paste0("cint_fun_",model)
  }
  
  min_integral <- 0
  max_integral <- Inf
  out <- vector("numeric", length(errors))
  for (i in seq_along(out)) {
    out[i] <- tryCatch(
      integrate(match.fun(coreFunction), 
                min_integral, max_integral, 
                pars, 
                radian = errors[i], stop.on.error = FALSE)$value, 
      error = function(e) tryCatch(
        integrate(match.fun(coreFunction), 
                  min_integral, max_integral, 
                  pars,
                  radian = if (errors[i] == 0) {
                    circular(.Machine$double.xmin)  
                  } else errors[i], stop.on.error = FALSE)$value), 
      error = function(e) tryCatch(
        integrate(match.fun(coreFunction), 
                  min_integral, max_integral,   
                  pars,
                  radian = if (errors[i] == 0) {
                    circular(.Machine$double.eps^2) 
                  } else errors[i], stop.on.error = FALSE)$value), error = function(e) NA)
  }

  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
}
