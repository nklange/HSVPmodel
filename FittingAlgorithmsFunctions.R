
# Genetic Algorithm ----------------------------------------------------------------------------------

# Bounds 

lowerRNminus <- rep(.Machine$double.eps,3)
upperRNminus <- c(100,5,200)

lowerRNplus <- rep(.Machine$double.eps,4)
upperRNplus <- c(400,5,200,200) # kappa_r drives up J1

lower <- list(lowerRNminus,lowerRNminus,lowerRNplus,lowerRNplus)
upper <- list(upperRNminus,upperRNminus,upperRNplus,upperRNplus)

# Function

fit_one_vp_ga <- function(data, rep, objective, model, ll_fun, lower, upper, ..., parallel = FALSE) {
  dp <- prep_data_index(data)
  
  out_list <- vector("list", rep)
  
  for (i in seq_len(rep)) {
    tic <- Sys.time()
    tmp <- tryCatch(ga(type = "real-valued", fitness = objective, 
                       ll_fun = ll_fun, 
                       error_list = dp$datalist, set_sizes = dp$set_sizes,
                       model = model,
                       ...,
                       #monitor = FALSE,
                       lower = lower, upper = upper, 
                       run = 5, parallel = parallel), 
                    error = function(e) NA)
    if (inherits(tmp, "ga")) {
      pars <- as_tibble(tmp@solution)
      colnames(pars) <- names(get_start_vp(model))
      pars$objective <- tmp@fitnessValue
      pars$iter <- tmp@iter
      pars$type <- tmp@type
      pars$model <- model
      pars$time <- Sys.time() - tic
      pars$id <- dp$id
      pars$exp <- dp$exp
      pars$leftout <- dp$leftout
      pars$rep <- i
      pars$cvid <- dp$cvid
      out_list[[i]] <- pars
    }
  }
  bind_rows(out_list)
}

# Numerical Integration NLMINB ----------------------------------------------------------

fit_one_vp_nlminb <- function(model, data, startpar, rep, objective, ll_fun, lower, upper, type, ...) {
  
  dp <- prep_data(data)
  out_list <- vector("list", rep)

  #startpar <- data.frame(do.call(rbind, startpar))
  startpar <- startpar %>% filter(cvid == dp$cvid) %>% select(names(get_start_vp(model)))
  startpar <- as.data.frame(startpar)
  
  for (i in seq_len(rep)) {
    
    tic <- Sys.time()
  
    start <- as.vector(t(startpar[i,]))

    tmp <- tryCatch(nlminb(start, objective = objective, 
                           error_list = dp$datalist, set_sizes = dp$set_sizes, 
                           ll_fun = ll_fun, model = model,..., 
                           type = type,
                           lower = .Machine$double.eps, upper = Inf,
                           control = list(eval.max = 300, iter.max = 300, trace = 1)), error = function(e) NA)
    if (is.list(tmp)) {
      
      out_list[[i]] <- bind_cols(
        spread(enframe(tmp$par), name, value) %>%
          setNames(sort(names(get_start_vp(model))))
        , as_tibble(tmp[c(2, 3, 6)])
        , tibble(
          time = Sys.time() - tic
        )
        , tibble (
          model = model
        )
        , tibble (
          id = dp$id
        )
        ,tibble (
           leftout = dp$leftout
         )
        , tibble (
          rep = i
        )
        ,tibble (
          exp = dp$exp
        )
         , tibble (
          cvid = dp$cvid
         )
        , start %>% #spread(enframe(start, "start"), start, value) %>%
          setNames(paste0("s_",names(get_start_vp(model))))
      )
    } 
  }
  bind_rows(out_list)
}
