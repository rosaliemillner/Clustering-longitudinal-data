#functions for MAP criterion, BIC calculation, computing optimal number of groups


MAP_criterion = function(Z){
  MAP=matrix(0,nrow(Z),ncol(Z))
  max_indices = max.col(Z)
  MAP[cbind(1:nrow(Z), max_indices)] = 1
  return(MAP)
} #indicates which sequence belongs to which cluster


BIC_and_params = function(G,data_sequences,E,itermax,tolerance,print=T){
  G_list = list()
  if(print){cat("Calculating BIC for G =",G, '\n')}
  n = length(data_sequences)
  EM_results = EM_algo(data_sequences,G,E,itermax,tolerance,print)
  Z=EM_results$Z
  logliklhd = EM_results$loglikelihood
  params = EM_results$params
  time = EM_results$execution_time_in_minutes
  q = G*E*(E+1)-1
  BIC_val = q*log(n) - 2*logliklhd
  return(list(BIC_val=BIC_val, Z=Z, params=params, execution_times=time))
} #performs the EM algo for a certain number G of mixture components and returns the BIC value of the model, as well as other parameters


model_results = function(data_sequences, C, E, itermax=50, print = T, tolerance=(10^(-5))){
  start_time=proc.time()
  BIC_list = lapply(C, BIC_and_params, data_sequences, E, itermax, tolerance, print)
  BIC_vect = sapply(BIC_list, function(x) x$BIC_val)
  #clusters = sapply(BIC_list, function(x) MAP_criterion(x$Z))
  Z = lapply(BIC_list, function(x) x$Z)
  groups = lapply(Z, function(M){ if(is.null(dim(M))){M} else{MAP_criterion(M)} })
  params_list = lapply(BIC_list, function(x) x$params)
  times = sapply(BIC_list, function(x) x$execution_times)
  min_G = which.min(BIC_vect)
  min_BIC = min(BIC_vect)
  if(print){
    cat("Minimal BIC for G =",min_G,"mixture components. BIC =", min_BIC ,'\n')
    cat("BIC list :", BIC_vect,'\n')
  }
  end_time=proc.time()
  time_minutes = (end_time["elapsed"]-start_time["elapsed"]) /60
  return(list(optimal_G=min_G, min_BIC=min_BIC, BIC_vect=BIC_vect, groups = groups, params = params_list, execution_times_per_G = times, total_execution_time_in_minutes = round(time_minutes, 2))) #maybe also return the likelihood?
} #runs the EM algo and returns a list with the recap of the outputs (including the semi-Markov parameters characterising each cluster.) as well as the optimal number of clusters using the BIC criterion

