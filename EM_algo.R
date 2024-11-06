#EM algorithm

library(optimx)
library(numDeriv)
library(sgd)
library(foreach)
library(doParallel)
library(abind)

EM_algo = function(data_sequences, G, E, itermax=50, tolerance=10^(-5), print=T){ #G is the "current" number of clusters
  start_time = proc.time()
  n=length(data_sequences)
  difflklh=c()
  epsilon = tolerance
  
  current_params = convert_params_structures(init_EM_algo(data_sequences,G,E),G,E) #initialization of parameters
  
  N = sapply(data_sequences, function(x) length(x$states)) #vector with the number of visited states/transitions per individual
  mat_indicator_start = matrix_indicator(data_sequences,1,E)
  tensor_nb_transitions = tensor_transitions(data_sequences,E)
  
  for (iter in 1:itermax){
    if(print){cat("Iter:",iter,'\n')}
    mat_likelihood = matrix_likelihood(data_sequences, current_params, G)
    
    if(G==1){ 
      Z = rep(1,n)
      
      current_params$alpha = as.matrix(colMeans(mat_indicator_start))
      
      sum_axis1 = apply(tensor_nb_transitions, MARGIN = c(1, 2), sum)
      #sum_axis2 = colSums(matrix(sum_axis1, nrow = E * E, ncol = n))
      sum_axis2 = rowSums(sum_axis1)
      P = matrix(sum_axis1, nrow=E) / matrix(rep(sum_axis2, each=E), nrow=E, ncol=E, byrow = TRUE)
    }
    else {
      Z = t(current_params$pi * t(mat_likelihood)) #mat nxG 
      Z = t(apply(Z, 1, function(x) { 
        x_sum = sum(x, na.rm = TRUE)
        if(x_sum==0 || is.na(x_sum)){return(rep(0,G))} else {return(x / x_sum)}})) #"Z = Z / rowSums(Z)"
      Z[is.nan(Z)]=0
      
      current_params$pi = colMeans(Z)
  
      current_params$alpha = sweep( t(mat_indicator_start)%*%Z , 2, colSums(Z), "/")
      
      numerator = array(0, dim = c(E, E, G))
      numerator = lapply(1:G, function(g) {
        numerator[, , g] = apply(tensor_nb_transitions, c(1, 2), function(s) {
          sum(Z[, g] * s) }) })
      numerator = abind(numerator, along = 3)
      
      denominator = apply(numerator, c(1,3), sum)
      denominator_extended = array(rep(denominator, each = E), dim = c(E, E, G))
      denominator_extended = aperm(denominator_extended, c(2,1,3))
      
      current_params$P = numerator / denominator_extended
      
    }
    num_cores = detectCores() - 1

    cl = makeCluster(num_cores)
    registerDoParallel(cl)
    clusterExport(cl, c("expected_partial_pen_loglikelihood", "gradient_approx_expected_partial_pen_loglikelihood", "matrix_indicator"))
    clusterExport(cl, c("data_sequences", "Pen_lg", "E", "G", "Z", "N"), envir = environment())
    
    clusterEvalQ(cl,{library(numDeriv)})  #library used for grad 
    
    tryCatch({
    results=foreach(u = (1:(G * E)), .combine = 'list') %dopar% {
      l = u %% E + 1
      g = (u - 1) %/% E + 1

      if (is.na(current_params$gamma_params[l,,g][1])) {
        return(NULL) # state not visited at all in the group
      }
      maximisation = optim(
        par = current_params$gamma_params[l,,g],
        fn = expected_partial_pen_loglikelihood,
        gr = gradient_approx_expected_partial_pen_loglikelihood,
        data_sequences = data_sequences,
        Z = Z,
        N = N,
        l = l,
        g = g,
        G = G,
        E = E,
        method = "L-BFGS-B",
        lower = c(10^(-3), 10^(-3)),
        control = list(fnscale = -1)
      )

      list(l = l, g = g, par = maximisation$par)
    } #fin parallelisation
    }, finally = {
      stopCluster(cl)
    })

    for (res in results) {
      if (!is.null(res)) {
        current_params$gamma_params[res$l, , res$g] = res$par
      }
    }
  
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    # for (u in 1:(E*G)){
    #   l = u%%E + 1
    #   g = (u-1)%/%E + 1
    #   cat("g = ",g,'\n')
    #   cat("l = ",l,'\n')
    # 
    #   if (is.na(current_params$gamma_params[l,,g][1])){ next } #state not visited at all in the group
    #   maximisation = optim(
    #     par = current_params$gamma_params[l,,g],
    #     fn = expected_partial_pen_loglikelihood, #first arg of the function should be the params wrt which we optimise
    #     gr = gradient_approx_expected_partial_pen_loglikelihood,
    #     data_sequences = data_sequences,
    #     Z = Z,
    #     N = N,
    #     l=l,
    #     g=g,
    #     G=G,
    #     method = "L-BFGS-B", #it is BFGS but with bounds  #or test descente de gradient avec approx du gradient --> SGD
    #     #lower=c(.Machine$double.eps,.Machine$double.eps),
    #     lower=c(10^(-3), 10^(-3)),
    #     control = list(fnscale = -1) #to maximise instead of minimise
    #   )
    #   current_params$gamma_params[l,,g] = maximisation$par
    #   print(maximisation$par)
    # }
    
    # cat("current_params$pi = ",current_params$pi,'\n')
    # cat("mat_likelihood = ",mat_likelihood,'\n')
    # cat("Z = ",Z,'\n')
    
    current_loglklhd = loglikelihood(current_params,mat_likelihood, Z, G) #different whether G=1 or G>1
    if (iter > 1){
      difflklh=c(difflklh,abs(current_loglklhd - old_loglklhd))
    }
    #cat("current liklhd = ",current_loglklhd,'\n')
    if ((iter>1) && (abs(current_loglklhd - old_loglklhd) < epsilon)){
      if(print){cat("Algo converged after", iter, "iterations.",'\n')}
      break
    }
    old_loglklhd = current_loglklhd
  }
  end_time = proc.time()
  time_minutes = (end_time["elapsed"]-start_time["elapsed"]) /60
  return(list(Z = Z, loglikelihood = current_loglklhd, params = current_params, diff_likelihood=difflklh, execution_time_in_minutes = round(time_minutes, 2)))
}
