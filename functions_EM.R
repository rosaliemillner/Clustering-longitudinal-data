#functions used in the EM_algo function

likelihood_ig = function(seq, params, g) {
  state_sequence = seq$states
  sojourn_times = seq$sojourn_times
  
  start_state = state_sequence[1]
  l = params$alpha[start_state, g] * dgamma(sojourn_times[1], shape = params$gamma_params[start_state, 1, g], rate = params$gamma_params[start_state, 2, g])
  if(is.na(l)){l=0}
  if (length(state_sequence) == 1) {return(l)}
  
  state_1_vec = state_sequence[-length(state_sequence)]
  state_2_vec = state_sequence[-1]
  
  shape_params = params$gamma_params[state_2_vec, 1, g]
  rate_params = params$gamma_params[state_2_vec, 2, g]
  
  valid_indices = !is.na(shape_params) & !is.na(rate_params) & !is.nan(shape_params) & !is.nan(rate_params)

  valid_shape_params = shape_params[valid_indices]
  valid_rate_params = rate_params[valid_indices]
  valid_sojourn_times = sojourn_times[-1][valid_indices]
  valid_state_1_vec = state_1_vec[valid_indices]
  valid_state_2_vec = state_2_vec[valid_indices]
  
  dgamma_values = dgamma(valid_sojourn_times, shape = valid_shape_params, rate = valid_rate_params)
  
  transition_probs = mapply(function(s1, s2) params$P[s1, s2, g], valid_state_1_vec, valid_state_2_vec)
  
  if(all(valid_indices==FALSE)){l=0}
  else{l = l * prod(transition_probs * dgamma_values)}
  
  return(l)
} #computes the likelihood for one cluster g (the params), for one sequence i


likelihood_ig_vectorized = Vectorize(function(i, g, params,data_sequences) {
  likelihood_ig(data_sequences[[i]], params, g)
}, vectorize.args = c("i", "g"))


matrix_likelihood = function(data_sequences, params, G){
  n=length(data_sequences)
  ind_indices = rep(1:n, each = G)
  groups_indices = rep(1:G, times = n)
  result_vector = likelihood_ig_vectorized(ind_indices, groups_indices, params, data_sequences)
  matrix = matrix(result_vector, nrow = n, ncol = G, byrow = TRUE)
  return(matrix)
} #returns a matrix nxG where each element in (i,g) is the likelihood for indiv i with the parameters in group g


loglikelihood = function(params, mat_likelihood, Z, G){
  n=length(mat_likelihood[,1])
  if (G>1){
    pi_outer = outer(rep(1,n),params$pi)
    L = Z * log(pi_outer * mat_likelihood)
    L[is.nan(L)]=0
    L = sum(L)
  } else {
    L = sum(log(mat_likelihood[mat_likelihood!=0]))
  }
  return(L)
} #computes the whole loglikelihood


matrix_indicator = function(data_sequences, k, E){
  n=length(data_sequences)
  kth_states = sapply(data_sequences, function(seq) if (k>length(seq$states)){return(0)} else {return(seq$states[k])})
  mat = matrix(0, nrow = n, ncol = E)
  mat[cbind(1:n, kth_states)] = 1
  return(mat)
} #returns a nxE matrix where the element in (i,j) is equal to 1 if the k-th visited state in the i-th sequence is j, and 0 otherwise


tensor_transitions = function(data_sequences,E){
  n = length(data_sequences)
  transition_tensor = array(0, dim = c(E, E, n)) #init
  
  transitions_count = lapply(data_sequences, function(seq) {
    seq_states = seq$states
    transitions = which(seq_states[-length(seq_states)] != seq_states[-1])
    next_states = seq_states[transitions + 1]
    transition_table = table(factor(seq_states[transitions], levels = 1:E), factor(next_states, levels = 1:E))
    as.numeric(transition_table)
  })
  
  transition_matrix = do.call(cbind, transitions_count)
  transition_tensor[] = transition_matrix
  transition_tensor = array(transition_tensor, dim = c(E, E, n))
  
  return(transition_tensor)
} #creates a tensor ExExn where element (h,j,i) is the number of transitions going from state h to state j in the i-th sequence


Pen_lg = function(gamma_params,N){
  a_lg = gamma_params[1]
  return(-(a_lg+log(a_lg))/sqrt(sum(N)))
} #computes the penalty for fixed state l and group g for penalized likelihood (here, the argument gamma_params is of length 2 and specific for l and g)


expected_partial_pen_loglikelihood = function(gamma_params, data_sequences, Z, N, l, g, G, E){ #gamma_params for fixed l (state) and g (mixture component) (vector dim 2)
  n=length(data_sequences)
  maxN = max(N)
  penalty = Pen_lg(gamma_params,N)
  g2=as.numeric(g)
  shape_lg = gamma_params[1]
  rate_lg = gamma_params[2]
  
  I = sapply(1:maxN, function(k) {
    matrix_indicator(data_sequences, k, E)[, l] }) #matrix of dim nx(maxN) where each I(i,k) is equal to 1 if the k-th visited state in sequence i is l, and 0 otherwise
  
  X = matrix(0,n,maxN) #matrix with the sojourn_times
  #X[cbind(rep(1:n, N), sequence(N))] = unlist(lapply(data_sequences, function(seq) seq$sojourn_times))
  X = t(sapply(1:n, function(i){
    sapply(1:maxN, function(k){
      if (k<=N[i]){X[i,k] = data_sequences[[i]]$sojourn_times[k]}
      else {0}
    }) }))
  valid_indices = (X > 0)
  
  # log_gamma_cdf_valid = log(dgamma(X[valid_indices], shape=shape_lg, rate=rate_lg)) #log(0) returns inf
  # log_gamma_cdf = matrix(0, nrow=n, ncol=maxN)
  # log_gamma_cdf[valid_indices] = log_gamma_cdf_valid
  
  gamma_cdf_valid = dgamma(X[valid_indices], shape=shape_lg, rate=rate_lg) #log(0) returns inf
  gamma_cdf = matrix(0, nrow=n, ncol=maxN)
  gamma_cdf[valid_indices] = gamma_cdf_valid
  
  valid_indices2 = (gamma_cdf>0)
  log_gamma_cdf = matrix(0, nrow=n, ncol=maxN)
  log_gamma_cdf[valid_indices2] = log(gamma_cdf[valid_indices2]) #matrix with the log of gamma cdf for the associated sojourn_times
  
  if (G==1){
    result = sum(Z * I[, , drop = FALSE] * log_gamma_cdf[, , drop = FALSE]) + penalty
  }
  else {
    #result = sum(Z[,g]*diag(I%*%t(log_gamma_cdf))) + penalty
    valid_indices3 = !is.na(Z[, g2])
    result = sum(Z[valid_indices3, g2] * I[valid_indices3, , drop = FALSE] * log_gamma_cdf[valid_indices3, , drop = FALSE]) + penalty
    }
  
  if (!is.finite(result)) { #prints in case of error
    cat("Non-finite result in expected_partial_pen_loglikelihood:\n")
    cat("shape_lg:", shape_lg, "rate_lg:", rate_lg, "\n")
    cat("log_gamma_cdf:", log_gamma_cdf, "\n")
    cat("Z:", Z[, g2], "\n")
    cat("I:", I, "\n")
    #print(result)
  }
  
  return(result)
} #computes the expected partial penalized likelihood for a fixed state l and a fixed mixture component g


gradient_approx_expected_partial_pen_loglikelihood = function(gamma_params, data_sequences, Z, N, l, g, G, E) {
  grad_func = function(params) {
    #expected_partial_pen_loglikelihood(params, data_sequences, Z, N, l, g, G)
    val <- expected_partial_pen_loglikelihood(params, data_sequences, Z, N, l, g, G, E)
    if (!is.finite(val)) {
      cat("Non-finite value in grad_func for params:", params, "\n")
      cat("Returning large penalty value to force optimization away from invalid region.\n")
      return(-1e10) # Return a large penalty value to steer optimization away from this region
    }
    return(val)
  }
  grad = grad(func = grad_func, x = gamma_params)
  
  if (any(!is.finite(grad))) {
    cat("Non-finite gradient for params:", gamma_params, "\n")
    cat("Returning zero gradient to avoid invalid update.\n")
    grad[!is.finite(grad)] <- 0 #Replace non-finite gradient values with zero
  }
  return(grad)
} 


convert_params_structures = function(current_params,G,E){
  #alpha_mat = do.call(cbind, current_params$alpha)
  alpha_mat = matrix(unlist(current_params$alpha), nrow = E, ncol = G)
  
  P_tensor = array(0, dim = c(E, E, G))
  for (g in 1:G) {
    P_tensor[,,g] = current_params$P[[g]]
  }
  gamma_tensor = array(0, dim = c(E, 2, G))
  for (g in 1:G) {
    gamma_tensor[,,g] = current_params$gamma_params[[g]]
  }
  return(list(
    pi = current_params$pi,
    alpha = alpha_mat,
    P = P_tensor,
    gamma_params = gamma_tensor
  ))
}