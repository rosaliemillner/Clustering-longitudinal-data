#initialisation of the EM algorithm

Y_mean = function(data_sequences,E){
  S=1:E
  n_indiv=length(data_sequences)
  Y_means=matrix(0, nrow=n_indiv, ncol=E)
  colnames(Y_means)=S
  for (i in 1:n_indiv){
    indiv_data_states = data_sequences[[i]]$states
    indiv_data_times = data_sequences[[i]]$sojourn_times
    state_list = split(indiv_data_times, indiv_data_states)
    
    means = sapply(state_list, mean)
    
    Y_means[i, names(means)]=means
  }
  return(Y_means) 
} #returns matrix with means of sojourn times per state per individual

Y_kmeans = function(data_sequences,Ym,G){
  kmeans_result = kmeans(Ym, centers = G)
  cluster_assignments = kmeans_result$cluster
  Xg_list = vector("list",G)
  for (g in 1:G){
    indices = which(cluster_assignments == g) 
    Xg_list[[g]] <- lapply(indices, function(i) data_sequences[[i]])
  }
  return(Xg_list)
} #performs k-means algo and returns lists w sequences for each cluster

Y_cluster_var = function(data_seq, E){ #here, data_seq is the list of sequences within one cluster
  data_states = lapply(data_seq, function(x) x[[1]])
  data_sojourn_times = lapply(data_seq, function(x) x[[2]])
  combined_states = do.call(c, data_states)
  combined_sojourn_times = do.call(c, data_sojourn_times)
  
  Yg_vars = rep(NA, E)
  
  state_list = split(combined_sojourn_times, combined_states)
  state_vars = sapply(state_list, var)
  
  state_vars[is.na(state_vars)] = 0 #(qd il y a un seul ele, var=NA)
  
  state_indices = as.numeric(names(state_vars))
  Yg_vars[state_indices] = state_vars

  return(Yg_vars)
} #returns a vector with the vars of sojourn times per state in one cluster

Y_cluster_mean = function(data_seq,E){ 
  Yg = Y_mean(data_seq,E)
  Ygm = apply(Yg,2,mean) 
  return(Ygm) 
} #returns a vector with the means of the means of sojourn times per state in one same cluster

gamma_parametres = function(Ym,Yv,G,E){
  gamma_params = list()
  for (g in 1:G){
    gamma_matrix_g = matrix(0,E,2)
    for (d in 1:E){
      mu=Ym[[g]][d]
      var=Yv[[g]][d]
      if (var==0 && !is.na(var)){var=1} #valeur arbitraire
      
      if (is.na(var)){gamma_matrix_g[d,]=NA} #state is not visited
      else {
        gamma_matrix_g[d,1]=(mu^2)/var
        gamma_matrix_g[d,2]=mu/var 
        }
    }
    gamma_params[[g]] = gamma_matrix_g
  }
  return(gamma_params)
} #performs the moments method to identify the parameters of the gamma distributions of the sojourn times

proportion_p_init = function(seq_from_a_cluster,E){
  list_only_states = lapply(seq_from_a_cluster, function(x) x[[1]]) #list of vectors, only the states
  alpha = lapply(1:E, function(num) { mean(unlist(lapply(list_only_states, function(x) x[1] == num))) })
  return(alpha)
} 

proportion_transitionmatrix = function(seq_from_a_cluster,E){
  only_states = lapply(seq_from_a_cluster, function(x) x[[1]])
  all_states = unlist(only_states)
  
  num_transitions = matrix(0, nrow = E, ncol = E)
  from_states = all_states[-length(all_states)]
  to_states = all_states[-1]
  
  transition_counts = table(from_states, to_states)
  
  num_transitions[as.numeric(rownames(transition_counts)), as.numeric(colnames(transition_counts))] = transition_counts

  transition_matrix = num_transitions / rowSums(num_transitions)
  transition_matrix[is.na(transition_matrix)] = 0
  return(transition_matrix)
}

init_EM_algo = function(data_sequences, G, E){  #with G the number of clusters and E the total number of states
  Ym = Y_mean(data_sequences,E) 
  data_sequences_clustered_list = Y_kmeans(data_sequences, Ym, G)
  Yg_cluster_mean = lapply(data_sequences_clustered_list, Y_cluster_mean, E) 
  Yg_cluster_var = lapply(data_sequences_clustered_list, Y_cluster_var, E)
  
  init_gamma_params_list = gamma_parametres(Yg_cluster_mean,Yg_cluster_var,G,E)
  init_pi = rep(1/G,G)  #uniform distribution
  init_alpha = lapply(data_sequences_clustered_list, proportion_p_init,E)
  init_P = lapply(data_sequences_clustered_list, proportion_transitionmatrix,E)
  
  return(list(pi=init_pi, alpha=init_alpha, P=init_P, gamma_params = init_gamma_params_list))
}
