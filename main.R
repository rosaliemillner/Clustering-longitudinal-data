#main, refer to this code to have a template of how to use the algorithm

rm(list=ls())

source("data_processing.R")
source("init_EM.R")
source("functions_EM.R")
source("EM_algo.R")
source("functions_main.R")


#--------- part of code TO ADAPT -----------

data = read.csv("data.csv", header = TRUE, sep = ",")
#OR data = read.table("data.txt", header = TRUE, sep = "\t")

#select the columns of interest (only longitudinal data, get rid of columns such as 'ID'...)

#E = #number of states in the data, characterising the state space S=1:E

  
#----------- data processing -----------

data_sequences = transfo_data(data)


#---------------- model ----------------

G_max = floor(E/2)  # maximal number of clusters to test (or set a smaller value arbitrarily)
C = 1:G_max  # number of clusters to run on and test on

result = model_results(data_sequences, C, E) # EM ALGORITHM, additional optional arguments: tolerance (tol) and max number of iterations (itermax)


# -------------- results --------------

BIC_for_each_G = result$BIC_vect

optimal_G = result$optimal_G

cluster_affectations_optimalG = result$groups[[optimal_G]]

markov_parameters_optimalG = result$params[[optimal_G]]


