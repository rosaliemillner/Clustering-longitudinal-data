#conversion of the data to the right format for the algorithm (lists of 2 lists: states & sojourn_times per sequence

#to adapt for the dataset of interest:
data_all = read.table("donnees.txt", header = TRUE, sep = "\t")
data = data_all[, 22:ncol(data_all)]
data = data[, -ncol(data)]


#transfo of the data to have the 2 sequences: successive states and associated sojourn times
individual_seqs = function(row) {
  states = c()
  times = c()
  p=length(row)
  for (i in 1:(p-1)) {
    if (is.na(row[i])){
      break 
    }
    if (is.na(row[i + 1]) || row[i] != row[i + 1]) {
      states = c(states, row[i])
      times = c(times, i)
    }
  }
  if (!is.na(row[p])){
    states = c(states, row[p])
    times = c(times, p)
  }
  sojourn_times = diff(times)
  sojourn_times = c(times[1],sojourn_times)
  return(list(states = unname(unlist(states)), sojourn_times = sojourn_times))
}

transfo_data = function(table){ #before applying this function, replace (if there are any) all the non-existing states by NA (if for instance some have been added to make all the sequences have the same length)
  data_sequences = lapply(1:nrow(table), function(i) individual_seqs(table[i,]))
  return(data_sequences)
} #to transform the dataset into the expected format for the algorithm: lists of 2 lists per sequence


#data_sequences = lapply(1:nrow(data), function(i) individual_seqs(data[i,]))
#saveRDS(data_sequences, file="data_sequences.rds")
data_sequences = readRDS("data_sequences.rds")


