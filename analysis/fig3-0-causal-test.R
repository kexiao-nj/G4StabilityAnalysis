library(ggplot2)
library(bnlearn)
library(igraph)
library(ggraph)
library(tidygraph)
library(showtext)
setwd(here::here("analysis"))

discrt_columns_by_median <- function(in_df, cols) {
    for (c in cols){
        # print(c)
        med <- median(in_df[, c])
        in_df[,paste(c, '.strength', sep='')] <- 'High'
        in_df[in_df[, c] < med, paste(c, '.strength', sep='')] <- 'Low'
        in_df[,paste(c, '.strength', sep='')] <- factor(in_df[,paste(c, '.strength', sep='')], 
            levels = c('Low','High'))
    }
    return(in_df)
}

# read prepared datasets
vars <- c('eG4Sig', 'Stability', 'ATACSig', 'phyloP', 'TFnum')
target_vals <- c('eG4Sig.strength', 'Stability.strength', 'phyloP.strength', 'ATACSig.strength', 'ChromState', 'TFnum.strength')
renamed_vals <- c('eG4Sig', 'Stability', 'phyloP', 'ATACSig', 'ChromState', 'TFNo')

## K562
in_df1 <- data.frame(read.csv('../prepared_datasets/K562.tsv', sep = "\t"))
in_df1 <- discrt_columns_by_median(in_df1, vars)
data1 <- in_df1[,target_vals]
data1$ChromState <- factor(data1$ChromState, levels = c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk", "7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF/Rpts", "13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies"))
names(data1) <- renamed_vals

## HepG2
in_df2 <- data.frame(read.csv('../prepared_datasets/HepG2.tsv', sep = "\t"))
in_df2 <- discrt_columns_by_median(in_df2, vars)
data2 <- in_df2[,target_vals]
data2$ChromState <- factor(data2$ChromState, levels = c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk", "7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF/Rpts", "13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies"))
names(data2) <- renamed_vals

## 293T
in_df3 <- data.frame(read.csv('../prepared_datasets/293T.tsv', sep = "\t"))
in_df3 <- discrt_columns_by_median(in_df3, vars)
data3 <- in_df3[,target_vals]
data3$ChromState <- factor(data3$ChromState, levels = c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk", "7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF/Rpts", "13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies"))
names(data3) <- renamed_vals

data_list <- list(data1, data2, data3)


########### main workflow ##############
#divide data into N groups
gen_Ngroup_combination <- function(df, P=0.9, N=10) {
  df_shuffled <- df[sample(nrow(df)), ]
  n <- nrow(df_shuffled)
  df_list <- split(df_shuffled, rep(1:N, each = n %/% N, length.out = n))

  bootstrapped_datasets <- list()
  for (i in 1:N) {
    selected_indices <- setdiff(1:N, i)
    bootstrapped_datasets[[i]] <- do.call(rbind, df_list[selected_indices])
  }
  return(bootstrapped_datasets)
}

gen_Msize_Ngroup_combination <- function(df, M=2000, N=10) {
  bootstrapped_datasets <- list()
  for (i in 1:N) {
    bootstrapped_datasets[[i]] <- df[sample(nrow(df), M), ]
  }
  return(bootstrapped_datasets)
}

gen_Ngroup_selfadaptingSize_combination <- function(df, M=2000, N=10) {
  bootstrapped_datasets <- list()
  if(is.numeric(M) & (M>0) & (M<=1)){
    # M is the proportion for sampling
    samplsize <- floor(M*nrow(df))
    # print(paste("+++++++++++++++++++++++", samplsize))
  } else if (is.numeric(M) & (M %% 1 == 0)) {
    # M is the size of sampling
    samplsize <- M
  } else {
    stop("Error: The value of parameter M should be either a fraction smaller than 1, or an integer greater than 1.")
  }
  
  for (i in 1:N) {
    bootstrapped_datasets[[i]] <- df[sample(nrow(df), samplsize), ]
  }
  return(bootstrapped_datasets)
}

gen_Ngroup_combination_multisample <- function(df_list, gen_fun='gen_Ngroup_selfadaptingSize_combination', M=0.9, N=10) {
  bootstrapped_dataset_list <- list()
  for (i in seq_along(df_list)){
    bootstrapped_dataset_list[[i]] <- do.call(gen_fun, list(df_list[[i]], M, N))
  }
  bootstrapped_multisets <- do.call(Map, c(list(rbind), bootstrapped_dataset_list))
  return(bootstrapped_multisets)
}

#bootstraping based on nine-group combination datasets * 10 times
one_bootstrap_experiment <- function(data_list, gen_method = 'gen_Ngroup_selfadaptingSize_combination', sla_fun = 'pc.stable', test_fun = 'x2', M=0.9, N=10) {
  bootstrapped_multisets <- gen_Ngroup_combination_multisample(data_list, gen_method, M, N)
  print(paste("nrow(bootstrapped_multisets[[1]]) = ", nrow(bootstrapped_multisets[[1]])))
  netlist <- vector(length(bootstrapped_multisets), mode = "list")
  # merged_matrix <- matrix(rep(0, ncol(data1)^2), nrow = ncol(data1), dimnames = list(renamed_vals, renamed_vals))
  merged_matrix <- matrix(rep(0, ncol(data1)^2), nrow = ncol(data1))
  for (i in seq_along(bootstrapped_multisets)){
    # randomize the columns
    bootsample <- bootstrapped_multisets[[i]]
    bootsample <- bootsample[, sample(ncol(bootsample), ncol(bootsample)), drop = FALSE]
    netlist[[i]] <- do.call(sla_fun, list(x = bootsample, test = test_fun))
    # print(netlist)
    g<-as.igraph(netlist[[i]])
    my_matrix <- as_adjacency_matrix(g)
    # print(my_matrix)
    merged_matrix <- merged_matrix + my_matrix[match(renamed_vals, rownames(my_matrix)), match(renamed_vals, colnames(my_matrix))]
    # print(merged_matrix)
  }
  return(list(netlist, merged_matrix))
}

#bootstraping experiments, multiple times
multi_bootstrap_experiment <- function(data_list, gen_method='gen_Ngroup_selfadaptingSize_combination', sla_fun='pc.stable', test_fun='x2', step=10, suffix='bs1-9', tardir='robusttest', M=0.9, N=10) {
  matrix_list <- vector(step, mode = "list")
  netlist_list <- vector(step, mode = "list")
  for(i in 1:step){
    tmplist <- one_bootstrap_experiment(data_list, gen_method, sla_fun, test_fun, M, N)
    netlist_list[[i]] <- tmplist[[1]]
    matrix_list[[i]] <- tmplist[[2]]
    print(matrix_list[[i]])
  }
  
  tardir2 <- paste('..', tardir, sep='/')
  if (!dir.exists(tardir2)) {
    dir.create(tardir2)
  }

  prefix <- 'netlist'
  netobj <- paste(prefix, sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
  saveRDS(netlist_list, file=paste(tardir2, netobj, sep='/'))

  prefix <- 'matrix'
  netobj <- paste(prefix, sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
  saveRDS(matrix_list, file=paste(tardir2, netobj, sep='/'))
}

# Now test different strategies, it might take a long time. So I suggest for result demonstration, you can just use the data generated prevously in the robusttest folder to save time
sla_fun <- 'pc.stable'
test_fun <- 'mc-x2'
step <- 100
suffix <- 'bs1-9'
tardir <- 'robusttest'

for (proportion in seq(0.1,0.6,0.1)){
  print("##############################################################")
  print(paste("proportion is", proportion))
  multi_bootstrap_experiment(data_list=data_list, gen_method='gen_Ngroup_selfadaptingSize_combination', sla_fun=sla_fun, test_fun=test_fun, step=step, suffix=suffix, tardir=tardir, M=proportion)
}


for (size in seq(1000,7000,1000)){
  print("##############################################################")
  print(paste("size is", size))
  multi_bootstrap_experiment(data_list=data_list, gen_method='gen_Ngroup_selfadaptingSize_combination', sla_fun=sla_fun, test_fun=test_fun, step=step, suffix=suffix, tardir=tardir, M=size)
}