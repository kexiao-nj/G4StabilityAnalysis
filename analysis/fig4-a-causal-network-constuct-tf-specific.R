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

bn2igraph <- function(bn, bn_str) {
  bn_igraph <- as.igraph(bn)
  E(bn_igraph)$strength <- bn_str$strength[match(
    paste(arcs(bn)[, 1], arcs(bn)[, 2]),
    paste(bn_str$from, bn_str$to)
  )]
  E(bn_igraph)$dirct <- bn_str$direction[match(
    paste(arcs(bn)[, 1], arcs(bn)[, 2]),
    paste(bn_str$from, bn_str$to)
  )]
  return(bn_igraph)
}

gen_datalist <- function(datalist, namelist) {
  outlist <- vector(length(datalist), mode = "list")
  for (i in seq_along(datalist)){
    outlist[[i]] <- datalist[[i]][,namelist]
  }
  return(outlist)
}


common_tf <- c('SP2', 'SP1', 'YY1', 'CTCF', 'FOXA1', 'TARDBP')
# read prepared datasets
vars <- c('Stability', 'ATACSig', 'phyloP')
target_vals <- c('Stability.strength', 'phyloP.strength', 'ATACSig.strength', 'ChromState')
renamed_vals <- c('Stability', 'phyloP', 'ATACSig', 'ChromState')

## K562
in_df1 <- data.frame(read.csv('../prepared_datasets/K562.tsv', sep = "\t"))
in_df1 <- discrt_columns_by_median(in_df1, vars)
data1 <- in_df1[,target_vals]
data1$ChromState <- factor(data1$ChromState, levels = c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk", "7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF/Rpts", "13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies"))
names(data1) <- renamed_vals
factor_data1 <- lapply(in_df1[,common_tf], function(col){
  factor(col, levels=c(0,1))
})
data1 <- cbind(data1,as.data.frame(factor_data1))

## HepG2
in_df2 <- data.frame(read.csv('../prepared_datasets/HepG2.tsv', sep = "\t"))
in_df2 <- discrt_columns_by_median(in_df2, vars)
data2 <- in_df2[,target_vals]
data2$ChromState <- factor(data2$ChromState, levels = c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk", "7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF/Rpts", "13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies"))
names(data2) <- renamed_vals
factor_data2 <- lapply(in_df2[,common_tf], function(col){
  factor(col, levels=c(0,1))
})
data2 <- cbind(data2,as.data.frame(factor_data2))

## 293T
in_df3 <- data.frame(read.csv('../prepared_datasets/293T.tsv', sep = "\t"))
in_df3 <- discrt_columns_by_median(in_df3, vars)
data3 <- in_df3[,target_vals]
data3$ChromState <- factor(data3$ChromState, levels = c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk", "7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF/Rpts", "13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies"))
names(data3) <- renamed_vals
factor_data3 <- lapply(in_df3[,common_tf], function(col){
  factor(col, levels=c(0,1))
})
data3 <- cbind(data3,as.data.frame(factor_data3))

data_list <- list(data1, data2, data3)


############### workflow for the tf networks ######################
gen_Ngroup_selfadaptingSize_combination <- function(df, M=2000, N=10) {
  bootstrapped_datasets <- list()
  if(is.numeric(M) & (M>0) & (M<=1)){
    # M is the proportion for sampling
    samplsize <- floor(M*nrow(df))
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
  merged_matrix <- matrix(rep(0, ncol(data_list[[1]])^2), nrow = ncol(data_list[[1]]))
  for (i in seq_along(bootstrapped_multisets)){
    # randomize the columns
    bootsample <- bootstrapped_multisets[[i]]
    bootsample <- bootsample[, sample(ncol(bootsample), ncol(bootsample)), drop = FALSE]
    netlist[[i]] <- do.call(sla_fun, list(x = bootsample, test = test_fun))
    g<-as.igraph(netlist[[i]])
    my_matrix <- as_adjacency_matrix(g)
    merged_matrix <- merged_matrix + my_matrix[match(names(data_list[[1]]), rownames(my_matrix)), match(names(data_list[[1]]), colnames(my_matrix))]
  }
  return(list(netlist, merged_matrix))
}

#bootstraping experiments, multiple times
multi_bootstrap_experiment <- function(data_list, gen_method='gen_Ngroup_selfadaptingSize_combination', sla_fun='pc.stable', test_fun='x2', step=10, suffix='bs1-9', tardir='robusttest', M=0.9, N=10) {
  matrix_list <- vector(step, mode = "list")
  netlist_list <- vector(step, mode = "list")
  # for(i in 1:step){
  #   tmplist <- one_bootstrap_experiment(data_list, gen_method, sla_fun, test_fun, M, N)
  #   netlist_list[[i]] <- tmplist[[1]]
  #   matrix_list[[i]] <- tmplist[[2]]
  #   print(matrix_list[[i]])
  # }
  
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


# Now test different strategies
sla_fun <- 'pc.stable'
test_fun <- 'mc-x2'
step <- 100
tardir <- 'tfrobustnet'

common_tf <- c('SP2', 'SP1', 'YY1', 'CTCF', 'FOXA1', 'TARDBP')
for (tf in common_tf){
  print("##############################################################")
  print(paste("TF is", tf))
  namelist <- c("Stability", "ChromState", "phyloP", 'ATACSig', tf)
  print(namelist)
  sigletf_datalist <- gen_datalist(data_list, namelist)
  multi_bootstrap_experiment(data_list=sigletf_datalist, gen_method='gen_Ngroup_selfadaptingSize_combination', sla_fun=sla_fun, test_fun=test_fun, step=step, suffix=paste(tf,'atac',sep=''), tardir=tardir, M=7000)
}
