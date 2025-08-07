library(ggplot2)
library(bnlearn)
library(igraph)
library(ggraph)
library(tidygraph)
library(showtext)
setwd(here::here("analysis"))

########################## gen accuracy and coverage #######################################
calc_acc_covr <- function(sla_fun='pc.stable', test_fun='mc-x2', step=100, suffix='bs1-9', tardir='robusttest', M=0.6, N=10){
  prefix <- 'matrix'
  netobj <- paste(prefix, sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')

  tardir2 <- paste('..', tardir, sep='/')
  matrix_list <- readRDS(paste(tardir2, netobj, sep='/'))

  column_names <- c(1:10)
  acc_df <- data.frame(matrix(ncol = length(column_names), nrow = step))
  colnames(acc_df) <- column_names
  cov_df <- data.frame(matrix(ncol = length(column_names), nrow = step))
  colnames(cov_df) <- column_names

  for(i in seq_along(matrix_list)){
    sum_mat <- sum(matrix_list[[i]])
    acc_row <- c()
    cov_row <- c()
    for(j in 1:N){
      acc_row <- append(acc_row, sum(matrix_list[[i]][matrix_list[[i]] >= j])/sum_mat)
      cov_row <- append(cov_row, mean(matrix_list[[i]][matrix_list[[i]] >= j])/N)
    }
    acc_df[i,] <- acc_row
    cov_df[i,] <- cov_row
  }

  prefix <- 'acc'
  netobj <- paste(prefix, sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
  saveRDS(acc_df, file=paste(tardir2, netobj, sep='/'))

  prefix <- 'cov'
  netobj <- paste(prefix, sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
  saveRDS(cov_df, file=paste(tardir2, netobj, sep='/'))
}

# calc_acc_covr(sla_fun='pc.stable', test_fun='mc-x2', step=100, suffix='bs1-9', tardir='robusttest', N=10)
sla_fun <- "pc.stable"
test_fun <- "mc-x2"
step <- 100
suffix <- 'bs1-9'
tardir <- 'robusttest'

for (proportion in seq(0.1,0.6,0.1)){
  print("##############################################################")
  print(paste("proportion is", proportion))
  calc_acc_covr(sla_fun=sla_fun, test_fun=test_fun, step=step, suffix=suffix, tardir=tardir, M=proportion, N=10)
}

for (size in seq(1000,7000,1000)){
  print("##############################################################")
  print(paste("size is", size))
  calc_acc_covr(sla_fun=sla_fun, test_fun=test_fun, step=step, suffix=suffix, tardir=tardir, M=size, N=10)
}


########################## calculate mean and std ##############################
calc_mean_std <- function(sla_fun='pc.stable', test_fun='mc-x2', step=100, suffix='bs1-9', tardir='robusttest', N=10){
  final_df <- data.frame()
  tardir2 <- paste('..', tardir, sep='/')
  for(M in seq(0.1,0.6,0.1)){
    tmp_df <- data.frame()
    netobj <- paste('acc', sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
    acc_df <- readRDS(paste(tardir2, netobj, sep='/'))
    for(i in 1:N){
      tmp_df[i,'acc'] <- mean(acc_df[,i])
      tmp_df[i,'accstd'] <- sd(acc_df[,i])
    }

    netobj <- paste('cov', sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
    cov_df <- readRDS(paste(tardir2, netobj, sep='/'))
    for(i in 1:N){
      tmp_df[i,'cov'] <- mean(cov_df[,i])
      tmp_df[i,'covstd'] <- sd(cov_df[,i])
    }
    tmp_df['Model'] <- paste('Proportion:', M, sep='')
    final_df <- rbind(final_df, tmp_df)
  }

  for(M in seq(1000,7000,1000)){
    tmp_df <- data.frame()

    netobj <- paste('acc', sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
    acc_df <- readRDS(paste(tardir2, netobj, sep='/'))
    for(i in 1:N){
      tmp_df[i,'acc'] <- mean(acc_df[,i])
      tmp_df[i,'accstd'] <- sd(acc_df[,i])
    }

    netobj <- paste('cov', sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
    cov_df <- readRDS(paste(tardir2, netobj, sep='/'))
    for(i in 1:N){
      tmp_df[i,'cov'] <- mean(cov_df[,i])
      tmp_df[i,'covstd'] <- sd(cov_df[,i])
    }
    tmp_df['Model'] <- paste('Size:', M, sep='')
    final_df <- rbind(final_df, tmp_df)
  }
  return(final_df)
}

mean_std_df <- calc_mean_std(sla_fun='pc.stable', test_fun='mc-x2', step=100, suffix='bs1-9', tardir='robusttest', N=10)


g1 <- ggplot(mean_std_df[is.element(mean_std_df$Model, c('Proportion:0.5', 'Proportion:0.6', 'Size:6000', 'Size:7000')), ], 
      aes(color=Model, x=cov, y=acc)) +
  geom_line(size = 1.2) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = acc-accstd, ymax = acc+accstd), width=0.01, size=0.1) +
  geom_errorbar(aes(xmin = cov-covstd, xmax = cov+covstd), width=0.01, size=0.1) +
  labs(x = "Coverage",
       y = "Accuracy",
       color = "") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2))
g1
ggsave(here::here("output-fig/fig3-a-cov-acc-curve.pdf"), g1, width = 6, height = 6, device = cairo_pdf)

g2 <- ggplot(mean_std_df[grepl('Proportion', mean_std_df$Model), ], 
      aes(color=Model, x=cov, y=acc)) +
  geom_errorbar(aes(ymin = acc-accstd, ymax = acc+accstd), width=0.01, size=0.1) +
  geom_errorbar(aes(xmin = cov-covstd, xmax = cov+covstd), width=0.01, size=0.1) +
  geom_line(size = 1.2) + 
  geom_point(size = 3) + 
  labs(x = "Coverage",
       y = "Accuracy",
       color = "") +
  theme_bw()
g2
ggsave(here::here("output-fig/figs3-prop-cov-acc-curve.pdf"), g2, width = 7, height = 6, device = cairo_pdf)

g3 <- ggplot(mean_std_df[grepl('Size', mean_std_df$Model), ], 
      aes(color=Model, x=cov, y=acc)) +
  geom_errorbar(aes(ymin = acc-accstd, ymax = acc+accstd), width=0.01, size=0.1) +
  geom_errorbar(aes(xmin = cov-covstd, xmax = cov+covstd), width=0.01, size=0.1) +
  geom_line(size = 1.2) + 
  geom_point(size = 3) + 
  labs(x = "Coverage",
       y = "Accuracy",
       color = "") +
  theme_bw()
g3
ggsave(here::here("output-fig/figs3-equl-cov-acc-curve.pdf"), g3, width = 7, height = 6, device = cairo_pdf)
