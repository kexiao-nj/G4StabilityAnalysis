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


# for the robust network
sla_fun <- 'pc.stable'
test_fun <- 'mc-x2'
step <- 100
suffix <- 'bs1-9'
tardir <- 'robusttest'
M <- 7000
N <- 10

netobj <- paste('matrix', sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
matrix_list <- readRDS(paste('..', tardir, netobj, sep='/'))
matrix_list_filtered <- vector(step, mode = "list")
for (i in seq_along(matrix_list)){
  mat <- matrix_list[[i]]
  mat[mat <= 6] <- 0
  mat[mat > 0] <- 1
  matrix_list_filtered[[i]] <- mat
}

# graph_from_adjacency_matrix(matrix_list_filtered[[1]], mode = "directed")
# as.bn(graph_from_adjacency_matrix(matrix_list_filtered[[1]], mode = "directed"))
network_list_filtered <- vector(step, mode = "list")
for (i in seq_along(matrix_list_filtered)){
  network_list_filtered[[i]] <- as.bn(graph_from_adjacency_matrix(matrix_list_filtered[[i]], mode = "directed"))
}


str.net <- custom.strength(network_list_filtered, names(data1))
avg.dag <- averaged.network(str.net, threshold = 0.9)
bn_igraph <- bn2igraph(avg.dag, str.net)
g2 <- ggraph(bn_igraph, layout = "nicely") +
  geom_edge_fan(aes(color = dirct, label = round(dirct,2)), 
    width = 3, label_size = 6, 
    arrow = arrow(length = unit(4, "mm")), end_cap = circle(8, unit="mm")
  ) +
  geom_node_point(size=10,shape=21,fill='cadetblue1',color='cadetblue1')+ 
  geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
  scale_edge_color_gradient(low = "#FFEDA0", high = "#F03B20", name = "Direction Prob", limits = c(0, 1) ) +
  # guides(edge_width = "none", edge_alpha = "none" ) +
  theme_graph(base_family = "DejaVu Sans", background = "white") +
  theme(
    legend.title = element_text(size = 14), # Change legend title font size
    legend.text = element_text(size=14)
    # legend.key.size = unit(1.5, 'cm') # Change legend key size
  )
g2
netobj <- paste('network', sla_fun, test_fun, M, N, step, suffix, 'pdf', sep='.')
ggsave(here::here("output-fig/fig3-b-most-robust-network.pdf"), g2,
  width = 8, height = 6, device = cairo_pdf)



# for all the equal-allocation network
draw_all_equal_alloc <- function(sla_fun, test_fun, M_vec, N, step, suffix, tardir, keywords='all_equal_alloc') {
  tgraph_list <- vector(length(M_vec), mode = "list")
  j <- 1
  if (keywords=='all_equal_alloc'){
    labl <- 'Equal-allocation:'
  }
  else {
    labl <- 'Proportional-allocation:'
  }
  
  for (M in M_vec){
    netobj <- paste('matrix', sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
    matrix_list <- readRDS(paste('..', tardir, netobj, sep='/'))
    matrix_list_filtered <- vector(step, mode = "list")

    for (i in seq_along(matrix_list)){
      mat <- matrix_list[[i]]
      mat[mat <= 6] <- 0
      mat[mat > 0] <- 1
      matrix_list_filtered[[i]] <- mat
    }

    graph_from_adjacency_matrix(matrix_list_filtered[[1]], mode = "directed")
    as.bn(graph_from_adjacency_matrix(matrix_list_filtered[[1]], mode = "directed"))

    network_list_filtered <- vector(step, mode = "list")
    for (i in seq_along(matrix_list_filtered)){
      network_list_filtered[[i]] <- as.bn(graph_from_adjacency_matrix(matrix_list_filtered[[i]], mode = "directed"))
    }
    str.net <- custom.strength(network_list_filtered, names(data1))
    avg.dag <- averaged.network(str.net, threshold = 0.9)
    bn_igraph <- bn2igraph(avg.dag, str.net)
    igr <- bn_igraph %>% set_vertex_attr("type", value = paste(labl, M))
    tgraph_list[[j]] <- as_tbl_graph(igr)
    j <- j + 1
  }
  tg_list <- do.call(bind_graphs, tgraph_list) %>%
    mutate(type = factor(type, levels = paste(labl, M_vec)))

  g_all <-ggraph(tg_list, layout = "nicely") +
    geom_edge_fan(aes(color = strength, label = round(dirct,2)), 
      width = 3, label_size = 6, 
      arrow = arrow(length = unit(4, "mm")), end_cap = circle(8, unit="mm")
    ) +
    geom_node_point(size=10,shape=21,fill='cadetblue1',color='cadetblue1')+ 
    geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
    scale_edge_color_gradient(low = "#FFEDA0", high = "#F03B20", name = "Direction Prob", limits = c(0, 1) ) +
    # guides(edge_width = "none", edge_alpha = "none" ) +
    theme_graph(base_family = "DejaVu Sans", background = "white") + 
    theme(strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold", size = 12))+
    facet_nodes(~type, scales = "free", ncol = 3)

  g_all
  netobj <- paste('network', sla_fun, test_fun, keywords, N, step, suffix, 'pdf', sep='.')
  ggsave(here::here(paste("output-fig/figs5-",netobj,sep='')), g_all,
    width = 20, height = 6, device = cairo_pdf)
}

sla_fun <- 'pc.stable'
test_fun <- 'mc-x2'
step <- 100
suffix <- 'bs1-9'
tardir <- 'robusttest'
N <- 10

M_vec <- c(4000, 5000, 6000)
draw_all_equal_alloc(sla_fun, test_fun, M_vec, N, step, suffix, tardir, keywords='all_equal_alloc')

M_vec <- c(0.4, 0.5, 0.6)
draw_all_equal_alloc(sla_fun, test_fun, M_vec, N, step, suffix, tardir, keywords='all_prop_alloc')
