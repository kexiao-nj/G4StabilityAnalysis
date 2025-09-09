library(ggplot2)
library(bnlearn)
library(igraph)
library(ggraph)
library(tidygraph)
library(showtext)
showtext_auto()
library(Matrix)
setwd(here::here("analysis"))

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

# for the tf network
sla_fun <- 'pc.stable'
test_fun <- 'mc-x2'
step <- 50
suffix <- 'SP1atac'
tardir <- 'tfrobustnet'
M <- 7000
N <- 10
tflist <- c('SP1','SP2',  'YY1', 'TARDBP', 'CTCF', 'FOXA1')

draw_tfnetwork_batch <- function(tflist, sla_fun='pc.stable', test_fun='mc-x2', step, tardir='tfrobustnet', M, N=10) {
  tardir2 <- paste('..', tardir, sep='/')
  if (!dir.exists(tardir2)) {
    dir.create(tardir2)
  }

  tgraph_list <- vector(length(tflist), mode = "list")

  for (j in seq_along(tflist)) {
    suffix <- paste(tflist[[j]], 'atac', sep = "")
    # namelist <- c("Stability", "ChromState", "phyloP", "ATACSig", tflist[[j]])
    namelist <- c("Stability", "ChromState", "phyloP", "ATACSig", "eG4Sig", tflist[[j]])
    netobj <- paste('matrix', sla_fun, test_fun, M, N, step, suffix, 'Rds', sep='.')
    matrix_list <- readRDS(paste(tardir2, netobj, sep='/'))
    matrix_list_filtered <- vector(step, mode = "list")

    for (i in seq_along(matrix_list)){
      mat <- matrix_list[[i]]
      mat[mat <= 5] <- 0
      mat[mat > 0] <- 1
      matrix_list_filtered[[i]] <- mat
    }
    network_list_filtered <- vector(step, mode = "list")
    for (i in seq_along(matrix_list_filtered)){
      network_list_filtered[[i]] <- as.bn(graph_from_adjacency_matrix(matrix_list_filtered[[i]], mode = "directed"), check.cycles = FALSE)
    }
    str.net <- custom.strength(network_list_filtered, namelist)
    avg.dag <- averaged.network(str.net, threshold = 1)
    bn_igraph <- bn2igraph(avg.dag, str.net)

    Matr <- as_adjacency_matrix(bn_igraph, attr='dirct')
    write.table(as.matrix(Matr), file=paste('../output-fig/table.robust.net', tflist[[j]], '.csv', sep=''), sep=",", quote = FALSE)

    bn_igraph <- bn_igraph %>% set_vertex_attr("type", value = tflist[[j]])
    tgraph_list[[j]] <- as_tbl_graph(bn_igraph)
  }
  tg_list <- do.call(bind_graphs, tgraph_list) %>%
    mutate(type = factor(type, levels = tflist))
  g_all <-ggraph(tg_list, layout = "fr") +
    geom_edge_fan(aes(color = dirct),
      width = 3, label_size = 6, 
      arrow = arrow(length = unit(4, "mm")), end_cap = circle(8, unit="mm")
    ) +
    geom_node_point(size=10,shape=21,fill='cadetblue1',color='cadetblue1')+ 
    geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
    scale_edge_color_gradient(low = "#FFEDA0", high = "#F03B20", name = "Direction Prob", limits = c(0, 1) ) +
    # guides(edge_width = "none", edge_alpha = "none" ) +
    theme_graph(base_family = "DejaVu Sans", background = "grey98") + 
    theme(strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold", size = 12))+
    facet_nodes(~type, scales = "free", ncol = 3)
  netobj <- paste('tf', sla_fun, test_fun, M, N, step, 'all', 'allshared.pdf', sep='.')
  ggsave(here::here("output-fig/fig4-tf-specific-networks.pdf"), g_all,
    width = 15, height = 10, device = cairo_pdf)

}

draw_tfnetwork_batch(tflist, M=M, step=step, tardir=tardir)
