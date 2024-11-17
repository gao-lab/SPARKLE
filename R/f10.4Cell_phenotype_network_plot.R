#' Cell Phenotype Network Plot
#'
#' This function generates a network plot to visualize cell phenotype interactions based on the provided interaction data.
#' The plot highlights significant interactions with edges colored and weighted according to the p-values and odds ratios (OR).
#'
#' @param cwas.test Data frame containing CWAS test data. If `best.model` is provided, this parameter is ignored.
#' @param best.model List containing the best model information. It should have a component `Chosen_model_info` which is a data frame with columns `celltype`, `Pvalue`, and `Beta`. Default is `NULL`.
#' @param interaction.best.model List containing the best interaction model information. It should have a component `Chosen_model_info` which is a data frame with columns `celltype`, `Pvalue`, and `Beta`. Default is `NULL`.
#' @param cell.color Character vector specifying the colors for the cells. Default is `c("#cd3023", "#b91e45", "#f9bcbd", "gray", "#9BBBE1", "#90C4E9", "#8095CE")`.
#' @param cell.size Numeric vector specifying the sizes for the cells. Default is `c(25, 20, 15, 14)`.
#' @param edge.width Numeric vector specifying the widths of the edges for different significance levels. Default is `c(2, 3, 4)`.
#' @param arrow.size Numeric value specifying the size of the arrows. Default is `0.75`.
#' @param arrow.width Numeric value specifying the width of the arrows. Default is `1.5`.
#' @param edge.color Character vector specifying the colors for the edges based on p-values and OR. Default is `c("darkred", "#b91e45", "#f9bcbd", "#9BBBE1", "#4198b9", "#2d4e76")`.
#' @param Interation.effect Logical value indicating whether to use interaction effect or moderation effect. Default is `TRUE`.
#'
#' @return A plot of the cell phenotype network.
#' @export
#'
Cell_phenotype_network_plot <- function (cwas.test = NULL, best.model = NULL, interaction.best.model = NULL, 
                 cell.color = c("#cd3023", "#b91e45", "#f9bcbd", "gray", 
                                         "#9BBBE1", "#90C4E9", "#8095CE"), cell.size = c(25, 
                                         20, 15, 14), edge.width = c(2, 3, 4), arrow.size = 0.75, 
                 arrow.width = 1.5, edge.color = c("darkred", "#b91e45", 
                                                            "#f9bcbd", "#9BBBE1", "#4198b9", "#2d4e76"),Interation.effect=T) 
  {
  if (is.null(best.model)) {
    if (is.null(cwas.test)) {
      stop("Please input caws data  data!", call. = FALSE)
    }
    best.all <- cwas_allmodel_cal(cwas.test)
    best.model <- cwas_autoselected_model(best.all)
    best.info <- best.model[["Chosen_model_info"]]
  }
  else {
    best.info <- best.model[["Chosen_model_info"]]
  }
  if (is.null(interaction.best.model)) {
    interaction <- cwas_2celltype_allmodel_cal(cwas.test, 
                                               interaction = T)
    interaction.best.model <- cwas_autoselected_model(interaction)
    interaction.best.info <- interaction.best.model[["Chosen_model_info"]]
  }
  else {
    interaction.best.info <- interaction.best.model[["Chosen_model_info"]]
  }
  
  
  
  edges <- strsplit(as.character(interaction.best.info$celltype), 
                    " ")
  edges <- do.call(rbind, edges)
  edges <- as.data.frame(edges)
  
  interaction.best.info$TargetCell <- edges[, 1]
  interaction.best.info$ModerateCell <- edges[, 2]
  best.info2 <- best.info
  best.info2$TargetCell <- best.info2$celltype
  best.info2$Beta.Single <- best.info2$Beta 
  best.info2 <- best.info2[,c("TargetCell","Beta.Single")]
  
  interaction.best.info <- merge(interaction.best.info,best.info2,by="TargetCell")
  if(Interation.effect){interaction.best.info$OR <- interaction.best.info$Beta }else{
    interaction.best.info$OR <- interaction.best.info$Beta * interaction.best.info$Beta.Single
    
  }
  
  
  edges_df <- data.frame(from = edges[, 2], to = edges[, 1], 
                         Pvalue = interaction.best.info$Pvalue, OR =interaction.best.info$OR )
  edges_df$weight <- ifelse(edges_df$Pvalue < 0.001, edge.width[1], 
                            ifelse(edges_df$Pvalue < 0.01, edge.width[2], ifelse(edges_df$Pvalue < 
                                                                                   0.05, edge.width[3], 0.1)))
  edges_df$color <- ifelse((edges_df$Pvalue < 0.001) & (edges_df$OR > 
                                                          0), edge.color[1], ifelse((edges_df$Pvalue < 0.01) & 
                                                                                      (edges_df$OR > 0), edge.color[2], ifelse((edges_df$Pvalue < 
                                                                                                                                  0.05) & (edges_df$OR > 0), edge.color[3], ifelse((edges_df$Pvalue < 
                                                                                                                                                                                      0.001) & (edges_df$OR < 0), edge.color[4], ifelse((edges_df$Pvalue < 
                                                                                                                                                                                                                                           0.01) & (edges_df$OR < 0), edge.color[5], ifelse((edges_df$Pvalue < 
                                                                                                                                                                                                                                                                                               0.05) & (edges_df$OR < 0), edge.color[6], NA))))))
  cellall <- unique(edges_df[, 1], edges_df[, 2])
  g <- igraph::graph_from_data_frame(d = edges_df, directed = TRUE)
  igraph::E(g)$width <- edges_df$weight
  igraph::E(g)$color <- edges_df$color
  
  
  best.info$celltype <- factor(best.info$celltype, levels = names(V(g)))
  best.info <- best.info[order(best.info$celltype), ]
  
  
  igraph::V(g)$color <- ifelse((best.info$Pvalue < 0.001) & 
                                 (best.info$Beta > 0), cell.color[1], ifelse((best.info$Pvalue < 
                                                                                0.01) & (best.info$Beta > 0), cell.color[2], ifelse((best.info$Pvalue < 
                                                                                                                                       0.05) & (best.info$Beta > 0), cell.color[3], ifelse((best.info$Pvalue < 
                                                                                                                                                                                              0.001) & (best.info$Beta < 0), cell.color[7], ifelse((best.info$Pvalue < 
                                                                                                                                                                                                                                                      0.01) & (best.info$Beta < 0), cell.color[6], ifelse((best.info$Pvalue < 
                                                                                                                                                                                                                                                                                                             0.05) & (best.info$Beta < 0), cell.color[5], cell.color[4]))))))
  igraph::V(g)$size <- ifelse(best.info$Pvalue < 0.001, cell.size[1], 
                              ifelse(best.info$Pvalue < 0.01, cell.size[2], ifelse(best.info$Pvalue < 
                                                                                     0.05, cell.size[3], cell.size[4])))
  igraph::V(g)$shape <- "circle"
  igraph::V(g)$label.cex <- 1.2
  igraph::V(g)$label.color <- "black"
    igraph::V(g)$frame.color <- "white"
      p1 <- plot(g, edge.width = igraph::E(g)$width, edge.color = adjustcolor(igraph::E(g)$color, 
                                                                              alpha.f = 0.7), edge.arrow.size = arrow.size, edge.arrow.width = arrow.width, 
                 vertex.label.cex = igraph::V(g)$label.cex, vertex.label.color = igraph::V(g)$label.color, 
                 vertex.size = igraph::V(g)$size, vertex.color = adjustcolor(igraph::V(g)$color, 
                                                                             alpha.f = 0.9), vertex.shape = igraph::V(g)$shape, 
                 vertex.frame.color = igraph::V(g)$frame.color, main = "Cell phenotype network plot", 
                 layout = igraph::layout_with_fr)
      
      if(Interation.effect){
        
        graphics::legend("topright", legend = c("Positive Moderation Effect (p < 0.001)", 
                                                "Positive Moderation Effect (p < 0.01)", "Positive Moderation Effect (p < 0.05)", 
                                                "Negative Moderation Effect (p < 0.05)", "Negative Moderation Effect (p < 0.01)", 
                                                "Negative Moderation Effect (p < 0.001)"), col = edge.color, 
                         lty = 1, lwd = 2, cex = 0.8, bg = "#e8e9ec")
        
        
      }else{
        
        
        graphics::legend("topright", legend = c("Synergistic effect (p < 0.001)", 
                                                "Synergistic effect (p < 0.01)", "Synergistic effect (p < 0.05)", 
                                                "Antagonistic effect (p < 0.05)", "Antagonistic effect (p < 0.01)", 
                                                "Antagonistic effect (p < 0.001)"), col = edge.color, 
                         lty = 1, lwd = 2, cex = 0.8, bg = "#e8e9ec")
        
        
      }
      
      graphics::legend("bottomright", legend = c("Cell rate ↑↑↑ (p < 0.001)", 
                                                 "Cell rate ↑↑ (p < 0.01)", "Cell rate ↑ (p < 0.05)", 
                                                 "No significant change", "Cell rate ↓ (p < 0.05)", 
                                                 "Cell rate ↓↓ (p < 0.01)", "Cell rate ↓↓↓ (p < 0.001)"), 
                       fill = cell.color, title = "Cell proportion changes with phenotype", 
                       cex = 0.8, bg = "#e8e8ea")
      return(g)
}
