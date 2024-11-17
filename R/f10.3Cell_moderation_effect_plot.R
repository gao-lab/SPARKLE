#' Cell Moderation Effect Plot
#'
#' This function generates a plot to visualize the cell moderation effect based on the given interaction information.
#'
#' @param interaction.best.info A data frame containing information about cell interactions.
#' @param cell.color The color of the cells in the plot. Default is "gray".
#' @param edge.width The width of the edges in the plot. Default is 3.
#' @param arrow.size The size of the arrows indicating edge directions. Default is 0.75.
#' @param arrow.width The width of the arrows indicating edge directions. Default is 1.5.
#' @param edge.color A vector containing colors for different significance levels of Pvalues.
#'                   Default is c("darkred", "black", "gray").
#'
#' @return A plot visualizing the cell moderation effect.
#'
#' @examples
#' Cell_moderation_effect_plot(interaction.best.info)
#'
#' @export
#'
#' @import igraph
#'
Cell_moderation_effect_plot <- function(interaction.best,cell.color="gray",edge.width=3,arrow.size=0.75,arrow.width=1.5,edge.color=c("darkred", "black", "gray")){



  interaction.best.info <- interaction.best[["Chosen_model_info"]]

  # 分割 celltype 列以获取节点
  edges <- strsplit(as.character(interaction.best.info$celltype), " ")
  edges <- do.call(rbind, edges)

  # 创建一个边数据框
  edges_df <- data.frame(from = edges[,1], to = edges[,2], Pvalue = interaction.best.info$Pvalue)

  # 根据 Pvalue 分配权重和颜色
  edges_df$weight <- ifelse(edges_df$Pvalue < 0.001, edge.width,
                            ifelse(edges_df$Pvalue < 0.01, edge.width,
                                   ifelse(edges_df$Pvalue < 0.05, edge.width, 0)))

  edges_df$color <- ifelse(edges_df$Pvalue < 0.001, edge.color[1],
                           ifelse(edges_df$Pvalue < 0.01, edge.color[2],
                                  ifelse(edges_df$Pvalue < 0.05, edge.color[3], NA)))

  # 过滤掉权重为0的边
  edges_df <- edges_df[edges_df$weight > 0, ]

  # 创建有向图对象
  g <- igraph::graph_from_data_frame(d = edges_df, directed = TRUE)

  # 设置边的宽度和颜色
  igraph::E(g)$width <- edges_df$weight
  igraph::E(g)$color <- edges_df$color
  igraph::V(g)$color <- cell.color
  igraph::V(g)$size <- degree(g, mode = "all") * 2 +5  # 节点大小与度相关
  igraph::V(g)$shape <- "circle"
  igraph::V(g)$label.cex <- 1.2
  igraph::V(g)$label.color <- "black"
    igraph::V(g)$frame.color <- "white"

    # 绘制网络图
    p1 <- plot(g,
               edge.width = igraph::E(g)$width,
               edge.color = adjustcolor(igraph::E(g)$color, alpha.f = 0.7),
               edge.arrow.size =  arrow.size,   # 增大箭头的大小
               edge.arrow.width = arrow.width,  # 增大箭头的宽度
               vertex.label.cex = igraph::V(g)$label.cex,
               vertex.label.color = igraph::V(g)$label.color,
               vertex.size = igraph::V(g)$size,
               vertex.color = adjustcolor(igraph::V(g)$color, alpha.f = 0.9),
               vertex.shape = igraph::V(g)$shape,
               vertex.frame.color = igraph::V(g)$frame.color,
               main = "Cell moderation effect plot",
               layout = layout_with_fr)  # 使用 Fruchterman-Reingold 布局算法

    # 添加图例
    graphics::legend("topright",
                     legend = c("Pvalue < 0.001", "Pvalue < 0.01", " Pvalue < 0.05"),
                     col =edge.color ,
                     lty = 1,
                     lwd = 2,
                     cex = 0.8,
                     bg ="white" )

}
