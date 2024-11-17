draw_pca <- function (exp, group_list, color = c("#2874C5", "#f87669", "#e6b707",
                                              "#868686", "#92C5DE", "#F4A582", "#66C2A5", "#FC8D62", "#8DA0CB",
                                              "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
                                              addEllipses = TRUE, style = "default", color.label = "Group",
          title = "", ...)
{

  if (!requireNamespace("factoextra", quietly = TRUE)) {
    stop("Package \"factoextra\" needed for this function to work.\n         Please install it by install.packages('factoextra')",
         call. = FALSE)
  }
  p1 <- all(apply(exp, 2, is.numeric))
  if (!p1)
    stop("exp must be a numeric matrix")
  p2 <- (sum(!duplicated(group_list)) > 1)
  if (!p2)
    stop("group_list must more than 1")
  dat <- as.data.frame(t(exp))
  dat.pca <- FactoMineR::PCA(dat, graph = FALSE)
  col = color[1:length(levels(group_list))]
  if (style == "default") {
    factoextra::fviz_pca_ind(dat.pca, geom.ind = "point",
                             col.ind = group_list, addEllipses = addEllipses,
                             palette = col, legend.title = "Groups", title = title,
                             ...)
  }
  else if (style == "ggplot2") {
    pdat = data.frame(dat.pca[["ind"]][["coord"]], Group = group_list)
    p = ggplot(pdat, aes(Dim.1, Dim.2)) + geom_point(aes(Dim.1,
                                                         Dim.2, fill = Group), shape = 21, color = "black") +
      scale_color_manual(values = color[1:nlevels(group_list)]) +
      scale_fill_manual(values = color[1:nlevels(group_list)]) +
      theme_classic() + theme(legend.position = "top") +
      labs(color = color.label, fill = color.label, title = title)
    if (addEllipses)
      p = p + stat_ellipse(aes(color = Group, fill = Group),
                           geom = "polygon", alpha = 0.3, linetype = 2)
    return(p)
  }
  else if (style == "3D") {
    colors = color[as.numeric(group_list)]
    pdat = data.frame(dat.pca[["ind"]][["coord"]], Group = group_list)
    scatterplot3d::scatterplot3d(pdat[, 1:3], color = "black",
                                 pch = 21, bg = colors, main = title)
    graphics::legend("bottom", col = "black", legend = levels(group_list),
                     pt.bg = color[1:nlevels(group_list)], pch = 21,
                     inset = -0.2, xpd = TRUE, horiz = TRUE)
  }
}
