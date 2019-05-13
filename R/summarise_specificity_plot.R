#' Summarised specificity plot
#'
#' This function visually summarises the specificity values outputted by
#' \code{query_gene_ctd}, across each ctd study.
#'
#' @param filtered_specificity
#'
#' @return Boxplot of specificity values of queried gene list across cell types
#'   in each study.
#' @export
#'

summarise_specificity_plot <- function(filtered_specificity){

  plot <- ggplot(data = filtered_specificity, aes(x = CellType,
                                                  y = Specificity)
  ) +
    geom_boxplot(aes()) +
    facet_wrap(vars(Species, Study), scales = "free_x") +
    labs(x = "Cell type", y = "Specificity", title = "") +
    scale_y_continuous(limits = c(0,1)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
          axis.title = element_text(size = 6, face = "bold"),
          strip.text = element_text(size = 6, face = "bold"),
          legend.position = "top",
          legend.text = element_text(size = 6),
          legend.key.size = unit(1,"line"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(t = -1, r = 0.5, b = 0.5, l = 0.5), "lines")
    )

  return(plot)

}
