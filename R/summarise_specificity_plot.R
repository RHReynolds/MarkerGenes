#' Summarised specificity plot
#'
#' This function visually summarises the specificity values outputted by
#' \code{query_gene_ctd}, across each ctd study.
#'
#' @param filtered_specificity A dataframe that contains the columns 'CellType',
#'   'Specificity', and 'Study', as typically outputted by
#'   \code{query_gene_ctd}. Dataframe can contain additional columns, but must
#'   just contain aforementioned.
#' @param flag.variable Default is NULL. Fill this with a column name within
#'   your dataframe if you wish to colour by certain groups within your table.
#'
#' @return Boxplot of specificity values of queried gene list across cell types
#'   in each study.
#' @export
#'

summarise_specificity_plot <- function(filtered_specificity, fill.variable = NULL){

  if(is.null(fill.variable)){
    plot <- ggplot(data = filtered_specificity, aes(x = MarkerGenes::reorder_within(x = CellType,
                                                                                    by = Specificity,
                                                                                    within = Study,
                                                                                    fun = median,
                                                                                    desc = TRUE),
                                                    y = Specificity)
    ) +
      geom_boxplot() +
      MarkerGenes::scale_x_reordered() +
      facet_wrap(vars(Species, Study), scales = "free_x") +
      labs(x = "Cell type", y = "Specificity", title = "") +
      scale_y_continuous(limits = c(0,1)) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
            axis.title = element_text(size = 6, face = "bold"),
            strip.text = element_text(size = 6, face = "bold"),
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.margin = unit(c(t = -1, r = 0.5, b = 0.5, l = 0.5), "lines")
      )
  } else{
  plot <- ggplot(data = filtered_specificity, aes(x = MarkerGenes::reorder_within(x = CellType,
                                                                                  by = Specificity,
                                                                                  within = Study,
                                                                                  fun = median,
                                                                                  desc = TRUE),
                                                  y = Specificity,
                                                  fill = filtered_specificity %>% .[[fill.variable]])
  ) +
    geom_boxplot() +
    MarkerGenes::scale_x_reordered() +
    facet_wrap(vars(Species, Study), scales = "free_x") +
    labs(x = "Cell type", y = "Specificity", title = "", fill = fill.variable) +
    scale_y_continuous(limits = c(0,1)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
          axis.title = element_text(size = 10, face = "bold"),
          strip.text = element_text(size = 8, face = "bold"),
          legend.position = "top",
          legend.text = element_text(size = 6),
          legend.key.size = unit(1,"line"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(t = -1, r = 0.5, b = 0.5, l = 0.5), "lines")
    )

  }

  return(plot)

}
