#' Summarised specificity plot
#'
#' This function visually summarises the mean specificity values outputted by
#' \code{query_gene_ctd}, across each ctd study.
#'
#' @param filtered_specificity A dataframe that contains the columns 'CellType',
#'   'Specificity', and 'Study', as typically outputted by
#'   \code{query_gene_ctd}. Dataframe can contain additional columns, but must
#'   just contain aforementioned.
#' @param fill.variable Default is NULL. Fill this with a column name within
#'   your dataframe if you wish to colour by certain groups within your table.
#'
#' @return Boxplot of mean specificity values of queried gene list across cell
#'   types in each study.
#' @export
#'

summarise_specificity_plot <- function(filtered_specificity, fill.variable = NULL){

  if(is.null(fill.variable)){
    plot <-
      ggplot2::ggplot(
        filtered_specificity,
        ggplot2::aes(
          x = MarkerGenes::reorder_within(x = CellType,
                                          by = Specificity,
                                          within = Study,
                                          fun = median,
                                          desc = TRUE),
          y = Specificity
        )
      ) +
      ggplot2::geom_boxplot() +
      MarkerGenes::scale_x_reordered() +
      ggplot2::facet_wrap(ggplot2::vars(Species, Study), scales = "free_x") +
      ggplot2::labs(x = "Cell type", y = "Specificity", title = "") +
      ggplot2::scale_y_continuous(limits = c(0,1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6),
        axis.text.x = ggplot2::element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        axis.title = ggplot2::element_text(size = 6, face = "bold"),
        strip.text = ggplot2::element_text(size = 6, face = "bold"),
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(t = -1, r = 0.5, b = 0.5, l = 0.5), "lines")
      )

  } else{

    plot <-
      ggplot2::ggplot(
        filtered_specificity,
        ggplot2::aes(
          x = MarkerGenes::reorder_within(x = CellType,
                                          by = Specificity,
                                          within = Study,
                                          fun = median,
                                          desc = TRUE),
          y = Specificity,
          fill = filtered_specificity %>% .[[fill.variable]]
        )
      ) +
      ggplot2::geom_boxplot() +
      MarkerGenes::scale_x_reordered() +
      ggplot2::facet_wrap(ggplot2::vars(Species, Study), scales = "free_x") +
      ggplot2::labs(x = "Cell type", y = "Specificity", title = "", fill = fill.variable) +
      ggplot2::scale_y_continuous(limits = c(0,1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6),
        axis.text.x = ggplot2::element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        axis.title = ggplot2::element_text(size = 10, face = "bold"),
        strip.text = ggplot2::element_text(size = 8, face = "bold"),
        legend.position = "top",
        legend.text = ggplot2::element_text(size = 6),
        legend.key.size = ggplot2::unit(1,"line"),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(t = -1, r = 0.5, b = 0.5, l = 0.5), "lines")
      )

  }

  return(plot)

}
