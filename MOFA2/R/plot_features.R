###################################################################
## Functions to visualise features and their values or loadings  ##
###################################################################

#' @title Barchart of feature loadings
#' @name plot_weights_signature
#' @description Barchar of the feature loadings.
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plot_weights_signature <- function(object, factors = "all", n_top = 10, ...) {
  
  W <- get_weights(object, factors = factors, as.data.frame = TRUE)

  top_features <- W %>% group_by(factor) %>% top_n(n = n_top, wt = abs(value)) %>% pull(feature) %>% unique

  df <- W %>% filter(feature %in% top_features)
  df$factor <- factor(df$factor, levels = unique(df$factor))

  df <- df[order(-abs(df$value)),]
  df$feature <- factor(df$feature, levels = unique(df$feature))


  df %>% 
    ggplot(., aes(x = feature, y = value)) + 
      geom_bar(stat = 'identity') + 
      facet_grid(factor ~ view, ...) +
      ylab("Feature loading") +
      theme(
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.y = element_text(size=18),
          axis.text.y = element_text(size=12),
          legend.text=element_text(size=14),

          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(angle = 0)
      )
}


#' @title Barchart of feature values
#' @name plot_features
#' @description Barchar of the feature values
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plot_features <- function(object, views = "all", groups = "all", factors = "all", 
              n_top = 10, ...) {
  
  W <- get_weights(object, factors = factors, as.data.frame = TRUE)

  top_features <- W %>% group_by(factor) %>% top_n(n = n_top, wt = abs(value)) %>% pull(feature) %>% unique

  # Select relevant features in data
  data <- get_data(object, views = views, groups = groups, as.data.frame = TRUE)
  # data <- lapply(data, function(m) {
  #   lapply(m, function(g) {
  #     g[rownames(g) %in% top_features,]
  #   })
  # })
  data <- data[data$feature %in% top_features,]

  data <- data[order(-abs(data$value)),]
  data$sample <- factor(data$sample, levels = unique(data$sample))

  message(paste0("Rendering ", length(top_features), " features across the dataset, this might take some time..."))
  data %>% 
    ggplot(., aes(x = sample, y = value)) + 
      geom_bar(aes(fill = group), stat = 'identity') + 
      facet_grid(feature ~ group, ...) +
      theme(
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),

          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),

          # Remove x axis ticks and labels
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),

          axis.title.y = element_text(size=18),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size=14),

          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(angle = 0)
      )
}
