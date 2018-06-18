#' @title Plot weights and correlations bewteen factors and expressions of the features in each view for both groups. 
#' @name plot_features_summary

library(gridExtra)
library(grid)

plot_features_summary <- function(object, features) {

  groups = groups_names(object)
  views = views_names(object)
  factors = factors_names(object)
  
  N_groups = length(groups)
  N_features = length(features)

  list_plots = list()
  
  idx_feature=0
  
  for (feature in features){
    
    idx_feature= idx_feature+1
    idx_column = 1
  
    #first plot : weights
    factors_values = c()
    views_values = c()
    weights = c()
    for (factor in factors){
      for (view in views){
        if (feature %in% features_names(object)[[view]]){
          factors_values = c(factors_values,factor)
          views_values = c(views_values,view)
          weights=c(weights,object@expectations$W[[view]][feature,factor])
        }
      }
    }
    df_weights = data.frame(factors_values,views_values,weights)
    name_value = "weight" #paste(feature,"_weight",sep="")
    colnames(df_weights) =  c("factor", "view", name_value)
    plotW = ggplot(df_weights, aes_string(x="factor", y=name_value, fill="view")) +
      geom_bar(position="dodge", stat="identity")
    list_plots[[(idx_feature-1)*(N_groups+1)+idx_column]] = plotW
    
    #following plots : correlations 
    
    for (group in groups){
      factors_values = c()
      views_values = c()
      corrs = c()
      for (factor in factors){
        for (view in views){
          if (feature %in% features_names(object)[[view]]){
            factors_values = c(factors_values,factor)
            views_values = c(views_values,view)
            corrs=c(corrs,cor(object@expectations$Z[[group]][,factor],object@training_data[[view]][[group]][feature,]))
          }
        }
      }
      df_corrs = data.frame(factors_values,views_values,corrs)
      name_value = "corr" # paste("corr_",feature,"_expression_with_factor_within_",group,sep="")
      colnames(df_corrs) =  c("factor", "view", name_value)
      plotCor = ggplot(df_corrs, aes_string(x="factor", y=name_value, fill="view")) +
        geom_bar(position="dodge", stat="identity")
      
      idx_column = idx_column + 1
      list_plots[[(idx_feature-1)*(N_groups+1)+idx_column]] = plotCor
    }
    
  }
  
  # Create row and column titles
  row.titles = features
  col.titles = c("weight")
  for (group in groups){
    col.titles = c(col.titles,paste("corr_within_",group,sep=""))
  }
  list_plots[1:N_groups+1] = lapply(1:N_groups+1, function(i) arrangeGrob(list_plots[[i]], top=col.titles[[i]]))
  list_idx = c()
  for (k in seq(1,N_features)){
    list_idx = c(list_idx,1+(k-1)*(N_groups+1))
  }
  list_plots[list_idx] = lapply(1:N_features, function(i) arrangeGrob(list_plots[[(i-1)*(N_groups+1)+1]], left=row.titles[[i]]))

  do.call("grid.arrange", c(list_plots,ncol= N_groups+1,nrow=N_features))    
  
}