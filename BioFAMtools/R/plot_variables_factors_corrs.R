#' @title Plot weights and correlations bewteen factors and expressions of the variables in each view for both groups. 
#' @name plot_variables_factors_corrs
#' @param variables : list of features names of the biofam model, or named vector where names are names of variables (covariates) and values are named vectors storing the values of the variables for each sample
#' @param viz_weight : option relevant only if variables are features names of the biofam model. If true, plot weights and not correlation across all groups in the first figures

library(gridExtra)
library(grid)

plot_variables_factors_corrs <- function(object, variables, viz_weight=FALSE) {

  groups = groups_names(object)
  views = views_names(object)
  factors = factors_names(object)
  
  N_groups = length(groups)
  N_variables = length(variables)

  list_plots = list()
  
  if (class(variables[[1]])=="character"){
    names_variables = variables
  } else {
    names_variables = names(variables)
  }
  
  idx_variable=0
  
  for (variable in variables){
    
    idx_variable = idx_variable+1
    name_variable = names_variables[[idx_variable]]
    
    #removing sample with missing values from variable if it is a named vector
    variable <- variable[names(which(!(is.na(variable))))]
    
    #first plot 
    idx_column = 1
    
    #first plot : correlation across all groups (or weights if viz_weigth is True and if variables are biofam features)
    factors_values = c()
    views_values = c()
    corrs = c()
    for (factor in factors){
      for (view in views){
        if (class(variable)=="character"){ #biofam feature
          if (name_variable %in% features_names(object)[[view]]){
            factors_values = c(factors_values,factor)
            views_values = c(views_values,view)
            if (viz_weight){
              corrs = c(corrs,object@expectations$W[[view]][name_variable,factor]) 
            }
            else{
              Z <- Reduce(rbind, object@expectations$Z)[,factor]
              Y <- do.call(c,lapply(t(object@training_data[[view]]),function(Y_group) Y_group[name_variable,]))
              common_samples <- intersect(names(which(!(is.na(Y)))),names(Z)) #removing samples with missing values
              corrs = c(corrs,cor(Z[common_samples],Y[common_samples]))
            }
          }
        }
        else{ #covariate
          factors_values = c(factors_values,factor)
          views_values = c(views_values,view)
          Z <- Reduce(rbind, object@expectations$Z)[,factor]
          corrs=c(corrs,cor(Z[names(variable)],variable))
        }
      }
    }
    df_corrs = data.frame(factors_values,views_values,corrs)
    if ((class(variable)=="character")&(viz_weight)){
      name_value = "weight"  
    } else{
      name_value = "corr"
    }
    colnames(df_corrs) =  c("factor", "view", name_value)
    df_corrs$factor <- factor(df_corrs$factor,levels=seq(1,max(as.integer(levels(df_corrs$factor)))))
    plotCor = ggplot(df_corrs, aes_string(x="factor", y=name_value, fill="view")) +
      geom_bar(position="dodge", stat="identity")
    
    list_plots[[(idx_variable-1)*(N_groups+1)+idx_column]] = plotCor
    
    #following plots : correlations within each group
    
    for (group in groups){
      factors_values = c()
      views_values = c()
      corrs = c()
      for (factor in factors){
        for (view in views){
          if (class(variable)=="character"){ #biofam feature
            if (name_variable %in% features_names(object)[[view]]){
              factors_values = c(factors_values,factor)
              views_values = c(views_values,view)
              Z <- object@expectations$Z[[group]][,factor]
              Y <- object@training_data[[view]][[group]][name_variable,]
              common_samples <- intersect(names(which(!(is.na(Y)))),names(Z)) #removing samples with missing values
              corrs=c(corrs,cor(Z[common_samples],Y[common_samples]))
            }
          }
          else{ #covariate
            factors_values = c(factors_values,factor)
            views_values = c(views_values,view)
            Z <- object@expectations$Z[[group]][,factor]
            corrs=c(corrs,cor(Z[names(variable)],variable))
          }
        }
      }
      df_corrs = data.frame(factors_values,views_values,corrs)
      name_value = "corr" 
      colnames(df_corrs) =  c("factor", "view", name_value)
      df_corrs$factor <- factor(df_corrs$factor,levels=seq(1,max(as.integer(levels(df_corrs$factor)))))
      plotCor = ggplot(df_corrs, aes_string(x="factor", y=name_value, fill="view")) +
        geom_bar(position="dodge", stat="identity")
      
      idx_column = idx_column + 1
      list_plots[[(idx_variable-1)*(N_groups+1)+idx_column]] = plotCor
    }
    
  }
  
  # Create row and column titles
  row.titles = names_variables
  col.titles = c("weight")
  for (group in groups){
    col.titles = c(col.titles,paste("corr_within_",group,sep=""))
  }
  
  list_plots[1:N_groups+1] = lapply(1:N_groups+1, function(i) arrangeGrob(list_plots[[i]], top=col.titles[[i]]))
  list_idx = c()
  for (k in seq(1,N_variables)){
    list_idx = c(list_idx,1+(k-1)*(N_groups+1))
  }
  list_plots[list_idx] = lapply(1:N_variables, function(i) arrangeGrob(list_plots[[(i-1)*(N_groups+1)+1]], left=row.titles[[i]]))

  do.call("grid.arrange", c(list_plots,ncol= N_groups+1,nrow=N_variables))    
  
}