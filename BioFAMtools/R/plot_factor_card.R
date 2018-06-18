
###########################################
## Functions to visualise latent factors ##
###########################################

#' @title Plot summary of informative plots relative to a specific factor
#' @name plot_factor_card

library(gridExtra)
library(grid)

plot_factor_card <- function(object, factor, groupwise=FALSE, fseas=NULL, gene_sets=NULL, color_by="group") {
  
  #TODO plot2 : add var explained in groups and views (bar plor : x = group, y = var explained, hue = view)
  #TODO possible : add quantitative sparsity

  views = views_names(object)
  N_views = length(views)
  
  #plot 1 : variance explained
  #TODO : bar plor : x = group, y = var explained, hue = view
  
  var_res = calculate_variance_explained(object, views = "all", groups = "all", factors = factor, groupwise=groupwise)
  groups_values = c()
  views_values = c()
  var_explained = c()
  for (group in names(var_res[["r2_total"]])){
    for (view in names(var_res[["r2_total"]][[group]])){
      groups_values = c(groups_values,group)
      views_values = c(views_values,view)
      var_explained=c(var_explained,var_res[["r2_total"]][[group]][[view]])
    }
  }
  df_var_res = data.frame(groups_values,views_values,var_explained)
  colnames(df_var_res) =  c("group", "view", "variance_explained")
  plotVar = ggplot(df_var_res, aes(x=group, y=var_explained, fill=view)) +
    geom_bar(position="dodge", stat="identity")
  
  #plot 2 : expression across samples
  
  plotZ = plot_factor_beeswarm(object, factors = factor, color_by = color_by, superimpose_groups = FALSE)
  
  plotViews = list()
  
  if(!(is.null(gene_sets))){
    gene_in_gsea = colnames(gene_sets)
  }
  
  idx=1

  for (view in views){
    
    #plot 3 : top weights 
    tmp = list()
    tmp[[1]] = plot_weights(object, view, factor, scale = F)

    #plot 4 : top pathways
    
    if(!(is.null(fseas))){
      
      if (length(fseas)!=length(views)){ #no pathways found
        print("fsea should contain as many FSEA objects as there are views")
        tmp[[2]] = grid.rect(gp=gpar(col="white"))
      }
      else{
        tmp[[2]] = lineplot_FSEA(fseas[[idx]], factor, threshold=0.05, max.pathways=25)
        if (length(tmp)==1){
          tmp[[2]] = grid.rect(gp=gpar(col="white"))
        }
      }
      
    }
    
    else{
      if(!(is.null(gene_sets))){
        
        gene_in_biofam = features_names(object)[[view]] #if transpose: sampleNames(best_model)
        common_genes = intersect(gene_in_gsea,gene_in_biofam)
        gene_sets_subset = gene_sets[,common_genes]
        
        metric_expr = "loading"
        statistic = "mean.diff" # "rank.sum" 
        sig_test = "parametric" #"cor.adj.parametric" # "permutation" # 
        
        #fsea=FSEA(object, view, data.matrix(gene_sets_subset),factors=factor,local.statistic=metric_expr,global.statistic=statistic,statistical.test=sig_test)
        #does not work for some reason
        fsea=FSEA(object, view, data.matrix(gene_sets_subset), local.statistic=metric_expr,global.statistic=statistic,statistical.test=sig_test)
        tmp[[2]]  = lineplot_FSEA(fsea, factor, threshold=0.05, max.pathways=25)
        if (length(tmp)==1){ #no pathways found
          tmp[[2]] = grid.rect(gp=gpar(col="white"))
        }
        
      }
      else{ 
        tmp[[2]] = grid.rect(gp=gpar(col="white"))
      }
    }
    
    plotViews[[idx]] = tmp
    idx = idx+1

  }

  if (N_views==1){
    list_plots = list(plotVar,plotViews[[1]][[1]],plotZ,plotViews[[1]][[2]])
    do.call("grid.arrange", c(list_plots,ncol=2,nrow=2))
  }
  else{
    if  (N_views==2){
      list_plots = list(plotVar,plotViews[[1]][[1]],plotViews[[1]][[2]])
      list_plots = append(list_plots,plotZ,plotViews[[2]][[1]],plotViews[[2]][[2]])
      do.call("grid.arrange", c(list_plots,ncol=3,nrow=2))
    }
    else{
      list_plots = list(plotVar,plotViews[[1]][[1]],plotViews[[1]][[2]])
      list_plots = append(list_plots,plotZ,plotViews[[2]][[1]],plotViews[[2]][[2]])
      idx=3
      for (view in views[seq(3,N_views)]){
        blankPlot = grid.rect(gp=gpar(col="white"))
        list_plots = append(list_plots,blankPlot,plotViews[idx][[1]],plotViews[idx][[2]])
        idx= idx+1
      }
      do.call("grid.arrange", c(list_plots,ncol=3,nrow=N_views))    
    }
  }
  
}