library(ggplot2)
library(reshape2)
library(gridExtra)


# function to overlay factor with positions 
plotFactPosition = function(object, fact='all', position_file=NA, separator=' ', clusters=NA){
  
  if ("SigmaZ" %in% names(object@Expectations)) {
    
    # get factors from model 
    f_sfa = getFactors(object)
    colnames(f_sfa) = paste('fact_',colnames(f_sfa), sep='')
    f_sfa = data.frame(f_sfa)
    
    # add length scale and significance to the factor names 
    # TODO check first that Parameters$Sigma$l is an attribute (or try catch)
    
    l_scales = paste('l=', round(object@Parameters$SigmaZ$l, 1), sep='')
    sig = paste('sig=', round(object@Parameters$SigmaZ$sig, 1), sep='')
    colnames(f_sfa) = paste(colnames(f_sfa), l_scales, sep=',  ')
    colnames(f_sfa) = paste(colnames(f_sfa), sig, sep=',  ')

    # load and cbind positions
    # TODO, only do this if there is no position file as an input, otherwise use that (for MOFA res)
    if(is.na(position_file)){
      positions = object@Parameters$SigmaZ$X
    }
    else{
      # positions = read.table(position_file, sep = separator)
      positions = read.table(position_file, sep = separator)
    }
    
    colnames(positions) = c('x', 'y')
    f_sfa = cbind(f_sfa, positions)
    
    if(!is.na(clusters)){
      clust = read.table(clusters)
      colnames(clust) = c('cluster')
      f_sfa = cbind(f_sfa, clust)
    }
    
    # all factors case
    if(fact=='all'){
      if(is.na(clusters)){
        f_sfa = melt(f_sfa, id.vars=c('x', 'y'))
        colnames(f_sfa) = c('x', 'y','factor', 'value')
        p = ggplot(f_sfa, aes(x=x, y=y))+
          geom_point(aes(fill=value), size=4, colour='black', pch=21)+
          theme_bw(base_size=15)+
          scale_fill_gradient2(low = '#004D7F', high='#B51700', mid='white')+
          facet_wrap(~factor)
        print(p)
      }
      else{
        f_sfa = melt(f_sfa, id.vars=c('x', 'y', 'cluster'))
        colnames(f_sfa) = c('x', 'y','cluster','factor', 'value')
        p = ggplot(f_sfa, aes(x=x, y=y))+
          geom_point(aes(fill=value), size=4, colour='black', pch=21)+
          theme_bw(base_size=15)+
          scale_fill_gradient2(low = '#004D7F', high='#B51700', mid='white')+
          facet_grid(cluster~factor)
        print(p)
      }
    }
   
  }
  
  if ("SigmaAlphaW" %in% names(object@Parameters)) {
    
    # get factors from model 
    f_sfa = getWeights(object)
    
    # add length scale and significance to the factor names 
    # TODO check first that Parameters$Sigma$l is an attribute (or try catch)
    
    views <- viewNames(object)
    
    plots <- list()
    
    for (m in views){ 
      colnames(f_sfa[[m]]) = paste('factors_',colnames(f_sfa[[m]]), sep='')
      f_sfa[[m]]= data.frame(f_sfa[[m]])

      l_scales = paste('l=', round(object@Parameters$SigmaAlphaW[[m]]$l, 1), sep='')
      sig = paste('sig=', round(object@Parameters$SigmaAlphaW[[m]]$sig, 1), sep='')
      colnames(f_sfa[[m]]) = paste(colnames(f_sfa[[m]]), l_scales, sep=',  ')
      colnames(f_sfa[[m]]) = paste(colnames(f_sfa[[m]]), sig, sep=',  ')
      
      # load and cbind positions
      # TODO, only do this if there is no position file as an input, otherwise use that (for MOFA res)
      if(is.na(position_file)){
        positions = object@Parameters$SigmaAlphaW[[m]]$X
      }
      else{
        # positions = read.table(position_file, sep = separator)
        positions = read.table(position_file[[m]], sep = separator)
      }
      
      colnames(positions) = c('x', 'y')
      f_sfa[[m]] = cbind(f_sfa[[m]], positions)
      
      if(!is.na(clusters)){
        clust = read.table(clusters[[m]])
        colnames(clust) = c('cluster')
        f_sfa[[m]] = cbind(f_sfa[[m]], clust)
      }
    
      # all factors case
      if(fact=='all'){
        if(is.na(clusters)){
          f_sfa[[m]] = melt(f_sfa[[m]], id.vars=c('x', 'y'))
          colnames(f_sfa[[m]]) = c('x', 'y','factor', 'value')
          p = ggplot(f_sfa[[m]], aes(x=x, y=y))+
            geom_point(aes(fill=value), size=4, colour='black', pch=21)+
            theme_bw(base_size=15)+
            scale_fill_gradient2(low = '#004D7F', high='#B51700', mid='white')+
            facet_wrap(~factor)+
            ggtitle(paste(m))
          plots[[m]] <- p
        }
        else{
          f_sfa[[m]] = melt(f_sfa[[m]], id.vars=c('x', 'y', 'cluster'))
          colnames(f_sfa[[m]]) = c('x', 'y','cluster','factor', 'value')
          p = ggplot(f_sfa[[m]], aes(x=x, y=y))+
            geom_point(aes(fill=value), size=4, colour='black', pch=21)+
            theme_bw(base_size=15)+
            scale_fill_gradient2(low = '#004D7F', high='#B51700', mid='white')+
            facet_grid(cluster~factor)+
            ggtitle(paste(m))
          plots[[m]] <- p
        }
      }
      
    }
    
    nCol <- floor(sqrt(length(plots)))
    do.call("grid.arrange", c(plots, ncol=nCol))

  }

}