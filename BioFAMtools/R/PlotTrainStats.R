
###########################################
## Functions to plot training statistics ##
###########################################

#' @rdname trainCurveFactors
#' @name trainCurveFactors
#' @title Training curve for the number of active factors
#' @description the BioFAModel starts with an initial number of factors and inactive factors are dropped during training if they explain small amounts of variation. 
#' This allows the model to automatically infer the dimensionality of the latent space.
#' The corresponding hyperparameters are defined in \code{\link{prepareBioFAM}}. \cr
#' All training statistics, including the number of active factors, can be fetch from the TrainStats slot of \code{\link{BioFAModel}} .
#' @param object a \code{\link{BioFAModel}} object.
#' @import ggplot2 scales
#' @export

trainCurveFactors <- function(object) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") { stop("'object' has to be an instance of BioFAModel") }
  
  # Collect training statistics
  idx = seq(1+object@TrainOptions$startdrop,length(object@TrainStats$activeK),object@TrainOptions$freqdrop)
  stat = object@TrainStats$activeK[idx] 
  data <- data.frame(time=idx, value=stat)
  
  # Plot
  p <- ggplot(data, aes_string(x="time", y="value")) +
    geom_line() +
    labs(title="", x="Iteration", y="Number of active latent varaibles")
    # scale_y_discrete(limits=c(min(data$value)-1, max(data$value)+1)) +
    # scale_y_continuous(limits=c(min(data$value)-1, max(data$value)+1), breaks=scales::pretty_breaks()) +
    theme(
      axis.title.x=element_text(colour="black",size=rel(1.75), margin=margin(20,0,3,0)),
      axis.title.y=element_text(colour="black",size=rel(1.75), margin=margin(0,20,0,3)),
      axis.text.x=element_text(colour="black",size=rel(1.5)),
      axis.text.y=element_text(colour="black",size=rel(1.5)),
      axis.ticks.x = element_line(colour="black"),
      axis.ticks.y = element_line(colour="black"),
      axis.line.x = element_line(color="black"),
      axis.line.y = element_line(color="black"),
      legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
    
  return(p)
}


#' @title Training curve for Evidence Lower Bound (ELBO)
#' @name trainCurveELBO
#' @rdname trainCurveELBO
#' @param object a \code{\link{BioFAModel}} object.
#' @param log boolean indicating whether to apply log transform
#' @description BioFAM inference is done using the variational Bayes algorithm, which maximises a quantity called the Evidence Lower Bound (ELBO).
#' The ELBO is supposed to increase monotonically up to convergence, but it can decrease substantially when dropping inactive factors.
#' For more details read the supplementary methods
#' The frequency of ELBO computation as well as the convergence criteria are defined as hyperparameters in \code{\link{prepareBioFAM}}. \cr
#' All Training statistics, including the ELBO, can be fetch from the TrainStats slot of \code{\link{BioFAModel}} .
#' @import ggplot2 scales
#' @export

trainCurveELBO <- function(object, log = F) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") { stop("'object' has to be an instance of BioFAModel") }
  
  # Fetch ELBO from TrainStats  
  idx = seq(1,length(object@TrainStats$elbo),object@TrainOptions$elbofreq)
  stat = object@TrainStats$elbo[idx]
  
  # Apply log transform
  if (log==T) { stat <- -log(-stat) }
  
  data <- data.frame(time=idx, value=stat)
  p <- ggplot2::ggplot(data, aes_string(x="time", y="value")) +
    geom_line() +
    labs(title="", x="Iteration", y="ELBO")
    theme(
      axis.title.x=element_text(colour="black",size=rel(1.75), margin=margin(20,0,3,0)),
      axis.title.y=element_text(colour="black",size=rel(1.75), margin=margin(0,20,0,3)),
      axis.text.x=element_text(colour="black",size=rel(1.5)),
      axis.text.y=element_text(colour="black",size=rel(1.5)),
      axis.ticks.x = element_line(colour="black"),
      axis.ticks.y = element_line(colour="black"),
      axis.line.x = element_line(color="black"),
      axis.line.y = element_line(color="black"),
      legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  
  if (log==F) {
      scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x))) }
      p <- p + scale_y_continuous(labels=scientific_10)
  }
  
  return(p)
}

