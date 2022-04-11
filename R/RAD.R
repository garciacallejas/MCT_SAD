#' obtain rank-abundance distributions from numeric vectors
#'
#' @param abundances either a (possibly named) numeric vector or a list of them
#'
#' @return either a dataframe or a list of them
#' @export
#'
#' @examples 
#' RAD(seq(1:20))
RAD <- function(abundances = NULL){
    if(inherits(abundances,"list")){
        
        rank.abundances <- list()
        
        for(i in 1:length(abundances)){
            if(!is.null(names(abundances[[i]]))){
                abund.df <- data.frame(species = names(abundances[[i]]),
                                       abundance = abundances[[i]])
            }else{
                abund.df <- data.frame(species = paste("sp",sprintf("%02d", 1:length(abundances[[i]])),sep=""),
                                       abundance = abundances[[i]])
            }
            
            rank.abundances[[i]] <- abund.df %>%
                mutate(relative.abundance = abundance/sum(abundance)) %>%
                mutate(species.rank = rank(-relative.abundance,ties.method = "first"))
        }
        
    }else{
        if(!is.null(names(abundances))){
            abund.df <- data.frame(species = names(abundances),abundance = abundances)
        }else{
            abund.df <- data.frame(species = paste("sp",sprintf("%02d", 1:length(abundances)),sep=""),
                                   abundance = abundances)
        }
        
        rank.abundances <- abund.df %>%
            mutate(relative.abundance = abundance/sum(abundance)) %>%
            mutate(species.rank = rank(-relative.abundance,ties.method = "first"))
    }
    
    return(rank.abundances)
}