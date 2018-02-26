### but you have to use #'


internal.frechet.coerce2matrix <- function(A)
{
    ### various ways of extracting a matrix of coordinates
    if (startsWith(class(A), "Spatial"))
    {
        M = coordinates(A); 
    }else{
        tryCatch(
            expr={M = as.matrix(A)},
            error=function(e){
                print(e)
                stop("Cannot extract coordinates: ")
            }
        )
        
    }
    

}


#' Calculate the Frechet Decision Using a Free Space Diagram
#'
#' @param T1 a valid trajectory
#' @param T2 a valid trajectory
#' @param eps the distance
#' @return A boolean (true if distance is smaller than epsilon
#' @family frechet.decision

frechet.decide <- function(t1,t2,eps, algorithm="duetschvahrenhold")
{
    M1 =internal.frechet.coerce2matrix(t1)
    M2 =internal.frechet.coerce2matrix(t2)
    if (algorithm == "duetschvahrenhold")
    {
        return (internal_dv_frechet_decide(M1,M2,eps));
    }
    if (algorithm == "bringmanbaldus")
    {
        return (internal_bb_frechet_decide(M1,M2,eps));
    }
    stop("Invalid algorithm. Only duetschvahrenhold and bringmanbaldus are valid decider implementations for now.");
}

