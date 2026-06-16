mkNNIndx <- function(coords, m, n_omp_threads = 1){
    
    if(!is.matrix(coords))
        stop("coords must be a matrix")
    
    storage.mode(coords) <- "double"
    
    n <- nrow(coords)
    r <- ncol(coords)

    if(m > n)
        stop("m must be less than or equal to the number of rows in coords")
    
    .Call(
        "mkNNIndx",
        coords,
        as.integer(n),
        as.integer(m),
        as.integer(r),
        as.integer(n_omp_threads),
        PACKAGE = "stLMM"
    )
}

mkNNIndxBrute <- function(coords, m, n_omp_threads = 1){

    if(!is.matrix(coords))
        stop("coords must be a matrix")

    storage.mode(coords) <- "double"

    n <- nrow(coords)
    r <- ncol(coords)

    if(m > n)
        stop("m must be less than or equal to the number of rows in coords")

    .Call(
        "mkNNIndxBrute",
        coords,
        as.integer(n),
        as.integer(m),
        as.integer(r),
        as.integer(n_omp_threads),
        PACKAGE = "stLMM"
    )
}
