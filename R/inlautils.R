## Find the indices into inla output config structures corresponding
## to a specific predictor, effect name, or inla.stack tag.
##
## result : an inla object
inla.output.indices = function(result, name=NULL, stack=NULL, tag=NULL, ...)
{
    if (!is.null(name) && !is.null(tag)) {
        stop("At most one of 'name' and 'tag' may be non-null.")
    }
    if (!is.null(tag)) {
        if (is.null(stack)) {
            stop("'tag' specified but 'stack' is NULL.")
        }
        tags <- names(stack$data$index)
        if (!(tag %in% tags)) {
            stop("'tag' not found in 'stack'.")
        }
    } else if (is.null(name) || (!is.null(name) && (name == ""))) {
        if ("APredictor" %in% result$misc$configs$contents$tag) {
            name <- "APredictor"
        } else {
            name <- "Predictor"
        }
    }

    ## Find variables
    if (!is.null(name)) {
        nameindex <- which(result$misc$configs$contents$tag == name)
        index <- (result$misc$configs$contents$start[nameindex] - 1L +
                  seq_len(result$misc$configs$contents$lengthstart[nameindex]))
    } else { ## Have tag
        index <- which(stack$data$index[[tag]])
    }

    index
}
