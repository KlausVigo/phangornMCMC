#' Dirichlet distribution
#'
#' Density and random generation for the Dirichlet distribution with mean equal
#' to mean and standard deviation equal to sd.
#'
#' The Dirichlet is the multidimensional of the bet distribution and the
#' canonical Bayesian Distribution for the parameters of a multinomial
#' distribution. \code{rdirichlet1(alpha)} is a shortcut for
#' \code{rdirichlet(1, alpha)}
#'
#' @param x A vector containing a single random deviate or matrix containg one
#' random deviate per row.
#' @param n Number of random vectors to generate.
#' @param alpha Vector or (for ddirichlet) matrix containing shape parameters.
#' @author Code original posted by Ben Bolker to R-News on Fri Dec 15 2000. See
#' https://stat.ethz.ch/pipermail/r-help/2000-December/009561.html. Ben
#' attributed the code to Ian Wilson i.wilson@maths.abdn.ac.uk. Subsequent
#' modifications by Gregory R. Warnes greg@warnes.net.
#'
#' @seealso \code{\link{dbeta}}, \code{\link{rbeta}}
#' @keywords cluster ~kwd2
#' @examples
#'
#' set.seed(42)
#' rdirichlet(5, c(1,1,1) )
#' set.seed(42)
#' rdirichlet1(c(1,1,1))
#'
#' @rdname rdirichlet
#' @export
rdirichlet <- function (n, alpha)
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
}


#' @rdname rdirichlet
#' @export
rdirichlet1 <- function (alpha)
{
    l <- length(alpha)
    x <- rgamma(l, alpha)
    x/sum(x)
}


#' @rdname rdirichlet
#' @export
ddirichlet <- function (x, alpha)
{
    dirichlet1 <- function(x, alpha) {
        logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
        s <- (alpha - 1) * log(x)
        s <- ifelse(alpha == 1 & x == 0, -Inf, s)
        exp(sum(s) - logD)
    }
    if (!is.matrix(x))
        if (is.data.frame(x))
            x <- as.matrix(x)
        else x <- t(x)
        if (!is.matrix(alpha))
            alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x),
                            byrow = TRUE)
        if (any(dim(x) != dim(alpha)))
            stop("Mismatch between dimensions of 'x' and 'alpha'.")
        pd <- vector(length = nrow(x))
        for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i,
                                                               ])
        pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
        pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
        pd
}
