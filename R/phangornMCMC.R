## based on coalescentMCMC from Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.

#' phangornMCMC
#' This function runs a Markov chain Monte Carlo (MCMC) algorithm to generate a
#' set of trees which is returned with their likelihoods.
#'
#' @param x a set of DNA sequences, typically an object of class "DNAbin".
#' @param ntrees the number of trees to output.
#' @param burnin the number of trees to discard as “burn-in”.
#' @param frequency the frequency at which trees are sampled.
#' @param tree0 the initial tree of the chain; by default, a fastME / UPGMA tree
#' from a JC69 distance is generated.
#' @param model the transition model
#' @param printevery an integer specifying the frequency at which to print the
#' numbers of trees proposed and accepted; set to 0 to cancel all printings.
#' @param bf base frequencies
#' @param Q rate matrix
#'
#' @examples
#' \dontrun{
#' require(phangorn)
#' data(yeast)
#' out <- phangornMCMC(yeast)
#' plot(out)
#' getMCMCtrees()
#' }
#'
phangornMCMC <- function(x, ntrees = 3000, burnin = 1000, frequency = 1,
    tree0 = NULL, model = NULL, printevery = 100, bf=baseFreq(x), Q=rep(1, 6),
    optBf = FALSE, optQ = FALSE, rooted=TRUE)
{
    on.exit({
        pml.free()
        if (k < nOut) {
            if (!k) stop("burn-in period not yet finished")
            TREES <- TREES[seq_len(k)]
            LL <- LL[seq_len(i)]
            if(optBf) BF <- BF[seq_len(i),]
            if(optQ) QQ <- QQ[seq_len(i),]
            params <- params[seq_len(i), ]
            warning(paste("MCMC interrupted after", i, "generations"))
        }
        ## compress the list of trees:
        attr(TREES, "TipLabel") <- TREES[[1L]]$tip.label
        lapply(TREES, function(x) x$tip.label <- NULL)
        class(TREES) <- "multiPhylo"

        suffix <- 1
        list.trees <- .get.list.trees()
        if (l <- length(list.trees))
            suffix <- 1 + as.numeric(sub("TREES_", "", list.trees[l]))
        suffix <- sprintf("%03d", suffix)
        assign(paste("TREES", suffix, sep = "_"), TREES,
               envir = .phangornMCMCenv)

        i <- i - 1L


        if(optBf) assign("BF", envir = .phangornMCMCenv)

        MCMCstats <- get("MCMCstats", envir = .phangornMCMCenv)
        MCMCstats[[suffix]] <- c(k, burnin, frequency, i, j)
        assign("MCMCstats", MCMCstats, envir = .phangornMCMCenv)

        LL <- cbind(LL, params)
        colnames(LL) <- c("logLik", para.nms)
        LL <- mcmc(LL, start = 1, end = i)
        return(LL)
    })

    verbose <- as.logical(printevery)

    if (is.null(tree0)) {
        dm <- dist.ml(x)
        if(rooted)tree0 <- upgma(dm)
        else tree0 <- fastme.bal(dm, nni = TRUE, spr = FALSE, tbr = FALSE)
    }
    # ~ 2 * minimal length
    n <- nTip <- length(tree0$tip.label)
    lamda_edge <- fitch(tree0, x) / (sum(attr(x, "weight")) * Ntip(tree0))

    nstates <- length(bf)
    nodeMax <- 2L * nTip - 2L + rooted
    nOut <- ntrees
    nOut2 <- ntrees * frequency + burnin

    INV <- Matrix(lli(x, tree0), sparse = TRUE)
    ll.0 <- numeric(attr(x, "nr"))

    x <- subset(x, tree0$tip.label)
    ##
    eig <- edQt(Q, bf)
    pml.init(x)
    getlogLik <- function(phy, x) {
        pml.fit(phy, x, bf = bf, eig = eig, INV = INV, ll.0 = ll.0)
    }

    TREES <- vector("list", nOut)
    if(optBf){
        BF <- matrix(NA_real_, nOut, length(bf))
        BF[1,] <- bf
    }
    if(optQ){
        QQ <- matrix(NA_real_, nOut, length(Q))
        QQ[1,] <- Q
    }

    LL <- numeric(nOut2)
    TREES[[1L]] <- tree0
    lnL0 <- getlogLik(tree0, x)
    LL[1L] <- lnL0

    params <- matrix(0, nOut2, np)

    i <- 2L
    j <- 0L # number of accepted moves
    k <- 0L # number of sampled trees

    if (verbose) {
        cat("Running the Markov chain:\n")
        cat("  Number of trees to output:", ntrees, "\n")
        cat("  Burn-in period:", burnin, "\n")
        cat("  Sampling frequency:", frequency, "\n")
        cat("  Number of generations to run:", ntrees * frequency + burnin, "\n")
        cat("Generation    Nb of accepted moves\n")
    }

    bt0 <- branching.times(tree0)
    params[1L, ] <- para0 <- getparams(tree0, bt0)

    nodesToSample <- (n + 2):nodeMax

    while (k < nOut) {
        if (verbose) if (! i %% printevery)
            cat("\r  ", i, "                ", j, "           ")

        # tree rearrangements rooted trees

        ## select one internal node excluding the root:
        target <- sample(nodesToSample, 1L) # target node for rearrangement
        THETA <- f.theta(bt0[target - n], para0) # the value of THETA at this node

        tr.b <- NeighborhoodRearrangement(tree0, n, nodeMax, target, THETA, bt0)
        ## do TipInterchange() every 10 steps:
        ## tr.b <-
        ##     if (! i %% 10) TipInterchange(tree0, n)
        ##     else NeighborhoodRearrangement(tree0, n, nodeMax, target, THETA, bt0)

        if (!(i %% frequency) && i > burnin) {
            k <- k + 1L
            TREES[[k]] <- tr.b
        }
        lnL.b <- getlogLik(tr.b, x)
        LL[i] <- lnL.b
        ## calculate theta for the proposed tree:
        bt <- branching.times(tr.b)
        params[i, ] <- para <- getparams(tr.b, bt)
        i <- i + 1L
        ACCEPT <- if (is.na(lnL.b)) FALSE else {
            if (lnL.b >= lnL0) TRUE
            else rbinom(1, 1, exp(lnL.b - lnL0))
        }
        if (ACCEPT) {
            j <- j + 1L
            lnL0 <- lnL.b
            tree0 <- tr.b
            para0 <- para
            bt0 <- bt
        }

        if(optBf){
#            bftmp <- rdirichlet1(nstates)
#            lltmp <- pml.fit(phy, x, bf = bf, eig = eig, INV = INV, ll.0 = ll.0)

        }


    }
    if (verbose) cat("\nDone.\n")
}

.get.list.trees <- function()
    ls(envir = .phangornMCMCenv, pattern = "^TREES_")

getMCMCtrees <- function(chain = NULL)
{
    list.trees <- .get.list.trees()
    l <- length(list.trees)
    if (is.null(chain)) {
        if (!l) return(NULL)
        if (l == 1)
            return(get(list.trees, envir = .phangornMCMCenv))
        ## l > 1:
        cat("Several lists of MCMC trees are stored:\n\n")
        for (i in 1:l) cat(i, ":", list.trees[i], "\n")
        cat("\nReturn which number? ")
        chain <- as.numeric(readLines(n = 1))
    } else {
        if (!l) {
            warning("no list of MCMC trees stored")
            return(NULL)
        }
        if (l < chain) {
            warning("no enough lists of MCMC trees stored")
            return(NULL)
        }
    }
    get(paste("TREES", sprintf("%03d", chain), sep = "_"),
        envir = .phangornMCMCenv)
}

saveMCMCtrees <- function(destdir = ".", format = "RDS", ...)
{
    format <- match.arg(toupper(format), c("RDS", "NEWICK", "NEXUS"))
    switch(format, RDS = {
        FUN <- saveRDS
        suffix <- ".rds"
    }, NEWICK = {
        FUN <- write.tree
        suffix <- ".tre"
    }, NEXUS = {
        FUN <- write.nexus
        suffix <- ".nex"
    })
    list.trees <- .get.list.trees()
    l <- length(list.trees)
    if (!l) warning("no list of trees to save") else {
        for (i in 1:l) {
            f <- list.trees[i]
            outfile <- paste(destdir, "/", f, suffix, sep = "")
            FUN(get(f, envir = .phangornMCMCenv), outfile, ...)
        }
    }
}

cleanMCMCtrees <- function()
    rm(list = .get.list.trees(), envir = .phangornMCMCenv)

getLastTree <- function(X) X[[length(X)]]

getMCMCstats <- function()
{
    cat("MCMC chain summaries (chains as columns):\n\n")
    get("MCMCstats", envir = .phangornMCMCenv)
}
