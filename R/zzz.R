.phangornMCMCenv <- new.env()

assign("MCMCstats",
       data.frame(row.names = c("Number of trees output",
                  "Burn-in period", "Sampling frequency",
                  "Number of generations", "Nb of accepted moves")),
       envir = .phangornMCMCenv)
