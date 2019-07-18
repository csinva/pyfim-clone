#!/usr/bin/Rscript

# simple example transactions
data <- list(c("a", "b", "c"),
             c("a", "d", "e"),
             c("b", "c", "d"),
             c("a", "b", "c", "d"),
             c("b", "c"),
             c("a", "b", "d"),
             c("d", "e"),
             c("a", "b", "c", "d"),
             c("c", "d", "e"),
             c("a", "b", "c"))

# write transaction data to a (temporary) file
cat(paste(sapply(data, function (d) paste(d, collapse="\t")),
          collapse="\n"), file="test.txt")

# set the mining parameters
zmin <- 2                       # minimum size = 2 items
zmax <- 4                       # maximum size = 4 items
                                # (default: no limit)
supp <- 10                      # minimum support = 10%
# supp <- -2                    # minimum support = 2 transactions
sep <- "\t"                     # item separator for output

# build the command to execute
cmd <- paste("fpgrowth",
             paste0("-s", supp),
             paste0("-m", zmin),
             paste0("-n", zmax),
             paste0("-k\"", sep, "\""),
             "-f\"\t\"", "-v\"\t%a\"",
             "test.txt", "test.out", sep=" ")

# execute the command
system(cmd)

# read the found frequent item sets (last element is support)
sets <- strsplit(scan("test.out", what=" ", sep="\n"),
                 split="\t", fixed=TRUE)
print(sets)
