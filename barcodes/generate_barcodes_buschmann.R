# Input n and dist as args
args <- commandArgs(TRUE)
if (length(args) != 2) {
    stop('Usage: generate_barcodes.R <n> <dist>')
}
n <- as.integer(args[1])
dist <- as.integer(args[2])
if (n < 10) {
    heuristic <- 'ashlock'
} else {
    heuristic <- 'conway'
}

# Load DNABarcodes library as well as "stringr" for string searching
library("DNABarcodes")
library("stringr")
# Create a preliminary pool
pool <- create.pool(n=n)
# Filter out candidate barcodes that contain GGC
pool.new <- pool[!str_detect(pool, "GGC")]
# Create a barcode set that corrects one substitution (default values)
bc <- create.dnabarcodes(n=n,
                         dist=dist,
                         metric='seqlev',
                         pool=pool.new,
                         heuristic=heuristic,
                         cores=1)

# Output results
fname = sprintf('barcodes%d-%d.txt', n, dist)
f <- file(fname)
writeLines(bc, f)
close(f)
