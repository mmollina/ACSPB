## Loading MAPpoly
require(mappoly)
require(tidyverse)

## Gtahering file name
file.name <- system.file("extdata/potato_example.csv", package = "mappoly")

## Loading data
dat <- read_geno_csv(file.in = file.name, ploidy = 4)
print(dat, detailed = T)
plot(dat)

## Filter individuals based on their genomic relationship
dat <- filter_individuals(dat)

## Filter based on missing data
dat <- filter_missing(dat, type = "marker", filter.thres = .05)
dat <- filter_missing(dat, type = "individual", filter.thres = .05)

## Filter based on segregation
seq.filt <- filter_segregation(dat, chisq.pval.thres = 0.05/dat$n.mrk)
seq.filt <- make_seq_mappoly(seq.filt)

## Remove redundant markers (will be included in the final map)
seq.red  <- elim_redundant(seq.filt)
plot(seq.red)

## Initial sequence with markers filtered
seq.init <- make_seq_mappoly(seq.red)
plot(seq.init)

## Genome order
go <- get_genomic_order(input.seq = seq.init) ## get genomic order of the sequence
go
plot(go)

## Pairrwise recombination fraction
ncores <- parallel::detectCores() - 1
tpt <- est_pairwise_rf2(seq.init, ncpus = ncores)
m <- rf_list_to_matrix(tpt) ## converts rec. frac. list into a matrix
sgo <- make_seq_mappoly(go) ## creates a sequence of markers in the genome order
plot(m, ord = sgo, fact = 5) ## plots a rec. frac. matrix using the genome order, averaging neighbor cells in a 5 x 5 grid
#plot(m, ord = names(which(sgo$chrom == 7)))

## Building linkage groups
g <- group_mappoly(m, expected.groups = 12, comp.mat = TRUE)
plot(g)
g
  
## Selecting Chromosome 1 
ch1 <- make_seq_mappoly(g, 9, ## Select LG1 
                       genomic.info = 1) 
m1 <- make_mat_mappoly(m, ch1)


## Ordering using MDS algorithm
mds.o1 <- mds_mappoly(input.mat = m1)
s1.mds <- make_seq_mappoly(mds.o1)
plot(m1, ord = s1.mds)

## Genome order for chromosome 1
gen.o1 <- get_genomic_order(ch1)
s1.gen <- make_seq_mappoly(gen.o1)
plot(m1, ord = s1.gen)

## Comparing orders
x<-match(s1.gen$seq.mrk.names, s1.mds$seq.mrk.names)
cor(x, seq_along(x))
plot(x)

## Phasing
tpt1 <- est_pairwise_rf(s1.mds, ncpus = ncores)
lg1.map <- est_rf_hmm_sequential(input.seq = s1.mds,
                                 start.set = 3,
                                 thres.twopt = 10,
                                 thres.hmm = 20,
                                 extend.tail = 20,
                                 info.tail = TRUE,
                                 twopt = tpt1,
                                 sub.map.size.diff.limit = 10,
                                 phase.number.limit = 20,
                                 reestimate.single.ph.configuration = TRUE,
                                 tol = 10e-3,
                                 tol.final = 10e-4)
print(lg1.map)
plot(lg1.map, mrk.names = TRUE, cex = 0.7)

## Re-estimating map 
lg1.map.up <- est_full_hmm_with_global_error(input.map = lg1.map, 
                                             error = 0.05, 
                                             verbose = TRUE)
plot(lg1.map.up, mrk.names = TRUE, cex = 0.7)

## Load full map
in.file <- "https://github.com/mmollina/SCRI/raw/main/docs/tetra/maps_updated.rda"
map_file <- tempfile()
download.file(in.file, map_file)
load(map_file)
plot_map_list(MAPs.up, col = mp_pallet3(12), title = "Linkage groups")

## Compute conditional probabilities of the haplotypes
g <- lapply(MAPs.up, calc_genoprob_error, step = 1, error = 0.05)
h <- calc_homologprob(g)
plot(h, lg = "all", ind = 2)

## Preferential pairing
p = calc_prefpair_profiles(g)
plot(p, min.y.prof = 0.1, max.y.prof = 0.7, P1 = "Atlantic", P2 = "B1829.5")

## Summary
summary_maps(MAPs.up)

## MAP vs genome
plot_genome_vs_map(MAPs.up)

## Export to QTLpoly
to.qtlpoly <- export_qtlpoly(g)

## Export map (publication)
export_map_list(lg1.map.up, file = "output_file.csv")

  
