require(mappoly)
setwd("~/repos/training/ACSPB/")
#### Reading data ####
dat <- read_geno_csv("family1_mappoly.csv", ploidy = 2)
plot(dat)
#### Filtering data ####
dat <- filter_missing(input.data = dat, type = "marker",
                      filter.thres = 0.05, inter = FALSE)
dat <- filter_missing(input.data = dat, type = "individual",
                      filter.thres = 0.05, inter = FALSE)
#dat <- filter_individuals(dat)
s.all <- filter_segregation(dat, 0.05/dat$n.mrk, inter = FALSE)
s.all <- make_seq_mappoly(s.all)
plot(s.all)
#### Two-point analysis ####
system.time(tpt <- est_pairwise_rf2(s.all, ncpus = 16))
###save(tpt, file = "all_tpt.rda", compress = "xz", compression_level = 9)
mat <- rf_list_to_matrix(tpt)

#### Grouping into chormosomes ####
grs <- group_mappoly(input.mat = mat, expected.groups = 12)
#### Building map for one chromosome ####

s1 <- make_seq_mappoly(grs, 1)
tpt1 <- est_pairwise_rf(s1, ncpus = 8)
m1 <- rf_list_to_matrix(tpt1)
o1 <- mds_mappoly(m1)
so1 <- make_seq_mappoly(o1)
map1 <- est_rf_hmm_sequential(input.seq = so1,
                              twopt = tpt1,
                              extend.tail = 20,
                              sub.map.size.diff.limit = 5)
map1.up <- est_full_hmm_with_global_error(input.map = map1,
                                          error = 0.05,
                                          tol = 10e-4,
                                          verbose = FALSE)

plot(map1.up, 
     mrk.names = T, 
     cex = .4, 
     P = "Parent 202", 
     Q = "Parent 4")


#### Ordering and gathering information to phase ####
LGS<-vector("list", 12)
for(i in 1:12){
  a <- make_seq_mappoly(grs, i)
  # SNP "PitaSNP040670" removed using diagnostic graphics
  tpt2 <- est_pairwise_rf(a, ncpus = 8)
  if(i == 4) a <- make_seq_mappoly(dat, a$seq.mrk.names[a$seq.mrk.names!="PitaSNP040670"])
  m <- rf_list_to_matrix(tpt2)
  o <- mds_mappoly(m)
  LGS[[i]] <- list(seq = make_seq_mappoly(o), tpt = tpt2)
  plot(m, ord =LGS[[i]]$seq, index = TRUE)
  #plot(m, ord =LGS[[i]]$seq$seq.mrk.names[310:330], index = TRUE)
}
mds.ord.whole.genome<-unlist(sapply(LGS, function(x) x$seq$seq.mrk.names))
plot(mat, ord = mds.ord.whole.genome, fact = 5)

#### Performing parallel phasing computation ####
phasing_and_hmm_rf <- function(X){
  fl <- paste0("output_map_ch_", X$ch, ".txt")
  sink(fl)
  map <- est_rf_hmm_sequential(input.seq = X$seq,
                               twopt = X$tpt,
                               extend.tail = 20,
                               sub.map.size.diff.limit = 5)
  sink()
  return(map)
}
dir.create("map_phasing_output/", showWarnings = FALSE)
setwd("map_phasing_output/")
cl <- parallel::makeCluster(12)
parallel::clusterEvalQ(cl, require(mappoly))
parallel::clusterExport(cl, "dat")
MAPs <- parallel::parLapply(cl, LGS, phasing_and_hmm_rf)
parallel::stopCluster(cl)

my.error.func<-function(X){
  x<-est_full_hmm_with_global_error(input.map = X,
                                    error = 0.05,
                                    tol = 10e-4,
                                    verbose = FALSE)
  return(x)
}
cl <- parallel::makeCluster(12)
parallel::clusterEvalQ(cl, require(mappoly))
parallel::clusterExport(cl, "dat")
MAP.err <- parallel::parLapply(cl,MAPs,my.error.func)
parallel::stopCluster(cl)

plot(lg1.err, mrk.names = T, cex = .4, P = "Parent 202", Q = "Parent 4")

Dat2 <- read.csv(file  = "Family2_Updated.txt", stringsAsFactors = FALSE, header = FALSE)
Dat2[1:10, 1:10]
