library(matrixStats)
library(xtable)

everything <- lapply(0:3, function (run_num) {

    nodes <- read.table(sprintf("nodes_positions_%d.txt", run_num))
    colnames(nodes) <- c("id", "location")
    nodes$pop <- cut(nodes$location, 25)
    if (FALSE) plot(nodes$location, col=nodes$pop) # sanity check

    # has one row per genome and one column per variant
    genotypes <- matrix(scan(sprintf("genotype_matrix_%d.txt", run_num)), nrow=nrow(nodes))

    # initializeMutationType("m1", 0.5, "g", 0.5, 1.0);  //Positive Direction, Dominant 
    # initializeMutationType("m2", 0.5, "g", -0.5, 1.0); //Negative Direction, Dominant
    # //initializeMutationType("m3", 0.5, "f", 0.0);       //Neither Direction, Neutral
    # initializeMutationType("m4", 0.5, "g", 0.5, 1.0);  //Positive Direction, Recessive
    # initializeMutationType("m5", 0.5, "g", -0.5, 1.0); //Negative Direction, Recessive
    # initializeMutationType("m6", 0.5, "g", 0.6, 1.0);  //Positive Direction, Additive
    # initializeMutationType("m7", 0.5, "g", -0.6, 1.0); //Negative Direction, Additive

    variants <- read.table(sprintf("variants_%d.txt", run_num))
    colnames(variants) <- c("id", "size", "type", "pos")
    stopifnot(all(variants$type != 3))
    stopifnot(ncol(genotypes) == nrow(variants))
    variants$freq <- colMeans(genotypes)
    variants$region <- cut(variants$pos, 10)
    variants$direction <- ifelse(variants$type %in% c(1, 4, 6), +1, -1)
    stopifnot(all(sign(variants$size) == variants$direction))
    variants$dominance <- c(1, 1, NA, 0, 0, 0.5, 0.5)[variants$type]

    # now we compute net effects by effect region and genome
    effects <- genotypes * variants$size[col(genotypes)]
    region_effects <- t(apply(effects, 1, tapply, variants$region, sum))
    nodes$phenotype <- rowSums(region_effects)

    # distribution of per-genome phenotypes
    if (FALSE) hist(nodes$phenotype)

    pos_effects <- t(apply(ifelse(effects > 0, effects, 0), 1, tapply, variants$region, sum))
    neg_effects <- t(apply(ifelse(effects < 0, effects, 0), 1, tapply, variants$region, sum))

    return(list(genotypes=genotypes,
                variants=variants,
                nodes=nodes, 
                effects=effects,
                region_effects=region_effects,
                pos_effects=pos_effects,
                neg_effects=neg_effects))
  } )

pdf(file='frequencies.pdf', width=8, height=5, pointsize=10)
lapply(everything, function (x) with(list2env(x), {
    # compute per-lake frequencies: will have one row per population
    #   and one column per lake
    lake_freqs <- apply(genotypes, 2, tapply, nodes$pop, mean)

    layout(t(1:2))

    # this shows how frequences alleles vary across lakes (one line per allele)
    matplot(lake_freqs, type='l',
            xlab='lake position',
            ylab='allele frequency')

    signed_freq <- variants$freq * variants$direction
    # check the variants fall in the effect regions
    plot(jitter(variants$pos, factor=100), signed_freq, col=variants$region,
         ylim=c(-1,1), 
         pch=20, xlab='position along the genome', ylab='allele frequency')
  } ))
dev.off()


pdf(file='effects.pdf', width=12, height=8, pointsize=10)
lapply(everything, function (x) with(list2env(x), {
    layout(matrix(1:12, nrow=3))
    par(mar=c(2, 2, 0, 0)+.1)
    for (k in 1:ncol(pos_effects)) {
        plot(pos_effects[,k], neg_effects[,k], pch=20, cex=2,
             xlim=c(0, max(pos_effects)),
             ylim=c(min(neg_effects), 0),
             asp=1,
             xlab='', ylab='')
        abline(0, -1)
    }
  } ))
dev.off()

pdf(file='effect_hists.pdf', width=12, height=8, pointsize=10)
layout(t(1:2))
lapply(everything, function (x) with(list2env(x), {
    hist(nodes$phenotype, xlab='genome phenotype', main='phenotype by genome', xlim=c(-10,10))
    hist(pos_effects + neg_effects, main='net effect by region', xlab='net effect', xlim=c(-6,6))
  } ))
dev.off()

# Find the total frequency, in each region, of any haplotypes with total effect below thresh
thresh <- -1.0
sapply(everything, function (x) with(list2env(x), {
    colMeans(region_effects < thresh)
  }))

# sum of the most negative haplotypes in each region, by population
best_haps <- sapply(everything, function (x) with(list2env(x), {
    tapply(1:nrow(region_effects), nodes$pop, function (z) {
        sum(colMins(region_effects[z,])) })
  }))

# expected total sum of fixed haplotypes, per lake,
# with prob of fixation 2 * beta * z * u, where u is the effect size; z is the deviation from optimum,
# and beta is 1/450
expected_fix <- sapply(everything, function (x) with(list2env(x), {
    pfix <- 2 * 10 * abs(neg_effects) / 450
    sapply( tapply(1:nrow(pfix), nodes$pop, function (k) {
                          pp <- pfix[k,]
                          xx <- abs(neg_effects)[k,]
                          sapply(1:ncol(pp), function (j) {
                                ord <- order(pp[,j], decreasing=TRUE)
                                p <- pp[,j][ord]
                                x <- xx[,j][ord]
                                cfix <- p * cumprod(c(1, 1-p[-length(p)]))
                                return(sum(x*cfix))
              } )
          } ), sum)
     }))
colnames(expected_fix) <- paste("m=", c(0.01, 0.1, 1, 10))

## This is the table of the text
xtable(data.frame(t(colMeans(expected_fix)), check.names=FALSE))

# \begin{tabular}{rrrr}
#   \hline
#   m= 0.01 & m= 0.1 & m= 1 & m= 10 \\ 
#   \hline
#   0.12 & 4.70 & 12.93 & 12.27 \\ 
#   \hline
# \end{tabular}

# sum of haplotype effect size, weighted by fitness
# with fitness beta * z * u, where u is the effect size; z is the deviation from optimum,
# and beta is 1/450
haplotype_fitness <- sapply(everything, function (x) with(list2env(x), {
    fitnesses <- 10 * abs(neg_effects) / 450
    # prob of fixation (to account for small s); eq'n 1.64 in Ewens
    fitnesses <- ifelse(fitnesses == 0, 1/400, fitnesses / (1 - exp(- 2 * 400 * fitnesses)))
    sapply( tapply(1:nrow(fitnesses), nodes$pop, function (k) {
                       fit <- fitnesses[k,]
                       fit <- sweep(fit, 2, colSums(fit), "/")
                       colSums(fit * abs(neg_effects[k,])) } ), sum)
  } ))


