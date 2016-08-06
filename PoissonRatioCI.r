#!/usr/bin/env Rscript
require(ggplot2)
require(rateratio.test)

#### FIGURE 4A, different read counts affect confidence interval ####

nvals = 500;  # number of different ratios to test
start_val = 3  # number of reads to start with, 1 gives really large number so starting with 3 for visibility
end_val = start_val + nvals - 1
seq_depth = 1000  # number of reads in library (doesn't actually matter as long as they are the same)

upperCI = data.frame(ratio=c(rep(0, nvals)), CIbound=c(rep("Upper", nvals)), reads=c(seq(start_val, nvals+start_val-1)))
lowerCI = data.frame(ratio=c(rep(0, nvals)), CIbound=c(rep("Lower", nvals)), reads=c(seq(start_val, nvals+start_val-1)))
for (i in start_val:end_val){
  cur_val = rateratio.test(c(i,i), c(seq_depth, seq_depth))$conf.int
  use_index = i-start_val+1
  lowerCI[use_index,"ratio"] = cur_val[1]
  upperCI[use_index,"ratio"] = cur_val[2]
}

CI = rbind(upperCI, lowerCI)
ggplot(data=CI, aes(x=reads, y=ratio, color=CIbound)) +
  geom_line()

#### FIGURE 4B, fewer read counts in RNaseR with different depths affect confidence interval ####


# assume 5-fold enrichment in RNaseR
# rate is 20 per 100 million reads in control library
# lambda = 20 for 100 million reads, 30 for 150 million reads, 40 for 200 million reads
# RNaseR library has 100 million reads and rate is 100 per 100 million reads
ratios1 = rpois(10000,100)/rpois(10000,20)
CI1 = rateratio.test(c(100,20), c(1000,1000))$conf.int

ratios2 = rpois(10000,100)/rpois(10000,30)
CI2 = rateratio.test(c(100,30), c(1000,1000))$conf.int

ratios3 = rpois(10000,100)/rpois(10000,40)
CI3 = rateratio.test(c(100,40), c(1000,1000))$conf.int

df1 = data.frame(ratio=ratios1, depth=c(rep("1x", 10000)))
df2 = data.frame(ratio=ratios2, depth=c(rep("1.5x", 10000)))
df3 = data.frame(ratio=ratios3, depth=c(rep("2x", 10000)))
df = rbind(df1, df2, df3)

ggplot(data=df, aes(x=ratio, color=depth)) + 
  geom_density() + xlab("Ratio") + ylab("Density") + scale_x_continuous(breaks=seq(1,20)) +
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=12,face="bold"), 
        legend.title=element_text(size=14,face="bold"), plot.title=element_text(size=18,face="bold")) 

# confidence intervals
print(CI1)
print(CI2)
print(CI3)
