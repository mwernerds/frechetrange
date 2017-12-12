### A very simple benchmark on some San Francisco trajectories
###
##
## Before using this script, install frechet, ggplot2 and microbenchmark

library("frechet");
library("ggplot2")
library("microbenchmark")


data(








tm <- microbenchmark(
    internal_frechet_decide_dv(M1,M2,2.99),
    internal_frechet_decide_bb(M1,M2,2.99),
    frechet.decide(M1,M2,2.99,"duetschvahrenhold"),
    frechet.decide(M1,M2,2.99,"bringmanbaldus"),
    times=10000L)

autoplot(tm)
