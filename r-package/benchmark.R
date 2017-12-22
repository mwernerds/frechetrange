### A very simple benchmark on some San Francisco trajectories
###
##
## Before using this script, install frechet, ggplot2 and microbenchmark

library("frechet");
library("ggplot2")
library("microbenchmark")

library(sp);
### Some basic test data

input = list(t1 = matrix(c(1,0,2,0,4,0,5,0), ncol=2),
             t2 = matrix(c(1,0,3,3,5,0), ncol=2))

formats = list (
    t.matrix = input,
    t.data.frame = list(t1 = data.frame(x = input[["t1"]],y = input[["t1"]]),
                       t2 = data.frame(x = input[["t2"]],y = input[["t2"]]))
    )

formats[["spdf"]]  = lapply(formats[["t.data.frame"]], function(x) SpatialPointsDataFrame(x, data.frame(ID=1:nrow(x))));

M1 = formats[["t.matrix"]][[1]];
M2 = formats[["t.matrix"]][[2]];





tm <- microbenchmark(
    internal_frechet_decide_dv(M1,M2,2.99),
    internal_frechet_decide_bb(M1,M2,2.99),
    frechet.decide(M1,M2,2.99,"duetschvahrenhold"),
    frechet.decide(M1,M2,2.99,"bringmanbaldus"),
    times=10000L)

autoplot(tm)
