### A very basic test suite for the trajcomp package
library(sp);
library(frechet)

### Some basic test data

input = list(t1 = matrix(c(1,0,2,0,4,0,5,0), ncol=2),
             t2 = matrix(c(1,0,3,3,5,0), ncol=2))

formats = list (
    t.matrix = input,
    t.data.frame = list(t1 = data.frame(x = input[["t1"]],y = input[["t1"]]),
                       t2 = data.frame(x = input[["t2"]],y = input[["t2"]]))
    )

formats[["spdf"]]  = lapply(formats[["t.data.frame"]], function(x) SpatialPointsDataFrame(x, data.frame(ID=1:nrow(x))));



### Test frechet.decide
quiet=lapply(1:length(formats), function(i){
    print(sprintf("Testing Format:    %s", names(formats)[i]))
    x = formats[[i]];
    print(sprintf("True class is:     %s",class(x[["t1"]])))
    frechet.decide(x[["t1"]],x[["t2"]],1.0)
});

### A very basic test suite for the trajcomp package
library(sp);
library(frechet)

### Some basic test data

input = list(t1 = matrix(c(1,0,2,0,4,0,5,0), ncol=2),
             t2 = matrix(c(1,0,3,3,5,0), ncol=2))

formats = list (
    t.matrix = input,
    t.data.frame = list(t1 = data.frame(x = input[["t1"]],y = input[["t1"]]),
                       t2 = data.frame(x = input[["t2"]],y = input[["t2"]]))
    )

formats[["spdf"]]  = lapply(formats[["t.data.frame"]], function(x) SpatialPointsDataFrame(x, data.frame(ID=1:nrow(x))));



### Test frechet.decide
quiet=lapply(1:length(formats), function(i){
    print(sprintf("Testing Format:    %s", names(formats)[i]))
    x = formats[[i]];
    print(sprintf("True class is:     %s",class(x[["t1"]])))
    frechet.decide(x[["t1"]],x[["t2"]],1.0)
});

### A very basic test suite for the trajcomp package
library(sp);
library(frechet)

### Some basic test data

input = list(t1 = matrix(c(1,0,2,0,4,0,5,0), ncol=2),
             t2 = matrix(c(1,0,3,3,5,0), ncol=2))

formats = list (
    t.matrix = input,
    t.data.frame = list(t1 = data.frame(x = input[["t1"]],y = input[["t1"]]),
                       t2 = data.frame(x = input[["t2"]],y = input[["t2"]]))
    )

formats[["spdf"]]  = lapply(formats[["t.data.frame"]], function(x) SpatialPointsDataFrame(x, data.frame(ID=1:nrow(x))));



### Test frechet.decide
quiet=lapply(1:length(formats), function(i){
    print(sprintf("Testing Format:    %s", names(formats)[i]))
    x = formats[[i]];
    print(sprintf("True class is:     %s",class(x[["t1"]])))
    frechet.decide(x[["t1"]],x[["t2"]],1.0)
});

### A very basic test suite for the trajcomp package
library(sp);
library(frechet)

### Some basic test data

input = list(t1 = matrix(c(1,0,2,0,4,0,5,0), ncol=2),
             t2 = matrix(c(1,0,3,3,5,0), ncol=2))

formats = list (
    t.matrix = input,
    t.data.frame = list(t1 = data.frame(x = input[["t1"]],y = input[["t1"]]),
                       t2 = data.frame(x = input[["t2"]],y = input[["t2"]]))
    )

formats[["spdf"]]  = lapply(formats[["t.data.frame"]], function(x) SpatialPointsDataFrame(x, data.frame(ID=1:nrow(x))));



### Test frechet.decide
quiet=lapply(1:length(formats), function(i){
    print(sprintf("Testing Format:    %s", names(formats)[i]))
    x = formats[[i]];
    print(sprintf("True class is:     %s",class(x[["t1"]])))
    frechet.decide(x[["t1"]],x[["t2"]],1.0)
});

### A very basic test suite for the trajcomp package
library(sp);
library(frechet)

### Some basic test data

input = list(t1 = matrix(c(1,0,2,0,4,0,5,0), ncol=2),
             t2 = matrix(c(1,0,3,3,5,0), ncol=2))

formats = list (
    t.matrix = input,
    t.data.frame = list(t1 = data.frame(x = input[["t1"]],y = input[["t1"]]),
                       t2 = data.frame(x = input[["t2"]],y = input[["t2"]]))
    )

formats[["spdf"]]  = lapply(formats[["t.data.frame"]], function(x) SpatialPointsDataFrame(x, data.frame(ID=1:nrow(x))));



### Test frechet.decide
quiet=lapply(1:length(formats), function(i){
    print(sprintf("Testing Format:    %s", names(formats)[i]))
    x = formats[[i]];
    print(sprintf("True class is:     %s",class(x[["t1"]])))
    frechet.decide(x[["t1"]],x[["t2"]],1.0)
});


print(sprintf("DV Decider: %d",internal_frechet_decide_dv(formats[["t.matrix"]][[1]]+3,formats[["t.matrix"]][[2]]+5,1)))

### and calculate it (distance is 2)
M1 = formats[["t.matrix"]][[1]];
M2 = formats[["t.matrix"]][[2]];

for (eps in seq(2.5,3.5,0.1))
{
    print(sprintf("Decide %.2f=%d",eps,frechet.decide(M1,M2,eps)))
    print(sprintf("Decide %.2f=%d",eps,frechet.decide(M1,M2,eps,"bringmanbaldus")))
}

### and a benchmark
library("ggplot2")
library("microbenchmark")

tm <- microbenchmark(
    internal_frechet_decide_dv(M1,M2,2.99),
    internal_frechet_decide_bb(M1,M2,2.99),
    frechet.decide(M1,M2,2.99,"duetschvahrenhold"),
    frechet.decide(M1,M2,2.99,"bringmanbaldus"),
    times=10000L)

autoplot(tm)



