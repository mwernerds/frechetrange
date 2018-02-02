### The Duetsch Vahrenhold Grid implementation on smaller San Francisco Data
library(frechet)

### Two helper functions to split the NaN-separated dataset into
### a List. Taken from library trajcomp (https://github.com/mwernerds/trajcomp)
getTrajectoryIDs <- function (Q)
{
    ret = matrix(nrow = nrow(Q));
    startPos = 1;
    k=1;
	for (i in 1:nrow(Q))
	{
	    if (is.nan(Q[i,1])){
	        ret[startPos:(i-1)] = k;
	        ret[i] = NaN;
	        k = k+1;
	        startPos = i+1;
		}
	}
	ret[startPos:nrow(Q)] = k;
    return (ret);
}

tsplit<-function(x) 
{
    indizes = getTrajectoryIDs(x);
    l = split(x,f=indizes);
    l$'NaN'= NULL;
    return(l);
}


### The Main Functionality
data(sanfrancisco)
datalist = tsplit(sanfrancisco);
queryIndex = 5;

### Create a handle for holding the data
ds = internal_createGrid();
### Feed the data into C++
tmp = lapply(datalist, function(x) internal_addTrajectoryToGrid(ds,as.matrix(x)));
### Now, build the index.
meshSize = 1.0;
internal_buildIndex(ds, meshSize);

### Now, query
result = internal_gridRangeQuery(ds,as.matrix(datalist[[queryIndex]]),.02);
print(sprintf("Range query found %d",length(result)));


### And visualize (query is red and thick, dataset is gray, found is black and small)
png("grid.png")
plot(sanfrancisco,t="l", col="gray");
lines(datalist[[queryIndex]],col="red", lwd=5)
lapply(result, function(x) lines(x,col="black"))
dev.off();
