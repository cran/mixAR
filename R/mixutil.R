## Do not edit this file manually.
## It has been automatically generated from *.org sources.

### return the last n elements of a vector.
lastn <- function(x, n){
    m <- length(x) - n               # error if m < 0
    if(m == 0) 
        x
    else 
        x[-seq_len(m)] # 2020-03-27 was: x[-(1:m)]
}
