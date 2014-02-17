get_valid_matrix = function(cols, numstates, rate) {
    found = 0
    while (found < cols) {
        root = sample(numstates)[1]
        temp = rTraitDisc(tree, model="ER", 
                          k=numstates, states=1:numstates, 
                          rate=rate, root.value=root,
                          ancestor=T)
        temp = as.matrix(temp)
        if (min(temp) < max(temp)) {
            if (found == 0) {
                m = temp
                roots = c(root)
            } else {
                m = cbind(m, temp)
                roots = c(roots, root)
            }
            found = found + 1
        }
    }
    return(list(m, roots))
}

get_continuous_matrix = function(cols) {
    r = rgamma(1, shape=1000, rate=1)
    m = replicate(cols, rTraitCont(tree, ancestor=T, sigma=r*0.5, root.value=r))
    m = round(as.matrix(m))
    return(list(m, r))
}

tree =  rtree(8)
data = get_valid_matrix(1, 8, 1)
data1 = get_continuous_matrix(1)
x = data[[1]]
x1 = data1[[1]]
x1 = ifelse(x1<0, 0, x1)
labels = x
plot(tree, show.tip.label=F)
Y <- labels[1:8]
A <- labels[-(1:8)]
nodelabels(A)
tiplabels(Y)