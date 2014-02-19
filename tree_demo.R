get_matrix = function(cols, numstates, rate, model) {
    found = 0
    while (found < cols) {
        root = sample(numstates)[1]
        temp = rTraitDisc(tree, model=model, 
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
    }?
    return(list(m, roots))
}

get_step_model = function(num_states) {
    mat=matrix(seq(1:num_states**2),num_states)
    for (i in 1:num_states) {
        for (j in 1:num_states) {
            mat[i,j] = abs(i-j)
        }
    }
    return(mat)
}


get_continuous_matrix = function(cols) {
    r = round(rgamma(1, shape=1000, rate=1))
    m = replicate(cols, rTraitCont(tree, ancestor=T, sigma=r*0.5, root.value=r))
    m = round(as.matrix(m))
    return(list(m, r))
}

num_states = 8
cols = 1
rate = 1
tree =  rtree(num_states)
data_er = get_matrix(cols, num_states, rate, "ER")
data_cont = get_continuous_matrix(cols)
data_step = get_matrix(cols, num_states, rate, get_step_model(num_states))
x_er = data_er[[1]]
x_cont = data_cont[[1]]
x_cont = ifelse(x_cont<0, 0, x_cont)
x_step = data_step[[1]]

pdf("./tree_demo.pdf", height=8.5, width=20)

par(mfrow=c(1,3))
plot(tree, show.tip.label=F)
labels = x_er
title(paste("ER, root=", data_er[[2]], sep=""))
Y <- labels[1:num_states]
A <- labels[-(1:num_states)]
nodelabels(A)
tiplabels(Y)    

plot(tree, show.tip.label=F)
labels = x_cont
title(paste("Continuous, root=", data_cont[[2]], sep=""))
Y <- labels[1:num_states]
A <- labels[-(1:num_states)]
nodelabels(A)
tiplabels(Y)    

plot(tree, show.tip.label=F)
labels = x_step
title(paste("Step, root=", data_step[[2]], sep=""))
Y <- labels[1:num_states]
A <- labels[-(1:num_states)]
nodelabels(A)
tiplabels(Y)    
dev.off()
