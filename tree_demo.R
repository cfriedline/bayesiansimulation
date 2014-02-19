get_matrix = function(cols, numstates, rate, model) {
    found = 0
    while (found < cols) {
        root = sample(numstates, size=1)
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
    }
    return(list(m, roots))
}

get_state_model = function(num_states) {
    mat=matrix(seq(1:num_states**2),num_states)
    for (i in 1:num_states) {
        for (j in 1:num_states) {
            mat[i,j] = abs(i-j)
        }
    }
    #print(mat)
    return(mat)
}

step_func = function(x, l) {
    probs = rep(1, num_states)
    m = matrix(1:num_states,1)
    for (i in 1:ncol(m)) {
        diff=abs(x-i)
        if (diff==0) {
            m[1,i] = 0
        } else {
            m[1,i] = (1/num_states)/diff
        }        
    }
    m[1,x] = 1-sum(m)
    #print(paste("ancestor =",x))
    #print(m)
    return(sample(num_states, size=1, prob=probs))
}
    

get_continuous_matrix = function(cols) {
    r = round(rgamma(1, shape=1000, rate=1))
    r = sample(100, size=1)
    m = replicate(cols, rTraitCont(tree, ancestor=T, sigma=r*100, root.value=r))
    m = round(as.matrix(m))
    return(list(m, r))
}

num_states = 8
cols = 1
rate = 1
tree =  rtree(num_states)
tree$edge.length = rep(0.5, length(tree$edge.length))
data_er = get_matrix(cols, num_states, rate, "ER")
data_cont = get_continuous_matrix(cols)
data_state = get_matrix(cols, num_states, rate, get_state_model(num_states))
data_step = get_matrix(cols, num_states, rate, step_func)
x_er = data_er[[1]]
x_cont = data_cont[[1]]
x_cont = ifelse(x_cont<0, 0, x_cont)
x_state = data_state[[1]]
x_step = data_step[[1]]

#pdf("./tree_demo.pdf", height=8.5, width=20)
par(mfrow=c(2,2))
data = list(x_er, x_cont, x_state, x_step)
titles = c('ER', "Continuous", "State", "Step")
for (i in 1:length(data)) {
    plot(tree, show.tip.label=F)
    labels = data[i][[1]]
    title(titles[i])
    Y <- labels[1:num_states]
    A <- labels[-(1:num_states)]
    nodelabels(A)
    tiplabels(Y)        
}

#dev.off()
