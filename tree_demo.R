library(ape)
library(plotrix)
rTraitDisc2 <-
    function(phy, model = "ER", k = if (is.matrix(model)) ncol(model) else 2,
             rate = 0.1, states = LETTERS[1:k], freq = rep(1/k, k),
             ancestor = FALSE, root.value = 1, ...)
    {
        if (is.null(phy$edge.length))
            stop("tree has no branch length")
        if (any(phy$edge.length < 0))
            stop("at least one branch length negative")
        
        if (is.character(model)) {
            switch(toupper(model), "ER" = {
                if (length(rate) != 1)
                    stop("`rate' must have one element")
                Q <- matrix(rate, k, k)
            }, "ARD" = {
                if (length(rate) != k*(k - 1))
                    stop("`rate' must have k(k - 1) elements")
                Q <- matrix(0, k, k)
                Q[col(Q) != row(Q)] <- rate
            }, "SYM" = {
                if (length(rate) != k*(k - 1)/2)
                    stop("`rate' must have k(k - 1)/2 elements")
                Q <- matrix(0, k, k)
                sel <- col(Q) < row(Q)
                Q[sel] <- rate
                Q <- t(Q)
                Q[sel] <- rate
            })
        }
        if (is.matrix(model)) {
            Q <- model
            if (ncol(Q) != nrow(Q))
                stop("the matrix given as `model' must be square")
        }
        
        phy <- reorder(phy, "pruningwise")
        n <- length(phy$tip.label)
        N <- dim(phy$edge)[1]
        ROOT <- n + 1L
        x <- integer(n + phy$Nnode)
        x[ROOT] <- as.integer(root.value)
        
        anc <- phy$edge[, 1]
        des <- phy$edge[, 2]
        el <- phy$edge.length
        
        if (is.function(model)) {
            environment(model) <- environment() # to find 'k'
            for (i in N:1) {
                x[des[i]] <- model(x[anc[i]], el[i], ...)
            }
        } else {
            #print(freq)
            freq <- rep(freq, each = k)
            #print(freq)
            Q <- Q * freq
            diag(Q) <- 0
            diag(Q) <- -rowSums(Q)
            for (i in N:1) {
                p <- matexpo(Q * el[i])[x[anc[i]], ]
                x[des[i]] <- sample.int(k, size = 1, FALSE, prob = p)
            }
        }
        
        if (ancestor) {
            if (is.null(phy$node.label)) phy <- makeNodeLabel(phy)
            names(x) <- c(phy$tip.label, phy$node.label)
        } else {
            x <- x[1:n]
            names(x) <- phy$tip.label
        }
        class(x) <- "factor"
        levels(x) <- states
        x
    }

get_matrix = function(cols, model) {
    found = 0
    while (found < cols) {
        root = sample(nrow(model), size=1)
        temp = rTraitDisc2(tree, model=model, 
                          states=1:nrow(model), 
                          root.value=root,
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

get_state_model = function(num_states,rate) {
    mat=matrix(seq(1:num_states**2),num_states)
    for (i in 1:num_states) {
        for (j in 1:num_states) {
            if (i != j) {
                mat[i,j] = (rate/num_states)*(num_states-(abs(i-j)))
            } else {
                mat[i,j] = 0
            }
        }
    }
    print(mat)
    return(mat)
}

get_er_model = function(num_states, rate) {
    mat = matrix(rep(rate, num_states**2), num_states)
    print(mat)
    mat
}

get_restricted_er_model = function(num_states, rate) {
    mat = get_er_model(num_states, rate)
    for (i in 1:nrow(mat)) {
        for (j in 1:ncol(mat)) {
            if (abs(i-j) == 1) {
                mat[i,j] = rate    
            } else if (i != j) {
                mat[i,j] = rate*0.00001
            } else {
                mat[i,j] = 0
            }
        }
    }
    print(mat)
    mat
}

get_continuous_matrix = function(cols) {
    r = round(rgamma(1, shape=1000, rate=1))
    r = sample(100, size=1)
    m = replicate(cols, rTraitCont(tree, ancestor=T, sigma=r*100, root.value=r))
    m = round(as.matrix(m))
    m = ifelse(m<0, 0, m)
    return(list(m, r))
}

num_states = 8
cols = 1
rate = 1
tree =  rtree(num_states)
tree$edge.length = rep(0.5, length(tree$edge.length))
er_model = get_er_model(num_states, rate)
data_er = get_matrix(cols, er_model)
res_er_model = get_restricted_er_model(num_states, rate)
data_res_er = get_matrix(cols, res_er_model)
data_cont = get_continuous_matrix(cols)
state_model = get_state_model(num_states, rate)
data_state = get_matrix(cols, state_model)

#pdf(paste("tree_demo_", rate, ".pdf", sep=""), height=8.5, width=20)
l = matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = TRUE)
#layout(l)
data = list(data_er, data_res_er, data_cont, data_state)
titles = c('ER', "Restricted ER", "Continuous", "State")
tables = list(er_model, res_er_model, matrix(), state_model)
for (i in 1:length(data)) {
    par(mfrow=c(1,2))
    plot(tree, show.tip.label=F)
    labels = data[i][[1]][[1]]
    root = data[i][[1]][[2]]
    title(paste(titles[i], " (root=", root, ")", sep=""))
    add.scale.bar()
    Y <- labels[1:num_states]
    A <- labels[-(1:num_states)]
    nodelabels(A)
    tiplabels(Y)
    plot.new()
    if (i!=3) {
        addtable2plot(0,0.4,table=tables[i][[1]], xpad=1, title=paste("Q matrix, rate=", rate, sep=""), cex=0.6, hlines=T, vlines=T)
    }
}
#dev.off()
