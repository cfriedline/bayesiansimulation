library(phangorn)

generate_first_triplet = function(bits) {
        triplet = replicate(bits, rTraitDisc(tree, model="ER", k=2,states=0:1))
        triplet = t(apply(triplet, 1, as.numeric))
        sums = rowSums(triplet)
        if (length(which(sums==0)) > 0 && length(which(sums==3)) == 1) {
            return(triplet)
        }
        return(generate_triplet(bits))
        }

get_valid_triplets = function(numsamples, needed, bits) {
            tryCatch({
                m = generate_first_triplet(bits)
                found = 0
                while (found < needed) {
                    triplet = replicate(bits, rTraitDisc(tree, model="ER", k=2,states=0:1))
                    triplet = t(apply(triplet, 1, as.numeric))
                    sums = rowSums(triplet)
                    if (length(which(sums==0)) > 0 && length(which(sums==3)) == 1) {
                        m = cbind(m, triplet)
                        found = found + 1
                    }

                }
            return(m)
            }, error = function(e){print(message(e))}, warning = function(e){print(message(e))})
        }

tree = rtree(10)
test = get_valid_triplets(10, 100, 3)