options(scipen = 6) # bias against scientific notation

file = "/Users/chris/projects/asmw5/results/1000-0/out/out.txt"
data = read.table(file, header=T)
attach(data)

methods = c("cluster", "nj", "pcoa")
dists = c("symm", "topo", "path")

for (i in 1:length(methods)) {
	for (j in 1:length(dists)) {
		metric = paste(methods[i], "_", dists[j], sep="")
		print(metric)
		reg_u = lm(get(paste("u_", metric, sep=""))~get(paste("u_", metric, "_norm", sep="")))
		reg_w = lm(get(paste("w_", metric, sep=""))~get(paste("w_", metric, "_norm", sep="")))
		summ_reg_u = summary(reg_u)
		summ_reg_w = summary(reg_w)

		pdf(paste(metric, ".pdf", sep=""))
		par(mfrow=c(2, 1))

		plot(get(paste("u_", metric, sep="")), get(paste("u_", metric, "_norm", sep="")), 
		main=paste("R^2 = ", round(summ_reg_u$adj.r.squared, 2), " (n = ", nrow(data), ", p = ", round(summ_reg_u$coefficients[2,4], 4), ")", sep=""), 
		xlab=paste("u_", metric, sep=""), ylab=paste("u_", metric, "_norm", sep=""))
	
		abline(reg_u)

		plot(get(paste("u_", metric, sep="")), get(paste("u_", metric, "_norm", sep="")), 
		main=paste("R^2 = ", round(summ_reg_w$adj.r.squared, 2), " (n = ", nrow(data), ", p = ", round(summ_reg_w$coefficients[2,4], 4), ")", sep=""), 
		xlab=paste("w_", metric, sep=""), ylab=paste("w_", metric, "_norm", sep=""))
		abline(reg_w)
		dev.off()
	}
	
}

