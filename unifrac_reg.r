file = "/Users/chris/projects/asmw5/results/1000-0/out/out.txt"
data = read.table(file, header=T)
attach(data)

metrics = c("cluster_path", "nj_path", "pcoa_path")

for (i in 1:length(metrics)) {
	metric = metrics[i]
	reg_u = lm(get(paste("u_", metric, sep=""))~get(paste("u_", metric, "_norm", sep="")))
	reg_w = lm(get(paste("w_", metric, sep=""))~get(paste("w_", metric, "_norm", sep="")))

	pdf(paste(metric, ".pdf", sep=""))
	par(mfrow=c(2, 1))

	plot(get(paste("u_", metric, sep="")), get(paste("u_", metric, "_norm", sep="")), 
	main=paste("R^2 = ", round(summary(reg_u)$adj.r.squared, 2), " (n = ", nrow(data), ")", sep=""), 
	xlab=paste("u_", metric, sep=""), ylab=paste("u_", metric, "_norm", sep=""))
	
	abline(reg_u)

	plot(get(paste("u_", metric, sep="")), get(paste("u_", metric, "_norm", sep="")), 
	main=paste("R^2 = ", round(summary(reg_w)$adj.r.squared, 2), " (n = ", nrow(data), ")", sep=""), 
	xlab=paste("w_", metric, sep=""), ylab=paste("w_", metric, "_norm", sep=""))
	abline(reg_w)
	dev.off()
}

