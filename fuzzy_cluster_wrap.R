## Fuzzy cluster wrapper
#
#  Created by Xi Jiang, Jul.31, 2014

fuzzy_cluster_wrap <- function(data_dir,n_cluster,steps=100,fuzziness=2) { 
 # data_dir: string of file location, and the input file is an ascii-txt containing a matrix of time bins by features
 # n_cluster: vector of cluster center numbers
 # steps: number of clustering steps, default at 100
 # fuzziness: extent to which cluster memberships can deviate, default at 2

	library("e1071")	# requires this library for the fuzzy-cluster "cmeans" function

	kdata <- read.table(data_dir)
	c <- data.matrix(kdata)
	
#	center_list <- NULL
#	size_list <- NULL
#	cluster_list <- NULL
#	membership_list <- NULL
	for (i in 1:length(n_cluster)) {
		cl <- cmeans(c,n_cluster[i],steps,method="cmeans",m=fuzziness)
#		center_list[i] = cl$center
#		size_list[i] = cl$size
#		cluster_list[i] = cl$cluster
#		membership_list[i] = cl$membership
		write.table(cl$center,file="fuzz_center",append=T)
		write.table(cl$size,file="fuzz_size",append=T)
		write.table(cl$cluster,file="fuzz_cluster",append=T)
		write.table(cl$membership,file="fuzz_membership",append=T)
		
	}
#	write.table(center_list,file="fuzz_center")
#	write.table(size_list,file="fuzz_size")
#	write.table(cluster_list,file="fuzz_cluster")
#	write.table(membership_list,file="fuzz_membership")
}
