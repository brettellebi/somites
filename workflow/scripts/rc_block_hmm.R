# get a random dpAD file
get_random_dpAB <- function() {
	dp_file_path = "/Users/tomas/projects/2016/medaka/F2_crosses/recombination_blocks/dp_AB/170925/"
	dp_files = dir(dp_file_path)[grep(".txt", dir(dp_file_path))]
	dp_file = paste(dp_file_path, "/", dp_files[round(runif(1, 1, length(dp_files)))], sep="")
return(read.table(dp_file))
}

# example function of binnning chromosome producnig a plot
bin_and_plot <- function(d, chr=5, max_count=15, bin_len = 10000) {
	cleand = d[d$V4 < max_count & d$V6 < max_count,]
	cleand$bin = floor(cleand$V2 / bin_len)
	cleand$ratio = cleand$V4 / (cleand$V4 + cleand$V6)

	binned = as.data.frame(tapply(cleand[cleand$V1 == chr,]$V4,cleand[cleand$V1 == chr,]$bin,sum))
	binned$v6sum = tapply(cleand[cleand$V1 == chr,]$V6,cleand[cleand$V1 == chr,]$bin,sum)
	binned$bin = tapply(cleand[cleand$V1 == chr,]$bin,cleand[cleand$V1 == chr,]$bin,unique)
	colnames(binned) = c("mat", "pat", "bin")
	binned$ratio = binned$mat / (binned$mat + binned$pat)
	
	plot(binned$bin, binned$ratio)
invisible(binned)
}

# clean and bin all chrosmomes from dpAB file
bin_all_chrs <- function(d, max_count=15, bin_len = 10000) {
	d_list = split( d , f = d$V1 )
	binned_chrs = lapply(d_list, bin_single, max_count=max_count,  bin_len=bin_len)
invisible(binned_chrs)
}

# bin and clean a single chromsome
bin_single <- function(d, max_count=15, bin_len = 10000) {
	cleand = d[d$V4 < max_count & d$V6 < max_count,]
	cleand$bin = floor(cleand$V2 / bin_len)
	cleand$ratio = cleand$V4 / (cleand$V4 + cleand$V6)
	binned = as.data.frame(tapply(cleand$V4,cleand$bin,sum))
	binned$v6sum = tapply(cleand$V6,cleand$bin,sum)
	binned$bin = tapply(cleand$bin,cleand$bin,unique)
	colnames(binned) = c("mat", "pat", "bin")
	binned$ratio = binned$mat / (binned$mat + binned$pat)
invisible(binned)
}

# bin all chromosome and plot
bin_all_chrs_and_plot <- function(d, max_count=15, bin_len = 10000) {
	binned = bin_all_chrs(d, max_count, bin_len)
	#par(mfrow=c(4,6), mar=c(2,2,2,2))
	#lapply(binned, function(x) { plot(x$bin, x$ratio, main="", xlab="bin", ylab="ratio", cex=0.5) })
	bdf = do.call(rbind, binned)
	bdf$chr = unlist(strsplit(rownames(bdf), "\\."))[seq(1, nrow(bdf)*2, by=2)]
return(bdf)
}

run_viterbi_and_plot <- function(all_binned) {
	all_binned$ratio[is.na(all_binned$ratio)] = 0.5
	input = data.frame(1, all_binned$bin, all_binned$ratio)
	v = ViteRbi(input, active=F)
	all_binned$state = v[,4]
	lbin = split(all_binned, f=all_binned$chr)
	#par(mfrow=c(4,6), mar=c(2,2,2,2))
	#lapply(lbin, function(x) { plot(x$bin, x$ratio, main="", xlab="bin", ylab="ratio", cex=0.5, col=x$state+1) })
invisible(lbin)
}


get_breakpoint_poitions <- function(all_calls) {
	all_calls = do.call(rbind, call_list)
	all_calls[,1] = as.numeric(as.character(all_calls[,1]))
	all_calls[,2] = as.numeric(as.character(all_calls[,2]))
	all_calls[,3] = as.numeric(as.character(all_calls[,3]))
	all_calls= all_calls[order(all_calls[,1], all_calls[,2], all_calls[,3]),]

	index = paste(all_calls[,1], all_calls[,2], all_calls[,3], sep="-")
	u_events = all_calls[!duplicated(index),]

	chrs = unique(u_events[,1])
	breakpoints = list()
	for(x in 1:length(chrs)) {
		calls = u_events[u_events[,1]==chrs[x],]
		positions = c(calls[,2], calls[,3])
		positions = positions[order(positions)]
		breakpoints[[x]] = unique(positions)
	}
return(breakpoints)
}

get_breakpoint_genotypes <- function(breakpoints, chr, call_list) {
	print("processing")
	calls = c(1, breakpoints)
	sample_genos =lapply(call_list, function(sample_calls){
			li = list()
			for(y in 1:(length(calls)-1)) {
				t = table(sample_calls$state[sample_calls$chr==chr & sample_calls$stop>= calls[y] & sample_calls$start< calls[y+1]])
				major_state = names(t[which.max(t)])
				if(is.null(major_state)) {
					major_state = NA
				}
				li[[y]] = major_state
			}
			unlist(li)
		})
return(do.call(rbind, sample_genos))
}

realised_relationship_matrix <- function(genotypes) {
	genotypes = do.call(cbind, genotypes)
	genotypes = as.matrix(genotypes)
	genotypes = matrix(sapply(genotypes, as.numeric), ncol=ncol(genotypes))
	genotypes[genotypes==0] = -1
	genotypes[genotypes==1] = 0
	genotypes[genotypes==2] = 1
	genotypes[is.na(genotypes)] = 0
	rrm = genotypes%*%t(genotypes)
return(rrm)
}

get_positions <- function(breakpoints) {
	n = NULL
	for(x in 1:length(breakpoints)) {
		b = cbind(c(1, breakpoints[[x]][-(length(breakpoints[[x]]))]), breakpoints[[x]])
		n = rbind(n, cbind(x, b))
	}
return(n)
}



run_viterbi_and_plot <- function(all_binned) {
    all_binned$ratio[is.na(all_binned$ratio)] = 0.5
    input = data.frame(1, all_binned$bin, all_binned$ratio)
    v = ViteRbi(input, active=F)
    all_binned$state = v[,4]
    lbin = split(all_binned, f=all_binned$chr)
    #par(mfrow=c(4,6), mar=c(2,2,2,2))
    #lapply(lbin, function(x) { plot(x$bin, x$ratio, main="", xlab="bin", ylab="ratio", cex=0.5, col=x$state+1) })
    invisible(lbin)
}

collect_chunked_data <- function(dp_files, index1, index2, bin_len = 5000) {
    l = list()
    for(x in index1:index2) {
        d = read.table(unlist(dp_files)[x])
        d = d[d[,1]!="MT",]
        d[,1] = as.numeric(as.character(d[,1]))
        d = d[order(d[,1], d[,2], d[,3]),]
        all_binned = bin_all_chrs_and_plot(d, max_count=15, bin_len=bin_len)
        d = data.frame(unlist(dp_files)[x], all_binned)
        l[[x]] = d
        print(x)
    }
    all_binned = do.call(rbind, l)
    all_binned$ratio[is.na(all_binned$ratio)] = 0.5
    input = data.frame(1, all_binned$bin, all_binned$ratio)
    v = ViteRbi(input, active=F)
    all_binned$state = v[,4]
    #sbin = split(all_binned, f=all_binned$x)
    #for(s in 1:length(sbin)) {
        #lbin = split(sbin[[s]], f=sbin[[s]]$chr)
        #par(mfrow=c(4,6), mar=c(2,2,2,2))
        #lapply(lbin, function(x) { plot(x$bin, x$ratio, main="", xlab="bin", ylab="ratio", cex=0.5, col=x$state+1) })
     #   print(x)
    #}
return(all_binned)
}

