args <- commandArgs(trailingOnly = TRUE)
network <- toString(args[1])
chip_net <- toString(args[2])
pwm_net <- toString(args[3])
fn_out <- toString(args[4])
max_rank <- as.numeric(args[5])*1000
if (length(args) < 6) {
	num_bins <- 10
} else {
	num_bins <- as.integer(args[6])
}

net <- read.table(network,header=TRUE)
pdna.inter <- read.table(chip_net,header=TRUE)
bsinfo <- read.table(pwm_net,header=TRUE)

source("util.r")

chip.bp.np.setsizes <- c()
regulators_evaluated_cnt <- c()

# chip_desired_recall <- 0.33
# np_conf_cutoff <- 0.02036

orf_universe <- intersect(union(unique(bsinfo$TARGET),unique(pdna.inter$TARGET)),unique(net$TARGET))
reg_universe <- intersect(unique(as.character(pdna.inter$REGULATOR)),unique(as.character(bsinfo$REGULATOR)))
interaction_universe_cnt <- length(reg_universe) * length(orf_universe)

np_conf_cutoffs <- net$CONFIDENCE[order(net$CONFIDENCE,decreasing=TRUE)][seq(max_rank/num_bins,max_rank,by=max_rank/num_bins)]
for (np_conf_cutoff in np_conf_cutoffs)
{
	regulators_evaluated <- c()
	for (chip_desired_recall in c(0.1)) # c(0.1,0.33,0.5)
	{
		bp.cutoff <- c()
		pdna.target.count <- c()
		predicted.bp.support <- c()
		
		chip.bp.np.sets <- c()

		i <- 1
		for (regulator in intersect(unique(as.character(pdna.inter$REGULATOR)),unique(as.character(bsinfo$REGULATOR))) )
		{
			targets <- union(as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]),as.character(bsinfo$TARGET[which(bsinfo$REGULATOR==regulator)]))
			pdna.evid <- rep(0,times=length(targets))
			binding.evid.max <- rep(0,times=length(targets))
			binding.evid.sum <- rep(0,times=length(targets))
			pdna.evid[match(as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]),targets)] <- 1
			binding.evid.max[match(as.character(bsinfo$TARGET[which(bsinfo$REGULATOR==regulator)]),targets)] <- bsinfo$MAXP[which(bsinfo$REGULATOR==regulator)]
			binding.evid.sum[match(as.character(bsinfo$TARGET[which(bsinfo$REGULATOR==regulator)]),targets)] <- bsinfo$SUMP[which(bsinfo$REGULATOR==regulator)]
			# cat(regulator,length(targets),"\n")
			pdna.target.count[i] <- sum(pdna.evid)
			binding.evid.comb <- apply(rbind(binding.evid.max,binding.evid.sum),2,max)	
			prc.max <- compute.prc(binding.evid.max,pdna.evid)
			prc.sum <- compute.prc(binding.evid.sum,pdna.evid)
			prc.comb <- compute.prc(binding.evid.comb,pdna.evid)
		
			if 	(length(which(prc.comb[,1]<=chip_desired_recall))>0)
			{
				bp.cutoff[i] <- binding.evid.comb[sort.list(binding.evid.comb,decreasing=TRUE)][length(which(prc.comb[,1]<=chip_desired_recall))]
			}
			else
			{
				bp.cutoff[i] <- 1
			}
		
			predicted.targets <- as.character(net$TARGET[which(net$REGULATOR==regulator)])[which(net$CONFIDENCE[which(net$REGULATOR==regulator)]>np_conf_cutoff)]
			predicted.targets <- intersect(predicted.targets,orf_universe) # Narrow to space under consideration

			if (length(predicted.targets) > 0) {
                regulators_evaluated <- c(regulators_evaluated, regulator)
            }

			predicted.targets.chip.supported <- intersect(predicted.targets,as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]))
			predicted.targets.chip.unsupported <- setdiff(predicted.targets,as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]))

			bp.targets <- targets[which(binding.evid.comb>bp.cutoff[i])]
			pdna.targets <- targets[which(pdna.evid>0)]
			
			chip.bp.np.sets <- rbind(chip.bp.np.sets,c(length(bp.targets),sum(pdna.evid),length(predicted.targets),length(intersect(bp.targets,pdna.targets)),length(intersect(bp.targets,predicted.targets)),length(intersect(predicted.targets,pdna.targets)),length(intersect(bp.targets,intersect(predicted.targets,pdna.targets)))))

			i <- i + 1
		}		
		chip.bp.np.setsize <- apply(chip.bp.np.sets,2,sum)
		chip.bp.np.setsizes <- rbind(chip.bp.np.setsizes,chip.bp.np.setsize)
	}
	regulators_evaluated_cnt <- rbind(regulators_evaluated_cnt, length(regulators_evaluated))
}

write.table(cbind(chip.bp.np.setsizes, rep(interaction_universe_cnt, length=nrow(chip.bp.np.setsizes)), regulators_evaluated_cnt),file=fn_out,col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
