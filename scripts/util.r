compute.prc <- function(q.adj,t.adj)
{
	q.nzi <- which(q.adj!=0)
	q.nz <- q.adj[q.nzi]
	t.nz <- t.adj[q.nzi]
	q.rnk <- t.nz[sort.list(q.nz,decreasing=TRUE)]
	n <- length(q.rnk)
	pr <- matrix(0,ncol=2,nrow=1) 
	if (length(q.rnk)>0)
	{
		pr <- matrix(0,ncol=2,nrow=n)
		tp <- 0
		p <- sum(t.adj)
		for (i in 1:n)
		{
			tp <- tp + q.rnk[i]
			pr[i,2] <- tp/i
			pr[i,1] <- tp/p
		}
	}
	pr
}

compute.auc.pr <- function(q.adj,t.adj)
{
	auc <- 0
	pr <- compute.prc(q.adj,t.adj)
	lastr <- 0;
	lastp <- 1;
	for (i in 1:dim(pr)[1])
	{
		if (pr[i,1]>lastr)
		{
				auc <- auc + min(pr[i,2],lastp) * (pr[i,1] - lastr)
				lastr <- pr[i,1]
		}
		lastp <- pr[i,2]
	}
	if (lastr < 1)
	{
		guess.remain <- length(t.adj) - dim(pr)[1]
		t.remain <- sum(t.adj) * (1-lastr)
		auc <- auc + (t.remain / guess.remain) * (1-lastr)
	}
auc
}


Pause <- function () {
    cat("Hit <enter> to continue...")
    readline()
    invisible()
}
