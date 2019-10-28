library(coda)
library(magrittr)
library(BAMMtools)

rN <- function(x) round(x, 3)

if(!exists('mcmcList')) {
	message('reading mcmc')
	mcmcList <- list(
	pr01 = read.csv('../DATA/cx.2019-04-24_mcmc_out.txt', header= T),
	pr50 = read.csv('../DATA/cx.2019-04-25.prior50_mcmc_out.txt', header = T)
	)
	burnstart <- floor(0.1 * nrow(mcmcList[[1]]))
	for(i in names(mcmcList)) mcmcList[[i]] <- mcmcList[[i]][burnstart:nrow(mcmcList[[i]]), ]
	}

if(!exists('tr.rates')) {
	message('making tree rates trees')
	tr.rates <- list(
		speciation =
			list(pr01 = getMeanBranchLengthTree(events$pr01, rate = 'speciation'),
				 pr50 = getMeanBranchLengthTree(events$pr50, rate = 'speciation')),
		extinction =
			list(pr01 = getMeanBranchLengthTree(events$pr01, rate = 'extinction'),
				 pr50 = getMeanBranchLengthTree(events$pr50, rate = 'extinction')),
		netDiv =
			list(pr01 = getMeanBranchLengthTree(events$pr01, rate = 'ndr'),
				 pr50 = getMeanBranchLengthTree(events$pr50, rate = 'ndr'))
		)# close list
	}# close if!exists

out <- timestamp(quiet = TRUE)

for(i in c('pr01', 'pr50')) {
	out <- c(out, paste('\nSTATS', i, '---------'))
	mL <- mcmcList[[i]]
	out <- c(out, paste('N_shifts ESS --', effectiveSize(mL$N_shifts) %>% rN))
	out <- c(out, paste('lnL ESS --', effectiveSize(mL$logLik) %>% rN))
	out <- c(out, paste('N shifts mean', mean(mL$N_shifts) %>% rN, '+/-', sd(mL$N_shifts) %>% rN))
	out <- c(out, paste('\nNet diversification rate ----'))
	out <- c(out, paste('\tRange:', paste(rN(range(tr.rates$netDiv[[i]]$phy$edge.length)), collapse = '--'), 'which is a',
		rN(max(tr.rates$netDiv[[i]]$phy$edge.length) / min(tr.rates$netDiv[[i]]$phy$edge.length)), '- fold range'))
	out <- c(out, paste('\tMean:', rN(mean(tr.rates$netDiv[[i]]$phy$edge.length)), '+/-',
							rN(sd(tr.rates$netDiv[[i]]$phy$edge.length)), '(s.d.)'
							))
	out <- c(out, paste('Speciation rate ----'))
	out <- c(out, paste('\tRange:', paste(rN(range(tr.rates$speciation[[i]]$phy$edge.length)), collapse = '--'), 'which is a',
		rN(max(tr.rates$speciation[[i]]$phy$edge.length) / min(tr.rates$speciation[[i]]$phy$edge.length)), '- fold range'))
	out <- c(out, paste('\tMean:', rN(mean(tr.rates$speciation[[i]]$phy$edge.length)), '+/-',
							rN(sd(tr.rates$speciation[[i]]$phy$edge.length)), '(s.d.)'
							))
	out <- c(out, paste('Extinction rate ----'))
	out <- c(out, paste('\tRange:', paste(rN(range(tr.rates$extinction[[i]]$phy$edge.length)), collapse = '--'), 'which is a',
		rN(max(tr.rates$extinction[[i]]$phy$edge.length) / min(tr.rates$extinction[[i]]$phy$edge.length)), '- fold range'))
	out <- c(out, paste('\tMean:', rN(mean(tr.rates$extinction[[i]]$phy$edge.length)), '+/-',
							rN(sd(tr.rates$extinction[[i]]$phy$edge.length)), '(s.d.)'
							))
	}

pdf('../OUT/SUPPLEMENT.rates.plot.pdf')

plot(tr.rates$netDiv$pr01$phy$edge.length,
	 tr.rates$netDiv$pr50$phy$edge.length,
	 xlab = 'Model averaged edge L, shifts prior = 1',
	 ylab = 'Model averaged edge L, shifts prior = 50')
dev.off()

out <- c(out, '---- STATS FOR BOTH ----')
out.linMod <- cor.test(tr.rates$netDiv$pr01$phy$edge.length,
						tr.rates$netDiv$pr50$phy$edge.length)
out <- c(out, paste('Prior 1 to prior 50 model-averaged edge lengths:',
					'\n\tr =', out.linMod$estimate %>% rN, '-- p =', out.linMod$p.value %>% rN))

message('Done with stats... writing to disk')
writeLines(out, format(Sys.time(), '../OUT/mcmcStats.%Y-%m-%d.%H.%M.txt'))
