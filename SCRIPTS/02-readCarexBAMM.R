library(ape)
library(BAMMtools)
library(ggplot2)
library(ggtree)
library(phyloch)
data(strat2012) # A data frame containing the stratigraphic chart by issued in 2012 by the International Commission on Stratigraphy.

plotBAMM = TRUE

tr <- read.tree('../DATA/Carex10Fossil.divTree.tre')
if(!exists('events')) {
	message('... reading in data sets ...')
	events <- list(pr01 = getEventData(tr, dir('../DATA/', patt = '24_event_data.txt', full = T), burnin = 0.1),
				   pr50 = getEventData(tr, dir('../DATA/', patt = 'prior50_event_data.txt', full = T), burnin = 0.1)
				   )

	events.table <- events.p <- list(pr01 = NA, pr50 = NA)
	for(i in c('pr01', 'pr50')) {
		message(paste('... starting', i, '...'))
		message('... doing data summaries...')
		events.table[[i]] <- do.call('rbind', events[[i]]$eventData)
		events.table[[i]] <- events.table[[i]][!events.table[[i]]$time == 0, ] # gets rid of basalmost node
		events.p[[i]] <-
		  (table(events.table[[i]][, 'node']) / length(events[[i]]$eventData)) %>%
		  sort(decreasing = TRUE) # probably of changes by node
		} # close i
	} # close making events

if(plotBAMM) {
	for(i in c('pr01', 'pr50')) {
		message(paste('... plotting tree', i, '...'))
		pdf(paste('../OUT/cx.phylo_10-fossils-BAMM.', i, '.polar.pdf', sep = ''), 8.5, 11)
		#eplot = plot(events, labels = F)
		eplot = plot(events[[i]], labels = F, method = 'polar')
		nodelabels(node = as.numeric(names(events.p[[i]])[which(events.p[[i]] > 0.3)]), pch = 21, bg = 'black', cex = 3)
		nodelabels(node = as.numeric(names(events.p[[i]])[which(events.p[[i]] > 0.3)]),
				   text = as.character(round(events.p[[i]][which(events.p[[i]] > 0.3)], 3)),
				   cex = 0.5,
				   frame = 'n',
				   col = 'white')

		#addBAMMshifts(events, par.reset=FALSE, cex=2)
		addBAMMlegend(eplot, location = "bottomleft", cex.axis = 0.4)
		#text(-1, mean(c(30, 250)), "Net diversification rate", cex = 0.55, srt = 90)
		#add.geoscale(tr)
		#axisGeo(GTS = strat2012, ages=T, cex = 0.6)
		dev.off()

		message('... plotting RTT...')
		pdf(paste('../OUT/cx.phylo_10-fossils-netdivTT.', i, '.pdf', sep = ''), 8.5, 11)
		plotRateThroughTime(events[[i]], ratetype = 'netdiv')
		dev.off()
		} # close i
	} # close plotBAMM
