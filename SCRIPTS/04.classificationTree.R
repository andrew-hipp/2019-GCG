## make a nice tree
library(ape)
library(phytools)
library(openxlsx)
library(magrittr)
library(phyloch) # if necessary, install from source: http://www.christophheibl.de/phyloch_1.5-3.tar.gz
data(strat2012) # A data frame containing the stratigraphic chart by issued in 2012 by the International Commission on Stratigraphy.

## format tree
tr.plot <- list(
  edge = events$pr01$edge,
  edge.length = events$pr01$edge.length,
  Nnode = events$pr01$Nnode,
  tip.label = events$pr01$tip.label
)
class(tr.plot) <- 'phylo'

## read tables
dat.class <- read.csv('../DATA/dat.classification.edited.csv', row.names=1, as.is = T)
dat.class <- dat.class[tr.plot$tip.label, ]

dat.calib <- read.xlsx('../DATA/calibrations.xlsx', 1)
dat.calib$node <- apply(dat.calib[, c('tip1', 'tip2')], 1,
            function(x) getMRCA(tr.plot, as.character(x)))
## get descendents from each of the key nodes
if(!exists('changeDesc')) {
  changeDesc <-
  lapply(as.numeric(names(events.p$pr01)[which(events.p$pr01 > 0.3)]),
         function(changeNode) {
           cbind(sp = tr.plot$tip.label[descendants(tr.plot, changeNode)],
                 sect = dat.class[tr.plot$tip.label[descendants(tr.plot, changeNode)], 'SECTION.EDITED']
               )
         })
  names(changeDesc) <- events.p$pr01[which(events.p$pr01 > 0.3)]
} # close changeDesk

pdf('../OUT/FIG.trClassn.pdf', 12, 9, useDingbats = FALSE)
plot(tr.plot, show.tip.label = F, edge.width = 0.5)
nodelabels( node = dat.calib$node,
            text = paste(dat.calib$clade, '\n', dat.calib$min, '-', dat.calib$max, sep= ''),
            cex = 0.7, bg = 'black', col = 'white'
          )
tiplabels(pch = '-',
          col = dat.class$SUBGENUS %>% as.factor %>% as.numeric,
          offset = 0.2)
axisGeo(GTS = strat2012, ages=T, cex = 0.6)
legend(0, 400, legend = rev(unique(dat.class$SUBGENUS)), pch = 19,
        col = rev(unique(dat.class$SUBGENUS) %>% as.factor %>% as.numeric),
        bty = 'n',
        title = 'SECTION', title.adj = 0)
dev.off()
