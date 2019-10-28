library(magrittr)
#library(ggtree)
library(ape)
library(phytools)
library(BAMMtools)
library(openxlsx)
source('../SCRIPTS/arc.cladelabels2.R')
# Components desired:
# 1. BAMM output
# 2. clade labels for all subgenera
# 3. clade labels for focus diversification clades
# 4. fossils at internal nodes
# 5 -- optional,

### TO MAKE THIS WORK:
### You have to use trace(plot.bammdata, edit = T) and change
###    plot.window(xlim = c(-1, 1) + c(-rb, rb) + c(-ofs,
### to
###    plot.window(xlim = c(-1.5, 1.4) + c(-rb, rb) + c(-ofs,

## Global variables
useFullFossilLabels <- FALSE # use full fossil labels or abbreviations
showFossils <- TRUE # show fossil calibrations at all?
showPP <- FALSE # show posteriors on shift nodes
doPDF <- TRUE # output to pdf?
roundTo <- 5
ppMin <- 0.1
sectOffset <- 1.21
bioOffset <- seq(from = 0.015, to = 0.105, by = 0.018)
names(bioOffset) <- c("Africa",
                     "South America",
                     "Pacific",
                     "Europe",
                     "Asia",
                     "North America")

bioCol = rev(c('North America' = '#ec2124',
            'Asia' = '#8bc63e',
            'Europe' = '#2aaae3',
            'Pacific' = '#119247',
            'South America' = '#f9931f',
            'Africa' = '#ffce52'))

## read and format tree
if(!(exists('events') & exists('tr.plot'))) {
  stop('Run scripts 02-readCarexBAMM first')
}

## read and format tree data
dat.bio <- read.nexus.data('../DATA/carex.biogeo6.nex.edit2')
dat.bio <- dat.bio %>% as.data.frame %>% t
dat.bio <- apply(dat.bio, 1:2, as.numeric)
row.names(dat.bio) <- gsub('xxxx', '|', row.names(dat.bio), fixed = T)
dat.bio.key <- read.csv('../DATA/carex.biogeo6.key.csv', as.is = T)
dimnames(dat.bio)[[2]] <- dat.bio.key$area
dat.bio <- dat.bio[, order(colSums(dat.bio))]

dat.bio.df <- as.data.frame(dat.bio)[tr.plot$tip.label, ]
for(i in names(dat.bio.df)) {
  dat.bio.df[[i]] <- ifelse(dat.bio.df[[i]] == 1, i, NA)
}

dat.calib <- read.xlsx('../DATA/calibrations.xlsx', 1)
dat.calib$node <- apply(dat.calib[, c('tip1', 'tip2')], 1,
            function(x) getMRCA(tr.plot, as.character(x)))

dat.class <- read.csv('../DATA/dat.classification.edited.csv', row.names=1, as.is = T)
dat.class <- dat.class[tr.plot$tip.label, ]

dat.changeSects <- read.xlsx('../DATA/changeSections.xlsx', 1)
dat.changeSects$node <- names(events.p$pr01)[match(dat.changeSects$pp, round(events.p$pr01, 5))]

dat.class.mrca <- data.frame(
  Subgenus = unique(dat.class$SUBGENUS),
  node = sapply(unique(dat.class$SUBGENUS), function(i) {
    findMRCA(tr.plot, row.names(dat.class)[dat.class$SUBGENUS == i])
  })
  )
dat.class.mrca <-
  dat.class.mrca[dat.class.mrca$Subgenus != 'unplaced', ]

## 0. format output
if(doPDF)
  pdf('../OUT/FIG02-v2.JSE.pdf', 6.69, 8.86)
  par(mar = c(0.1, 0.1, 0.1, 0.1))

## 1. plot BAMM
eplot = plot.bammdata(events$pr01, labels = F, spex= 'netdiv',
            method = 'polar', xlim = c(-10, 10),
            par.reset = FALSE)
if(showPP) {
  nodelabels(node = as.numeric(names(events.p$pr01)[which(events.p$pr01 > ppMin)]), pch = 21, bg = 'black', cex = 3)
  nodelabels(node = as.numeric(names(events.p$pr01)[which(events.p$pr01 > ppMin)]),
       text = as.character(round(events.p$pr01[which(events.p$pr01 > ppMin)], roundTo)),
       cex = 0.3,
       frame = 'n',
       col = 'white')
     }

addBAMMlegend(eplot, location = c(-1.48, -1.45, -1.75, -1.48), cex.axis = 0.4)
text(-1.54, -1.38, "Net diversification\nrate",
      pos = 4, cex = 0.6)

## 2. Add subgenera
for(i in 1:dim(dat.class.mrca)[1])
    arc.cladelabels2(tree = tr.plot,
        dat.class.mrca$Subgenus[i] %>% as.character,
        dat.class.mrca$node[i] %>% as.character,
        orientation= if(dat.class.mrca$Subgenus[i] %in% c("Siderosticta","Psyllophora"))
        "horizontal" else "curved",
        ln.offset = sectOffset - 0.07,
        lab.offset = sectOffset - 0.04,
        cex = 0.6,
        lwd = 2.5,
        arcGap = 1,
        mark.node=FALSE)

## 3. Add sections
for(i in 1:dim(dat.changeSects)[1])
    arc.cladelabels2(tree = tr.plot,
        dat.changeSects$label[i] %>% as.character,
        dat.changeSects$node[i] %>% as.character,
        orientation = "horizontal",
        lwd = 1.5,
        col = 'gray',
        ln.offset = sectOffset,
        lab.offset = sectOffset + 0.02,
        arcGap = 1,
        cex = 0.5,
        mark.node=TRUE)

## 4. Add fossils
if(showFossils) {
    if(useFullFossilLabels) {
    fossilLabs <- paste(dat.calib$clade, '\n', dat.calib$min, '-', dat.calib$max, sep= '')
  } else fossilLabs <- dat.calib$abbreviation
  nodelabels( node = dat.calib$node,
            pch = 19,
              cex = 2, col = 'black'
            )
  nodelabels( node = dat.calib$node,
              text = fossilLabs,
              cex = 0.4, frame = 'n', col = 'white'
            )
  }

## 5. add biogeography
assign("last_plot.phylo",
        c(get("last_plot.phylo", envir = .PlotPhyloEnv), align.tip.label = FALSE),
        envir = .PlotPhyloEnv)
for(i in names(dat.bio.df)) {
  bioWhich <- which(!is.na(dat.bio.df[[i]]))
  tiplabels(
            # text = '-',
            pch = 19,
            cex = 0.27,
            col = bioCol[i],
            tip = bioWhich,
            offset = bioOffset[i],
            frame = 'none'
            )
}
## 6. add legend
legend(-1, -1.8, yjust = 0, names(bioCol), pch = 19, col = bioCol,
      bty = 'n', cex = 0.5, pt.cex = 0.6)
text(-1, -1.38, "Species\ndistributions",
    pos = 4, cex = 0.6)
## 99. close graphics device
if(doPDF) dev.off()
