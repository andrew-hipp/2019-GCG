library(magrittr)
library(ggtree)
library(ape)
library(phytools)
# taxi to airport: 6:45 p.m., order conf #78983

# set out what to do
doBars = T
doBioTree = T
doClassTree = F

# read and format data
dat.bio <- read.nexus.data('../DATA/carex.biogeo6.nex.edit2')
dat.bio <- dat.bio %>% as.data.frame %>% t
dat.bio <- apply(dat.bio, 1:2, as.numeric)
row.names(dat.bio) <- gsub('xxxx', '|', row.names(dat.bio), fixed = T)
dat.bio.key <- read.csv('../DATA/carex.biogeo6.key.csv', as.is = T)
dimnames(dat.bio)[[2]] <- dat.bio.key$area
dat.bio <- dat.bio[, order(colSums(dat.bio))]

tre.bio <- read.nexus('../DATA/carex.biogeo.tre.nex.edit2')
tre.bio$tip.label <- gsub('xxxx', '|', tre.bio$tip.label, fixed = T)

dat.bio.df <- as.data.frame(dat.bio)
for(i in names(dat.bio.df)) {
  dat.bio.df[[i]] <- ifelse(dat.bio.df[[i]] == 1, i, NA)
}

dat.class <- read.csv('../DATA/dat.classification.edited.csv', row.names=1, as.is = T)
dat.class <- dat.class[tr.plot$tip.label, ]

dat.class.mrca <- data.frame(
  Subgenus = unique(dat.class$SUBGENUS),
  node = sapply(unique(dat.class$SUBGENUS), function(i) {
    findMRCA(tr.plot, row.names(dat.class)[dat.class$SUBGENUS == i])
  })
  )
dat.class.mrca <-
  dat.class.mrca[dat.class.mrca$Subgenus != 'unplaced', ]
# make areas barplot
jpeg('../OUT/areaBars.jpg', 800, 600)
par(mar = c(5, 8, 4, 2) + 0.1)
barplot(colSums(dat.bio),
        horiz = TRUE, las = 2,
        col = rev(c('North America' = '#ec2124',
                    'Asia' = '#8bc63e',
                    'Europe' = '#2aaae3',
                    'Pacific' = '#119247',
                    'South America' = '#f9931f',
                    'Africa' = '#ffce52')), # matched to Santi fig 2019-07-28
        xlab = 'Species sampled'
      ) # close barplot
dev.off()

## plot tree
if(doBioTree) {
  p <- ggtree(tr.plot, layout = 'fan', open.angle = 12)
  p <- gheatmap(p, dat.bio.df,
                width = 0.2,
                colnames_angle = 270,
                font.size = 7,
                hjust = 0
              )
  p <- p + scale_fill_manual(
        values = c('North America' = '#ec2124',
                    'Asia' = '#8bc63e',
                    'Europe' = '#2aaae3',
                    'Pacific' = '#119247',
                    'South America' = '#f9931f',
                    'Africa' = '#ffce52')
                  )
  p <- p + theme(legend.position = 'none')
  for(i in 1:dim(dat.class.mrca)[1]) {
    message(paste('adding geom_cladelabel for', i))
    p <- p + geom_cladelabel(node=dat.class.mrca[i, 'node'],
                               label=dat.class.mrca[i, 'Subgenus'],
                               align=T,
                               color = c(paste('gray', i+1, '0', sep = ''),
                                        'black'),
                               offset = 8,
                               #angle=270,
                               hjust='center', offset.text = 4, barsize=4,
                               fontsize = 9
                              )
  }
  jpeg('../OUT/biogeogTree.jpg', 3200, 2400)
  print(p)
  dev.off()
  # print(p)
} # close doBioTree
