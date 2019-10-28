arc.cladelabels2 <- function (tree = NULL, text, node = NULL, ln.offset = 1.02, lab.offset = 1.06,
    cex = 1, orientation = "curved", arcGap = 0, ...)
{
    obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (obj$type != "fan")
        stop("method works only for type=\"fan\"")
    h <- max(sqrt(obj$xx^2 + obj$yy^2))
    if (hasArg(mark.node))
        mark.node <- list(...)$mark.node
    else mark.node <- TRUE
    if (hasArg(interactive))
        interactive <- list(...)$interactive
    else {
        if (is.null(node))
            interactive <- TRUE
        else interactive <- FALSE
    }
    if (interactive)
        node <- getnode()
    if (hasArg(lwd))
        lwd <- list(...)$lwd
    else lwd <- par()$lwd
    if (hasArg(col))
        col <- list(...)$col
    else col <- par()$col
    if (hasArg(lend))
        lend <- list(...)$lend
    else lend <- par()$lend
    if (hasArg(clockwise))
        clockwise <- list(...)$clockwise
    else clockwise <- TRUE
    if (hasArg(n))
        n <- list(...)$n
    else n <- 0.05
    if (mark.node)
        points(obj$xx[node], obj$yy[node], pch = 21, bg = "red")
    if (is.null(tree)) {
        tree <- list(edge = obj$edge, tip.label = 1:obj$Ntip,
            Nnode = obj$Nnode)
        class(tree) <- "phylo"
    }
    d <- getDescendants(tree, node)
    d <- sort(d[d <= Ntip(tree)])
    deg <- atan(obj$yy[d]/obj$xx[d]) * 180/pi
    ii <- intersect(which(obj$yy[d] >= 0), which(obj$xx[d] <
        0))
    deg[ii] <- 180 + deg[ii]
    ii <- intersect(which(obj$yy[d] < 0), which(obj$xx[d] < 0))
    deg[ii] <- 180 + deg[ii]
    ii <- intersect(which(obj$yy[d] < 0), which(obj$xx[d] >=
        0))
    deg[ii] <- 360 + deg[ii]
    plotrix:::draw.arc(x = 0, y = 0, radius = ln.offset * h, deg1 = min(deg) + (arcGap/2),
        deg2 = max(deg)-(arcGap/2), lwd = lwd, col = col, lend = lend, n = n)
    if (orientation == "curved")
        plotrix:::arctext(text, radius = lab.offset * h, middle = mean(range(deg *
            pi/180)), cex = cex, clockwise = clockwise)
    else if (orientation == "horizontal") {
        x0 <- lab.offset * cos(median(deg) * pi/180) * h
        y0 <- lab.offset * sin(median(deg) * pi/180) * h
        text(x = x0, y = y0, label = text, adj = c(if (x0 >=
            0) 0 else 1, if (y0 >= 0) 0 else 1), offset = 0,
            cex = cex)
    }
}
