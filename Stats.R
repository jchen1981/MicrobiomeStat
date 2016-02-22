
vkey2 <- function (map, title = NA, side = 2, stretch = 1.4, x, y, 
		wh)
{
	opar <- par(xpd = NA)
	on.exit(par(opar))
	n <- length(map$breaks) + 1
	dy <- strheight("A")
	aspect <- diff(grconvertX(1:2, from = "inches"))/diff(grconvertY(1:2,
					from = "inches"))
	dx <- dy * aspect
	if (missing(wh)) {
		
		wh <- 1:(n-1)
		
	}
	labs <- format(map$breaks[wh])
	maxlabwidth <- max(strwidth(labs))
	if (missing(x)) {
		x <- grconvertX(1, from = "nfc") - (2 * dx)
		if (side == 4)
			x <- x - maxlabwidth - dx
	}
	else {
		if (is.list(x)) {
			y <- x$y
			x <- x$x
		}
	}
	if (missing(y))
		y <- par("usr")[3] + dy
	ybord <- y + ((0:(n - 1)) * dy * stretch)
	rect(x, ybord[-n], x + dx, ybord[-1], col = map$colors, border = NA)
	if (side == 4) {
		xtext <- x + dx
		text(x = x, y = ybord[n] + (1.5 * dy), title, adj = c(0,
						0))
	}
	if (side == 2) {
		xtext <- x
		text(x = x + dx, y = ybord[n] + (1.5 * dy), title, adj = c(1,
						0))
	}
	text(x = xtext, y = ybord[wh] + 0.5 * dy, labels = labs, pos = side)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}


heatmap.3 <- function(x,
		Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
		distfun = dist,
		hclustfun = hclust,
		dendrogram = c("both","row", "column", "none"),
		symm = FALSE,
		scale = c("none","row", "column"),
		na.rm = TRUE,
		revC = identical(Colv,"Rowv"),
		add.expr,
		breaks,
		symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
		col = "heat.colors",
		colsep,
		rowsep,
		sepcolor = "white",
		sepwidth = c(0.05, 0.05),
		cellnote,
		notecex = 1,
		notecol = "cyan",
		na.color = par("bg"),
		trace = c("none", "column","row", "both"),
		tracecol = "cyan",
		hline = median(breaks),
		vline = median(breaks),
		linecol = tracecol,
		margins = c(5,5),
		ColSideColors = NULL,
		RowSideColors = NULL,
		side.height.fraction=0.3,
		cexRow = 0.2 + 1/log10(nr),
		cexCol = 0.2 + 1/log10(nc),
		labRow = NULL,
		labCol = NULL,
		key = TRUE,
		keysize = 1.5,
		density.info = c("none", "histogram", "density"),
		denscol = tracecol,
		symkey = max(x < 0, na.rm = TRUE) || symbreaks,
		densadj = 0.25,
		main = NULL,
		xlab = NULL,
		ylab = NULL,
		lmat = NULL,
		lhei = NULL,
		lwid = NULL,
		NumColSideColors = 1,
		NumRowSideColors = 1,
		KeyValueName="Value",...){
	
	invalid <- function (x) {
		if (missing(x) || is.null(x) || length(x) == 0)
			return(TRUE)
		if (is.list(x))
			return(all(sapply(x, invalid)))
		else if (is.vector(x))
			return(all(is.na(x)))
		else return(FALSE)
	}
	
	x <- as.matrix(x)
	scale01 <- function(x, low = min(x), high = max(x)) {
		x <- (x - low)/(high - low)
		x
	}
	retval <- list()
	scale <- if (symm && missing(scale))
				"none"
			else match.arg(scale)
	dendrogram <- match.arg(dendrogram)
	trace <- match.arg(trace)
	density.info <- match.arg(density.info)
	if (length(col) == 1 && is.character(col))
		col <- get(col, mode = "function")
	if (!missing(breaks) && (scale != "none"))
		warning("Using scale=\"row\" or scale=\"column\" when breaks are",
				"specified can produce unpredictable results.", "Please consider using only one or the other.")
	if (is.null(Rowv) || is.na(Rowv))
		Rowv <- FALSE
	if (is.null(Colv) || is.na(Colv))
		Colv <- FALSE
	else if (Colv == "Rowv" && !isTRUE(Rowv))
		Colv <- FALSE
	if (length(di <- dim(x)) != 2 || !is.numeric(x))
		stop("`x' must be a numeric matrix")
	nr <- di[1]
	nc <- di[2]
	if (nr <= 1 || nc <= 1)
		stop("`x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2)
		stop("`margins' must be a numeric vector of length 2")
	if (missing(cellnote))
		cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
	if (!inherits(Rowv, "dendrogram")) {
		if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
					c("both", "row"))) {
			if (is.logical(Colv) && (Colv))
				dendrogram <- "column"
			else dedrogram <- "none"
			warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
					dendrogram, "'. Omitting row dendogram.")
		}
	}
	if (!inherits(Colv, "dendrogram")) {
		if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
					c("both", "column"))) {
			if (is.logical(Rowv) && (Rowv))
				dendrogram <- "row"
			else dendrogram <- "none"
			warning("Discrepancy: Colv is FALSE, while dendrogram is `",
					dendrogram, "'. Omitting column dendogram.")
		}
	}
	if (inherits(Rowv, "dendrogram")) {
		ddr <- Rowv
		rowInd <- order.dendrogram(ddr)
	}
	else if (is.integer(Rowv)) {
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd))
			stop("row dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Rowv)) {
		Rowv <- rowMeans(x, na.rm = na.rm)
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd))
			stop("row dendrogram ordering gave index of wrong length")
	}
	else {
		rowInd <- nr:1
	}
	if (inherits(Colv, "dendrogram")) {
		ddc <- Colv
		colInd <- order.dendrogram(ddc)
	}
	else if (identical(Colv, "Rowv")) {
		if (nr != nc)
			stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
		if (exists("ddr")) {
			ddc <- ddr
			colInd <- order.dendrogram(ddc)
		}
		else colInd <- rowInd
	}
	else if (is.integer(Colv)) {
		hcc <- hclustfun(distfun(if (symm)
									x
								else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd))
			stop("column dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Colv)) {
		Colv <- colMeans(x, na.rm = na.rm)
		hcc <- hclustfun(distfun(if (symm)
									x
								else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd))
			stop("column dendrogram ordering gave index of wrong length")
	}
	else {
		colInd <- 1:nc
	}
	retval$rowInd <- rowInd
	retval$colInd <- colInd
	retval$call <- match.call()
	x <- x[rowInd, colInd]
	x.unscaled <- x
	cellnote <- cellnote[rowInd, colInd]
	if (is.null(labRow))
		labRow <- if (is.null(rownames(x)))
					(1:nr)[rowInd]
				else rownames(x)
	else labRow <- labRow[rowInd]
	if (is.null(labCol))
		labCol <- if (is.null(colnames(x)))
					(1:nc)[colInd]
				else colnames(x)
	else labCol <- labCol[colInd]
	if (scale == "row") {
		retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
		x <- sweep(x, 1, rm)
		retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
		x <- sweep(x, 1, sx, "/")
	}
	else if (scale == "column") {
		retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
		x <- sweep(x, 2, rm)
		retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
		x <- sweep(x, 2, sx, "/")
	}
	if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
		if (missing(col) || is.function(col))
			breaks <- 16
		else breaks <- length(col) + 1
	}
	if (length(breaks) == 1) {
		if (!symbreaks)
			breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
					length = breaks)
		else {
			extreme <- max(abs(x), na.rm = TRUE)
			breaks <- seq(-extreme, extreme, length = breaks)
		}
	}
	nbr <- length(breaks)
	ncol <- length(breaks) - 1
	if (class(col) == "function")
		col <- col(ncol)
	min.breaks <- min(breaks)
	max.breaks <- max(breaks)
	x[x < min.breaks] <- min.breaks
	x[x > max.breaks] <- max.breaks
	if (missing(lhei) || is.null(lhei))
		lhei <- c(keysize, 4)
	if (missing(lwid) || is.null(lwid))
		lwid <- c(keysize, 4)
	if (missing(lmat) || is.null(lmat)) {
		lmat <- rbind(4:3, 2:1)
		
		if (!is.null(ColSideColors)) {
			#if (!is.matrix(ColSideColors))
			#stop("'ColSideColors' must be a matrix")
			if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
				stop("'ColSideColors' must be a matrix of nrow(x) rows")
			lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
			#lhei <- c(lhei[1], 0.2, lhei[2])
			lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
		}
		
		if (!is.null(RowSideColors)) {
			#if (!is.matrix(RowSideColors))
			#stop("'RowSideColors' must be a matrix")
			if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
				stop("'RowSideColors' must be a matrix of ncol(x) columns")
			lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
			#lwid <- c(lwid[1], 0.2, lwid[2])
			lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
		}
		lmat[is.na(lmat)] <- 0
	}
	
	if (length(lhei) != nrow(lmat))
		stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
	if (length(lwid) != ncol(lmat))
		stop("lwid must have length = ncol(lmat) =", ncol(lmat))
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	
	if (!is.null(RowSideColors)) {
		if (!is.matrix(RowSideColors)){
			par(mar = c(margins[1], 0, 0, 0.5))
			image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
			box()
		} else {
			par(mar = c(margins[1], 0, 0, 0.5))
			rsc = t(RowSideColors[,rowInd, drop=F])
			rsc.colors = matrix()
			rsc.names = names(table(rsc))
			rsc.i = 1
			for (rsc.name in rsc.names) {
				rsc.colors[rsc.i] = rsc.name
				rsc[rsc == rsc.name] = rsc.i
				rsc.i = rsc.i + 1
			}
			rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
			image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
			
			if (length(rownames(RowSideColors)) > 0) {
				axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
			}
			box()
		}
	}
	
	if (!is.null(ColSideColors)) {
		
		if (!is.matrix(ColSideColors)){
			par(mar = c(0.5, 0, 0, margins[2]))
			image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
			box()
		} else {
			par(mar = c(0.5, 0, 0, margins[2]))
			csc = ColSideColors[colInd, , drop=F]
			csc.colors = matrix()
			csc.names = names(table(csc))
			csc.i = 1
			for (csc.name in csc.names) {
				csc.colors[csc.i] = csc.name
				csc[csc == csc.name] = csc.i
				csc.i = csc.i + 1
			}
			csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
			image(csc, col = as.vector(csc.colors), axes = FALSE)
			if (length(colnames(ColSideColors)) > 0) {
				axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
			}
			box()
		}
	}
	
	par(mar = c(margins[1], 0, 0, margins[2]))
	x <- t(x)
	cellnote <- t(cellnote)
	if (revC) {
		iy <- nr:1
		if (exists("ddr"))
			ddr <- rev(ddr)
		x <- x[, iy]
		cellnote <- cellnote[, iy]
	}
	else iy <- 1:nr
	image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
	retval$carpet <- x
	if (exists("ddr"))
		retval$rowDendrogram <- ddr
	if (exists("ddc"))
		retval$colDendrogram <- ddc
	retval$breaks <- breaks
	retval$col <- col
	if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
		mmat <- ifelse(is.na(x), 1, NA)
		image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
				col = na.color, add = TRUE)
	}
	axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
			cex.axis = cexCol)
	if (!is.null(xlab))
		mtext(xlab, side = 1, line = margins[1] - 1.25)
	axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
			cex.axis = cexRow)
	if (!is.null(ylab))
		mtext(ylab, side = 4, line = margins[2] - 1.25)
	if (!missing(add.expr))
		eval(substitute(add.expr))
	if (!missing(colsep))
		for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
	if (!missing(rowsep))
		for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
	min.scale <- min(breaks)
	max.scale <- max(breaks)
	x.scaled <- scale01(t(x), min.scale, max.scale)
	if (trace %in% c("both", "column")) {
		retval$vline <- vline
		vline.vals <- scale01(vline, min.scale, max.scale)
		for (i in colInd) {
			if (!is.null(vline)) {
				abline(v = i - 0.5 + vline.vals, col = linecol,
						lty = 2)
			}
			xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
			xv <- c(xv[1], xv)
			yv <- 1:length(xv) - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (trace %in% c("both", "row")) {
		retval$hline <- hline
		hline.vals <- scale01(hline, min.scale, max.scale)
		for (i in rowInd) {
			if (!is.null(hline)) {
				abline(h = i + hline, col = linecol, lty = 2)
			}
			yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
			yv <- rev(c(yv[1], yv))
			xv <- length(yv):1 - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (!missing(cellnote))
		text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
				col = notecol, cex = notecex)
	par(mar = c(margins[1], 0, 0, 0))
	if (dendrogram %in% c("both", "row")) {
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	}
	else plot.new()
	par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
	if (dendrogram %in% c("both", "column")) {
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	}
	else plot.new()
	if (!is.null(main))
		title(main, cex.main = 1.5 * op[["cex.main"]])
	if (key) {
		par(mar = c(5, 4, 2, 1), cex = 0.75)
		tmpbreaks <- breaks
		if (symkey) {
			max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
			min.raw <- -max.raw
			tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
			tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
		}
		else {
			min.raw <- min(x, na.rm = TRUE)
			max.raw <- max(x, na.rm = TRUE)
		}
		
		z <- seq(min.raw, max.raw, length = length(col))
		image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
				xaxt = "n", yaxt = "n")
		par(usr = c(0, 1, 0, 1))
		lv <- pretty(breaks)
		xv <- scale01(as.numeric(lv), min.raw, max.raw)
		axis(1, at = xv, labels = lv)
		if (scale == "row")
			mtext(side = 1, "Row Z-Score", line = 2)
		else if (scale == "column")
			mtext(side = 1, "Column Z-Score", line = 2)
		else mtext(side = 1, KeyValueName, line = 2)
		if (density.info == "density") {
			dens <- density(x, adjust = densadj, na.rm = TRUE)
			omit <- dens$x < min(breaks) | dens$x > max(breaks)
			dens$x <- dens$x[-omit]
			dens$y <- dens$y[-omit]
			dens$x <- scale01(dens$x, min.raw, max.raw)
			lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
					lwd = 1)
			axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
			title("Color Key\nand Density Plot")
			par(cex = 0.5)
			mtext(side = 2, "Density", line = 2)
		}
		else if (density.info == "histogram") {
			h <- hist(x, plot = FALSE, breaks = breaks)
			hx <- scale01(breaks, min.raw, max.raw)
			hy <- c(h$counts, h$counts[length(h$counts)])
			lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
					col = denscol)
			axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
			title("Color Key\nand Histogram")
			par(cex = 0.5)
			mtext(side = 2, "Count", line = 2)
		}
		else title("Color Key")
	}
	else plot.new()
	retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
			high = retval$breaks[-1], color = retval$col)
	invisible(retval)
}

load_package <- function (package.names=c('MASS', 'ggbiplot', 'vegan', 'pscl', 'glmmADMB', 'aod', 'nlme', 'MiRKAT', 'matrixStats', 'gplots', 'scales', 'ggplot2', 'GUniFrac', 'rpart', 'qvalue', 'DESeq2',
				'phangorn', 'phyloseq', 'RColorBrewer', 'squash', 'rhdf5', 'biom', 'reshape', 'randomForest', 'Boruta', 'ade4')) {
	for (package.name in package.names) {
		require(package.name, character.only=T)
	}
}

getPermuteMatrixBlock <- function(permutations, strata) {
	strata <- factor(strata)
	strata <- factor(as.numeric(strata))
	res <- sapply(1:permutations, function(i) {
				strata1 <- strata
				levels(strata1) <- sample(levels(strata1))
				order(as.numeric(as.character(strata1)))
			})
	t(res)
}

# adonis2: permutate the covariate instead of data matrix - validated
adonis2 <- function (formula, data = NULL, permutations = 999, method = "bray", 
		strata = NULL, block.perm = TRUE, contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
		...) 
{
	TOL <- 1e-07
	Terms <- terms(formula, data = data, keep.order=TRUE)
	lhs <- formula[[2]]
	lhs <- eval(lhs, data, parent.frame())
	formula[[2]] <- NULL
	rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
	op.c <- options()$contrasts
	options(contrasts = c(contr.unordered, contr.ordered))
	rhs <- model.matrix(formula, rhs.frame)
	options(contrasts = op.c)
	grps <- attr(rhs, "assign")
	qrhs <- qr(rhs)
	rhs <- rhs[, qrhs$pivot, drop = FALSE]
	rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
	grps <- grps[qrhs$pivot][1:qrhs$rank]
	u.grps <- unique(grps)
	nterms <- length(u.grps) - 1
	st <- ifelse(nterms == 1, 2, nterms)
	H.s <- lapply(st:(nterms+1), function(j) {
				Xj <- rhs[, grps %in% u.grps[1:j]]
				qrX <- qr(Xj, tol = TOL)
				Q <- qr.Q(qrX)
				tcrossprod(Q[, 1:qrX$rank])
			})
	if (inherits(lhs, "dist")) {
		if (any(lhs < -TOL)) 
			stop("dissimilarities must be non-negative")
		dmat <- as.matrix(lhs^2)
	}
	else {
		dist.lhs <- as.matrix(vegdist(lhs, method = method, ...))
		dmat <- dist.lhs^2
	}
	n <- nrow(dmat)
	G <- -sweep(dmat, 1, rowMeans(dmat))/2
	SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
	if (nterms == 1) {
		SS.Exp.each <- SS.Exp.comb
	} else {
		SS.Exp.each <- SS.Exp.comb[2] - SS.Exp.comb[1]
	}
	
	H.snterm <- H.s[[length(H.s)]]
	tIH.snterm <- t(diag(n) - H.snterm)
	
#	if (length(H.s) > 1) {
#		H.s[[2]] <- H.s[[2]] - H.s[[1]]
#	}
	
	SS.Res <- sum(G * tIH.snterm)
	df.Exp <- sum(grps == u.grps[length(u.grps)])
	df.Res <- n - qrhs$rank
	
	if (inherits(lhs, "dist")) {
		beta.sites <- qr.coef(qrhs, as.matrix(lhs))
		beta.spp <- NULL
	}
	else {
		beta.sites <- qr.coef(qrhs, dist.lhs)
		beta.spp <- qr.coef(qrhs, as.matrix(lhs))
	}
	colnames(beta.spp) <- colnames(lhs)
	colnames(beta.sites) <- rownames(lhs)
	
	F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
	
	f.test <- function(tH, G, df.Exp, df.Res, tIH.snterm) {
		(sum(G * tH)/df.Exp)/(sum(G * tIH.snterm)/df.Res)
	}
	
	rhs.1 <- rhs[, grps %in% u.grps[1:nterms], drop=F]
	rhs.2 <- rhs[, grps %in% u.grps[(nterms+1)], drop=F]	
	
	if (missing(strata)) strata <- NULL
	
	if (block.perm == FALSE) {
		permute.ind <- vegan:::getPermuteMatrix(permutations, n, strata = strata)
	} else {
		if (is.null(strata)) stop('Block permutation requires strata!\n')
		strata.u <- unique(strata)
		reorder.ind <- unlist(lapply(strata.u, function(x) which(strata == x)))
		expand.ind <- rep(1:length(strata.u), sapply(strata.u, function(x) sum(strata == x)))
		rhs.2.u <- rhs.2[sapply(strata.u, function(x) which(strata == x)[1]), , drop=F]
	}
	
	# Not right
	f.perms <- 
			sapply(1:permutations, function(i) {
						if (block.perm == FALSE) {
							rhs.2 <- rhs[permute.ind[i, ], grps %in% u.grps[(nterms+1)], drop=F]	
						} else {
							rhs.2[reorder.ind, ] <- rhs.2.u[sample(nrow(rhs.2.u)), , drop=F][expand.ind, , drop=F]
						}		
						Xj <- cbind(rhs.1, rhs.2)
						qrX <- qr(Xj, tol = TOL)
						Q <- qr.Q(qrX)
						tH.snterm <- t(tcrossprod(Q[, 1:qrX$rank]))
						tIH.snterm <- diag(n) - tH.snterm
						
						if (nterms > 1) {
							tH.snterm <- tH.snterm - t(H.s[[1]])
						}							
						f.test(tH.snterm, G, df.Exp, df.Res, 
								tIH.snterm)
					}
			)
	
	f.perms <- round(f.perms, 12)
	F.Mod <- round(F.Mod, 12)
	SumsOfSqs = c(SS.Exp.each, SS.Res, SS.Exp.comb[length(SS.Exp.comb)] + SS.Res)
	tab <- data.frame(Df = c(df.Exp, df.Res, n - 1), SumsOfSqs = SumsOfSqs, 
			MeanSqs = c(SS.Exp.each/df.Exp, SS.Res/df.Res, NA), F.Model = c(F.Mod, 
					NA, NA), R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)], 
			P = c((sum(f.perms >= F.Mod) + 1)/(permutations + 
								1), NA, NA))
	rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps[length(u.grps)]], 
			"Residuals", "Total")
	colnames(tab)[ncol(tab)] <- "Pr(>F)"
	attr(tab, "heading") <- "Terms added sequentially (first to last)\n"
	class(tab) <- c("anova", class(tab))
	out <- list(aov.tab = tab, call = match.call(), coefficients = beta.spp, 
			coef.sites = beta.sites, f.perms = as.matrix(f.perms), model.matrix = rhs, 
			terms = Terms)
	class(out) <- "adonis"
	out
}

# Test 1 - validated
# pv1 <- pv2 <- numeric(100)
# for ( i in 1:100) {
#   x1 <- rnorm(20); x2 <- x1 + rnorm(20); y <- 0.5*x1 + x2 + rnorm(20);
#   pv1[i] <- adonis(dist(y) ~ x1 + x2)$aov.tab[2, 6]
#   pv2[i] <- adonis2(dist(y) ~ x1 + x2, block.perm=F)$aov.tab[1, 6]
# }
# hist(pv1)
# hist(pv2)
# plot(-log10(pv1), -log10(pv2))
# 
# # Test 2 - validated
# pv1 <- pv2 <- numeric(100)
# for ( i in 1:100) {
#   x2 <- factor(sample(1:2, 20, repl=T)); x1 <- rnorm(20) + 0.25*as.numeric(x2); y <- 0.5*x1 + rnorm(20);
#   pv1[i] <- adonis(dist(y) ~ x1 + x2)$aov.tab[2, 6]
#   ID <- c(1:20, 1:20)
#   x1 <- c(x1, x1); x2 <- factor(c(x2, x2)); y <- c(y, y)
#   pv2[i] <- adonis2(dist(y) ~ x1 + x2, strata=ID)$aov.tab[1, 6]
# }
# hist(pv1)
# hist(pv2)
# plot(pv1, pv2)
# 
# pv1 <- pv2 <- numeric(100)
# for ( i in 1:100) {
#   x2 <- factor(sample(1:2, 20, repl=T)); x1 <- rnorm(20) + 0.25*as.numeric(x2); y <- 0.5*x1 + 1*as.numeric(x2) + rnorm(20);
#   pv1[i] <- adonis(dist(y) ~ x1 + x2)$aov.tab[2, 6]
#   ID <- c(1:20, 1:10)
#   x1 <- c(x1, x1[1:10]); x2 <- c(x2, x2[1:10]); z <- 0.2*rnorm(10);
#   y <- c(y[1:10] - z, y[11:20], y[1:10] + z)
#   ind <- sample(30)
#   y <- y[ind]; x1 <- x1[ind]; x2 <- x2[ind]; ID <- ID[ind]
#   pv2[i] <- adonis2(dist(y) ~ x1 + x2, strata=ID)$aov.tab[1, 6]
# }
# hist(pv1)
# hist(pv2)
# plot(-log10(pv1), -log10(pv2))
# abline(0, 1, col='red')
########################################
# This generates the matrix columns-wise
# From JnPaulson
generate_matrix <- function(x){
	indptr  = x$sample$matrix$indptr+1
	indices = x$sample$matrix$indices+1
	data    = x$sample$matrix$data
	nr = length(x$observation$ids)
	
	counts = sapply(2:length(indptr),function(i){
				x = rep(0,nr)
				seq = indptr[i-1]:(indptr[i]-1)
				x[indices[seq]] = data[seq]
				x
			})
	rownames(counts) = x$observation$ids
	colnames(counts) = x$sample$ids
	# I wish this next line wasn't necessary
	lapply(1:nrow(counts),function(i){
				counts[i,]
			})
}
generate_metadata <- function(x){
	metadata = x$metadata
	metadata = lapply(1:length(x$ids),function(i){
				id_metadata = lapply(metadata,function(j){
							if(length(dim(j))>1){ as.vector(j[,i,drop=FALSE]) }
							else{ j[i] }
						})
				list(id = x$ids[i],metadata=id_metadata)
			})
	return(metadata)
}
namedList <- function(...) {
	L <- list(...)
	snm <- sapply(substitute(list(...)),deparse)[-1]
	if (is.null(nm <- names(L))) nm <- snm
	if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
	setNames(L,nm)
}

read_hdf5_biom<-function(file_input){
	x = h5read(file_input,"/",read.attributes = TRUE)
	data = generate_matrix(x)
	rows = generate_metadata(x$observation)
	columns = generate_metadata(x$sample)
	shape = c(length(data),length(data[[1]])) # dim(data)
	
	# Experimental -- need to actually load these from file
	id = attr(x,"id")
	vs = attr(x,"format-version")
	format = sprintf("Biological Observation Matrix %s.%s",vs[1],vs[2])
	format_url = attr(x,"format-url")
	type = "OTU table"
	#type=attr(x,"type")
	generated_by = attr(x,"generated-by")
	date = attr(x,"creation-date")
	matrix_type = "dense"
	matrix_element_type = "int"
	
	namedList(id,format,format_url,type,generated_by,date,matrix_type,matrix_element_type,
			rows,columns,shape,data)
}
############################################
# comm are otu counts, row: otus, column samples
# intersect.no: Pairwise ratio calculated on pairs with at least 'intersect.no' common taxa
GMPR <- function (comm, intersect.no=4) {
	sapply(1:ncol(comm),  function(i) {
				x <- comm[, i]
				pr <- sapply(1:ncol(comm),  function(j) {
							y <- comm[, j]
							ind <- x != 0 & y != 0
							if (sum(ind) >= intersect.no) {
								res <- median(x[ind] / y[ind])
							} else {
								res <- NA
							}
						})
				if (sum(is.na(pr)) != 0) warning(paste0('Sample '), i, ' does not have any common taxa with some other samples!\n')
				exp(mean(log(pr[!is.na(pr)])))
			}
	)
}

uniquefy_taxa_names <- function (data.obj) {
	for (level in names(data.obj$abund.list)) {
		obj <- data.obj$abund.list[[level]]
		rownames(obj) <- gsub('unclassified', paste0('Unclassified', substr(level, 1, 1)), rownames(obj))
		data.obj$abund.list[[level]] <- obj
	}
	return(data.obj)
}


load_data <- function (otu.file, map.file, tree.file=NULL, parseFunction=parse_taxonomy_greengenes, version='Old', species=TRUE, filter.no=1, rep.seq=NULL,
	     norm='TSS', intersect.no=4, winsor=FALSE, winsor.qt=0.97,
		 ko.file=NULL, cog.file=NULL, ko.ann.file=NULL,
		 meta.sep='', quote="", comment="", read.gg=FALSE, rff=FALSE, dep=NULL, seed=1234) {
	# ko and cog file are not rarefied	
	# filter.no: filter the OTUs with read support less than filter.no (default is filtering singleton); singleton will not be filtered after rarefaction
	# winsorization and GMPR should be further studied. Current default is false and GMPR is on the genus level
	
	set.seed(seed)
	cat("Load meta file...\n")
	if (grepl("csv$", map.file)) {
		meta.dat <- read.csv(map.file, header=T, check.names=F, row.names=1, comment=comment, quote=quote)
	} else {
		meta.dat <- read.table(map.file, header=T, check.names=F, row.names=1, comment=comment, sep=meta.sep, quote=quote)
	}
	
	# Load Tree
	if (!is.null(tree.file)) {
		cat("Load tree file ...\n")
		if (read.gg == F) {
			tree.12 <- read.tree(tree.file)
		} else {
			tree.12 <- read_tree_greengenes(tree.file)
		}
		
		if (is.rooted(tree.12) == F) {
			tree.12 <- midpoint(tree.12)
		}
		
	} else {
		tree.12 <- NULL
	}

	cat("Load OTU file...\n")
	if (version != 'New') {
		biom.obj <-  import_biom(otu.file, parseFunction = parseFunction)  
	} else {
		temp <-  read_hdf5_biom(otu.file)
		otu.file <- paste0(otu.file, '.old')
		write_biom(temp, otu.file) 
		biom.obj <- import_biom(otu.file, parseFunction = parseFunction)  
	}

	otu.tab.12 <- otu_table(biom.obj)@.Data
	
	otu.ind <- rowSums(otu.tab.12 != 0) > filter.no
	otu.tab.12 <- otu.tab.12[otu.ind, ]
	# OTU names
	otu.name.full <- as.matrix(biom.obj@tax_table[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')])
	otu.name.full <- otu.name.full[otu.ind, ]
	
	otu.name.12 <- otu.name.full[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
	otu.name.12[is.na(otu.name.12)] <- 'unclassified'
	otu.name.12 <- otu.name.12@.Data
	
	otu.name.12[, ] <- gsub('\\[', '', otu.name.12)
	otu.name.12[, ] <- gsub('\\]', '', otu.name.12)
	
	if (rff == TRUE) {
		if (is.null(dep)) {
			otu.tab.12 <- t(Rarefy(t(otu.tab.12))$otu.tab.rff)
		} else {
			otu.tab.12 <- t(Rarefy(t(otu.tab.12), dep)$otu.tab.rff)
		}	
	} 

	samIDs <- intersect(rownames(meta.dat), colnames(otu.tab.12))
	if (length(samIDs) == 0)  stop('Sample names in the meta file and biom file differ?\n')
	
	meta.dat <- meta.dat[samIDs, ]
	otu.tab.12 <- otu.tab.12[, samIDs]
	

   # Create abundance list
    cat("Create taxa abundance list ...\n")
	abund.list.12 <- list()
	hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
	for (hierach in hierachs) {	
		if (hierach != 'Phylum') {
			single.names <- otu.name.12[, hierach]
		#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name.12[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name.12[, 'Phylum']
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
		}
		family <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list.12[[hierach]] <- family
	}
	
	if (species) {
		abund.list.12[['Species']] <- otu.tab.12
		rownames(abund.list.12[['Species']]) <- paste0("OTU", rownames(otu.tab.12), ":", otu.name.12[, 'Phylum'], ";", otu.name.12[, 'Genus'])
	}
	

	if (norm == 'GMPR') {
		sf <- GMPR(abund.list.12[['Genus']], intersect.no)
		names(sf) <- samIDs
		warning('GMPR is only suitable for samples from a same body location!\n')
	} else {
		if (norm == 'TSS') {
			sf <- colSums(otu.tab.12)
		} else {
			sf <- NULL
		}
	}
	
	if (winsor == TRUE) {
		# Addressing the outlier (97% percent) or at least one outlier
		abund.list.12 <- sapply(abund.list.12, function(genus) {
					genus.p <- t(t(genus) / sf)
					genus.p <- apply(genus.p, 1, function(x) {
								cutoff <- quantile(x, winsor.qt)
								x[x >= cutoff] <- cutoff
								x
							}
					)
					# column/row switch
					genus.w <- t(round(genus.p * sf))
				})
		# OTU table
		otu.tab.12.p <- t(t(otu.tab.12) / sf)
		otu.tab.12.p <- apply(otu.tab.12.p, 1, function(x) {
					cutoff <- quantile(x, winsor.qt)
					x[x >= cutoff] <- cutoff
					x
				}
		)
		# column/row switch
		otu.tab.12 <- t(round(otu.tab.12.p * sf))
	}

		
	# Rarefaction/Normalizing factors are not calculated for functional data
	if (!is.null(ko.file)) {
		cat("Load kegg file...\n")
		ko <- read_biom(ko.file)
		ko.dat <- as.matrix(biom_data(ko))
		ko.dat <- ko.dat[, intersect(colnames(ko.dat), samIDs)]
		# Rarefaction?
	  
	    if (is.null(ko.ann.file)) {
			# Old - back compatability
			ko.ann <- observation_metadata(ko)
			ko.ann <- cbind(KEGG_Pathways1=sapply(ko.ann, function(x) x['KEGG_Pathways1']), 
					KEGG_Pathways2=sapply(ko.ann, function(x) x['KEGG_Pathways2']), 
					KEGG_Pathways3=sapply(ko.ann, function(x) x['KEGG_Pathways3']))
			rownames(ko.ann) <- rownames(ko.dat)		
			ko.ann[is.na(ko.ann)] <- 'Unclassified'
			
			hierachs <- c("KEGG_Pathways1", "KEGG_Pathways2", "KEGG_Pathways3")
			for (hierach in hierachs) {	
				tax.family <- ko.ann[, hierach]
				family <- aggregate(ko.dat, by=list(tax.family), FUN=sum)
				rownames(family) <- family[, 1]
				family <- as.matrix(family[, -1])
				abund.list.12[[hierach]] <- family
			}
		} else {
			# New
			load(ko.ann.file)
			#
			kos <- rownames(ko.dat)
			abund.list.12[["KEGG_Pathways3"]] <- NULL
			kos.id <- NULL
			for (ko.item in names(kegg.map)) {
				kos.common <- intersect(kos, kegg.map[[ko.item]])
				if (length(kos.common) != 0) {
					abund.list.12[["KEGG_Pathways3"]] <- rbind(abund.list.12[["KEGG_Pathways3"]], colSums(ko.dat[kos.common, , drop=F]))
					kos.id <- c(kos.id, ko.item)
				}
			}
			rownames(abund.list.12[["KEGG_Pathways3"]]) <- kos.id
			
			abund.list.12[["KEGG_Metabolism"]] <- abund.list.12[["KEGG_Pathways3"]][intersect(kos.id, unlist(kegg.ann[['Metabolism']])), ]
			rownames(abund.list.12[["KEGG_Metabolism"]]) <- paste0('M', rownames(abund.list.12[["KEGG_Metabolism"]])) 
			
			abund.list.12[["KEGG_Defense"]] <- NULL
			kos.id <- NULL
			for (ko.item in names(defense.map)) {
				kos.common <- intersect(kos, defense.map[[ko.item]])
				if (length(kos.common) != 0) {
					abund.list.12[["KEGG_Defense"]] <- rbind(abund.list.12[["KEGG_Defense"]], colSums(ko.dat[kos.common, , drop=F]))
					kos.id <- c(kos.id, ko.item)
				}
			}
			rownames(abund.list.12[["KEGG_Defense"]]) <- kos.id
			
			abund.list.12[["KEGG_Toxin"]] <- NULL
			kos.id <- NULL
			for (ko.item in names(toxin.map)) {
				kos.common <- intersect(kos, toxin.map[[ko.item]])
				if (length(kos.common) != 0) {
					abund.list.12[["KEGG_Toxin"]] <- rbind(abund.list.12[["KEGG_Toxin"]], colSums(ko.dat[kos.common, , drop=F]))
					kos.id <- c(kos.id, ko.item)
				}
			}
			rownames(abund.list.12[["KEGG_Toxin"]]) <- kos.id
		}
	
	}
	
	if (!is.null(cog.file)) {
		cat("Load cog file...\n")
		cog <- read_biom(cog.file)
		cog.dat <- as.matrix(biom_data(cog))
		cog.dat <- cog.dat[, intersect(colnames(ko.dat), samIDs)]
		# rarefaction?
		cog.ann <- observation_metadata(cog)
		hierachs <- c("COG_Category1", "COG_Category2")
		for (hierach in hierachs) {	
			tax.family <- sapply(cog.ann, function(x) x[hierach])
			family <- aggregate(cog.dat, by=list(tax.family), FUN=sum)
			rownames(family) <- family[, 1]
			family <- as.matrix(family[, -1])
			abund.list.12[[hierach]] <- family
		}
	}
		
	data.obj <- list(otu.tab=otu.tab.12, otu.name=otu.name.12, abund.list=abund.list.12, meta.dat=meta.dat, tree=tree.12,
			otu.name.full=otu.name.full, size.factor=sf, call=match.call())
}

# THe GMPR was performed on the genus level, it may not be a good idea
winsor_data <- function (data.obj, norm='GMPR', intersect.no=4, winsor.qt=0.97) {
	
	otu.tab.12 <- data.obj$otu.tab
	abund.list.12 <- data.obj$abund.list
	
	if (norm == 'GMPR') {
		warning('GMPR is only suitable for samples from the same body location!\n')
		sf <- GMPR(abund.list.12[['Genus']], intersect.no)
		names(sf) <- colnames(otu.tab.12)
	} 
	if (norm == 'TSS') {
		sf <- colSums(otu.tab.12)
	}
	
	# Addressing the outlier (97% percent) or at least one outlier
	abund.list.12 <- sapply(abund.list.12, function(genus) {
				genus.p <- t(t(genus) / sf)
				genus.p <- apply(genus.p, 1, function(x) {
							cutoff <- quantile(x, winsor.qt)
							x[x >= cutoff] <- cutoff
							x
						}
				)
				# column/row switch
				genus.w <- t(round(genus.p * sf))
			})
	# OTU table
	otu.tab.12.p <- t(t(otu.tab.12) / sf)
	otu.tab.12.p <- apply(otu.tab.12.p, 1, function(x) {
				cutoff <- quantile(x, winsor.qt)
				x[x >= cutoff] <- cutoff
				x
			}
	)
	# column/row switch
	otu.tab.12 <- t(round(otu.tab.12.p * sf))
    data.obj$otu.tab <- otu.tab.12
	data.obj$abund.list <- abund.list.12
	data.obj$size.factor <- sf
	data.obj
}

construct_distance <- function (data.obj, unifrac.file=NULL,  Phylum='All', dist.RData=NULL, save.RData=NULL, filter.no=0, rff=FALSE, dep=NULL, seed=1234) {
	set.seed(seed)
	
	if (!is.null(dist.RData)) {
		load(dist.RData, envir=.GlobalEnv )
	} else {
		dist.list.12 <- list()
		cat("Generalized UniFrac ...\n")
		
		otu.tab <- t(data.obj$otu.tab)
		
		if (rff == TRUE) {
			if (is.null(dep)) {
				otu.tab <- Rarefy(otu.tab)$otu.tab.rff
			} else {
				otu.tab <- Rarefy(otu.tab, dep)$otu.tab.rff
			}	
		}
				
		if (Phylum != 'All') {
			ind <- data.obj$otu.name[, 'Phylum'] == Phylum
			otu.tab <- otu.tab[, ind]
		}
		
		# Filter otus with reads <= filter.no
		otu.tab <- otu.tab[, colSums(otu.tab) > filter.no]
        
		# Remove samples with no reads
		if (sum(rowSums(otu.tab) == 0) >= 1) {
			otu.tab <- otu.tab[rowSums(otu.tab) != 0, ]
			warnings('Some samples do not have reads after rarefaction! Please be careful!\n')
		}
		
		unifrac12 <- GUniFrac(otu.tab, data.obj$tree)$unifracs
	
		dist.list.12[['WUniFrac']] <- unifrac12[, , 'd_1']
		dist.list.12[['GUniFrac']] <- unifrac12[, , 'd_0.5']
		if (is.null(unifrac.file)) {
			dist.list.12[['UniFrac']] <- unifrac12[, , 'd_UW']
		} else {
			# The orders may be different
			dist.list.12[['UniFrac']] <- as.matrix(read.table(unifrac.file, row.names=1, header=T)) # Rarefaction
		}

		# Need speed up
		# Suggest using rarefied counts 
	    # If case/control has different sequencing depth, it will result in false clustering! 
		dist.list.12[['BC']] <-as.matrix(vegdist(otu.tab))
		
		genus <- t(data.obj$abund.list[['Genus']])
		genus <- genus / rowSums(genus)
		dist.list.12[['Euc']] <- as.matrix(dist(genus))
		
		genus <- sqrt(genus)
		dist.list.12[['Hel']] <-as.matrix(dist(genus))
#		dist.list.12[['JS']] <- as.matrix(distance(otu_table(data.obj$abund.list[['Genus']], taxa_are_rows=T), method='jsd'))
		
		if (!is.null(save.RData)) {
			save(dist.list.12, file=save.RData)
		}
	}

    return(dist.list.12)
}


outlier_detect <- function (data.obj, dist.obj, min.dep=2000) {
	# Future development
	samIDs <- colnames(data.obj$otu.tab)[colSums(data.obj$otu.tab) >= 2000]
	return(samIDs)
}

subset_data <- function (data.obj, samIDs) {
	data.obj$meta.dat <- data.obj$meta.dat[samIDs, ]
	data.obj$otu.tab <- data.obj$otu.tab[, samIDs]
	data.obj$abund.list <- lapply(data.obj$abund.list, function(x) x[, samIDs])
	if (length(data.obj$size.factor)) data.obj$size.factor <- data.obj$size.factor[samIDs]
	if (length(data.obj$ko.list) != 0) {
		data.obj$ko.list <- lapply(data.obj$ko.list, function(x) x[, samIDs])
	}
	if (length(data.obj$cog.list) != 0) {
		data.obj$cog.list <- lapply(data.obj$cog.list, function(x) x[, samIDs])
	}
	data.obj
}

subset_dist <- function (dist.obj, samIDs) {
	lapply(dist.obj, function(x) x[samIDs, samIDs])
}


perform_sequence_stat_analysis <- function (data.obj, ann='') {
	sink(paste0('Sequence_Analysis_Statistics_', ann, '.txt'))
	otu.tab <- data.obj$otu.tab
	
	# Sequencing depth
	otu.abund <- rowSums(otu.tab)
	sam.abund <- colSums(otu.tab)
	otu.prev <- rowSums(otu.tab!=0)/ncol(otu.tab)
	
	otu.abund <- otu.abund[otu.abund >= 1]
	sam.abund <- sam.abund[sam.abund >= 1]
	cat('16S rDNA targeted sequencing yields ', mean(sam.abund), 'reads/sample on average (range:', min(sam.abund), '-', max(sam.abund), ').')
	cat('Clustering of these 16S sequence tags produces ', sum(otu.abund > 0), ' OTUs at 97% similarity level.')
	
	pdf(paste0('Sequence_Analysis_Statistics_', ann, '.pdf'), height=5, width=5)
	obj <- ggplot2::ggplot(data=data.frame(x=otu.abund), aes(x=x)) + geom_histogram(col='black', fill='gray') + ylab('Frequency') + xlab('Abundance(Total counts)') +
			scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000, 100000))
	print(obj)
	obj <- ggplot2::ggplot(data=data.frame(x=sam.abund), aes(x=x)) + geom_histogram(col='black', fill='gray')  + ylab('Frequency') + xlab('Sequencing depth')
	print(obj)
	obj <- ggplot2::ggplot(data=data.frame(x=otu.prev), aes(x=x))  + ylab('Frequency') + xlab('Prevalence(Occurence frequency)') + geom_histogram(col='black', fill='gray')
	print(obj)
	dev.off()
	
	phy.abund <- data.obj$abund.list[['Phylum']]
	fam.abund <- data.obj$abund.list[['Family']]
	gen.abund <- data.obj$abund.list[['Genus']]
	
	phy.prev <- rowSums(phy.abund != 0) / ncol(phy.abund)
	fam.prev <- rowSums(fam.abund != 0) / ncol(phy.abund)
	gen.prev <- rowSums(gen.abund != 0) / ncol(phy.abund)
	
	phy.abund <- rowMeans(t(t(phy.abund) / sam.abund))
	fam.abund <- rowMeans(t(t(fam.abund) / sam.abund))
	gen.abund <- rowMeans(t(t(gen.abund) / sam.abund))
	
	cat('These OTUs belong to ', sum(phy.abund > 0), ' phyla,', sum(fam.abund > 0), ' families and ', sum(gen.abund > 0), 'genera.')
	
	phy.prev <- sort(phy.prev, decr=T)
	phy.prev <- round(phy.prev[phy.prev >= 0.90] * 100, 2)
	
	fam.prev <- sort(fam.prev, decr=T)
	fam.prev <- round(fam.prev[fam.prev >= 0.90] * 100, 2)
	
	gen.prev <- sort(gen.prev, decr=T)
	gen.prev <- round(gen.prev[gen.prev >= 0.90] * 100, 2)
	
	cat('\nThe most prevalent phyla are ', paste(paste0(names(phy.prev), '(', phy.prev, '%)'), collapse=' '), ';')
	cat('\nThe most prevalent families are ', paste(paste0(names(fam.prev), '(', fam.prev, '%)'), collapse=' '), ';')
	cat('\nand the most prevalent genera are ', paste(paste0(names(gen.prev), '(', gen.prev, '%)'), collapse=' '), '.')
	
	phy.abund <- sort(phy.abund, decr=T)
	phy.abund <- round(phy.abund[phy.abund >= 0.05] * 100, 2)
	
	fam.abund <- sort(fam.abund, decr=T)
	fam.abund <- round(fam.abund[fam.abund >= 0.05] * 100, 2)
	
	gen.abund <- sort(gen.abund, decr=T)
	gen.abund <- round(gen.abund[gen.abund >= 0.05] * 100, 2)
	
	cat('\nThe most abundant phyla are ', paste(paste0(names(phy.abund), '(', phy.abund, '%)'), collapse=' '), ';')
	cat('\nThe most abundant families are ', paste(paste0(names(fam.abund), '(', fam.abund, '%)'), collapse=' '), ';')
	cat('\nand the most abundant genera are ', paste(paste0(names(gen.abund), '(', gen.abund, '%)'), collapse=' '), '.')
	sink()
}

perform_demograph_analysis <- function (data.obj, grp.name) {
	obj <- summary(data.obj$meta.dat)
	write.csv(obj, "meta.data.summary.csv", quote=F)

	grp <- data.obj$meta.dat[, grp.name]
	res <- NULL
	pv.vec <- NULL
	if (is.factor(grp)) {
		for (var1 in setdiff(colnames(data.obj$meta.dat), grp.name)) {
			temp <- data.obj$meta.dat[, var1]
			res <- rbind(res, c("", var1, rep("", nlevels(grp) - 1)))
			res <- rbind(res, c("", levels(grp)))
			if (is.factor(temp)) {
				res <- rbind(res, cbind(levels(temp), table(temp, grp)))
				if (nlevels(temp) == 1) {
					pv <- NA
				} else {
					err <- try(
							pv <- formatC(fisher.test(table(temp, grp))$p.value, digit=3)
					)
					if (inherits(err, "try-error")) {
						pv <- NA
					}
				}
				res <- rbind(res, c('Fisher p', pv, rep("", nlevels(grp) - 1)))

			} else {
				res <- rbind(res, c('mean', aggregate(temp, by=list(grp), FUN='mean')[, 2]))
				res <- rbind(res, c('sd', aggregate(temp, by=list(grp), FUN='sd')[, 2]))
				
				err <- try(
						pv <- formatC(summary(aov(temp ~ grp))[[1]][1, 'Pr(>F)'], digit=3)
				)
				if (!inherits(err, "try-error")) {
					res <- rbind(res, c('ANOVA p', pv, rep("", nlevels(grp) - 1)))
				} else {
					res <- rbind(res, c('ANOVA p', 'NA', rep("", nlevels(grp) - 1)))
					pv <- NA
				}

			}
			res <- rbind(res, rep("", nlevels(grp)+1))
			pv.vec <- c(pv.vec, pv)
		}
		
	} else {

		
	}

	write.csv(res, "meta.data.by.grp.csv", row.names=F, quote=F)
	names(pv.vec) <- setdiff(colnames(data.obj$meta.dat), grp.name)
	return(pv.vec)
}


generate_rarefy_curve <- function (data.obj, phylo.obj, grp.name, depth=NULL, npoint=10, iter.no=5,
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), ann='') {
	cat("Create rarefaction curves!\n")
	if (is.null(depth)) {
		depth <- min(sample_sums(phylo.obj))
		phylo.even <- rarefy_even_depth(phylo.obj, rngseed=12345)
	} else {
		if (depth > min(sample_sums(phylo.obj))) {
			ind <- sample_sums(phylo.obj) >= depth
			cat(sum(!ind), " samples do not have sufficient number of reads!\n")
			sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
			data.obj <- subset_data(data.obj, ind)
		}
		phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
	}
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]

	res <- NULL
	incr <- depth %/% npoint
	sink('temp.txt')
	for (dep in c(10, incr*(1:npoint))) {
		x <- 0
		for (i in 1:iter.no) {
			phylo.even <- rarefy_even_depth(phylo.obj, dep, rngseed=12345+i)
			x <- x + estimate_richness(phylo.even, measures=measures)
		}
		
		res <- rbind(res, t(x[, measures, drop=F]/iter.no))
	}
	colnames(res) <- rownames(df)
	sink()
	
	pdf(paste0("Alpha_diversity_Rarefaction_Curve_", ann, ".pdf"), width=5, height=5)
	for (i in 1:length(measures)) {
		measure <- measures[i]
		cat("Measure: ", measure, "\n")
		res2 <- res[(0:(npoint))*length(measures)+i, , drop=F]
		m <- t(apply(res2, 1, function(x) tapply(x, grp, mean)))
		se <- t(apply(res2, 1, function(x) tapply(x, grp, function(y) sd(y)/sqrt(length(y)))))
		uci <- m+se
		lci <- m-se
		
		m <- melt(m)
		uci <- melt(uci)
		lci <- melt(lci)
		
		res2 <- cbind(c(10, incr*(1:npoint)), m[, 2:3], uci[, 3], lci[, 3])
		colnames(res2) <- c('Depth', 'Group', 'mean', 'max', 'min')
		
		res2 <- as.data.frame(res2)
		res2$Group <- factor(res2$Group, levels=levels(grp))
		
		#write.table(res2, paste0("Alpha_diversity_Rarefaction_", ann, "_", measure, ".txt"))
		
		obj <- ggplot(res2, aes(x=Depth, y=mean, color=Group, group=Group)) +
				geom_errorbar(aes(ymin=min, ymax=max), alpha=0.5, width=.25, position=position_dodge(.2)) + 
				geom_line() + 
				geom_point(size=3, shape=21, fill="white") +
				labs(y=measure) +
				theme(legend.justification=c(1,0), legend.position=c(1,0))
		print(obj)		
	}
	dev.off()
}

perform_alpha_test <- function (data.obj, phylo.obj, rarefy=TRUE, depth=NULL, iter.no=10, 
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'),  model='lm', 
		formula=NULL, grp.name=NULL, adj.name=NULL, ann='', ...) {
	# Implement future random effects model
	
	if (rarefy == TRUE) {
		if (is.null(depth)) {
			depth <- min(sample_sums(phylo.obj))
		} else {
			if (depth > min(sample_sums(phylo.obj))) {
				ind <- sample_sums(phylo.obj) >= depth
				cat(sum(!ind), " samples do not have sufficient number of reads!\n")
				sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
				data.obj <- subset_data(data.obj, ind)
			}
		}
			
		x <- 0
		sink('temp.txt')
		for (i in 1:iter.no) {
			phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345+i)
			x <- x + estimate_richness(phylo.even, measures=measures)
		}
		sink()
		x <- x / iter.no
	} else {
		x <- estimate_richness(phylo.obj, measures=measures)
	}

	df <- data.obj$meta.dat
	
	if (rarefy == T) {
		sink(paste0('Alpha_diversity_test_results_rarefied_', ann, '.txt'))
	} else {
		sink(paste0('Alpha_diversity_test_results_unrarefied_', ann, '.txt'))
	}
	
	date()
	if (is.null(formula)) {
		if (is.null(adj.name)) {
			formula <- paste('~', grp.name)
		} else {
			formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
		}
	}
	
	for (measure in measures) {
		cat("Alpha diversity:", measure, "\n")
		xx <- x[, measure]
		if (model == 'lm') {
			cat('Linear model:\n')
			prmatrix(summary(lm(as.formula(paste('xx ', formula)), df, ...))$coefficients)
		}
		if (model == 'lme') {
			df$xx <- xx
			cat('Linear mixed effects model:\n')
			prmatrix(summary(lme(as.formula(paste('xx ', formula)), df, method='ML', ...))$tTable)
		}
		cat("\n")
	}
	sink()

}

generate_alpha_boxplot <- function (data.obj, phylo.obj, rarefy=TRUE, depth=NULL, grp.name, strata=NULL, 
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson')) {	
	# To be completed - jetter when strata is not null
	if (rarefy == TRUE) {
		if (is.null(depth)) {
			depth <- min(sample_sums(phylo.obj))
		} else {
			if (depth > min(sample_sums(phylo.obj))) {
				ind <- (sample_sums(phylo.obj) >= depth)
				cat(sum(!ind), " samples do not have sufficient number of reads!\n")
				
				sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
				data.obj <- subset_data(data.obj, ind)
			}
		}

		phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
		x <- estimate_richness(phylo.even, measures=measures)
	} else {
		x <- estimate_richness(phylo.obj, measures=measures)
	}
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	hei <- 5
	if (is.null(strata)) {
		wid <- 5
	} else {
		wid <- 5.5
	}
	if (rarefy == T) {
		pdf(paste0('Alpha_diversity_boxplot_rarefied.pdf'), height=hei, width=wid)
	} else {
		pdf(paste0('Alpha_diversity_boxplot_unrarefied.pdf'), height=hei, width=wid)
	}
	
	if (is.null(strata)) {
		for (measure in measures) {
			cat(measure, '\n')
			xx <- x[, measure]		
			df2 <- data.frame(Value=xx, Group=grp)
			dodge <- position_dodge(width=0.75)
			obj <- ggplot(df2, aes(x=Group, y=Value, col=Group)) +
					geom_boxplot(position=dodge,  outlier.colour = NA) + 
					geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
					labs(y=measure) +
					theme(legend.position="none")
			print(obj)
		}	
	} else {
		for (measure in measures) {
			cat(measure, '\n')
			xx <- x[, measure]		
			grp2 <- df[, strata]
			df2 <- data.frame(Value=xx, Group=grp, Strata=grp2)
			
			dodge <- position_dodge(width=0.95)
			obj <- ggplot(df2, aes(x=Strata, y=Value, col=Group)) +
					geom_boxplot(position=dodge,  outlier.colour = NA) + 
					geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
					labs(y=measure, x=strata) 
			print(obj)
		}
		
	}
	dev.off()
}

perform_alpha_test_otu <- function (data.obj, formula, grp.name, ...) {
	# Assume two groups
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	dep <- colSums(data.obj$otu.tab)
	obs <- colSums(data.obj$otu.tab !=0)
	obj <- lm(as.formula(paste('log(obs)', formula, '+ log(dep)')), df, ...)
	sink('Alpha_diversity_test_results_observed_OTU.txt')
	date()
	cat('Linear model correcting for differential sequencing depth:\n')
	prmatrix(summary(obj)$coefficients)
	sink()
	
	pdf("Alpha_diversity_test_results_observed_OTU.pdf")
	plot(dep, obs, log='xy', xlab="Sequencing Depth", ylab="Observed OTU Number", 
			bg=c('blue', 'red')[grp], pch=21)
	grp.levels <- levels(grp)
	abline(lm(obs[grp==grp.levels[1]] ~ dep[grp==grp.levels[1]]), col='red')
	abline(lm(obs[grp==grp.levels[2]] ~ dep[grp==grp.levels[2]]), col='blue')
	legend('left', levels(grp), fill=c('blue', 'red'))
	dev.off()
}


generate_ordination <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name, adj.name=NULL, emp.lev=NULL, strata=NULL, pc.separate, pca.method='cmd', ann=NULL, sub=NULL,
		clab=1.0, cex.pt=1.25, ellipse=T, cstar= 1, wid=5, hei=5, pdf=TRUE, ...) {
	# Implment strata
	# To be completed, add continuous case
	strata0 <- strata
	
	if (pdf) {
		if (is.null(ann)) {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', grp.name, '.pdf'), width=wid, height=hei)
		} else {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', ann, '.pdf'), width=wid, height=hei)
		}
	}

	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	
	if (is.null(emp.lev)) {
		grp <- factor(grp, levels=c(setdiff(levels(grp), emp.lev), emp.lev))
	}
	
	if (!is.null(strata)) {
		strata <- factor(df[, strata])
	} else {
		strata <- factor(grp)
	}
	
	darkcols <- hue_pal(l=40)(nlevels(grp))
	lightcols <- hue_pal(c=45, l=80)(nlevels(grp))
	pchs <- rep(c(21, 22, 23, 24, 25), ceiling(nlevels(strata) / 5))[1:nlevels(strata)]
	
	for (dist.name in dist.names) {
		dist.temp <- dist.obj[[dist.name]]
		if (!is.null(adj.name)) {
			adj <- as.data.frame(df[, adj.name])
			obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
			dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
			dist.temp <- dist(dat2)
		} 
		if (pca.method == 'cmd') {
			obj <- cmdscale(as.dist(dist.temp), k=2, eig=T)
			pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
			y <- cbind(obj$points[, 1], obj$points[, 2])
			
			xlab <- paste0('PC1(', pve[1], '%)')
			ylab <- paste0('PC2(', pve[2], '%)')
		} 

		if (pca.method == 'nmds') {
				obj <- metaMDS(as.dist(dist.temp), k=2)
				y <- cbind(obj$points[, 1], obj$points[, 2])
				xlab <- 'NMDS1'
                ylab <- 'NMDS2'
		} 
		
		if (pca.method == 'pls') {
			require(mixOmics)
			# Test
			obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
			obj <- plsda(obj, grp, ncomp=2)
			y <- cbind(obj$variates$X[, 1], obj$variates$X[, 2])
			xlab <- 'PLS1'
			ylab <- 'PLS2'
		}
		
		colnames(y) <- c("PC1", "PC2")		
		plot(y[, 1], y[, 2], type='n', xlim=range(y[, 1]) * 1.2, ylim=range(y[, 2])*1.2, 
				xlab=xlab, ylab=ylab)
#		points(y[, 1], y[, 2], bg='yellow', col=darkcols[grp], 
#				type='p', pch = pchs[grp], cex=cex.pt)
		col1 <- lightcols[1:nlevels(grp)]
		col2 <- darkcols[1:nlevels(grp)]
		col3 <- darkcols[grp]
		if (!is.null(emp.lev)) {
			col1 <- col2
			col1[which(levels(grp) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
			col1[which(levels(grp) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
			col2[which(levels(grp) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
			col2[which(levels(grp) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
            col3[grp !=  emp.lev] <- rgb(0.5, 0.5, 0.5, 0.5)
			col3[grp ==  emp.lev] <- rgb(1.0, 0, 0, 1.0)

		}
		if (ellipse == T) {
			s.class(y, 
					fac = grp,
					cstar = cstar,
					clab = 0,
					cpoint = 0,
					axesell = F,
					col = col1,
					grid = TRUE,
					add.plot=T,
					...
			)
		}
		points(y[, 1], y[, 2], bg=col3, col='black',
				type='p', pch = pchs[strata], cex=cex.pt)
		s.class(y, 
				fac = grp,
				cstar =0,
				cellipse = 0,
				clab = clab,
				cpoint = 0,
				axesell = F,
				col = col2,
				grid = TRUE,
				add.plot=T
		)
		if (!is.null(strata0)) {
			legend('topleft', legend=(levels(strata)), pch=pchs[1:nlevels(strata)])	
		}

		title(main=paste(dist.name, "distance"), sub=sub)
#		text(-0.25, 0.3, "PERMANOVA p=0.016")
	}
	
	if (pdf) {
		dev.off()
	}
}

generate_distance_barplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name, strata=NULL, within=T, between=T) {
	strata.name <- strata
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	grp.btws <- outer(grp.levels, grp.levels, paste, sep="_")
	grp.btws <- grp.btws[lower.tri(grp.btws)]
	
	if (!is.null(strata.name)) {
		strata <- df[, strata.name]
	} else {
		strata <- factor(rep(1, nrow(df))) # pseudo strata
	}
	res.df <- NULL
	for (dist.name in dist.names) {
		for (stratum in levels(strata)) {
			ind <- strata %in% stratum
			dist.sub <- dist.obj[[dist.name]][ind, ind]
			df2 <- df[ind, ]
			if (between) {
				for (grp.btw in grp.btws) {
					ind1 <- df2[, grp.name] == unlist(strsplit(grp.btw, "_"))[1]
					ind2 <- df2[, grp.name] == unlist(strsplit(grp.btw, "_"))[2]	
					temp <- as.vector(dist.sub[ind1, ind2])
					res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Between', 
									DistanceType=grp.btw, Distance=mean(temp), sd=sd(temp)))
				}	
			}

			if (within) {
				for (grp.wth in grp.levels) {
					ind1 <- df2[, grp.name] == grp.wth
					temp <- dist.sub[ind1, ind1]			
					temp <- temp[lower.tri(temp)]
					res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Within',
									DistanceType=grp.wth, Distance=mean(temp), sd=sd(temp)))
				}	
			}
		}
		
	}
	if (between & within) {
		res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
		levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "_"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
	} else {
		if (between) {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "_"), paste, collapse=' vs\n'))
		} else {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
		}
	}
	
	if (is.null(strata.name)) {
		pdf(paste("Beta_diversity", "btw", between, "wth", within, "no_strata_barplot.pdf", sep="_"),
				width=5, height=5)
		
		limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
		dodge <- position_dodge(width=0.9)
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			obj1 <- ggplot(temp, aes(x=DistanceType, y=Distance, fill=DistanceType)) + 
					geom_bar(position=dodge, stat="identity", width=0.75) + 
					geom_bar(position=dodge, stat="identity", width=0.75, colour="black", show_guide=FALSE, size=0.25) +
					geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
					labs(y=paste(dist.name, "Distance"), x='') +
					theme(legend.position="none")
					
			print(obj1)
			
		}
		dev.off()
	} else {
		pdf(paste("Beta_diversity", "btw", between, "wth", within, "strata", strata.name, "barplot.pdf", sep="_"),
				width=2.5*(nlevels(strata)-1) + 5, height=5)
		
		limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
		dodge <- position_dodge(width=0.9)
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			obj1 <- ggplot(temp, aes(x=Strata, y=Distance, fill=DistanceType)) + 
					geom_bar(position=dodge, stat="identity") + 
					geom_bar(position=dodge, stat="identity", colour="black", show_guide=FALSE, size=0.25) +
					geom_errorbar(limits, position=dodge, size=0.25, width=0.5) +
					labs(y=paste(dist.name, "Distance"), x=strata.name)
			print(obj1)
			
		}
		dev.off()
		
	}
	
}

generate_distance_boxplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name, strata=NULL, within=F, between=T) {
	strata.name <- strata
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	grp.btws <- outer(grp.levels, grp.levels, paste, sep="_")
	grp.btws <- grp.btws[lower.tri(grp.btws)]
	
	if (!is.null(strata.name)) {
		strata <- df[, strata.name]
	} else {
		strata <- factor(rep(1, nrow(df))) # pseudo strata
	}
	res.df <- NULL
	for (dist.name in dist.names) {
		for (stratum in levels(strata)) {
			ind <- strata %in% stratum
			dist.sub <- dist.obj[[dist.name]][ind, ind]
			df2 <- df[ind, ]
			if (between) {
				for (grp.btw in grp.btws) {
					ind1 <- df2[, grp.name] == unlist(strsplit(grp.btw, "_"))[1]
					ind2 <- df2[, grp.name] == unlist(strsplit(grp.btw, "_"))[2]	
					temp <- as.vector(dist.sub[ind1, ind2])
					n <- length(temp)
					res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Between', n),
									DistanceType=rep(grp.btw, n), Distance=temp))
				}	
			}
			
			if (within) {
				for (grp.wth in grp.levels) {
					ind1 <- df2[, grp.name] == grp.wth
					temp <- dist.sub[ind1, ind1]			
					temp <- temp[lower.tri(temp)]
					n <- length(temp)
					res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Within', n),
									DistanceType=rep(grp.wth, n), Distance=temp))
				}	
			}
		}
		
	}
	if (between & within) {
		res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
		levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "_"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
	} else {
		if (between) {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "_"), paste, collapse=' vs\n'))
		} else {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
		}
	}
	
	if (is.null(strata.name)) {
		pdf(paste("Beta_diversity", "btw", between, "wth", within, "no_strata_boxplot.pdf", sep="_"),
				width=5, height=5)
		
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			
			dodge <- position_dodge(width=0.95)		
			obj1 <- ggplot(temp, aes(x=DistanceType, y=Distance, col=DistanceType)) + 
					geom_boxplot(position=dodge, outlier.colour = NA) + 
#					geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
					labs(y=paste(dist.name, "Distance"), x='') +
					theme(legend.position="none")
			
			print(obj1)
			
		}
		dev.off()
	} else {
		pdf(paste("Beta_diversity", "btw", between, "wth", within, "strata", strata.name, "boxplot.pdf", sep="_"),
				width=2.5*(nlevels(strata)-1) + 5, height=5)
		
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			dodge <- position_dodge(width=0.95)		
			obj1 <- ggplot(temp, aes(x=Strata, y=Distance, col=DistanceType)) + 
					geom_boxplot(position=dodge, outlier.colour = NA) + 
		#			geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
					labs(y=paste(dist.name, "Distance"), x=strata.name)
			print(obj1)
			
		}
		dev.off()
		
	}
	
}


generate_clustering <- function(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), meta.info, cluster.method='average', 
		is.labRow=F, cex.lab=NULL, ann="", wid=10, hei=6) {
	
	pdf(paste0("Beta_diversity_Hierachical_clustering_", ann, ".pdf"), width=wid, height=hei)
	for (dist.name in dist.names) {

		df <- data.obj$meta.dat
		dist.sub <- dist.obj[[dist.name]]
		dend <- hclust(as.dist(dist.sub), cluster.method)
		
		rainbow3 <- function(x) {
			rainbow(x + 2)[1:x]
		}
		jet3 <- function(x) {
			jet(x + 2)[1:x]
		}
		colFnsC <- c(colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue')), colorRampPalette(c('black', 'red')))
		colFnsF <- c(rainbow3, jet3)
		key.list <- list()
		mat <- NULL
		for (keyID in meta.info) {
			x <- df[, keyID]
			i <- 0
			j <- 0
			if (is.factor(x)) {
				key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)	
				i <- (i + 1) %% length(colFnsF)
				mat <- cbind(mat, key.list[[keyID]]$colors[x])
			} else {
				key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
				j <- (j + 1) %% length(colFnsC)
				mat <- cbind(mat, cmap(x, key.list[[keyID]]))
			}
		}
		colnames(mat) <- meta.info

		par(oma = c(1, 2, 1, 2))
		par(mar = c(5,4,4,3)+0.1)  # make space for color keys
		if (is.labRow == F) {
			labRow <- ""
		} else {
			labRow <- rownames(df)
		}
		if (is.null(cex.lab)) {
			cex.lab <- 25 / nrow(df) 
		}
		
		dendromat(dend, mat, labRow=labRow,
				ylab = 'Distance', main = paste(dist.name, "distance"), cex.lab=cex.lab)
		
		par(oma=c(0, 0, 0, 0))
		y.cord <- (1/length(meta.info)) * (0:(length(meta.info) - 1))
		k <- 1
		for (keyID in meta.info) {
			x <- df[, keyID]
			if (is.factor(x)) {
				vkey2(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
			} else {
				vkey(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
			}
			k <- k + 1
		}		
	}
	dev.off()
}

PermanovaG2 <- function (formula, dat = NULL, ...) 
{
	save.seed <- get(".Random.seed", .GlobalEnv)
	lhs <- formula[[2]]
	lhs <- eval(lhs, dat, parent.frame())
	rhs <- as.character(formula)[3]
	f.perms <- -Inf
	f.stat <- -Inf
	for (i in 1:(dim(lhs)[3])) {
		assign(".Random.seed", save.seed, .GlobalEnv)
		Y <- as.dist(lhs[, , i])
		formula2 <- as.formula(paste("Y", "~", rhs))
		obj <- adonis2(formula2, dat, ...)
		f.perms <- ifelse(f.perms > obj$f.perms, f.perms, obj$f.perms)
		temp <- obj$aov.tab[1:ncol(f.perms), "F.Model"]
		f.stat <- ifelse(f.stat > temp, f.stat, temp)
	}
	pv <- (rowSums(t(f.perms) >= f.stat) + 1)/(nrow(f.perms) + 
				1)
	aov.tab <- data.frame(F.Model = f.stat, p.value = pv)
	rownames(aov.tab) <- rownames(obj$aov.tab)[1:ncol(f.perms)]
	list(aov.tab = aov.tab)
}

perform_permanova_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		PermanovaG.dist=c('UniFrac',  'WUniFrac', 'BC'),
		formula=NULL,  grp.name=NULL, adj.name=NULL, pairwise=F, block.perm=F, strata=NULL, ann='', ...) {
	# PermanovaG not implemented for block permutation
	df <- data.obj$meta.dat
	if (!is.null(strata)) {
		if (is.character(strata)) {
			strata <- df[, strata]
		}
	}
	if (!is.null(formula)) {
		ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
		df <- df[ind, ]
		if (!is.null(strata)) {
			strata <- strata[ind]
		}
		sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '.txt'))
		date()
		cat('\nPERMANOVA test: \n')
		for (dist.name in dist.names) {
			cat(dist.name, " distance: \n")
			dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
			if (block.perm == F) {
				obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
			} else {
				obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
			}
			
			prmatrix(obj$aov.tab)
			cat("\n")
		}
		if (block.perm == F) {
			cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
			response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
			for (dist.name in PermanovaG.dist) {
				response[, , dist.name] <- as.dist(dist.obj[[dist.name]][ind, ind])
			}
			obj <- PermanovaG(as.formula(paste("response", formula)), df,  strata=strata, ...)
			prmatrix(obj$aov.tab)
			cat("\n")
		}
		cat("\n")
		sink()
	} else {
		if (pairwise == F) {
			if (is.null(adj.name)) {
				formula <- paste('~', grp.name)
			} else {
				formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
			}
			
			ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
			df <- df[ind, ]
			if (!is.null(strata)) {
				strata <- strata[ind]
			}
			sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '.txt'))
			date()
			cat('\nPERMANOVA test: \n')
			for (dist.name in dist.names) {
				cat(dist.name, " distance: \n")
				dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
				if (block.perm == F) {
					obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
				} else {
					obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
				}
				prmatrix(obj$aov.tab)
				cat("\n")
			}
			
			if (block.perm == F) {
				cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
				response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
				for (dist.name in PermanovaG.dist) {
					response[, , dist.name] <- dist.obj[[dist.name]][ind, ind]
				}
				obj <- PermanovaG(as.formula(paste("response", formula)), df,  strata=strata, ...)
				prmatrix(obj$aov.tab)
				cat("\n")
			}
			
			cat("\n")
			sink()
			
		} else {
			sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwise.txt'))
			date()
			cat('\nPairwise PERMANOVA test: \n')
			grp <- factor(df[, grp.name])
			grp.levels <- levels(grp)
			grp.nlevels <- nlevels(grp)
			pmat.all <- NULL
			rmat.all <- NULL
			pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
			colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
			for (dist.name in dist.names) {
				cat(dist.name, " distance: \n")
				pmat <- matrix(NA, grp.nlevels, grp.nlevels)
				colnames(pmat) <- rownames(pmat) <- grp.levels
				rmat <- matrix(NA, grp.nlevels, grp.nlevels)
				colnames(rmat) <- rownames(rmat) <- grp.levels
				for (i in 1:(grp.nlevels-1)) {
					grp.level1 <- grp.levels[i]
					for (j in (i+1):grp.nlevels) {
						
						grp.level2 <- grp.levels[j]
						cat(grp.level1, ' vs ', grp.level2, '\n')
						ind <- grp %in% c(grp.level1, grp.level2)
						df2 <- subset(df, ind)
						df2[, grp.name] <- factor(df2[, grp.name])
						dist.mat <- dist.obj[[dist.name]][ind, ind]
						strata2 <- strata[ind]
						
						if (is.null(adj.name)) {
							formula <- paste('~', grp.name)
						} else {
							formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
						}
						
						ind2 <- apply(df2[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
						df2 <- df2[ind2, ]
						dist.mat2 <- as.dist(dist.mat[ind2, ind2])
						strata2 <- strata2[ind2]
						if (block.perm == F) {
							obj <- adonis(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2, ...)
						} else {
							obj <- adonis2(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2,  ...)
						}
						prmatrix(obj$aov.tab)
						cat("\n")
	
						if (block.perm == F) {
							pmat[i, j] <- pmat[j, i] <- obj$aov.tab[length(adj.name)+1, 6]
							rmat[i, j] <- rmat[j, i] <- obj$aov.tab[length(adj.name)+1, 5]
						} else {
							pmat[i, j] <- pmat[j, i] <- obj$aov.tab[1, 6]
							rmat[i, j] <- rmat[j, i] <- obj$aov.tab[1, 5]
						}
						
						# PERMANOVA G after last distance
						if (block.perm == F) {
							if (dist.name == dist.names[length(dist.names)]) {
								cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
								response <- array(NA, c(sum(ind2), sum(ind2), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
								for (dist.name in PermanovaG.dist) {
									response[, , dist.name] <- dist.mat[ind2, ind2]
								}
								obj <- PermanovaG(as.formula(paste("response", formula)), df2,  strata=strata2, ...)
								prmatrix(obj$aov.tab)
								cat("\n")
								pmat.G[i, j] <- pmat.G[j, i] <- obj$aov.tab[length(adj.name)+1, 2]
							}
						}
					}
				}
				cat("\n")
				pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
				rmat.all <- rbind(rmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(rmat), rep("", grp.nlevels))
			}
			cat("\n")
			sink()
			write.csv(pmat.all, paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwiseP.csv'))
			write.csv(rmat.all, paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwiseR.csv'))
			write.csv(pmat.G, paste0('Beta_diversity_PERMANOVA_G_test_', ann, '_pairwiseP.csv'))
		}
	}
}

randomForestTest <- function(x, y, perm.no=999, ...) {
	iris.rf <- randomForest(x=x, y=y, importance=FALSE)
	to1 <- mean(abs((as.numeric(y) - 1) - predict(iris.rf,  type='prob')[, 2])^2) #mean squarred error
	to2 <- mean(y !=  predict(iris.rf,   type='response'))  #mean prediction error
	tp <- sapply(1:perm.no, function(i) {
				if (i %% 10 == 0) cat('.')
				y.p <- sample(y)
				iris.rf <- randomForest(x=x, y=y.p, importance=FALSE)
				
				t1 <- mean(abs((as.numeric(y.p) - 1) - predict(iris.rf,  type='prob')[, 2])^2) 
				t2 <- mean(y.p !=  predict(iris.rf,   type='response'))
				c(t1, t2)
			}
	
	)
	pv1 <- (sum(tp[1, ] <= to1) + 1) / (perm.no + 1)
	pv2 <- (sum(tp[2, ] <= to2) + 1) / (perm.no + 1)
	c(pv1=pv1, pv2=pv2)
}

perform_rf_test <- function (data.obj, grp.name, taxa.level='Genus', perm.no=999, prev=0.1, minp=0.000, ann='',...) {
	if (taxa.level == 'Species') {
		if (taxa.level %in% names(data.obj$abund.list)) {
			ct <- data.obj$abund.list[[taxa.level]]
		} else {
			# Accomodate different version
			ct <- data.obj$otu.tab
			rownames(ct) <- paste0("OTU", rownames(ct), ":", data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
			data.obj$abund.list[['Species']] <- ct
		}
	} else {
		ct <- data.obj$abund.list[[taxa.level]]
	}
	
	prop <- t(t(ct) / colSums(ct))
	prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
	prop <- t(prop)
	
	grp <- data.obj$meta.dat[, grp.name]
	if (!is.factor(grp)) stop('Current random Forest test is designed to deal with binary factor data!\n')

	cat('RF test P values  ...\n')
	obj <- randomForestTest(prop, grp, perm.no, ...)
	sink(paste0('OverallAssoc_RF_test_', taxa.level, '_', ann, '.txt'))
	cat ('Prediction Probability:', obj[1], '\n')
	cat ('Binary Response:', obj[2], '\n')
	sink()
}

perform_mirkat_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'WUniFrac', 'BC'), 
		grp.name=NULL, adj.name=NULL, pairwise=F,  ann='', ...) {
	
	# MiRKAT not implemented for correlated data
	df <- data.obj$meta.dat
	
	if (pairwise == F) {
		
		ind <- apply(df[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
		df <- df[ind, ]
		
		grp <- df[, grp.name]
		if (is.factor(grp)) {
			grp <- as.numeric(grp) - 1
		}
		if (!is.null(adj.name)) {
			# No intercept
			adj <- model.matrix(~ ., data.frame(df[, adj.name]))
			# Remove collinear terms
			qadj <- qr(adj, tol = 1e-07)
			adj <- adj[, qadj$pivot, drop = FALSE]
			adj <- adj[, 1:qadj$rank, drop = FALSE]
			# Remove intercept
			adj <- adj[, colSums(adj==1) != nrow(adj)]
		} else {
			adj <- NULL
		}
		
		sink(paste0('Beta_diversity_MiRKAT_test_', ann, '.txt'))
		date()
		cat('\nMiRKAT test combining ', paste(dist.names, collapse=','), '\n')
		Ks <- list()
		for (dist.name in dist.names) {
			Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][ind, ind])
		}
		obj <- MiRKAT(grp, X=adj, Ks)
		cat('Individual P value: ')
		prmatrix(t(obj$indivP))
		cat('\nOmnibus P value: ')
		cat(obj$omnibus_p)
		cat("\n")
		sink()
		
	} else {
		sink(paste0('Beta_diversity_MiRKAT_test_', ann, '_pairwise.txt'))
		date()
		cat('\nPairwise MiRKAT test: \n')
		grp <- factor(df[, grp.name])
		grp.levels <- levels(grp)
		grp.nlevels <- nlevels(grp)
		pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
		colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
		parr <- array(NA, c(grp.nlevels, grp.nlevels, length(dist.names)), dimnames=list(grp.levels, grp.levels, dist.names))
		
		for (i in 1:(grp.nlevels-1)) {
			grp.level1 <- grp.levels[i]
			for (j in (i+1):grp.nlevels) {
				
				grp.level2 <- grp.levels[j]
				cat(grp.level1, ' vs ', grp.level2, '\n')
				ind <- grp %in% c(grp.level1, grp.level2)
				df2 <- subset(df, ind)
				df2[, grp.name] <- factor(df2[, grp.name])
				ind2 <- apply(df2[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
				df2 <- df[ind2, ]
				
				grp2 <- as.numeric(df2[, grp.name]) - 1
				if (!is.null(adj.name)) {
					# No intercept
					adj <- model.matrix(~ ., data.frame(df2[, adj.name]))
					# Remove collinear terms
					qadj <- qr(adj, tol = 1e-07)
					adj <- adj[, qadj$pivot, drop = FALSE]
					adj <- adj[, 1:qadj$rank, drop = FALSE]
					# Remove intercept
					adj <- adj[, colSums(adj==1) != nrow(adj)]
				} else {
					adj <- NULL
				}
				
				Ks <- list()
				for (dist.name in dist.names) {
					Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][rownames(df2), rownames(df2)])
				}
				obj <- MiRKAT(grp2, X=adj, Ks)
				pmat.G[i, j] <- pmat.G[j, i] <- obj$omnibus_p
				parr[i, j, ] <- parr[j, i, ] <- obj$indivP
				cat('Individual P value: ')
				prmatrix(t(obj$indivP))
				cat('\nOmnibus P value: ')
				cat(obj$omnibus_p)
				cat("\n")
				
			}
		}
		cat("\n")
		sink()
		
		pmat.all <- NULL
		for (dist.name in dist.names) {
			pmat <- parr[, , dist.name]
			pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
		}
		
		write.csv(pmat.all, paste0('Beta_diversity_MiRKAT_test_', ann, '_pairwiseP.csv'))
		write.csv(pmat.G, paste0('Beta_diversity_MiRKAT_O_test_', ann, '_pairwiseP.csv'))
	}
}


perform_betadisper_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), grp.name) {
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	sink('Beta_diversity_BETADISPER_test.txt')
	date()
	cat('\nBetadisper test: \n')
	for (dist.name in dist.names) {
		cat(dist.name, " distance: \n")
		dist.mat <- as.dist(dist.obj[[dist.name]])
		obj <- betadisper(dist.mat, grp)
		prmatrix(anova(obj))
		cat("\n")
	}
	cat("\n")
	sink()
}

distance_compare_test <- function (dist.mat, ind1, ind2, ind3, ID2=NULL, ID3=NULL, alternative='greater', nperm=999) {
	
	ind23 <- c(ind2, ind3)
	if (!is.null(ID2) & !is.null(ID3)) {
		ID23 <- c(ID2, ID3)
		ID2u <- unique(ID2)
		ID3u <- unique(ID3)
		ID23u <- c(ID2u, ID3u)
		n2u <- length(ID2u)
		n3u <- length(ID3u)
	}
	n2 <- length(ind2)
	n3 <- length(ind3)
	
	dist12 <- dist.mat[ind1, ind2]
	dist13 <- dist.mat[ind1, ind3]
	
	stat.obs <- mean(dist12) - mean(dist13)
	stat.perm <- sapply(1:nperm, function(i) {
				if (is.null(ID2) & is.null(ID3)) {
					ind23.p <- sample(ind23)
					ind2.p <- ind23.p[1:n2]
					ind3.p <- ind23.p[(n2+1):(n2+n3)]
				} else {
					ID23u.p <- sample(ID23u)
					ID2.p <- ID23u.p[1:n2u]
					ID3.p <- ID23u.p[(n2u+1):(n2u+n3u)]
					ind2.p <- ind23[ID23 %in% ID2.p]
					ind3.p <- ind23[ID23 %in% ID3.p]
				}

				dist12.p <- dist.mat[ind1, ind2.p]
				dist13.p <- dist.mat[ind1, ind3.p]
				mean(dist12.p) - mean(dist13.p)
			})
	if (alternative == 'greater') {
		pv <- mean(c(stat.perm >= stat.obs, TRUE))
	} else {
		pv <- mean(c(stat.perm <= stat.obs, TRUE))
	}  
	pv
}

perform_distance_compare_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name, level1, level2, level3, subject=NULL, alternative='greater', nperm=999, seed=123, ann='') {
	sink(paste0('Beta_diversity_test_', level1, '-', level2, '_', alternative, '_', level1, '-', level3, '_', ann, '.txt'))
	set.seed(seed)
	cat('Testing the distance ', level1, '-', level2, ' is ', alternative, ' than ', level1, '-', level3, '\n')
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	IDs <- df[, subject]
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	ind1 <- which(grp == level1)
	ind2 <- which(grp == level2)
	ind3 <- which(grp == level3)
	if (is.null(subject)) {
		ID2 <- NULL
		ID3 <- NULL
	} else {
		ID2 <- as.character(IDs[grp == level2])
		ID3 <- as.character(IDs[grp == level3])
	}
	for (dist.name in dist.names) {
		cat(dist.name, ' Distance:p=') 
		pv <- distance_compare_test(dist.obj[[dist.name]], ind1, ind2, ind3, ID2, ID3, alternative, nperm)
		cat(pv, '\n')
	}
	sink()
	
}


perform_taxa_compare_test <- function (data.obj, grp.name, level1, level2, level3, alternative='greater', nperm=1000, seed=123,
		taxa.levels=c('Phylum', 'Family', 'Genus'), taxa.name='All', prev=0.1, minp=0.002) {
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	ind1 <- which(grp == level1)
	ind2 <- which(grp == level2)
	ind3 <- which(grp == level3)
	
	pv.list <- list()
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		prop <- t(t(prop) / colSums(prop))
		pv.vec <- NULL
		if (taxa.name == 'All') {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
		
		for (taxon in rownames(prop)) {
			cat(".")
			dist.mat <- as.matrix(dist(rank(prop[taxon, ])))
			pv <- distance_compare_test(dist.mat, ind1, ind2, ind3, alternative, nperm)
			pv.vec <- c(pv.vec, pv)
		}
		names(pv.vec) <- rownames(prop)
		pv.list[[LOI]] <- pv.vec
		cat("\n")
	}
	
	for (LOI in taxa.levels) {
		pv.vec <- pv.list[[LOI]]

		qv.vec <- p.adjust(pv.vec, 'fdr')

		
		write.csv(data.frame("p value"=pv.vec, "q value"=qv.vec), paste0("Taxa_TrendAnalysis_", LOI, ".csv"))
	}
	
}


generate_taxa_boxplot <- function (data.obj,  grp.name, strata=NULL, scale='P', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'), 
		taxa.name='All', rm.outlier=T, prev=0.1, minp=0.002, ann='All') {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		if (scale == 'logP') {
			prop <- prop + 0.5
		} 
		
		prop <- t(t(prop) / colSums(prop))
		
		headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
		names(headnames) <- rownames(prop)
		
		if (taxa.name == 'All') {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
		
		if (scale == 'logP') {
			prop <- log10(prop)
		}

		if (scale == 'sqrtP') {
			prop <- sqrt(prop)
		}
		
		if (scale == 'binary') {
			temp <- prop != 0
			prop[temp] <- 'Presence'
			prop[!temp] <- 'Absence'
		}
		hei <- 5
		if (is.null(strata)) {
			wid <- 5
		} else {
			wid <- 5.5
		}
		
		if (scale == 'P') {
			ylab <- 'Proportion'
		} 
		if (scale == 'logP') {
			ylab <- 'log10(Proportion)'
		}
		if (scale == 'sqrtP') {
			ylab <- 'sqrt(Proportion)'
		}
		if (scale == 'binary') {
			ylab <- 'Count'
		} 
		
		pdf(paste("Taxa_Boxplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
		if (is.null(strata)) {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund), quantile(taxon.abund, 0.95) * 1.25)
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					df2 <- data.frame(Value=taxon.abund, Group=grp)	
					dodge <- position_dodge(width=0.9)
					obj <- ggplot(df2, aes(x=Group, y=Value, col=Group)) +
							geom_boxplot(position=dodge, outlier.colour = NA) + 
							geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
							labs(y=ylab, title=headnames[taxon]) +
							ylim(ylims) +
							theme(legend.position="none")
					print(obj)
				} else {
					df2 <- data.frame(Value=taxon.abund, Group=grp)	

					obj <- ggplot(df2, aes(x=Group, fill=Value)) +
							geom_bar(width=.5) +
							labs(y=ylab, title=headnames[taxon]) +
							theme(legend.title=element_blank())

					print(obj)
				}

			}	
		} else {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund), quantile(taxon.abund, 0.95) * 1.25)
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					grp2 <- df[, strata]
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
					dodge <- position_dodge(width=0.9)
					obj <- ggplot(df2, aes(x=Strata, y=Value, col=Group, fill=Group)) +
							geom_point(position=position_jitterdodge(dodge.width=0.9), size=3.0, alpha=0.6) +
							geom_boxplot(fill='white', position=dodge, outlier.colour = NA) + 
							labs(y=ylab, x=strata, title=headnames[taxon]) +
							ylim(ylims)
					print(obj)
				} else {
					grp2 <- df[, strata]
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
					obj <- ggplot(df2, aes(x=Group, fill=Value)) +
							geom_bar(width=.5) +
							labs(y=ylab, title=headnames[taxon]) +
							facet_wrap(~ Strata) + 
							theme(legend.title=element_blank())
					
					print(obj)
				}
			}
			
		}
		dev.off()
	}		
}


generate_taxa_scatterplot <- function (data.obj,  grp.name, strata=NULL, scale='P', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'), 
		taxa.name='All', rm.outlier=T, prev=0.1, minp=0.002, ann='All') {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		if (scale == 'logP') {
			prop <- prop + 0.5
		} 
		
		prop <- t(t(prop) / colSums(prop))
		
		headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
		names(headnames) <- rownames(prop)
		
		if (taxa.name == 'All') {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
		
		if (scale == 'logP') {
			prop <- log10(prop)
		}
		
		if (scale == 'sqrtP') {
			prop <- sqrt(prop)
		}
		
		if (scale == 'binary') {
			temp <- prop != 0
			prop[temp] <- 'Presence'
			prop[!temp] <- 'Absence'
		}
		hei <- 5
		if (is.null(strata)) {
			wid <- 5
		} else {
			wid <- 5.5
		}
		
		if (scale == 'P') {
			ylab <- 'Proportion'
		} 
		if (scale == 'logP') {
			ylab <- 'log10(Proportion)'
		}
		if (scale == 'sqrtP') {
			ylab <- 'sqrt(Proportion)'
		}
		if (scale == 'binary') {
			ylab <- 'Count'
		} 
		
		pdf(paste("Taxa_Scatterplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
		if (is.null(strata)) {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund), quantile(taxon.abund, 0.95) * 1.25)
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					df2 <- data.frame(Value=taxon.abund, Group=grp)	
					dodge <- position_dodge(width=0.9)
					obj <- ggplot(df2, aes(x=Group, y=Value)) +
							geom_point() +
							geom_smooth() +
							labs(y=ylab, x=grp.name, title=headnames[taxon]) +
					#		ylim(ylims) +
							theme(legend.position="none")
					print(obj)
				} else {
					df2 <- data.frame(Value=taxon.abund, Group=grp)	
					
					obj <- ggplot(df2, aes(y=Group, x=Value)) +
							geom_boxplot() +
							geom_jitter() +
							labs(y=grp.name, x='Presence/Absence', title=headnames[taxon]) +
							theme(legend.title=element_blank())
					
					print(obj)
				}
				
			}	
		} else {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund), quantile(taxon.abund, 0.95) * 1.25)
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					grp2 <- df[, strata]
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
					obj <- ggplot(df2, aes(x=Group, y=Value)) +
							geom_point() +
							geom_smooth() +
							labs(y=ylab, x=grp.name, title=headnames[taxon]) +
							facet_wrap(~ Strata) + 
							theme(legend.title=element_blank())			
					print(obj)
					
				} else {				
					grp2 <- factor(df[, strata])
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
					dodge <- position_dodge(width=0.9)
					obj <- ggplot(df2, aes(x=Strata, y=Group, col=Value, fill=Value)) +
							geom_point(position=position_jitterdodge(dodge.width=0.9), size=3.0, alpha=0.6) +
							geom_boxplot(fill='white', position=dodge, outlier.colour = NA) + 
							labs(y=grp.name, x=strata, title=headnames[taxon]) 
					print(obj)
				}
			}
			
		}
		dev.off()
	}		
}

taxa_barplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='sqrt',  xsize=10, ylab='Proportion', error='se', cutoff=0.00001) {
	
	if (scale == 'log') {
		prop[prop <= cutoff] <- cutoff
	}
	
    grp <- factor(df[, grp.name])
	
	if (is.null(strata)) {
		df2 <- data.frame(Group=grp, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Taxa', 'Value')
		
		# Could be revised
		temp1 <- aggregate(Value ~ Group + Taxa, df2, mean)

		if (error == 'se') {
			temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) sd(x) / sqrt(length(x)))
		}
		if (error == 'sd') {
			temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) sd(x))
		}
		
		
		df2 <- cbind(temp1, temp2[, 3])
		colnames(df2) <- c('Group', 'Taxa', 'Mean', 'SE')
		
		limits <- aes(ymax = Mean + SE, ymin = Mean - SE)
		dodge <- position_dodge(width=0.95)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) + 
				geom_bar(position=dodge, stat="identity") + 
				geom_bar(position=dodge, stat="identity", colour="black", show_guide=FALSE, size=0.25) +
				geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		
		# To be revised
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	} else {
		grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
		df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')
		
		# Could be revised
		temp1 <- aggregate(Value ~ Group + Strata + Taxa, df2, mean)
		if (error == 'se') {
			temp2 <- aggregate(Value ~ Group + Strata + Taxa, df2, function(x) sd(x) / sqrt(length(x)))
		}
		if (error == 'sd') {
			temp2 <- aggregate(Value ~ Group + Strata + Taxa, df2, function(x) sd(x))
		}
		
		
		df2 <- cbind(temp1, temp2[, 4])
		colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Mean', 'SE')
		
		limits <- aes(ymax = Mean + SE, ymin = Mean - SE)
		dodge <- position_dodge(width=0.95)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) + 
				geom_bar(position=dodge, stat="identity") + 
				geom_bar(position=dodge, stat="identity", colour="black", show_guide=FALSE, size=0.25) +
				geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
				facet_wrap(~Strata, ncol=1) +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	}
	return(obj1)
	
}

taxa_boxplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='none', xsize=10, ylab='Proportion', cutoff=0.00001) {
	
	
	grp <- factor(df[, grp.name])
	
	if (scale == 'log') {
		prop[prop <= cutoff] <- cutoff
	}
	
	if (is.null(strata)) {
		df2 <- data.frame(Group=grp, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Taxa', 'Value')
		

		dodge <- position_dodge(width=0.95)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) + 
				geom_boxplot(position=dodge) + 
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	} else {
		grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
		df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')
		
		dodge <- position_dodge(width=0.95)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) + 
				geom_boxplot(position=dodge) + 
				facet_wrap(~Strata, ncol=1) +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		# To be revised
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	}
	return(obj1)
	
}

generate_taxa_biplot <- function (data.obj, taxa, trans='sqrt', grp.name, ann='', ...) {
	
    grp <- data.obj$meta.dat[, grp.name]

	prop <- NULL
	for (LOI2 in names(data.obj$abund.list)) {
		ct <- data.obj$abund.list[[LOI2]]
		ct <- ct + 1
		prop0 <- t(t(ct) / colSums(ct))
		prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])	
	}
	colnames(prop) <- colnames(prop0)
	if (nrow(prop) != length(taxa)) {
		warnings('Some taxa not found in abundance lists! Please check the names!\n')
	}
	
	if (trans == 'normal') 	prop <- t(apply(prop, 1, function(x) qqnorm(x, plot=F)$x))
	if (trans == 'log') prop <- log(prop)
	if (trans == 'sqrt') prop <- sqrt(prop)
	if (trans == 'rank') prop <- t(apply(prop, 1, rank))
	
	wine.pca <- prcomp(t(prop), scale. = TRUE)
	g <- ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, 
			groups = grp, ellipse = TRUE, circle = FALSE, ...) 
	g <- g + scale_color_discrete(name = '') + theme_bw()
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top') 
	pdf(paste0("Taxa_Biplot_", ann, ".pdf"), height=6, width=6)
	print(g)
	dev.off()
}

#  Ad presence and absence bar
generate_taxa_barplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
		taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8)) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	for (i in 1:length(taxa.levels)) {
		LOI <- taxa.levels[i]
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		prop <- t(t(prop) / colSums(prop))
				
		if (taxa.name == 'All') {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
			
		# decreasing
		prop <- prop[rev(order(rowMeans(prop))), ]
		
		if (is.null(wids) | is.null(heis)) {
			wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
	        hei <- 7
		} else {
			wid <- wids[i]
			hei <- heis[i]
		}
			
		pdf(paste("Taxa_Barplot_Aggregate", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)		
		obj1 <- taxa_barplot_aggregate (prop, df, grp.name, strata, scale, xsize[i]) 
		print(obj1)
		dev.off()		
	}
}

generate_taxa_boxplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
		taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8)) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	
	for (i in 1:length(taxa.levels)) {
		LOI <- taxa.levels[i]
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		prop <- t(t(prop) / colSums(prop))
		
		if (taxa.name == 'All') {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
		
		# decreasing
		prop <- prop[rev(order(rowMeans(prop))), ]
		
		if (is.null(wids) | is.null(heis)) {
			wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
			hei <- 7
		} else {
			wid <- wids[i]
			hei <- heis[i]
		}
		
		pdf(paste("Taxa_Boxplot_Aggregate", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)		
		obj1 <- taxa_boxplot_aggregate (prop, df, grp.name, strata, scale, xsize[i]) 
		print(obj1)
		dev.off()		
	}
}

# Add combined barplot with error bar and presence/absence bar (currently presence/absence bar is in generate_taxa_boxplot
generate_taxa_barplot <- function (data.obj,  grp.name, strata=NULL, scale='P', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
		taxa.name='All', rm.outlier=T, prev=0.1, minp=0.002, ann='All') {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		if (scale == 'logP') {
			prop <- prop + 0.5
		} 
		prop <- t(t(prop) / colSums(prop))
		
		headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
		names(headnames) <- rownames(prop)
		
		if (taxa.name == 'All') {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
		
		if (scale == 'logP') {
			prop <- log10(prop)
		}
		
		
		hei <- 5
		if (is.null(strata)) {
			wid <- 5
		} else {
			wid <- 4 * nlevels(df[, strata])
		}
		
		if (scale == 'P') {
			ylab <- 'Proportion'
		} else {
			ylab <- 'log10(Proportion)'
		}
		
		pdf(paste("Taxa_Barplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
		if (is.null(strata)) {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				taxon.abund2 <- taxon.abund
				if (scale == 'P') {
					if (rm.outlier == T) {
						ylims <- c(0, quantile(taxon.abund, 0.95) * 1.25)
						taxon.abund[taxon.abund > quantile(taxon.abund, 0.95) * 1.25] <- quantile(taxon.abund, 0.95) * 1.25
					} else {
						ylims <- c(0, max(taxon.abund))
					}
					
				} else {
					ylims <- range(taxon.abund)
				}
				df2 <- data.frame(x=factor(1:length(taxon.abund)), Value=taxon.abund, Group=grp)			
				dodge <- position_dodge(width=0.99)
				obj <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
						geom_bar(position=dodge, stat='identity' ) + 
						facet_grid(. ~ Group, scales='free_x', space="free") +
						labs(y=ylab, title=headnames[taxon]) +
						ylim(ylims) +
						xlab('') +
						theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
						theme(legend.position="none")
				mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, mean))
				obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint))		
				mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, median))
				obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)
				
				print(obj)
			}	
		} else {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				taxon.abund2 <- taxon.abund
				if (scale == 'P') {
					if (rm.outlier == T) {
						ylims <- c(0, quantile(taxon.abund, 0.95) * 1.25)
						taxon.abund[taxon.abund > quantile(taxon.abund, 0.95) * 1.25] <- quantile(taxon.abund, 0.95) * 1.25
					} else {
						ylims <- c(0, max(taxon.abund))
					}
					
				} else {
					ylims <- range(taxon.abund)
				}
				
				grp2 <- df[, strata]
				obj.list <- list()
				for (level in levels(grp2)) {
					ind <- grp2 %in% level
					df2 <- data.frame(x=factor(1:length(taxon.abund[ind])), Value=taxon.abund[ind], Group=grp[ind])	
	
					dodge <- position_dodge(width=0.99)
					obj0 <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
							geom_bar(position=dodge, stat='identity' ) + 
							facet_grid(. ~ Group, scales='free_x', space="free") +
							labs(y=ylab, title=headnames[taxon]) +
							ylim(ylims) +
							xlab(paste(strata, level)) +
							theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
							theme(legend.position="none")
					mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], mean))
					obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint))		
					mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], median))
					obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)
					obj.list[[level]] <- obj0		
				}
				multiplot(plotlist=obj.list, cols=nlevels(grp2))
			}	
			
		}
		dev.off()
	}		
}


twopart.test <- function(x1, x2, zero.p=0.2) {
	
	n1 <- length(x1)
	n2 <- length(x2)
	p1 <- mean(x1 != 0)
	p2 <- mean(x2 != 0)
	m1 <- sum(x1 != 0)
	m2 <- sum(x2 != 0)
	p12 <- (m1 + m2) / (n1 + n2)
	q12 <- 1 - p12
	
	if (q12 >= zero.p) {
		Z <- (abs(p1 - p2) - (1/(2*n1) + 1/(2*n2))) / sqrt(p12 * q12 * (1/n1 + 1/n2))
		x1 <- x1[x1!=0]
		x2 <- x2[x2!=0]
		R1 <- sum(rank(c(x1, x2))[1:length(x1)])
		ti <- as.vector(table(c(x1, x2)))
		W <- (abs(R1 - m1*(m1+m2+1)/2) - 1/2) / sqrt((m1*m2/12)*(m1+m2+1-sum(ti*(ti^2-1))/(m1+m2)/(m1+m2-1)))
		X2 <- Z^2 + W^2
		res <- list()
		res$stat <- X2
		res$p.value <- 1 - pchisq(X2, 2)
		res$Z <- Z
		res$W <- W
		res$test <- 'TwoPart'
	} else {
		res <- wilcox.test(x1, x2)
		res$test <- 'Wilcox'
	}
	res
}

getPermuteMatrix <- function (perm, N, strata = NULL) 
{
	if (length(perm) == 1) {
		perm <- how(nperm = perm)
	}
	if (!missing(strata) && !is.null(strata)) {
		if (inherits(perm, "how") && is.null(getBlocks(perm))) 
			setBlocks(perm) <- strata
	}
	if (inherits(perm, "how")) 
		perm <- shuffleSet(N, control = perm)
	if (is.null(attr(perm, "control"))) 
		attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"), 
						nperm = nrow(perm)), class = "how")
	perm
}

# Need to be investigated
# a. Speed up 
# b. Permutation method （response, covariate or residual permutation) - Use PERMANOVA procedure - Response Permutation (Has some problem)
# c. Revise - add permutation stratified by subject
permute_differential_analysis <- function (meta.dat, prop, grp.name, adj.name=NULL, strata=NULL, sqrt.trans=TRUE, resid.perm=TRUE, perm.no=999) {
	# Square root transformation
	if (sqrt.trans) {
		Y <- sqrt(prop)
	} else {
		Y <- prop
	}

	if (resid.perm) {
		# Permute the residual
		if (is.null(adj.name)) {
			Y <- t(resid(lm(as.formula(paste('t(Y) ~ 1')), meta.dat)))
		} else {
			Y <- t(resid(lm(as.formula(paste('t(Y) ~ ', adj.name)), meta.dat)))
		}

	}
	
	n <- ncol(prop)
	I <- diag(n)
	if (is.null(adj.name)) {
		M0 <- model.matrix(~ 1, meta.dat)
	} else {
		df0 <- meta.dat[, c(adj.name), drop=F]
		M0 <- model.matrix( ~., df0)
	}
	
	qrX0 <- qr(M0, tol = 1e-07)
	Q0 <- qr.Q(qrX0)
	Q0 <- Q0[, 1:qrX0$rank, drop=FALSE]
	
	
#	P0 <- M0 %*% solve(t(M0) %*% M0) %*% t(M0)
		
	df1 <- meta.dat[, c(adj.name, grp.name), drop=F]
	M1 <- model.matrix( ~., df1)
	
	qrX1 <- qr(M1, tol = 1e-07)
	Q1 <- qr.Q(qrX1)
	Q1<- Q1[, 1:qrX1$rank, drop=FALSE]
	
	TSS <- rowSums(Y^2)
	MSS1 <- rowSums((Y %*% Q1)^2)
	MSS0 <- rowSums((Y %*% Q0)^2)
	F0 <- (MSS1 - MSS0) /  (TSS - MSS1) 
	
#	P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
#	F0 <- diag(Y %*% (P1 - P0) %*% t(Y)) / diag(Y %*% (I - P1) %*% t(Y))
#	df3 <- df1
	

	p <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
	perm.no <- nrow(p)
	
	Fp <- sapply(1:perm.no, function(i) {
				if (i %% 100 == 0) cat('.')
				Yp <- Y[, p[i, ]]
#				df3[, grp.name] <- sample(df1[, grp.name])
#				M1 <- model.matrix( ~., df3)
#				P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
				MSS1p <- rowSums((Yp %*% Q1)^2)
				MSS0p <- rowSums((Yp %*% Q0)^2)
				(MSS1p - MSS0p) /  (TSS - MSS1p) 
			})
	res <- cbind(Fp >= F0, 1)
	res <- rowMeans(res)
	res[is.na(res)] <- 1
	res
}


perform_differential_analysis <- function (data.obj, method=NULL, grp.name, adj.name=NULL, subject=NULL, 
		taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus'), 
		prev=0.1, minp=0.002, medianp=NULL, mt.method='fdr', cutoff=0.15, ann='', seed=123, ...) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	ind <- !is.na(grp)
	data.obj <- subset_data(data.obj, ind)
	grp <- grp[ind]
	df <- df[ind, ]
	
	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}
	
	dep <- colSums(data.obj$otu.tab)
	diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
	
	if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
		cat("Signficant sequencing depth confounding!\n For nonparametric test, rarefaction should be performed first!\n")
		cat('Without rarefaction, there will be many false postives (rare taxa)!\n')
		cat('Automatic rarefaction to the minimum seqence depth is performed!\n')
#		cat("For parametric test with sequence depth adjustment (DESeq2), please be cautious about the results!\n There may be potential residual sequence depth confounding!\n")
		# potential revising
		otu.tab <- t(Rarefy(t(data.obj$otu.tab))$otu.tab.rff)
		data.obj <- subset_data(data.obj, colnames(otu.tab))
		data.obj$otu.tab <- otu.tab
		otu.name <- data.obj$otu.name
		# Re-summarize the taxa data
		for (hierach in taxa.levels) {
			if (hierach != 'Species') {
				if (hierach != 'Phylum') {
					tax.family <- paste(otu.name[, 'Phylum'], otu.name[, hierach], sep=";")
				} else {
					tax.family <- otu.name[, 'Phylum']
				}
				family <- aggregate(otu.tab, by=list(tax.family), FUN=sum)
				rownames(family) <- family[, 1]
				family <- as.matrix(family[, -1])
				data.obj$abund.list[[hierach]] <- family
			} else {
				data.obj$abund.list[['Species']][, ] <- data.obj$otu.tab
			}
		}
	}

	
	if (is.factor(grp)) {
		if (is.null(method)) {
			if (nlevels(grp) == 2) {
				method <- 'wilcox'
			} else {
				method <- 'kruskal'
			}
		}
		if (method == 'wilcox.pair') {
			if (nlevels(grp) != 2) stop("Wilcox test requires two groups!\n")
			if (is.null(subject)) stop("Paired wilcox needs subject information!\n")
			subject <- factor(df[, subject])
			ind1 <- ind2 <- NULL
			for(sub in levels(subject)) {
				temp1 <- which(as.numeric(grp) == 1 & subject == sub)
				temp2 <- which(as.numeric(grp) == 2 & subject == sub)
				if (length(temp1) != 0 & length(temp2) != 0) {
					ind1 <- c(ind1, temp1)
					ind2 <- c(ind2, temp2)
				}
				
			}
		#	if (length(ind1) != length(ind2))  warning('Some subjects are not paired!\n')
		}
		
		if (method == 'perm.pair') {
			subject <- factor(df[, subject])
		}
		
		
		pv.list <- qv.list <-  fc.list <- pc.list <- m.list <- nzm.list  <- prv.list <- list()
		res.final <- NULL
		for (LOI in taxa.levels) {
			cat(LOI, "\n")
			ct <- data.obj$abund.list[[LOI]]
			prop <- t(t(ct) / colSums(ct))
			if (!is.null(minp)) {
				prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
				ct <- ct[rownames(prop), , drop=FALSE]
			}
			
			if (!is.null(medianp)) {
				nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
				prop <- prop[nz.mean > medianp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
				ct <- ct[rownames(prop), , drop=FALSE]
			}
			
#			if (method == 'DESeq2') {
#				if (is.null(adj.name)) {
#					formula <- as.formula(paste('~', grp.name))
#				} else {
#					formula <- as.formula(paste('~', paste(paste(adj.name, collapse='+'), '+', grp.name)))
#				}
#				dds <- DESeqDataSetFromMatrix(countData = ct + 1,
#						colData = df,
#						design = as.formula(formula))
#				dds <- DESeq(dds)
#				res <- results(dds)
#				pv.de2 <- res[, 'pvalue']
#				names(pv.de2) <- rownames(res)
#			}
			 		
			if (method == 'perm') {
				set.seed(seed)
				pv.de2 <- permute_differential_analysis(df, prop, grp.name, adj.name, ...)
				names(pv.de2) <- rownames(prop)
			}
			
			if (method == 'perm.pair') {
				set.seed(seed)
				pv.de2 <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject, ...)
				names(pv.de2) <- rownames(prop)
			}
			
			pv.vec <- m.vec <- nzm.vec <- prv.vec <-  fc.vec <- pc.vec <-  NULL
			
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				
				pv <- fc <- pc <- m <- nzm <- prv <- NULL
				if (method == 'wilcox') {
					if (nlevels(grp) != 2) stop("Wilcox test requires two groups!\n")
					pv <- wilcox.test(taxon.abund ~ grp)$p.value
					
				}
				if (method == 'wilcox.pair') {
					pv <- wilcox.test(taxon.abund[ind1], taxon.abund[ind2], paired=T)$p.value
				}
				if (method == 'twopart') {
					if (nlevels(grp) != 2) stop("Two part test requires two groups!\n")
					grp1 <- taxon.abund[as.numeric(grp)==1]
					grp2 <- taxon.abund[as.numeric(grp)==2]
					pv <- twopart.test(grp1, grp2)$p.value
				}	
				if (method == 'kruskal') {
					if (nlevels(grp) <= 2) warning("Kruskal-wallis test requires three or more groups!\n")
					pv <- kruskal.test(taxon.abund ~ grp)$p.value
				}
#				if (method == 'DESeq2') {
#					pv <- pv.de2[taxon]
#				}
				
				if (method == 'perm') {
					pv <- pv.de2[taxon]
				}
				
				if (method == 'perm.pair') {
					pv <- pv.de2[taxon]
				}
				m <- tapply(taxon.abund, grp, function(x) mean(x))
				nzm <- tapply(taxon.abund, grp, function(x) mean(x[x != 0]))
				prv <- tapply(taxon.abund, grp, function(x) sum(x != 0))				
				
				if (nlevels(grp) == 2) {
					grp.no <- table(grp)
					fc <- log2(m[1] / m[2])
					pc <- prv[1] / grp.no[1] / prv[2] * grp.no[2]
				} else {
					pc <- fc <- NA
				}
				
				pv.vec <- rbind(pv.vec, pv)
				fc.vec <- rbind(fc.vec, fc)
				m.vec <- rbind(m.vec, m)
				nzm.vec <- rbind(nzm.vec, nzm)
				pc.vec <- rbind(pc.vec, pc)
				prv.vec <- rbind(prv.vec, prv / table(grp))	
			}
			
			
			temp <- p.adjust(pv.vec[, 1], 'fdr')
			
			qv.vec <- matrix(temp, ncol=1)
			
			rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(pc.vec) <- rownames(m.vec) <- rownames(nzm.vec) <- rownames(prv.vec) <- rownames(prop)
			colnames(pv.vec) <- 'Pvalue'
			colnames(qv.vec) <- 'Qvalue'
			colnames(fc.vec) <- 'logFoldChange'
			colnames(pc.vec) <- 'PrevalChange'
			colnames(m.vec) <- paste(levels(grp), 'Mean')
			colnames(nzm.vec) <- paste(levels(grp), 'nzMean')
			colnames(prv.vec) <- paste(levels(grp), 'preval')
			
			pv.list[[LOI]] <- pv.vec
			qv.list[[LOI]] <- qv.vec
			fc.list[[LOI]] <- fc.vec
			pc.list[[LOI]] <- pc.vec
			m.list[[LOI]] <- m.vec
			nzm.list[[LOI]] <- nzm.vec
			prv.list[[LOI]] <- prv.vec
			
			res <- cbind(pv.vec, qv.vec, m.vec, nzm.vec, fc.vec, prv.vec, pc.vec)
			rownames(res) <- rownames(prop)
			write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
			
			if (mt.method == 'fdr') {
				res.final <- rbind(res.final, res[res[, 'Qvalue'] <= cutoff, , drop=F])
			}
			if (mt.method == 'raw') {
				res.final <- rbind(res.final, res[res[, 'Pvalue'] <= cutoff, , drop=F])
			}
			
		}
		
		if (!is.null(res.final)) {
			colnames(res.final) <- colnames(res)
			write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
		}
		return(list(pv.list=pv.list, fc.list=fc.list, pc.list=pc.list, qv.list=qv.list, m.list=m.list))
	} else {
		if (is.null(method)) {
			method <- 'Spearman'
		}
		# Continuous case - currently only has DESeq2 
		pv.list <- qv.list <-  fc.list <-  m.list <- list()
		res.final <- NULL
		for (LOI in taxa.levels) {
			cat(LOI, "\n")
			ct <- data.obj$abund.list[[LOI]]
			prop <- t(t(ct) / colSums(ct))
			if (!is.null(minp)) {
				prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
				ct <- ct[rownames(prop), , drop=FALSE]
			}
			
			if (!is.null(medianp)) {
				nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
				prop <- prop[nz.mean > medianp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
				ct <- ct[rownames(prop), , drop=FALSE]
			}
			
#			if (method == 'DESeq2') {
#				if (is.null(adj.name)) {
#					formula <- as.formula(paste('~', grp.name))
#				} else {
#					formula <- as.formula(paste('~', paste(paste(adj.name, collapse='+'), '+', grp.name)))
#				}
#				dds <- DESeqDataSetFromMatrix(countData = ct + 1,
#						colData = df,
#						design = as.formula(formula))
#				dds <- DESeq(dds)
#				res <- results(dds)
#				pv.de2 <- res[, 'pvalue']
##				fc.de2 <- res[, 'log2FoldChange']
#			# May need to revise
#				fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
#				names(pv.de2) <- rownames(res)
#			}
			
			if (method == 'perm') {
				set.seed(seed)
				pv.de2 <- permute_differential_analysis(df, prop, grp.name, adj.name, ...)
				names(pv.de2) <- rownames(prop)
				# Place holder
				fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
			}
			
			if (method == 'Spearman') {
				if (!is.null(adj.name)) {
					stop("Spearman test can't adjust covariates!")
				}
				pv.de2 <- apply(prop, 1, function(x) {
							cor.test(x, grp, method='spearman')$p.value
						})
				names(pv.de2) <- rownames(prop)
				# Place holder
				fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
			}
			pv.vec <- matrix(pv.de2, ncol=1)	
			qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)
			fc.vec <- matrix(fc.de2, ncol=1)
			m.vec <- matrix(rowMeans(prop), ncol=1)
			
			rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(m.vec)  <- rownames(prop)
			colnames(pv.vec) <- 'Pvalue'
			colnames(qv.vec) <- 'Qvalue'
			colnames(fc.vec) <- 'SpearmanCorr'
			colnames(m.vec) <- 'Mean'
			
			
			pv.list[[LOI]] <- pv.vec
			qv.list[[LOI]] <- qv.vec
			fc.list[[LOI]] <- fc.vec
			m.list[[LOI]] <- m.vec
			
			
			res <- cbind(m.vec, fc.vec, pv.vec, qv.vec)
			rownames(res) <- rownames(prop)
			write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
			
			if (mt.method == 'fdr') {
				res.final <- rbind(res.final, res[res[, 'Qvalue'] <= cutoff, , drop=F])
			}
			if (mt.method == 'raw') {
				res.final <- rbind(res.final, res[res[, 'Pvalue'] <= cutoff, , drop=F])
			}
		}
		
		if (!is.null(res.final)) {
			colnames(res.final) <- colnames(res)
			write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
		}
		return(list(pv.list=pv.list, fc.list=fc.list, qv.list=qv.list, m.list=m.list))
	}
	
}

#Owls <- transform(Owls, 
#		Nest=reorder(Nest,NegPerChick), 
#		logBroodSize=log(BroodSize), 
#		NCalls=SiblingNegotiation)
#
#fit_zinbinom <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+ 
#				offset(logBroodSize)+(1|Nest), 
#		data=Owls, 
#		zeroInflation=TRUE, 
#		family="nbinom")

perform_differential_analysis_para_single_RE <- function (taxon, ldep, grp.name, adj.name=NULL, subject=NULL, df, method='NB', LRT=FALSE) {
	# ldep: log depth (size factor)
	if (!is.null(adj.name) || !is.null(subject)) {
		if (sum(grepl(grp.name, c(adj.name, subject, 'ldep')))) {
			stop('grp.name could not be part of adj.name or subject, or there will be problem!\n')
		}
	}
	
	df$ldep <- ldep
	df$taxon <- taxon
	if (LRT) warning('Currently, only Wald test is implemented!\n')
	
	if (is.null(adj.name) & is.null(subject)) {
		grp.name.adj.name.subject <- grp.name
	} else {
		if (is.null(subject)) {
			grp.name.adj.name.subject <- paste(grp.name, '+', adj.name)
		} else {
			if (grepl('\\|', subject)) {
				grp.name.adj.name.subject <- paste(grp.name, '+', adj.name, '+ (', subject, ')')
			} else {
				grp.name.adj.name.subject <- paste(grp.name, '+', adj.name, '+ (1|', subject, ')')
			}
		}
	}
	if (method == 'NB') {
		m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df, 
				zeroInflation=FALSE, family='nbinom')
	}
	
	if (method == 'B') {
		taxon2 <- as.numeric(taxon != 0)
		df$taxon2 <- taxon2
		m1.nb <- glmmadmb(as.formula(paste('taxon2 ~', grp.name.adj.name.subject, '+ ldep')), data = df, 
				zeroInflation=FALSE, family='binomial')
	}
	
	if (method == 'ZINB1') {
		m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df, 
				zeroInflation=TRUE, family='nbinom')
	}
	
	if (method == 'ZINB0') {
		m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df, 
				zeroInflation=TRUE, family='truncnbinom')
	}
	
	code <- list(m1.conv=m1.nb$convmsg)	
	pv.nb <- wald.test(b = coef(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
	method <- paste(method, 'Wald')
	
	
	aic.nb <- -2 * m1.nb$loglik + 2 * m1.nb$npar
	
	coef.nb <- coef(m1.nb)
	fc.nb <- coef.nb[grep(grp.name, names(coef.nb))]
	
	ci.nb <- confint.default(m1.nb)
	obj <- ci.nb[grep(grp.name, rownames(ci.nb)), ]
	
	if (is.vector(obj)) {
		fc.lc.nb <- obj[1]
		fc.uc.nb <- obj[2]
	} else {
		fc.lc.nb <- obj[, 1]
		fc.uc.nb <- obj[, 2]
		names(fc.lc.nb) <- paste(names(fc.lc.nb), '2.5%')
		names(fc.uc.nb) <- paste(names(fc.uc.nb), '97.5%')
	}
	return(list(method=method, pv=pv.nb, lfc=fc.nb, lfc.lci=fc.lc.nb, lfc.uci=fc.uc.nb, aic=aic.nb, code=code))
}



perform_differential_analysis_para_single_FE <- function (taxon.abund, ldep, grp.name, adj.name=NULL, subject=NULL, df, method='NB', LRT=FALSE) {
	# ldep: log depth (size factor)
	if (!is.null(adj.name) || !is.null(subject)) {
		if (sum(grepl(grp.name, c(adj.name, subject, 'ldep')))) {
			stop('grp.name could not be part of adj.name or subject, or there will be problem!\n')
		}
	}
	if (LRT & method == 'OP') warning('Overdispersed Poisson does not support LRT! Wald test used!\n')
	
	if (is.null(adj.name)) {
		grp.name.adj.name <- grp.name
	} else {
		grp.name.adj.name <- paste(grp.name, '+', adj.name)
	}
	if (method == 'NB') {
		m1.nb <- glm.nb(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')), data = df)
		if (LRT) {
			m0.nb <- update(m1.nb, as.formula(paste('. ~ . -',  grp.name)))
			code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary, m0.conv=m0.nb$converged, m0.bound=m0.nb$boundary)
			pv.nb <- anova(m1.nb, m0.nb)['Pr(Chi)'][2, ]
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary)	
			pv.nb <- wald.test(b = coef(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}

		aic.nb <- summary(m1.nb)$aic
		
		coef.nb <- coef(m1.nb)
		fc.nb <- coef.nb[grep(grp.name, names(coef.nb))]
		
		ci.nb <- confint.default(m1.nb)
		obj <- ci.nb[grep(grp.name, rownames(ci.nb)), ]
		
		if (is.vector(obj)) {
			fc.lc.nb <- obj[1]
			fc.uc.nb <- obj[2]
		} else {
			fc.lc.nb <- obj[, 1]
			fc.uc.nb <- obj[, 2]
			names(fc.lc.nb) <- paste(names(fc.lc.nb), '2.5%')
			names(fc.uc.nb) <- paste(names(fc.uc.nb), '97.5%')
		}
		return(list(method=method, pv=pv.nb, lfc=fc.nb, lfc.lci=fc.lc.nb, lfc.uci=fc.uc.nb, aic=aic.nb, code=code))
	}
	if (method == 'B') {
		taxon.abund2 <- as.numeric(taxon.abund != 0)
		m1.b <- glm(as.formula(paste('taxon.abund2 ~', grp.name.adj.name, '+ ldep')), data = df, family=binomial)
		if (LRT) {
			m0.b <- update(m1.b, as.formula(paste('. ~ . -',  grp.name)))
			code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary, m0.conv=m0.b$converged, m0.bound=m0.b$boundary)
			pv.b <- pchisq(2 * (logLik(m1.b) - logLik(m0.b)), df = df.residual(m0.b) - df.residual(m1.b), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary)	
			pv.b <- wald.test(b = coef(m1.b), Sigma = vcov(m1.b), Terms = grep(grp.name, names(coef(m1.b))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}

		aic.b <- summary(m1.b)$aic
		coef.b <- coef(m1.b)
		fc.b <- coef.b[grep(grp.name, names(coef.b))]
		
		ci.b <- confint.default(m1.b)
		obj <- ci.b[grep(grp.name, rownames(ci.b)), ]
		
		if (is.vector(obj)) {
			fc.lc.b <- obj[1]
			fc.uc.b <- obj[2]
		} else {
			fc.lc.b <- obj[, 1]
			fc.uc.b <- obj[, 2]
			names(fc.lc.b) <- paste(names(fc.lc.b), '2.5%')
			names(fc.uc.b) <- paste(names(fc.uc.b), '97.5%')
		}
		return(list(method=method, pv=pv.b, lfc=fc.b, lfc.lci=fc.lc.b, lfc.uci=fc.uc.b, aic=aic.b, code=code))
	}
	
	if (method == 'OP') {
		# No LRT
		m1.op <- glm(as.formula(paste('taxon.abund ~', grp.name.adj.name)), offset=ldep, data = df, family=quasipoisson)
		code <- list(m1.conv=m1.op$converged, m1.bound=m1.op$boundary)

		# pv.op <- pchisq(2 * (logLik(m1.op) - logLik(m0.op)), df = df.residual(m0.op) - df.residual(m1.op), lower.tail=FALSE) # LRT not applicable
		coef.op <- coef(m1.op)		
		pv.op <- wald.test(b = coef.op, Sigma = vcov(m1.op), Terms = grep(grp.name, names(coef.op)))$result$chi2['P']
		method <- paste(method, 'Wald')
		fc.op <- coef.op[grep(grp.name, names(coef.op))]
		
		ci.op <- confint.default(m1.op)
		obj <- ci.op[grep(grp.name, rownames(ci.op)), ]
		
		if (is.vector(obj)) {
			fc.lc.op <- obj[1]
			fc.uc.op <- obj[2]
		} else {
			fc.lc.op <- obj[, 1]
			fc.uc.op <- obj[, 2]
			names(fc.lc.op) <- paste(names(fc.lc.op), '2.5%')
			names(fc.uc.op) <- paste(names(fc.uc.op), '97.5%')
		}
		return(list(method=method, pv=pv.op, lfc=fc.op, lfc.lci=fc.lc.op, lfc.uci=fc.uc.op, aic=NULL, code=code))
	}
	
	if (method == 'ZINB0') {
		m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')),
				data = df, dist = "negbin", EM = TRUE)
		if (LRT) {
			if (is.null(adj.name)) {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep)')),
						data = df, dist = "negbin", EM = TRUE)
			} else {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep)')),
						data = df, dist = "negbin", EM = TRUE)
			}
			code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
			# LRT
			pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.zinb$converged)	
			pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}
		
		aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)
		
		coef.zinb <- coef(m1.zinb)
		fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
		
		ci.zinb <- confint.default(m1.zinb)
		obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
		
		if (is.vector(obj)) {
			fc.lc.zinb <- obj[1]
			fc.uc.zinb <- obj[2]
		} else {
			fc.lc.zinb <- obj[, 1]
			fc.uc.zinb <- obj[, 2]
			names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
			names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
		}
		return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
	}
	
	
	if (method == 'ZINB1') {
		m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) | ldep')),
				data = df, dist = "negbin", EM = TRUE)
		if (LRT) {
			if (is.null(adj.name)) {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
						data = df, dist = "negbin", EM = TRUE)
			} else {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) | ldep')),
						data = df, dist = "negbin", EM = TRUE)
			}
			code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
			# LRT
			pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.zinb$converged)	
			pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}
		
		aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)
		
		coef.zinb <- coef(m1.zinb)
		fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
		
		ci.zinb <- confint.default(m1.zinb)
		obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
		
		if (is.vector(obj)) {
			fc.lc.zinb <- obj[1]
			fc.uc.zinb <- obj[2]
		} else {
			fc.lc.zinb <- obj[, 1]
			fc.uc.zinb <- obj[, 2]
			names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
			names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
		}
		return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
	}
	
	
	if (method == 'ZINB2') {
		m2.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) |', grp.name.adj.name, '+ ldep')),
				data = df, dist = "negbin", EM = TRUE)
		if (LRT) {
			if (is.null(adj.name)) {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
						data = df, dist = "negbin", EM = TRUE)
			} else {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) |', adj.name, ' + ldep')),
						data = df, dist = "negbin", EM = TRUE)
			}
			
			code <- list(m1.conv=m2.zinb$converged, m0.conv=m0.zinb$converged)
			# LRT
			pv2.zinb <-  pchisq(2 * (logLik(m2.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m2.zinb), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m2.conv=m2.zinb$converged)	
			pv2.zinb <- wald.test(b = coef(m2.zinb), Sigma = vcov(m2.zinb), Terms = grep(grp.name, names(coef(m2.zinb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}

		aic2.zinb <- -2 * logLik(m2.zinb) + 2 * (m2.zinb$n - m2.zinb$df.residual)
		
		coef.zinb <- coef(m2.zinb)
		fc2.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
		
		ci.zinb <- confint.default(m2.zinb)
		obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
		
		if (is.vector(obj)) {
			fc2.lc.zinb <- obj[1]
			fc2.uc.zinb <- obj[2]
		} else {
			fc2.lc.zinb <- obj[, 1]
			fc2.uc.zinb <- obj[, 2]
			names(fc2.lc.zinb) <- paste(names(fc2.lc.zinb), '2.5%')
			names(fc2.uc.zinb) <- paste(names(fc2.uc.zinb), '97.5%')
		}
		return(list(method=method, pv=pv2.zinb, lfc=fc2.zinb, lfc.lci=fc2.lc.zinb, lfc.uci=fc2.uc.zinb, aic=aic2.zinb, code=code))
	}

}


perform_differential_analysis_para <- function (data.obj,  grp.name, adj.name=NULL, subject=NULL, RE=FALSE, method='Adaptive1', ZINB='ZINB1', LRT=FALSE, 
		taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus'), winsor=TRUE, winsor.qt=0.95, norm='GMPR', intersect.no=4,
		prev=0.1, minp=0.002, medianp=NULL, mt.method='fdr', cutoff=0.15, ann='', ...) {
	# To be completed
	# subject holds the random effects formula
	if (!RE) {
		if (!(method %in% c('ZINB', 'B', 'NB', 'OP', 'Adaptive1', 'Adaptive2'))) stop('The speficied model is not supported!\n')
		perform_differential_analysis_para_single <- perform_differential_analysis_para_single_FE
		if (!is.null(subject)) warning('subject will not be used. Are you sure you want to run fixed effects model? ')
	} else {
		if (!(method %in% c('ZINB', 'B', 'NB', 'Adaptive1', 'Adaptive2'))) stop('The speficied model does not have random effects implementation!\n')
		if (ZINB != 'ZINB1') stop('Currently only ZINB1 is supported!\n')
		if (is.null(subject)) warning('subject is not supplied. Fixed effects model will be used instead!\n')
		perform_differential_analysis_para_single <- perform_differential_analysis_para_single_RE
	}

	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	ind <- !is.na(grp)
	data.obj <- subset_data(data.obj, ind)
	grp <- grp[ind]
	df <- df[ind, ]
	
	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}
	
	dep <- colSums(data.obj$otu.tab)
	diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
	
	if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
		cat("Signficant sequencing depth confounding!\n")
		cat("For parametric test with sequence depth adjustment, please be cautious about the results!\n")
		cat("There may be potential residual sequence depth confounding!\n")
	}
	
	pv.list <- qv.list <-  fc.list <- fc.lc.list <- fc.uc.list <- met.list <- list()
	res.final <- NULL
	
	if (norm == 'GMPR') {
		dep <- GMPR(data.obj$abund.list[['Genus']], intersect.no)
	}
	
	if (norm == 'TSS') {
		dep <- colSums(data.obj$abund.list[['Genus']])
	}
	
	ldep <- log(dep)
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")

		taxon.ct <- data.obj$abund.list[[LOI]]
		

		if (winsor == TRUE) {
			# Addressing the outlier (97% percent) or at least one outlier
			
			taxon.ct.p <- t(t(taxon.ct) / dep)
			taxon.ct.p <- apply(taxon.ct.p, 1, function(x) {
						cutoff <- quantile(x, winsor.qt)
						x[x >= cutoff] <- cutoff
						x
					}
			)
			# column/row switch
			taxon.ct <- t(round(taxon.ct.p * dep))
			
		}
		
		prop <- t(t(taxon.ct) / colSums(taxon.ct))
		if (!is.null(minp)) {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
		
		if (!is.null(medianp)) {
			nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
			prop <- prop[nz.mean > medianp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
			
		pv.vec <-  fc.vec <- fc.lc.vec <- fc.uc.vec <- met.vec <- conv.vec <- NULL
		obj <- NULL
		for (taxon in rownames(taxon.ct)) {
			cat('.')
			taxon.abund <- taxon.ct[taxon, ]
			
			######## Logistic regression ###############
			if (method == 'B') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='B', LRT)) 
			######## Overdispersed Poisson regression #########	
			if (method == 'OP') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='OP', LRT)) 
			######## Negative binomial regression #########
			if (method == 'NB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT))
			######## Zeroinflated negbinomial regression 1 ########
			if (method == 'ZINB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT)) 
			
			if (method == 'Adaptive1') {
				error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT))
				error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT))
				if (class(error1) != 'try-error' & class(error2) != 'try-error') {
					if (obj1$aic < obj2$aic) {
						obj <- obj1
					} else {
						obj <- obj2
					}
					error <- error1
				} else {
					# pv == 0 indicates some problems in fitting
					if (class(error1) != 'try-error' & obj1$pv != 0) {
						obj <- obj1
						error <- error1
					} else {
						if (class(error2) != 'try-error' & obj1$pv != 0) {
							obj <- obj2
							error <- error2
						} else {
							error <- error2
						}
					}
				}	
			}
			
			if (method == 'Adaptive2') {
				error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT))
				
				if (class(error1) == 'try-error' | obj1$pv == 0) {
					error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT))
					if (class(error2) != 'try-error') {
						obj <- obj2
						error <- error2
					} else {
						error <- error2
					}
				} else {
					error <- error1
					obj <- obj1
				}
				
			}
			
			# Random Effects model 
			# ZINB, B, NB, Adpative1 is implemented based on glmmADMB
			
			
			if (class(error) == "try-error") {
				obj$pv <- obj$lfc <- obj$lfc.lci <- obj$lfc.uci <- obj$method <- NA
			}
			pv.vec <- rbind(pv.vec, obj$pv)
			fc.vec <- rbind(fc.vec, obj$lfc)
			fc.lc.vec <- rbind(fc.lc.vec, obj$lfc.lci)
			fc.uc.vec <- rbind(fc.uc.vec, obj$lfc.uci)
			met.vec <- rbind(met.vec, obj$method)
			
		}
		cat('\n')

		qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)
		
		rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(fc.uc.vec) <- rownames(fc.lc.vec) <- rownames(met.vec) <- rownames(prop)
		colnames(pv.vec) <- 'Pvalue'
		colnames(qv.vec) <- 'Qvalue'
        colnames(met.vec) <- 'Method'

		
		pv.list[[LOI]] <- pv.vec
		qv.list[[LOI]] <- qv.vec
		fc.list[[LOI]] <- fc.vec
		fc.lc.list[[LOI]] <- fc.lc.vec
		fc.uc.list[[LOI]] <- fc.uc.vec
		met.list[[LOI]] <- met.vec
 
		res <- cbind(pv.vec, qv.vec, fc.vec, fc.lc.vec, fc.uc.vec, met.vec)
		rownames(res) <- rownames(prop)
		write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
		
		if (mt.method == 'fdr') {
			res.final <- rbind(res.final, res[as.numeric(res[, 'Qvalue']) <= cutoff, , drop=F])
		}
		if (mt.method == 'raw') {
			res.final <- rbind(res.final, res[ as.numeric(res[, 'Pvalue']) <= cutoff, , drop=F])
		}
	}
	
	if (!is.null(res.final)) {
		colnames(res.final) <- colnames(res)
		res.final <- res.final[rowSums(is.na(res.final)) == 0, , drop=F]
		write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
	}
	return(list(pv.list=pv.list, qv.list=qv.list, fc.list=fc.list, fc.uc.list=fc.uc.list, fc.lc.list=fc.lc.list, met.list=met.list))
}


plot_effect_size <- function (month,  value, pos.lab, neg.lab, ylab, hjust1=1.3, hjust2=-0.3, lab.size=3, xsize=10) {
	month <- factor(month, levels=month)
	dtm <- data.frame(month=month, value=value)
	dtm$colour <- factor(ifelse(dtm$value < 0, neg.lab, pos.lab), levels=c(pos.lab, neg.lab))
	dtm$hjust <- ifelse(dtm$value > 0, hjust1, hjust2)
	obj <- ggplot(dtm, aes(month, value, label = month, hjust = hjust)) + 
			geom_text(aes(y = 0, colour = colour), size=lab.size) + 
			geom_bar(stat = "identity", aes(fill = colour)) +
			theme(axis.text.x = element_text(size=xsize)) +
			ylim(c(-max(abs(value))*1.1, max(abs(value))*1.1)) +
			coord_flip() + 
			scale_x_discrete(breaks = NULL) +
			labs(x = "", y = ylab) +
			theme(legend.position="top", legend.title=element_blank())
}


plot_effect_size2 <- function (fold.dat.plot1, ylabel='log2(Fold change)', is.ln=TRUE, ord=TRUE) {
	
	if (is.ln) {
		fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] <- fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] * log2(exp(1))
	}
	Alphas <- seq(1, 99, 2) / 100
#Multiplier <- qnorm(1 - Alphas / 2)
	Multiplier <- seq(1, 0.01, len=50)
	zzTransparency <<- 1/(length(Multiplier)/4)
	
#	fold.dat.plot1 <- data.frame(IV=rownames(fold.dat.plot), Estimate=fold.dat.plot[, 1], LCI=fold.dat.plot[, 2], UCI=fold.dat.plot[, 3])
	fold.dat.plot1 <- data.frame(cbind(fold.dat.plot1, Scalar=rep(Multiplier, each = nrow(fold.dat.plot1))))
#	fold.dat.plot1[,  -1] <- apply(fold.dat.plot1[, -1], 2, function(x){as.numeric(as.character(x))})
	fold.dat.plot1$Emphasis <- by(1 - seq(0, 1, length = length(Multiplier) + 1)[-(length(Multiplier) + 1)],
			as.character(round(Multiplier, 5)), mean)[as.character(round(fold.dat.plot1$Scalar, 5))]
	
	fold.dat.plot1$IV <- factor(fold.dat.plot1$IV, unique(fold.dat.plot1$IV))
	
	if (ord) {
		fold.dat.plot1 <- fold.dat.plot1[order(as.character(fold.dat.plot1$IV)), ]
	}
	
	OutputPlot <- ggplot2::qplot(data = fold.dat.plot1, x = IV, y = Estimate,
			ymin = Estimate - (Estimate -LCI)*Scalar, ymax = Estimate + (UCI - Estimate)*Scalar,
			ylab = NULL, xlab = NULL, alpha = I(zzTransparency), colour = I(gray(0)), geom = "blank")
	
	OutputPlot <- OutputPlot + geom_hline(yintercept = 0, lwd = I(7/12), colour = I(hsv(0/12, 7/12, 7/12)), alpha = I(5/12))
	OutputPlot <- OutputPlot + geom_linerange(data = fold.dat.plot1, aes(size = 1/Emphasis), alpha = I(zzTransparency), colour = I(gray(0)))
	OutputPlot <- OutputPlot + scale_size_continuous() + guides(size=FALSE)
#OutputPlot <- OutputPlot + facet_grid(~ ModelName)
	OutputPlot <- OutputPlot + coord_flip() + geom_point(aes(x = IV, y = Estimate), colour = I(gray(0))) + theme_bw() + ylab(ylabel)
}


visualize_differential_analysis <- function (data.obj, diff.obj,  grp.name=NULL, strata=NULL, test='Nonpara', mt.method='fdr', scale='sqrt', cutoff=0.15,
		taxa.levels=c('Phylum', 'Family', 'Genus'), ord=TRUE, eff.type='logP',
		xsize=10, ann='', hei1=NULL, wid1=NULL, hei2=NULL, wid2=NULL) {
	
	# uniquefy names
	# For backward compatibility. Newer version will not need this and below. The old version has 'unclassified'. Newer version has 'Unclassified'. Case difference.
	data.obj <- uniquefy_taxa_names(data.obj)
	
	for (name1 in names(diff.obj)) {
		obj0 <- diff.obj[[name1]]
		for (level in names(obj0)) {
			obj <- obj0[[level]]
			rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
			obj0[[level]] <- obj
		}
		diff.obj[[name1]] <- obj0
	}
	#pv.list <- diff.obj$pv.list
	fc.list <- diff.obj$fc.list
	qv.list <- diff.obj$qv.list
	pv.list <- diff.obj$pv.list
	if (test == 'Para') {
		fc.lc.list <- diff.obj$fc.lc.list
		fc.uc.list <- diff.obj$fc.uc.list
	}
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	ind <- !is.na(grp)
	data.obj <- subset_data(data.obj, ind)
	grp <- grp[ind]
	df <- df[ind, ]
	
	prop <- NULL
	eff <- eff.lc <- eff.uc <- NULL
	taxa.names <- NULL
	if (is.null(taxa.levels)) {
		LOIs <- names(qv.list)
	} else {
		LOIs <- taxa.levels
		if (sum(!(taxa.levels %in% names(qv.list)))) {
			stop('Taxa levels are not contained in differential abundance analysis results!\n')
		}
	}
	for (LOI in LOIs) {
		pv.vec <- pv.list[[LOI]]
		fc.vec <- fc.list[[LOI]]
		#qv.vec <- qvalue(pv.vec[, 1])$qvalues
		qv.vec <- qv.list[[LOI]]
		
		if (test == 'Para') {
			fc.lc.vec <- fc.lc.list[[LOI]]
			fc.uc.vec <- fc.uc.list[[LOI]]
		}
		
		if (mt.method == 'fdr') {
			taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
			taxa.name <- taxa.name[!is.na(taxa.name)]
		}

		if (mt.method == 'raw') {
			taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
			taxa.name <- taxa.name[!is.na(taxa.name)]
		}
		
		if (length(taxa.name) != 0) {
			prop0 <- data.obj$abund.list[[LOI]]
			prop0 <- t(t(prop0) / colSums(prop0))
			prop0 <-  prop0[taxa.name, , drop=F]
			if (ord == TRUE) {
				prop0 <- prop0[rev(order(rowMeans(prop0))), , drop=F]
			}
			prop <- rbind(prop, prop0)
			# currently using fold change
		    if (test == 'Para') {
				eff <- rbind(eff, fc.vec[taxa.name, , drop=F])
				eff.lc <- rbind(eff.lc, fc.lc.vec[taxa.name, , drop=F])
				eff.uc <- rbind(eff.uc, fc.uc.vec[taxa.name, , drop=F])
			} else {
				if (eff.type == 'LFC') {
					eff <- c(eff, fc.vec[taxa.name, ])
				}
				if (eff.type == 'logP') {
					eff <- c(eff, sign(fc.vec[taxa.name, ]) * (-log10(pv.vec[taxa.name, ])))
				}
				
			}
			taxa.names <- c(taxa.names, taxa.name)
		}
	}
	
	if (length(taxa.names) == 0) {
		cat('No differential taxa! \n')
	} else {
		if (length(taxa.names) >= 2) {
			if (is.null(wid1) | is.null(hei1)) {
				wid1 <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
				hei1 <- 7
			}
	
			if (!is.null(grp.name)) {
				obj <- taxa_barplot_aggregate(prop, df, grp.name, strata, scale, xsize)
				pdf(paste("Taxa_DifferentialAbundance_AbundanceBarplot", scale, mt.method, cutoff, ann, ".pdf", sep="_"), height=hei1, width=wid1)
				print(obj)
				dev.off()
				
				obj1 <- taxa_boxplot_aggregate(prop, df, grp.name, strata, scale, xsize) 
				pdf(paste("Taxa_DifferentialAbundance_AbundanceBoxplot", scale, mt.method, cutoff, ann, ".pdf", sep="_"), height=hei1, width=wid1)
				print(obj1)
				dev.off()

			}

			# currently fold change
		    if (test == 'Para') {
				rownames(eff) <- rownames(eff.lc) <- rownames(eff.uc) <- taxa.names
				for (k in 1:ncol(eff)) {
					fold.dat.plot1 <- data.frame(Estimate=eff[, k], LCI=eff.lc[, k], UCI=eff.uc[, k], IV=taxa.names)
					obj <- plot_effect_size2(fold.dat.plot1)
					pdf(paste("Taxa_DifferentialAbundance_EffectBarplot", colnames(eff)[k], mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
					print(obj)
					dev.off()
				}
			} else {
				if (!is.na(eff[1])) {
					names(eff) <- taxa.names
					eff <- eff[!is.na(eff) & is.finite(eff)]
					eff <- sort(eff)
					taxa.names2 <- names(eff)
					if (is.null(wid2) | is.null(hei2)) {
						hei2 <- 4 + length(taxa.names2) / 20 * 3 
						wid2 <- 6
					}
					if (eff.type == 'LFC') {
						obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Log2 fold change')
						pdf(paste("Taxa_DifferentialAbundance_LFCBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
						print(obj)
						dev.off()	
					}

					if (eff.type == 'logP') {
						obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='-log10(P)')
						pdf(paste("Taxa_DifferentialAbundance_logPBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
						print(obj)
						dev.off()	
					}					
				}
			}

			# create heatmp
#			taxa.names2 <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
			
			if (length(taxa.names) >= 2) {
				generate_taxa_heatmap(data.obj, taxa.levels='All', sam.ord=order(grp), taxa=taxa.names, meta.info=grp.name, Colv=F, dendrogram='row', 
						ann=paste0(mt.method, '_', cutoff, '_', ann, '_Unclustered'))
				generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=grp.name,  ann=paste0(mt.method, '_', cutoff, '_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=grp.name,  data.type='R', ann=paste0(mt.method, '_', cutoff, '_Rank_', ann))
				try(
				generate_taxa_biplot(data.obj, taxa=taxa.names, trans='sqrt', grp.name, ann=paste0(mt.method, '_', cutoff, '_', ann), varname.size = 1.5)	
		        )
			}

		}
	}
}



create_lefse_format <- function(data.obj, diff.obj, grp.name, mt.method='fdr', cutoff=0.15, prev=0.1, minp=0.002, 
		lefse.dir="/data2/microbiome/jeff/tools/nsegata-lefse/", ann="") {
	
	
	# To be improved - currently no effect size
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	levels(grp) <- paste0(1:nlevels(grp), levels(grp))

	qv.list <- diff.obj$qv.list
	pv.list <- diff.obj$pv.list
	m.list <- diff.obj$m.list
	
	otu.name.12 <- data.obj$otu.name
	otu.tab.12 <- data.obj$otu.tab
	tax.family.a <- NULL
	for (i in 1:6) {
		tax.family <- apply(otu.name.12[, 1:i, drop=F], 1, paste, collapse=".")
		phlan.tab <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
		rownames(phlan.tab) <- phlan.tab[, 1]
		phlan.tab <- as.matrix(phlan.tab [, -1, drop=F])
		phlan.tab <- t(t(phlan.tab) / colSums(phlan.tab))	
		phlan.tab <- phlan.tab[rowSums(phlan.tab!=0) > prev*ncol(phlan.tab) & rowMaxs(phlan.tab) > minp, , drop=F]
		tax.family.a <- c(tax.family.a, rownames(phlan.tab))
		#tax.family.a <- c(tax.family.a, unique(tax.family))
	}
	alias.a <- sapply(strsplit(tax.family.a, "\\."), function(x) {
				if (length(x) >= 2) {
					if (length(x) == 2) {
						return(x[2])
					} else {
						return(paste0(x[2], ";", x[length(x)]))
					}
				} else {
					return(x)
				}
			})
	
	taxa.names <- NULL
	abundant.grp.names <- NULL
	
	for (LOI in setdiff(names(qv.list), 'Species')) {
		qv.vec <- qv.list[[LOI]]
		pv.vec <- pv.list[[LOI]]
		#qv.vec <- qvalue(pv.vec[, 1])$qvalues
		m.vec <- m.list[[LOI]]
		if (mt.method == 'fdr') {
			taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
			abundant.grp.name <- apply(m.vec, 1, function(x) levels(grp)[which.max(x)])[qv.vec <= cutoff]
		}

		if (mt.method == 'raw') {
			taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
			abundant.grp.name <- apply(m.vec, 1, function(x) levels(grp)[which.max(x)])[pv.vec <= cutoff]
		}
		
		if (length(taxa.name) != 0) {
			taxa.names <- c(taxa.names, taxa.name)
			abundant.grp.names <- c(abundant.grp.names, abundant.grp.name)
		}
	}
	
	# remove 'unclassified'
    abundant.grp.names <- abundant.grp.names[!grepl('unclassified', taxa.names, ignore.case=T)]
	taxa.names <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
	
	taxa.names <- intersect(taxa.names, alias.a)
	
	ind <- match(taxa.names, alias.a)
	phlan <- cbind(tax.family.a, '1.5', '\t', '\t', '-')
	phlan[ind, 2:5] <- cbind('3.5',  abundant.grp.names, '3.0', '-')
	write.table(phlan,  'Lefse.LDA.txt', row.names=F, col.names=F, quote=F, sep='\t')
#	cmd1 <- paste0('python ', lefse.dir,  'plot_cladogram.py Lefse.LDA.txt Taxa_DifferentialAbundance_', mt.method, cutoff, '_Cladogram_', ann, '.pdf', ' --format pdf')
#	cmd2 <- paste0("python ", lefse.dir, 'plot_res.py Lefse.LDA.txt Taxa_DifferentialAbundance_', mt.method, cutoff, '_LDA_', ann, '.pdf --format pdf')
#	system(cmd1)
#	system(cmd2)
}

perform_lefse_analysis <- function (data.obj,  grp.name, sub.grp.name=NULL, prev=0.1, minp=0.002,
		lefse.dir="/data2/microbiome/jeff/tools/nsegata-lefse/", ann="") {
	otu.name.12 <- data.obj$otu.name
	otu.tab.12 <- data.obj$otu.tab
	meta.dat <- data.obj$meta.dat
	
	levels(meta.dat[, grp.name]) <- paste0(1:nlevels(meta.dat[, grp.name]), levels(meta.dat[, grp.name]))
	
	tax.family <- apply(otu.name.12, 1, paste, collapse="|")
	lefse.tab <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
	rownames(lefse.tab) <- lefse.tab[, 1]
	lefse.tab <- as.matrix(lefse.tab [, -1])

	lefse.tab <- t(t(lefse.tab) / colSums(lefse.tab))
	
	lefse.tab <- lefse.tab[rowMaxs(lefse.tab) > minp & rowSums(lefse.tab!=0) > prev*ncol(lefse.tab), , drop=FALSE]

	if (is.null(sub.grp.name)) {
		header <- rbind(class=as.character(meta.dat[colnames(lefse.tab), grp.name]), 
				id=colnames(lefse.tab))
	} else {
		header <- rbind(class=as.character(meta.dat[colnames(lefse.tab), grp.name]), 
				subclass=as.character(meta.dat[colnames(lefse.tab), sub.grp.name]), 
				id=colnames(lefse.tab))
	}

	
	lefse.tab <- rbind(header, lefse.tab)
	write.table(lefse.tab, "lefse.txt", sep="\t", col.names=F, quote=F)
	
#	system(paste0('mkdir LefSe_', ann))
#	output <- paste0('LefSe_', ann, '/')
#	if (is.null(sub.grp.name)) {
#		cmd1 <- paste0("python ", lefse.dir, "format_input.py lefse.txt ", output, "temp.in ", "-c 1 -u 2 -o 1000000")
#	} else {
#		cmd1 <- paste0("python ", lefse.dir, "format_input.py lefse.txt ", output, "temp.in ", "-c 1 -s 2 -u 3 -o 1000000")
#	}
#
#	cmd2 <- paste0("python ", lefse.dir, "run_lefse.py ", output, "temp.in ",  output, "lda.res")
#	cmd3 <- paste0("python ", lefse.dir, "plot_res.py ",  output, "lda.res ", output, "lda.pdf --format pdf")
#	cmd4 <- paste0("python ", lefse.dir, "plot_cladogram.py ",  output, "lda.res ", output, "cladogram.pdf --format pdf --labeled_stop_lev 6 --abrv_stop_lev  6")
#	
#	system(cmd1)
#	system(cmd2)
#	system(cmd3)
#	system(cmd4)
}

plot.Boruta2 <- function (x, colCode = c("green", "yellow", "red", "blue"), sort = TRUE,
		whichShadow = c(TRUE, TRUE, TRUE), col = NULL, xlab = "Attributes", ids=NULL,
		ylab = "Importance", ...)
{
	if (class(x) != "Boruta")
		stop("This function needs Boruta object as an argument.")
	lz <- lapply(1:ncol(x$ImpHistory), function(i) x$ImpHistory[is.finite(x$ImpHistory[, i]), i])
	names(lz) <- colnames(x$ImpHistory)
	numShadow <- sum(whichShadow)
	lz <- lz[c(rep(TRUE, length(x$finalDecision)), whichShadow)]
	col <- Boruta:::generateCol(x, colCode, col, numShadow)
	if (sort) {
		ii <- order(sapply(lz, median))
		lz <- lz[ii]
		col <- col[ii]
	}
	names(lz) <- gsub("^X",  "", names(lz))
	if (is.null(ids)) {
		len <- sum(x$finalDecision %in% c('Confirmed', 'Tentative'))
		ind <- (length(lz) - len) : length(lz)
	} else {
		ind <- match(ids, names(lz))
	}
	
	boxplot(lz[ind], xlab = xlab, ylab = ylab, col = col[ind], ...)
	invisible(x)
	names(lz[ind])
}

createROC <- function (pv.list, lab.list, pos.lab='1', file.name='ROC.pdf') {
	require(ROCR)
	n <- length(pv.list)
	aucs <- numeric(n)
	names(aucs) <- names(pv.list)

	cols <- scales::hue_pal()(n)
	ltys <- rep(c(1, 2), ceiling(n/2))[1:n]
	pdf(file.name)
	for (i in 1:n) {
		
  		    cat("*")
			pv.mat <- pv.list[[i]]
			lab.mat <- lab.list[[i]]
			
			pred <- prediction(pv.mat, lab.mat==pos.lab)
			perf <- performance(pred, "tpr", "fpr")
			aucs[i] <- mean(unlist(performance(pred, 'auc')@y.values))
			plot(perf, avg="threshold", col=cols[i], lty=ltys[i], lwd=2,  add=ifelse(i==1, FALSE, TRUE),  main='ROC curve')
				
		}

	legend("right", legend=paste0(names(pv.list), "(AUC:", round(aucs, 3), ")"), col=cols, lty=ltys, lwd=2,  bty="n")
	dev.off()
}

#invlogit <- function(x) {
#	rbinom(length(x), 1, exp(x) / (1 + exp(x)))
#}
#pv.list <- list(x1=rnorm(100), x2=rnorm(100))
#lab.list <- list(x1=invlogit(pv.list[['x1']]), x2=invlogit(pv.list[['x2']]))
#createROC(pv.list, lab.list)

predictionRF <- function (data.obj,  resp.name, formula=NULL, taxa.level='Species', binary=FALSE, prev=0.1, minp=0.002, B=50, seed=123, 
		boruta.level=c('Confirmed', 'Tentative'), ann='',...) {

	#sink(paste0("Taxa_RandomForest_", taxa.level, ".txt"))
	date()
	response <- data.obj$meta.dat[, resp.name]
	
	if (!is.null(formula)) {
		adj.var <- unlist(strsplit(unlist(strsplit(formula, '\\s*~\\s*'))[2], '\\s*\\+\\s*'))
		adj <- data.obj$meta.dat[, adj.var, drop=F]
	} 

	if (taxa.level == 'Species') {
		if (taxa.level %in% names(data.obj$abund.list)) {
			ct <- data.obj$abund.list[[taxa.level]]
		} else {
			# Accomodate different version
			ct <- data.obj$otu.tab
			rownames(ct) <- paste0("OTU", rownames(ct), ":", data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
			data.obj$abund.list[['Species']] <- ct
		}
	} else {
		ct <- data.obj$abund.list[[taxa.level]]
	}

	prop <- t(t(ct) / colSums(ct))
	prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
	prop <- t(prop)
	
	if (binary == TRUE) {
		prop <- (prop != 0)
	}
	
    original.names <- colnames(prop)
	set.seed(seed)
	
	if (is.null(formula)) {
		performance <- matrix(0, nrow=B, ncol=2)
		colnames(performance) <- c("RF_M","Guess")
		roc.list <- list(RF_M=NULL)
		lab.list <- roc.list
	} else {
		performance <- matrix(0, nrow=B, ncol=4)
		colnames(performance) <- c("RF_M", "RF_CF", "RF_M+CF", "Guess")
		roc.list <- list("RF_M"=NULL, "RF_CF"=NULL, "RF_M+CF"=NULL)
		lab.list <- roc.list
	}

	
	colnames(prop) <- gsub(";", "_", colnames(prop))
	colnames(prop) <- gsub(":", "_", colnames(prop))
	colnames(prop) <- gsub("-", "_", colnames(prop))
	colnames(prop) <- gsub("\\.", "_", colnames(prop))
	names(original.names) <- colnames(prop)
	if (!is.null(formula)) {
		names(adj.var) <- adj.var
		original.names <- c(original.names, adj.var)
	}
	
	
	if (!is.null(formula)) {
		padj <- cbind(prop, adj)
	} else {
		padj <- prop
	}
	
	I <- nrow(prop)
	cat('Begin to bootstrap ...\n')
	for(b in 1:B){
		if (b %% 50 == 0) cat (".")
		err <- try({
			bsample <- sample(1:I, I, replace=T)
			if (is.factor(response)) {
				rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
				performance[b, "RF_M"] <- mean(response[-bsample]!= rf1$test$predicted)
				performance[b, "Guess"]<- mean(response[-bsample] != levels(response)[which.max(tabulate(response[bsample]))])
				roc.list[['RF_M']] <- c(roc.list[['RF_M']], rf1$test$vote[, 1])
				lab.list[['RF_M']] <- c(lab.list[['RF_M']], response[-bsample])
				if (!is.null(formula)) {
					rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
					rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
					performance[b, "RF_CF"] <- mean(response[-bsample]!= rf2$test$predicted)
					performance[b, "RF_M+CF"] <- mean(response[-bsample]!= rf3$test$predicted)
					roc.list[['RF_CF']] <- c(roc.list[['RF_CF']], rf2$test$vote[, 1])
					roc.list[['RF_M+CF']] <- c(roc.list[['RF_M+CF']], rf3$test$vote[, 1])
					lab.list[['RF_CF']] <- c(lab.list[['RF_CF']], response[-bsample])
					lab.list[['RF_M+CF']] <- c(lab.list[['RF_M+CF']], response[-bsample])
				}
				
			} else {
				rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], ...)
				performance[b, "RF_M"] <- mean((response[-bsample] - rf1$test$predicted)^2)
				performance[b, "Guess"]<- mean((response[-bsample] - mean(response[bsample]))^2)
				
				if (!is.null(formula)) {
					
					rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample],  ...)
					rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample],  ...)
					performance[b, "RF_CF"] <- mean((response[-bsample] - rf2$test$predicted)^2)
					performance[b, "RF_M+CF"] <- mean((response[-bsample]- rf3$test$predicted)^2)
				}
			}
		})
		if (inherits(err, "try-error")) {
			next
		}

	}
	cat("\n")
	if (is.null(formula)) {
		cat("Fridman.test p value: ", friedman.test(performance)$p.value, "\n")
	} else {
		cat("Fridman.test p value (M+CF vs CF): ", friedman.test(performance[, c('RF_M+CF', 'RF_CF')])$p.value, "\n")
	}
	
	pdf(paste0("Taxa_Random_forest_misclassification_barplot_", taxa.level, '_', ann, ".pdf"), height=5, width=4)
	if (!is.null(formula)) {
		performance2 <- performance[, c('Guess', 'RF_CF', 'RF_M', 'RF_M+CF')]
	} else {
		performance2 <- performance
	}
	if (is.factor(response)) {
		boxplot(performance2, col="#4DAF4A", ylab='Classification error', las=2)
	} else {
		boxplot(performance2, col="#4DAF4A", ylab='PMSE', las=2)
	}	
	dev.off()
	
	# ROC curve - only compare to the first level of the factor 
    if (is.factor(response)) {
		lab.list <- lapply(lab.list, function(x) {
					y <- factor(x)
					levels(y) <- levels(response)
					y
				}
		)
		createROC(roc.list, lab.list, pos.lab=levels(response)[1], file.name=paste0("Taxa_Random_forest_ROC_", taxa.level, '_', ann, ".pdf"))
	}
	
	if (is.factor(response)) {
		rf <- randomForest(x=padj, y=response, sampsize=table(response), importance=T, ...)
	} else {
		rf <- randomForest(x=padj, y=response,  importance=T, ...)
	}
	
	importance <- as.data.frame(rf$importance)
	
	if (is.factor(response)) {
		write.csv(importance[rev(order(importance$MeanDecreaseAccuracy)), ],
				paste0("Taxa_RandomForest_Ranking_MeanDecreaseAccuracy_", taxa.level,'_', ann, ".csv"))
		write.csv(importance[rev(order(importance$MeanDecreaseGini)), ], 
				paste0("Taxa_RandomForest_Ranking_MeanDecreaseGini_", taxa.level, '_', ann, ".csv"))
	} else {
		write.csv(importance[rev(order(importance[, '%IncMSE'])), ],
				paste0("Taxa_RandomForest_Ranking_IncMSE_", taxa.level, '_', ann, ".csv"))
		write.csv(importance[rev(order(importance$IncNodePurity)), ], 
				paste0("Taxa_RandomForest_Ranking_IncNodePurity_", taxa.level, '_', ann, ".csv"))
	}

	
# Boruta Feature Selection - All subset selection (can't solve the confounding problem)
	obj.Boruta <- Boruta(padj, response, doTrace = 2)	
	write.csv(obj.Boruta$finalDecision, paste0("Taxa_Random_forest_Boruta_Feature_Selection_", taxa.level, '_', ann, ".csv"))
	
	pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_", taxa.level, '_', ann,  ".pdf"), height=6, width=10)
	par(mar=par('mar') + c(3, 0, 0, 0))
	plot(obj.Boruta, main = "Feature selection by Boruta", ylab="Importance z-score", lwd = 0.5, las = 3, xlab = "",
			cex=1 / (ncol(prop)/50), cex.axis=0.25*200/ncol(prop), yaxt='n')
	axis(2, cex.axis=1)
	dev.off()
	#sink()
	
	pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_Significant_Only_", taxa.level, '_', ann, ".pdf"), height=6, width=6)
	par(mar=par('mar') + c(3, 0, 0, 0))
	otu.ids <- plot.Boruta2(obj.Boruta, main = "Feature selection by Boruta", lwd = 0.5, las = 3, 
			ylab="Importance z-score", xlab = "", cex.axis=1, yaxt='n')
	axis(2, cex.axis=1)
	dev.off()
	
	if ('Confirmed' %in% boruta.level) {
		taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed')]] 
		if (!is.null(formula)) {
			taxa.names <- setdiff(taxa.names, adj.var)
		}
		
		if (length(taxa.names) > 1) {
			generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12,
					margins=c(5, 15), ann=paste0('BorutaFeatures_P_Confirmed_', ann))
			generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='R',
					margins=c(5, 15), ann=paste0('BorutaFeatures_R_Comfirmed_', ann))
			generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, Colv=F, dendrogram='row', sam.ord=order(response),
					margins=c(5, 15), ann=paste0('BorutaFeatures_P_Confirmed_Unclusterded_', ann))
			generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='R', Colv=F, dendrogram='row', sam.ord=order(response),
					margins=c(5, 15), ann=paste0('BorutaFeatures_R_Comfirmed_Unclusterded_', ann))
			
			if (is.factor(response)) {
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				build.decision.tree(data.obj,  resp.name=resp.name, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann)) 
			} else {
				data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
				levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
			}
		}
	}

	if ('Tentative' %in% boruta.level) {
		taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed',  'Tentative')]] 
		if (!is.null(formula)) {
			taxa.names <- setdiff(taxa.names, adj.var)
		}
		
		if (length(taxa.names) > 1) {
			generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12,
					margins=c(5, 15), ann=paste0('BorutaFeatures_P_Tentative_', ann))
			generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='R',
					margins=c(5, 15), ann=paste0('BorutaFeatures_R_Tentative_', ann))
			
			generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, Colv=F, dendrogram='row', sam.ord=order(response),
					margins=c(5, 15), ann=paste0('BorutaFeatures_P_Tentative_Unclustered_', ann))
			generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='R', Colv=F, dendrogram='row', sam.ord=order(response),
					margins=c(5, 15), ann=paste0('BorutaFeatures_R_Tentative_Unclustered_', ann))
			
			if (is.factor(response)) {
				
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				build.decision.tree(data.obj,  resp.name=resp.name, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann)) 
			} else {
				data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
				levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
			}
		}
	}	
}

col.func <- colorRampPalette(c("blue", "cyan", "red", "yellow"))

generate_taxa_heatmap <- function (data.obj, taxa.levels='Genus', taxa='All', meta.info, sam.ord=NULL, data.type='P',  prev=0.1, minp=0.002, 
		row.col.dat='Phyla', phy.no=4, sepwidth=0.01, colsep=NULL, rowsep=NULL,
		white='white', colFunc=jet, Rowv=T, Colv=T, dendrogram='both', margins=c(5, 15), in.grid=F, sepcolor='black', is.labCol=T, cexCol=1, cexRow=NULL,
		omas=c(1, 1, 1, 8), width=12, height=6, ann='All') {
	
	df <- data.obj$meta.dat

	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		
		if (LOI == 'All') {
			if (taxa[1] == 'All') {
				stop("Please specify the taxa names that will be included in the heatmap!\n")
			} 
			prop <- NULL
			for (LOI2 in names(data.obj$abund.list)) {
				ct <- data.obj$abund.list[[LOI2]]
				prop0 <- t(t(ct) / colSums(ct))
				prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])	
	
			}
			colnames(prop) <- colnames(prop0)
			if (nrow(prop) != length(taxa)) {
				warnings('Some taxa not found in abundance lists! Please check the names!\n')
			}
			
  		} else {
			ct <- data.obj$abund.list[[LOI]]
			prop <- t(t(ct) / colSums(ct))
			
			if (taxa != 'All') {
				prop <- prop[taxa, , drop=FALSE]
			} else {
				prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
			}			
		}

		# Sort row and column
		if (is.null(sam.ord)) {
			prop <- prop[order(rownames(prop)), , drop=FALSE]
		} else {
			prop <- prop[order(rownames(prop)), sam.ord, drop=FALSE]
		}
		
		if (is.labCol) {
			labCol <- colnames(prop)
		}  else {
			labCol <- ''
		}
		# Deal with zeros
		if (data.type == 'B') {
			col.scheme <- c("lightyellow", "red")
			prop[, ] <- as.numeric(prop != 0)
			breaks <- c(-0.01, 0.01, 1.1)
		} 
		if (data.type == 'P'){
			#col.scheme <- c(white, colFunc(50))		
			col.scheme = c(white, brewer.pal(11, "Spectral"))
			ind.temp <- prop != 0
			minp <- min(prop[prop!=0])/1.1
			prop[prop==0] <- minp
			prop <- log10(prop)
		#	breaks <- c(log10(minp)-0.01, seq(log10(minp)+0.01, 0, len=51))	
	        breaks <- c(log10(minp)-0.01, seq(log10(minp)+0.01, 0, len=12))
		}
		if (data.type == 'R'){
			col.scheme <- c('white', colorRampPalette(c("green", "black", "red"))(ncol(prop)-1))
			prop <- t(apply(prop, 1, function(x) {
								temp <- rank(x[x!=0])
								s <- (ncol(prop) - 1) / (max(temp) - min(temp))
								temp <- 1 + (temp - min(temp)) * s
								x[x!=0] <- temp
								x
								}))
			breaks <- seq(0, ncol(prop), len=ncol(prop)+1)
		}
		
		if (nrow(prop) > 150) {
			labRow <- ''
			cexRow <- 1
		} else {
			labRow <- rownames(prop)
			if (is.null(cexRow)) {
				cexRow <- ifelse(0.5 * 60 / nrow(prop) > 1, 1, 0.5 * 60 / nrow(prop))
			}

		}
		
		if (is.null(colsep)) {
			if (in.grid == T) {
				colsep <- c(0:ncol(prop))  
			} else {
				colsep <- c(0, ncol(prop))  
			}
		}
		
		if (is.null(rowsep)) {
			if (in.grid == T) {
				rowsep <- c(0:nrow(prop))
			} else {
				rowsep <- c(0, nrow(prop))
			}
		}
		
		# Determine the col/rowside colors
		rainbow3 <- function(x) {
			rainbow(x + 2)[1:x]
		}
		jet3 <- function(x) {
			jet(x + 2)[1:x]
		}
		colFnsC <- c(colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue')), colorRampPalette(c('black', 'red')))
		colFnsF <- c(rainbow3, jet3)
		key.list <- list()
		colsidecol <- NULL
		for (keyID in meta.info) {
			if (is.null(sam.ord)) {
				x <- df[, keyID]
			} else {
				x <- df[sam.ord, keyID]
			}		
			i <- 0
			j <- 0
			if (is.factor(x)) {
				key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)	
				i <- (i + 1) %% length(colFnsF)
				colsidecol <- cbind(colsidecol, key.list[[keyID]]$colors[x])
			} else {
				key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
				j <- (j + 1) %% length(colFnsC)
				colsidecol <- cbind(colsidecol, cmap(x, key.list[[keyID]]))
			}
		}
		colnames(colsidecol) <- meta.info
		
		# add aunbdance key 
		if (data.type == 'B') {
			prop.cmap <- list(breaks=c("Absence", "Presence"), colors=c("lightyellow", "red"), base=NA, col.na=NA, right=F, include.lowest=F) 
			KeyName <- 'Abundance'
		} 
#		if (data.type == 'P'){
#			prop.cmap <- makecmap(prop[ind.temp], n = 5, colFn = colFunc)
#		}
#		if (data.type == 'R'){
#			prop.cmap <- makecmap(as.vector(prop), n = 5, colFn = colFunc)
#		}
		
		if (data.type == 'P'){
			KeyName <- 'log(Proportion)'	
		}
		if (data.type == 'R'){
			KeyName <- 'Rank'
		}
			
		# add row col key
        if (row.col.dat == 'Phyla') {
			if (LOI %in% c('Class', 'Order', 'Family', 'Genus', 'Species')) {
				if (LOI == 'Species') {
					phy <- sapply(strsplit(rownames(prop), ":"), function(x) x[2])
					phy <- sapply(strsplit(phy, ";"), function(x) x[1])
				} else {
					phy <- sapply(strsplit(rownames(prop), ";"), function(x) x[1])
				}
				temp <- sort(table(phy), decr=T)
				if (length(temp) > phy.no) {
					rare.phy <- names(temp)[-(1:phy.no)]
					phy[phy %in% rare.phy] <- 'Other'
					phy <- factor(phy, levels=c(names(temp)[(1:phy.no)], 'Other'))
				} else {
					phy <- factor(phy)
				}
				rowsidecol<- rainbow(nlevels(phy))[phy]
				rowsidecol <- rbind(rowsidecol, rowsidecol)
				rownames(rowsidecol) <- c('', '')
				phy.cmap <- list(breaks=levels(phy), colors=rainbow(nlevels(phy)), base=NA, col.na=NA, right=F, include.lowest=F) 
			} else {
				rowsidecol <- NULL
			}
		} else {
			rowsidecol <- NULL
		}
		
		pdf(paste0('Taxa_Heatmap_', LOI, '_', ann, '.pdf'), width=width, height=height)
		par(oma = omas)

		if (data.type == 'R' | data.type == 'B') {
			dist2 <- dist
		}
		if (data.type == 'P') {
			# Better clustering of taxa
			dist2 <- function(x) as.dist((1-cor(t(x)))/2)
		}

		# Pearson correlation distance
		heatmap.3(prop, 
				Rowv=Rowv, 
				Colv=Colv, 
				distfun = dist2,
				dendrogram=dendrogram,
				scale='none',
				col=col.scheme, 
				breaks=breaks, 
				symbreaks=F,
				trace='none',
				margins= margins, 
				colsep = colsep,  
				rowsep = rowsep,  
				sepcolor= sepcolor, 
				sepwidth=c(sepwidth, sepwidth),
				ColSideColors=colsidecol,
				RowSideColors=rowsidecol,
				cexRow=cexRow,
				labRow=labRow,
				labCol=labCol,
				cexCol=cexCol,
				key=(data.type != 'B'), density.info='none', symkey=F, KeyValueName=KeyName,	
				NumColSideColors= 0.5 *length(meta.info),
				NumRowSideColors= 0.5
		)

		par(cex=0.75)
		par(oma=c(0, 0, 1, 0)) 
		
		if (!is.null(rowsidecol) & row.col.dat == 'Phyla') {
			y.cord <- (1/(length(meta.info)+2)) * (0:((length(meta.info)+2) - 1))
			vkey2(phy.cmap, 'Phylum', x=0, y=-0.2, stretch=1.2)

		}

		y.cord <- (1/(length(meta.info))) * (0:(length(meta.info) - 1))
		k <- 1
		for (keyID in meta.info) {
			x <- df[, keyID]
			if (is.factor(x)) {
				vkey2(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
			} else {
				vkey(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
			}
			k <- k + 1
		}	
		
		if (data.type == 'B') {
			vkey2(prop.cmap, KeyName, y=1, stretch=1.2)		
		} 
#		if (data.type == 'P'){
#			vkey(prop.cmap, 'log10(Proportion)', y=1, stretch=1.2)
#		}	
#		if (data.type == 'R'){
#			vkey(prop.cmap, 'Rank', y=1, stretch=1.2)
#		}	
		dev.off()	
	}
}

generate_stacked_barplot <- function(data.obj, grp.name=NULL, taxa.levels=c('Phylum', 'Family', 'Genus'), agg.cutoff=0.005, 
		hei1=6, wid1=9, hei2=6, wid2=9, margin=10, ann='') {
	if (is.null(grp.name)) {
		grp <- 1:nrow(data.obj$meta.dat)
	} else {
		grp <- data.obj$meta.dat[, grp.name]
	}
	

	pdf(paste0("Taxa_Stacked_Barplot_Overall_Compo_", ann, ".pdf"), height=hei1, width=wid1)
	
	par(mar=par('mar') + c(0, margin, 0, 0))
	
	name.list <- list()
	col.list <- list()
	for (taxa.level in taxa.levels) {
		abund0 <- data.obj$abund.list[[taxa.level]]
		abund0 <- t(t(abund0) / colSums(abund0))
		
		abund1 <- abund0[rowMeans(abund0) >= agg.cutoff, , drop=F]
		abund2 <- abund0[rowMeans(abund0) < agg.cutoff, , drop=F]
		
		prop <- rbind(abund1, Other=colSums(abund2))
		colnames(prop) <- colnames(abund0)
		name.list[[taxa.level]] <- rownames(abund1)
		
		#rand.col <- c(sample(rainbow(nrow(prop)*2), nrow(prop)-1), 'gray')
		rand.col <- c(rep_len(brewer.pal(12, "Paired"), nrow(prop)-1), 'gray')
		col.list[[taxa.level]] <- rand.col

		cex.legend = ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
		prop <- prop[, order(grp)]
		barplot(prop, col=rand.col, ylab='Proportion', las=2, legend.text=rownames(prop), cex.names=0.5,
				args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), main=taxa.level)
		
	}
	dev.off()
	
	# Generate averaged over stack barplot
	cex.legend <- 10
	if (!is.null(grp.name)) {
		pdf(paste0("Taxa_Stacked_Barplot_Grouped_Compo_", ann, ".pdf"), height=hei2, width=wid2)
		
		par(mar=par('mar') + c(0, margin/length(taxa.levels), 0, 0))
		par(mfrow=c(1, length(taxa.levels)))
		for (taxa.level in taxa.levels) {
			abund0 <- data.obj$abund.list[[taxa.level]]
			abund0 <- t(t(abund0) / colSums(abund0))
			
			abund0 <- t(apply(abund0, 1, function(x) {
								tapply(x, grp, mean)
							}))
			
			abund1 <- abund0[name.list[[taxa.level]], ]
			abund2 <- 1 - colSums(abund1)
			
			prop <- rbind(abund1, Other=abund2)
			colnames(prop) <- colnames(abund0)
			
			newsize <- ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
			cex.legend <- ifelse(cex.legend < newsize, cex.legend, newsize) 
			#prop <- prop[, order(grp)]
			barplot(prop, col=col.list[[taxa.level]], ylab='Proportion', las=2, cex.names=1, main=taxa.level)
			#		legend.text=rownames(prop), args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-2.2, 0)))
			
		}
		
		par(mar=c(0, 0, 0, 0))
		par(oma=c(0, 0, 0, 0))
		for (taxa.level in taxa.levels) {
			abund0 <- data.obj$abund.list[[taxa.level]]
			abund0 <- t(t(abund0) / colSums(abund0))
			
			abund0 <- t(apply(abund0, 1, function(x) {
								tapply(x, grp, mean)
							}))
			
			abund1 <- abund0[name.list[[taxa.level]], ]
			abund2 <- 1 - colSums(abund1)
			
			prop <- rbind(abund1, Other=abund2)
			colnames(prop) <- colnames(abund0)
			
		#	cex.legend = ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
			#prop <- prop[, order(grp)]
	        plot(1, type="n", axes=FALSE, xlab="", ylab="")
			legend('left', legend=rownames(prop), bty='n', fill=col.list[[taxa.level]], cex=cex.legend)
			
		}
		dev.off()
	}

}

build.decision.tree <- function(data.obj,  resp.name, taxa.level='Species', binary=FALSE, taxa, ann='All') {
	ann <- paste(taxa.level, ann, sep="_")
	response <- data.obj$meta.dat[, resp.name]
	
	ct <- data.obj$abund.list[[taxa.level]]
	prop <- t(t(ct) / colSums(ct))
	prop <- prop[taxa, , drop=F] 	
	if (binary == TRUE) {
		prop <- (prop != 0)
	}
	
	dat <- as.data.frame(t(prop))
	dat <- data.frame(dat, response)
	try(
		if (is.factor(response)) {
			fit <- rpart(response ~ ., method="class", data=dat)
			post(fit, file = paste0("Taxa_Unpruned_Classification_tree_", ann, ".ps"), title = "Unpruned Classification Tree")
			pfit<- prune(fit, cp= fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
			post(pfit, file = paste0("Taxa_Pruned_Classification_", ann, ".ps"), title = "Pruned Classification Tree")
			
		} else {
			fit <- rpart(response ~ ., method="anova", data=dat)
			post(fit, file = paste0("Taxa_Unpruned_Regression_tree_", ann, ".ps"), title = "Unpruned Regression Tree")
			pfit<- prune(fit, cp= fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
			post(pfit, file = paste0("Taxa_Pruned_Regression_tree_", ann, ".ps"), title = "Pruned RegressionTree")
		}
	)	
}

bootstrap.pwr <- function(pcs, formula, dat, ns=NULL, perm.no=199, iter.no=100) {
	if (is.null(ns)) {
		ns <- nrow(dat)*seq(1, 4, len=8)
	}
	pvs <- numeric(length(ns))
	names(pvs) <- paste(ns)
	for (n in ns) {
		cat('.')
		temp <- sapply(1:iter.no, function(i) {
			bt.ind <- sample(1:nrow(dat), n, repl=T)
			dat.bt <- dat[bt.ind, ]
			dist.bt <- dist(pcs[bt.ind, ])
			aov.tab <- adonis(as.formula(paste('dist.bt', formula)), dat=dat.bt, permutations=perm.no)$aov.tab
	        pv <- aov.tab[nrow(aov.tab)-2, ncol(aov.tab)]
		})
        pvs[paste(n)] <- mean(temp <= 0.05)
	}

	return(pvs)
}


perform_power_analysis <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), formula=NULL,  
		grp.name=NULL, adj.name=NULL, ann='', ...) {
	if (is.null(formula)) {
		if (is.null(adj.name)) {
			formula <- paste('~', grp.name)
		} else {
			formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
		}
	}	 
	df <- data.obj$meta.dat
    pvm <- NULL
	for (dist.name in dist.names) {
		cat('*')
		dist.mat <- dist.obj[[dist.name]]
		pcs <- cmdscale(dist.mat, k=nrow(dist.mat)-1)
		pvs <- bootstrap.pwr(pcs, formula, df, ...)
		pvm <- rbind(pvm, pvs)
	}
	colnames(pvm) <- names(pvs)
	rownames(pvm) <- dist.names
	
	pvdf <- melt(pvm)
	colnames(pvdf) <- c('Distance_type', 'Sample_size', 'Value')
	pdf(paste0("Power_curve_bootstrap_", ann, '.pdf')) 
	g.obj <- ggplot(pvdf, aes(x=Sample_size, y=Value)) +
				geom_point() +
				geom_line() +
				ylab('Power') +
				facet_wrap(~ Distance_type, ncol=2) +
				theme_bw()
	print(g.obj)
	dev.off()
	
	return(pvdf)
}


