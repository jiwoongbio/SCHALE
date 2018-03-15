schale.available.samples <- function(table, suffix.pattern) {
	numeric.colnames <- colnames(table[, unlist(lapply(table, is.numeric))])
	available.samples <- sub(suffix.pattern[1], "", numeric.colnames[grep(suffix.pattern[1], numeric.colnames)])
	if(length(suffix.pattern) >= 2) {
		for(index in 2:length(suffix.pattern)) {
			available.samples <- intersect(available.samples, sub(suffix.pattern[index], "", numeric.colnames[grep(suffix.pattern[index], numeric.colnames)]))
		}
	}
	available.samples
}

schale.variant.table <- function(variant.table, normal.sample, tumor.samples = NULL) {
	if(is.matrix(variant.table)) variant.table <- data.frame(variant.table)
	if(!is.data.frame(variant.table)) stop("variant.table must be data.frame or matrix.")

	available.samples <- schale.available.samples(variant.table, c("\\.refDepth$", "\\.altDepth$"))
	if(length(available.samples) == 0) stop("variant.table must have [sample].refDepth and [sample].altDepth columns.")
	if(!(all(normal.sample %in% available.samples) && length(normal.sample) == 1)) stop("normal.sample is not appropriate.")
	if(is.null(tumor.samples)) tumor.samples <- available.samples
	tumor.samples <- setdiff(tumor.samples, normal.sample)
	if(!(all(tumor.samples %in% available.samples) && length(tumor.samples) >= 1)) stop("tumor.samples is not appropriate.")

	samples <- c(normal.sample, tumor.samples)
	for(sample in samples) {
		variant.table[, paste(sample, "depth", sep = ".")] <- rowSums(variant.table[, paste(sample, c("refDepth", "altDepth"), sep = ".")])
		variant.table[, paste(sample, "BAF", sep = ".")] <- variant.table[, paste(sample, "altDepth", sep = ".")] / variant.table[, paste(sample, "depth", sep = ".")]
	}
	variant.table[, paste(samples, "depth", sep = ".")] <- edgeR::equalizeLibSizes(edgeR::DGEList(counts = variant.table[, paste(samples, "depth", sep = ".")]))$pseudo.counts
	for(tumorSample in tumor.samples) {
		variant.table[, paste(tumorSample, "LOH", sep = ".")] <- abs(variant.table[, paste(tumorSample, "BAF", sep = ".")] - 0.5) - abs(variant.table[, paste(normal.sample, "BAF", sep = ".")] - 0.5)
		variant.table[, paste(tumorSample, "CNA", sep = ".")] <- log(variant.table[, paste(tumorSample, "depth", sep = ".")] / variant.table[, paste(normal.sample, "depth", sep = ".")], base = 2)
	}
	variant.table
}

schale.segment.table <- function(variant.table, tumor.samples = NULL, chromosomes = NULL, centromere.table = NULL, alpha = 0.01, undo.alpha = 0.01) {
	if(is.matrix(variant.table)) variant.table <- data.frame(variant.table)
	if(!is.data.frame(variant.table)) stop("variant.table must be data.frame or matrix.")
	if(!all(c("chromosome", "position") %in% colnames(variant.table))) stop("variant.table must have chromosome and position columns.")

	available.samples <- schale.available.samples(variant.table, c("\\.LOH$", "\\.CNA$"))
	if(length(available.samples) == 0) stop("variant.table must have [sample].LOH and [sample].CNA columns.")
	if(is.null(tumor.samples)) tumor.samples <- available.samples
	if(!(all(tumor.samples %in% available.samples) && length(tumor.samples) >= 1)) stop("tumor.samples is not appropriate.")

	if(is.null(chromosomes)) {
		chromosomes <- as.character(unique(variant.table$chromosome))
	}
	if(!is.null(centromere.table)) {
		if(!all(c("chromosome", "start", "end") %in% colnames(centromere.table))) stop("centromere.table must have chromosome, start and end columns.")
	}
	segment.out <- data.frame()
	for(tumorSample in tumor.samples) {
		LOH.object <- DNAcopy::CNA(genomdat = variant.table[, paste(tumorSample, "LOH", sep = ".")], chrom = variant.table$chromosome, maploc = variant.table$position, presorted = FALSE, sampleid = paste(tumorSample, "LOH", sep = "."), data.type = "binary")
		CNA.object <- DNAcopy::CNA(genomdat = variant.table[, paste(tumorSample, "CNA", sep = ".")], chrom = variant.table$chromosome, maploc = variant.table$position, presorted = FALSE, sampleid = paste(tumorSample, "CNA", sep = "."), data.type = "logratio")
		CNA.object <- DNAcopy::smooth.CNA(CNA.object)
		segment.out <- rbind(segment.out, DNAcopy::segment(LOH.object, alpha = alpha, min.width = 2, verbose = 0)$out)
		segment.out <- rbind(segment.out, DNAcopy::segment(CNA.object, alpha = alpha, min.width = 2, verbose = 0)$out)
	}
	schale.segment.p.value <- function(head.index, tail.index) {
		head.segment.table <- chromosome.segment.table[head.index, ]
		tail.segment.table <- chromosome.segment.table[tail.index, ]
		head.variant.table <- chromosome.variant.table[head.segment.table$start <= chromosome.variant.table$position & chromosome.variant.table$position <= head.segment.table$end, ]
		tail.variant.table <- chromosome.variant.table[tail.segment.table$start <= chromosome.variant.table$position & chromosome.variant.table$position <= tail.segment.table$end, ]
		columns <- paste(rep(tumor.samples, each = 2), rep(c("LOH", "CNA"), times = length(tumor.samples)), sep = ".")
		standardDeviations <- unlist(lapply(columns, function(column) {
			sqrt(sum(c(head.variant.table[, column] - head.segment.table[, paste(column, "mean", sep = ".")], tail.variant.table[, column] - tail.segment.table[, paste(column, "mean", sep = ".")]) ^ 2) / (head.segment.table$number.of.variants + tail.segment.table$number.of.variants))
		}))
		if(head.segment.table$number.of.variants + tail.segment.table$number.of.variants >= length(columns) && sum(standardDeviations == 0) == 0) {
			centerValues <- as.numeric(head.segment.table[, paste(columns, "mean", sep = ".")] + tail.segment.table[, paste(columns, "mean", sep = ".")]) / 2
			values <- rbind(data.frame(group = "head", t((t(head.variant.table[, columns]) - centerValues) / standardDeviations), check.names = FALSE), data.frame(group = "tail", t((t(tail.variant.table[, columns]) - centerValues) / standardDeviations), check.names = FALSE))
			values <- cbind(values, princomp(values[, columns])$scores)
			t.test(x = values$Comp.1[values$group == "head"], y = values$Comp.1[values$group == "tail"], var.equal = TRUE)$p.value
		} else {
			0
		}
	}
	segment.table <- data.frame()
	for(chromosome in chromosomes) {
		chromosome.variant.table <- variant.table[variant.table$chromosome == chromosome, ]
		if(dim(chromosome.variant.table)[1] > 0) {
			chromosome.segment.table <- data.frame(chromosome = chromosome, start = sort(unique(segment.out[segment.out$chrom == chromosome, "loc.start"])), end = sort(unique(segment.out[segment.out$chrom == chromosome, "loc.end"])))
			number.of.segments <- dim(chromosome.segment.table)[1]
			if(!is.null(centromere.table)) {
				chromosome.centromere.table <- centromere.table[centromere.table$chromosome == chromosome, ]
				index <- which(chromosome.segment.table$start < chromosome.centromere.table$start & chromosome.centromere.table$end < chromosome.segment.table$end)
				if(length(index) == 1) {
					# p
					chromosome.segment.table.p <- data.frame(chromosome = chromosome, start = chromosome.segment.table$start[index], end = max(chromosome.variant.table$position[chromosome.variant.table$position < chromosome.centromere.table$start]))
					if(index > 1) {
						chromosome.segment.table.p <- rbind(chromosome.segment.table[1:(index - 1), ], chromosome.segment.table.p)
					}
					# q
					chromosome.segment.table.q <- data.frame(chromosome = chromosome, start = min(chromosome.variant.table$position[chromosome.centromere.table$end < chromosome.variant.table$position]), end = chromosome.segment.table$end[index])
					if(index < number.of.segments) {
						chromosome.segment.table.q <- rbind(chromosome.segment.table.q, chromosome.segment.table[(index + 1):number.of.segments, ])
					}
					# p + q
					chromosome.segment.table <- rbind(chromosome.segment.table.p, chromosome.segment.table.q)
					number.of.segments <- dim(chromosome.segment.table)[1]
				}
			}
			for(index in 1:number.of.segments) {
				selection <- chromosome.segment.table$start[index] <= chromosome.variant.table$position & chromosome.variant.table$position <= chromosome.segment.table$end[index]
				chromosome.segment.table$number.of.variants[index] <- sum(selection)
				for(tumorSample in tumor.samples) {
					chromosome.segment.table[index, paste(tumorSample, "LOH.mean", sep = ".")] <- mean(chromosome.variant.table[selection, paste(tumorSample, "LOH", sep = ".")])
					chromosome.segment.table[index, paste(tumorSample, "CNA.mean", sep = ".")] <- mean(chromosome.variant.table[selection, paste(tumorSample, "CNA", sep = ".")])
				}
			}
			if(number.of.segments > 1) {
				p.values <- unlist(lapply(1:(number.of.segments - 1), function(head.index) {schale.segment.p.value(head.index, head.index + 1)}))
				while(sum(p.values > undo.alpha) > 0) {
					head.index <- which.max(p.values)
					tail.index <- head.index + 1
					head.number.of.variants <- chromosome.segment.table$number.of.variants[head.index]
					tail.number.of.variants <- chromosome.segment.table$number.of.variants[tail.index]
					chromosome.segment.table$end[head.index] <- chromosome.segment.table$end[tail.index]
					chromosome.segment.table$number.of.variants[head.index] <- head.number.of.variants + tail.number.of.variants
					for(column in paste(rep(tumor.samples, each = 2), rep(c("LOH.mean", "CNA.mean"), times = length(tumor.samples)), sep = ".")) {
						chromosome.segment.table[head.index, column] <- (chromosome.segment.table[head.index, column] * head.number.of.variants + chromosome.segment.table[tail.index, column] * tail.number.of.variants) / (head.number.of.variants + tail.number.of.variants)
					}
					chromosome.segment.table <- chromosome.segment.table[-tail.index, ]
					number.of.segments <- number.of.segments - 1
					p.values <- p.values[-head.index]
					if(head.index > 1) {
						p.values[head.index - 1] <- schale.segment.p.value(head.index - 1, head.index)
					}
					if(head.index < number.of.segments) {
						p.values[head.index] <- schale.segment.p.value(head.index, head.index + 1)
					}
				}
			}
			segment.table <- rbind(segment.table, chromosome.segment.table)
		}
	}
	segment.table
}

schale.genome.plot <- function(variant.table, segment.table, tumor.sample, control.sample = NULL, chromosomes = NULL, chromosome.lengths = NULL, main = NULL, axis.lwd = 1, lwd = 2) {
	if(is.matrix(variant.table)) variant.table <- data.frame(variant.table)
	if(!is.data.frame(variant.table)) stop("variant.table must be data.frame or matrix.")
	if(!all(c("chromosome", "position") %in% colnames(variant.table))) stop("variant.table must have chromosome and position columns.")

	if(is.matrix(segment.table)) segment.table <- data.frame(segment.table)
	if(!is.data.frame(segment.table)) stop("segment.table must be data.frame or matrix.")
	if(!all(c("chromosome", "start", "end") %in% colnames(segment.table))) stop("segment.table must have chromosome, start and end columns.")

	available.samples <- intersect(schale.available.samples(variant.table, c("\\.LOH$", "\\.CNA$")), schale.available.samples(segment.table, c("\\.LOH\\.mean$", "\\.CNA\\.mean$")))
	if(length(available.samples) == 0) stop("variant.table must have [sample].LOH and [sample].CNA columns. segment.table must have [sample].LOH.mean and [sample].CNA.mean columns.")
	if(!(all(tumor.sample %in% available.samples) && length(tumor.sample) == 1)) stop("tumor.sample is not appropriate.")

	if(is.null(chromosomes)) {
		chromosomes <- as.character(unique(variant.table$chromosome))
		chromosomes <- chromosomes[grep("^(chr)?([0-9]*|X|Y)$", chromosomes)]
	}
	if(is.null(chromosome.lengths)) {
		chromosome.lengths <- unlist(lapply(chromosomes, function(chromosome) {
			max(variant.table$position[variant.table$chromosome == chromosome])
		}))
		names(chromosome.lengths) <- chromosomes
	}
	chromosome.table <- data.frame(chromosome = chromosomes, length = chromosome.lengths[chromosomes], row.names = chromosomes)
	chromosome.table$label <- sub("^chr", "", chromosome.table[, 1])
	chromosome.table$offset <- 0
	for(index in 2:(dim(chromosome.table)[1])) {
		chromosome.table$offset[index] <- chromosome.table$offset[index - 1] + chromosome.table$length[index - 1] + 1
	}
	xlim <- c(0, chromosome.table$offset[dim(chromosome.table)[1]] + chromosome.table$length[dim(chromosome.table)[1]] + 1)

	variant.table <- variant.table[variant.table$chromosome %in% chromosome.table[, 1], ]
	variant.table$x <- variant.table$position + chromosome.table[as.character(variant.table$chromosome), "offset"]

	segment.table <- segment.table[segment.table$chromosome %in% chromosome.table[, 1], ]
	segment.table$x.start <- segment.table$start + chromosome.table[as.character(segment.table$chromosome), "offset"]
	segment.table$x.end   <- segment.table$end   + chromosome.table[as.character(segment.table$chromosome), "offset"]

	if(is.null(control.sample)) {
		variant.table[, paste("plot", c("LOH", "CNA"), sep = ".")] <- variant.table[, paste(tumor.sample, c("LOH", "CNA"), sep = ".")]
		segment.table[, paste("plot", c("LOH.mean", "CNA.mean"), sep = ".")] <- segment.table[, paste(tumor.sample, c("LOH.mean", "CNA.mean"), sep = ".")]
	} else {
		variant.table[, paste("plot", c("LOH", "CNA"), sep = ".")] <- variant.table[, paste(tumor.sample, c("LOH", "CNA"), sep = ".")] - variant.table[, paste(control.sample, c("LOH", "CNA"), sep = ".")]
		segment.table[, paste("plot", c("LOH.mean", "CNA.mean"), sep = ".")] <- segment.table[, paste(tumor.sample, c("LOH.mean", "CNA.mean"), sep = ".")] - segment.table[, paste(control.sample, c("LOH.mean", "CNA.mean"), sep = ".")]
	}

	par(mar = c(5, 5, 4, 5) + 0.1)
	plot(x = NULL, y = NULL, xlim = xlim, ylim = c(-2.0, 2.0), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
	points(variant.table$x, variant.table$plot.CNA, pch = 20, cex = 0.1, col = "#8080FF")

	par(new = TRUE)
	plot(x = NULL, y = NULL, xlim = xlim, ylim = c(-0.5, 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
	points(variant.table$x, variant.table$plot.LOH, pch = 20, cex = 0.1, col = "#FF8080")

	abline(h = 0,                           lwd = axis.lwd, col = "#000000")
	abline(v = chromosome.table$offset[-1], lwd = axis.lwd, col = "#000000")
	axis(1, chromosome.table$offset + chromosome.table$length / 2, chromosome.table$label, lwd = 0)

	par(new = TRUE)
	plot(x = NULL, y = NULL, xlim = xlim, ylim = c(-2.0, 2.0), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
	segments(segment.table$x.start, segment.table$plot.CNA.mean, segment.table$x.end, segment.table$plot.CNA.mean, lwd = lwd, col = "#0000FF")
	axis(2, col = "#0000FF", col.axis = "#0000FF", yaxp = c(-2.0, 2.0, 4), lwd = axis.lwd)

	par(new = TRUE)
	plot(x = NULL, y = NULL, xlim = xlim, ylim = c(-0.5, 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
	segments(segment.table$x.start, segment.table$plot.LOH.mean, segment.table$x.end, segment.table$plot.LOH.mean, lwd = lwd, col = "#FF0000")
	axis(4, col = "#FF0000", col.axis = "#FF0000", yaxp = c(-0.5, 0.5, 2), lwd = axis.lwd)

	title(main = main, xlab = "Chromosome")
	mtext("CNA", side = 2, line = 3, col = "#0000FF")
	mtext("LOH", side = 4, line = 3, col = "#FF0000")
}

schale.segment.plot <- function(segment.table, tumor.sample, control.sample = NULL, ...) {
	if(is.matrix(segment.table)) segment.table <- data.frame(segment.table)
	if(!is.data.frame(segment.table)) stop("segment.table must be data.frame or matrix.")
	if(!("number.of.variants" %in% colnames(segment.table))) stop("segment.table must have number.of.variants column.")

	available.samples <- schale.available.samples(segment.table, c("\\.LOH\\.mean$", "\\.CNA\\.mean$"))
	if(length(available.samples) == 0) stop("segment.table must have [sample].LOH.mean and [sample].CNA.mean columns.")
	if(!(all(tumor.sample %in% available.samples) && length(tumor.sample) == 1)) stop("tumor.sample is not appropriate.")
	if(!is.null(control.sample)) {
		if(!(all(control.sample %in% available.samples) && length(control.sample) == 1)) stop("control.sample is not appropriate.")
	}

	segment.table$point.size <- log(segment.table$number.of.variants / mean(segment.table$number.of.variants) + 1, base = 2)
	if(is.null(control.sample)) {
		segment.table[, paste("plot", c("LOH.mean", "CNA.mean"), sep = ".")] <- segment.table[, paste(tumor.sample, c("LOH.mean", "CNA.mean"), sep = ".")]
	} else {
		segment.table[, paste("plot", c("LOH.mean", "CNA.mean"), sep = ".")] <- segment.table[, paste(tumor.sample, c("LOH.mean", "CNA.mean"), sep = ".")] - segment.table[, paste(control.sample, c("LOH.mean", "CNA.mean"), sep = ".")]
	}
	plot(segment.table[, paste("plot", c("LOH.mean", "CNA.mean"), sep = ".")], pch = 20, cex = segment.table$point.size, xlab = "LOH", ylab = "CNA", xlim = c(-0.5, 0.5), ylim = c(-2.0, 2.0), ...)
}

schale.calculate.CNA.offset <- function(segment.table, tumor.samples = NULL, cutoffs = c(0.01, 0.5, 100)) {
	if(is.matrix(segment.table)) segment.table <- data.frame(segment.table)
	if(!is.data.frame(segment.table)) stop("segment.table must be data.frame or matrix.")
	if(!("number.of.variants" %in% colnames(segment.table))) stop("segment.table must have number.of.variants column.")

	available.samples <- schale.available.samples(segment.table, c("\\.LOH\\.mean$", "\\.CNA\\.mean$"))
	if(length(available.samples) == 0) stop("segment.table must have [sample].LOH.mean and [sample].CNA.mean columns.")
	if(is.null(tumor.samples)) tumor.samples <- available.samples
	if(!(all(tumor.samples %in% available.samples) && length(tumor.samples) >= 1)) stop("tumor.samples is not appropriate.")

	if(!(is.numeric(cutoffs) && length(cutoffs) == 3)) stop("cutoffs is not appropriate.")

	CNA.offset <- unlist(lapply(tumor.samples, function(tumor.sample) {
		median(segment.table[abs(segment.table[, paste(tumor.sample, "LOH.mean", sep = ".")]) <= cutoffs[1] & segment.table[, paste(tumor.sample, "CNA.mean", sep = ".")] <= cutoffs[2] & segment.table[, "number.of.variants"] >= cutoffs[3], paste(tumor.sample, "CNA.mean", sep = ".")])
	}))
	names(CNA.offset) <- tumor.samples
	CNA.offset
}

schale.normalize.CNA.offset <- function(table, CNA.offset) {
	if(!(is.numeric(CNA.offset) && is.character(names(CNA.offset)))) stop("CNA.offset is not appropriate.")

	if(all(paste(names(CNA.offset), "CNA", sep = ".") %in% colnames(table))) {
		table[, paste(names(CNA.offset), "CNA", sep = ".")] <- t(t(table[, paste(names(CNA.offset), "CNA", sep = ".")]) - CNA.offset)
	}
	if(all(paste(names(CNA.offset), "CNA.mean", sep = ".") %in% colnames(table))) {
		table[, paste(names(CNA.offset), "CNA.mean", sep = ".")] <- t(t(table[, paste(names(CNA.offset), "CNA.mean", sep = ".")]) - CNA.offset)
	}
	table
}

schale.correlation.plot <- function(segment.table, sample.x, sample.y, type = "CNA", ...) {
	if(is.matrix(segment.table)) segment.table <- data.frame(segment.table)
	if(!is.data.frame(segment.table)) stop("segment.table must be data.frame or matrix.")
	if(!("number.of.variants" %in% colnames(segment.table))) stop("segment.table must have number.of.variants column.")

	available.samples <- schale.available.samples(segment.table, c("\\.LOH\\.mean$", "\\.CNA\\.mean$"))
	if(length(available.samples) == 0) stop("segment.table must have [sample].LOH.mean and [sample].CNA.mean columns.")
	if(!(all(sample.x %in% available.samples) && length(sample.x) == 1)) stop("sample.x is not appropriate.")
	if(!(all(sample.y %in% available.samples) && length(sample.y) == 1)) stop("sample.y is not appropriate.")

	segment.table[, c("x", "y")] <- segment.table[, paste(c(sample.x, sample.y), type, "mean", sep = ".")]
	segment.table$point.size <- log(segment.table$number.of.variants / mean(segment.table$number.of.variants) + 1, base = 2)
	lm.out <- lm(y ~ x, data = segment.table, weights = segment.table$point.size)

	minx <- min(segment.table$x)
	miny <- minx * lm.out$coefficients[2] + lm.out$coefficients[1]
	if(min(segment.table$y) < miny) {
		miny <- min(segment.table$y)
		minx <- (miny - lm.out$coefficients[1]) / lm.out$coefficients[2]
	}

	maxx <- max(segment.table$x)
	maxy <- maxx * lm.out$coefficients[2] + lm.out$coefficients[1]
	if(max(segment.table$y) > maxy) {
		maxy <- max(segment.table$y)
		maxx <- (maxy - lm.out$coefficients[1]) / lm.out$coefficients[2]
	}

	plot(segment.table$x, segment.table$y, pch = 20, cex = segment.table$point.size, xlim = c(minx, maxx), ylim = c(miny, maxy), xlab = sample.x, ylab = sample.y, pty = "s", ...)
	abline(lm.out)
}
