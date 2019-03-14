# An R file to be sourced.
#
# Rearrangement plot (aka chromothripsis paper type plot) for plotting
# rearrangement breakpoints and copy number.
#
# After sourcing workspace will have a variable 'chr_lens', and three functions:
# arc(), window_means() and plot_rearrangements().
#
# Author: Yilong Li
# Last update: 2015-03-02, added BAF plotting.
###############################################################################

# To read in this file:
# source("/Users/yl3/Documents/workspace/rg_ordering_pilot/rearrangement_plot.R")


library(quantsmooth)  # For ideogram

if (file.exists("~/work_related/programs/R_scripts/hg19_chr_sizes.txt")) {
    chr_lens = read.table("~/work_related/programs/R_scripts/hg19_chr_sizes.txt", header=F, sep="\t", row.names=1, colClasses = c("character", "numeric"))
} else {
    chr_lens = read.table("/nfs/users/nfs_y/yl3/programs/ucsc/hg19.chrom_sizes", header=F, sep="\t", row.names=1, colClasses = c("character", "numeric"))
}
temp = rownames(chr_lens)
chr_lens = chr_lens[,1]
names(chr_lens) = temp


arc = function(x0, x1, y, xr, yr, col, lwd) {
    x = (x0 + x1)/2  # Center of arc
    xr = x - x0 	 # x-radius of arc

    apply(
        cbind(x, y, xr, yr, col),
        1,
        function(z) {
            x   = as.numeric(z[1])
            y   = as.numeric(z[2])
            xr  = as.numeric(z[3])
            yr  = as.numeric(z[4])
            col = z[5]
            x_points = seq(x - xr, x + xr, length.out = 200)
            y_points = y + yr * sqrt( 1  -  ( (x_points-x) / xr )^2 )

            lines(
                x_points,
                y_points,
                col = col,
                lwd = lwd
            )
        }
    )

    return()
}

window_means = function (coords, cns, min_pos, max_pos, win_size) {
    cut_levels = cut(
        coords,
        seq(min_pos, max_pos + win_size, win_size),
        labels = F,
        include_lowest = T,
        right = F
    )

    cut_values = by(
        cns,
        cut_levels,
        function(x) mean(x, na.rm = T)
    )

    return(
        cut_values[ as.character(1:ceiling((max_pos-min_pos+1)/win_size)) ]
    )
}

plot_rearrangements = function(
    bedpe, chrs, cn_bedgraph = NULL, segments = NULL,
    yrange = NULL, ideogram=T, cn_cex=0.3, lwd = 0.75, cn_win_size = NULL,
    BFB_ids = c(), arrow_ln = 0.15, xlim = NULL, chr_lim = NULL, annot = NULL, main = NULL,
    muts = NULL, allele_freqs = NULL, allele_freq_segments = NULL, xaxt = T
) {
    chrs = as.character(chrs)
    if (!is.null(chr_lim)) {
        chr_lim = as.character(chr_lim)
    }
    chr_cum_lns = c(0, cumsum(chr_lens[chrs])[-length(chrs)])
    names(chr_cum_lns) = chrs
    xrange_size = cumsum(chr_lens[chrs])

    if (is.null(xlim)) {
        xlim = c(1, sum(chr_lens[chrs]))

        if (is.null(yrange)) {
            yrange = c(0, quantile(cn_bedgraph[cn_bedgraph[,1] %in% chrs, 4], p = 0.999, na.rm = T))
            yrange = c(floor(yrange[1]), ceiling(yrange[2]))
        }
    }
    else if (length(chrs) > 1) {
        if (!is.null(chr_lim)) {
            stop()
        }

        if (is.null(yrange)) {
            yrange = c(0, quantile(cn_bedgraph[cn_bedgraph[,1] %in% chr_lim & cn_bedgraph[,2] > xlim[1] & cn_bedgraph[,3] < xlim[2], 4], p = 0.999, na.rm = T))
            yrange = c(floor(yrange[1]), ceiling(yrange[2]))
        }

        xlim = chr_cum_lns[chr_lim] + xlim
    }
    else {
        if (is.null(yrange)) {
            yrange = c(0, quantile(cn_bedgraph[cn_bedgraph[,1] == chrs & cn_bedgraph[,2] > xlim[1] & cn_bedgraph[,3] < xlim[2], 4], p = 0.999, na.rm = T))
            yrange = c(floor(yrange[1]), ceiling(yrange[2]))
        }
    }

    if (!is.null(cn_bedgraph)) {
        cn = cn_bedgraph[cn_bedgraph[,1] %in% chrs, ]
    }
    else if (!is.null(segments)) {
    }
    else {
        stop("Either cn_bedgraph or segments must be provided")
    }

    # Set CN window size dynamically
    if (is.null(cn_win_size)) {
        if (is.null(xlim) || xlim[2] - xlim[1] >= 1e7) {
            cn_win_size = 1e4
        }
        else if (xlim[2] - xlim[1] >= 1e6) {
            cn_win_size = 1e3
        }
        else {
            cn_win_size = 5e2
        }
    }

    yrange_size = yrange[2] - yrange[1]

    td_col         = "darkorange4"
    del_col        = "darkslateblue"
    tail_tail_col  = "cadetblue4"
    head_head_col  = "chartreuse2"
    # td_col = "#E41A1C"
    # del_col = "#1E90FF"
    # head_head_col = "#9932CC"
    # tail_tail_col = "#4DAF4A"
    # inter_chrs_col = "darkorchid1"

    # Create the plot
    # par(mar = c(5, 4, 2, 2) + .1)
    layout(matrix(c(rep(1, 7), 2)))
    par(mar = c(0, 4, 4, 4) + .1, oma = c(5 + .1, 0, 0, 0))

    plot(
        c(),
        ylim = c(yrange[1] - 0.5, yrange[2] + 1.4*yrange_size),
        xlim = xlim,
        bty  = "n",
        yaxt = "n",
        xaxt = "n",
        xlab = "",
        ylab = "",
        yaxs = "i",
        xaxs = "i",
        main = main
    )


    # Shaded grid
    for (i in yrange[1]:yrange[2]) {
        polygon(
            c(1, sum(chr_lens[chrs]), sum(chr_lens[chrs]), 1),
            i - 0.5 + c(0, 0, 1, 1),
            col=rgb(.1, .1, .1, ifelse(i %% 2 == 0, 0.1, 0.05)),
            lty=0
        )
    }


    # Line to separate chromosomes
    if (length(chrs) > 1) {
        segments(
            x0 = cumsum(chr_lens[chrs])[-length(chrs)],
            y0 = yrange[1] - 0.5,
            y1 = yrange[2] + 0.5
        )
    }


    # Plot CN
    win_size = cn_win_size
    for (c in chrs) {
        if (!is.null(cn_bedgraph)) {
            # sel = cn[,1] == c
            # When plotting only one chromome with zoomed-in image, only plot
            # the data that will be visible in the graph.
            if (length(chrs) == 1 && !all(xlim == c(1, chr_lens[chrs]))) {
                sel = cn[,1] == c & cn[,2] > xlim[1] - 1e5 & cn[,3] < xlim[2] + 1e5
                x = win_size/2 + seq(xlim[1]-1e5, xlim[2]+1e5, win_size)
                y = window_means(rowMeans(cn[sel, 2:3]), cn[sel, 4], xlim[1]-1e5, xlim[2]+1e5, win_size)
                points(
                    x = x,
                    y = y,
                    pch = 16,
                    cex = cn_cex
                )
            }
            else {
                sel = cn[,1] == c
                points(
                    x = chr_cum_lns[c] + win_size/2 + seq(1, chr_lens[c], win_size),
                    y = window_means(rowMeans(cn[sel, 2:3]), cn[sel, 4], 1, chr_lens[c], win_size),
                    pch = 16,
                    cex = cn_cex
                )
            }
        }

        if (!is.null(segments)) {
            sel = segments[,1] == c

            segments(
                x0 = chr_cum_lns[c] + segments[sel, 2] + 1,
                x1 = chr_cum_lns[c] + segments[sel, 3],
                y0 = segments[sel, 4],
                lwd = 2,
                col = "blue"
            )
        }
    }
    axis(2, at = axisTicks(usr=yrange, log=F), las=2)
    title(ylab = "Copy number")


    # Plot rearrangements: First dotted lines
    segments(
        x0 = c(1, 1),
        x1 = c(sum(chr_lens[chrs]), sum(chr_lens[chrs])),
        y0 = c(yrange[2] + 0.3*yrange_size, yrange[2] + 0.75*yrange_size),
        lty = 3
    )
    abline(h = yrange[2] + 0.5)

    # Then 'intra-chromosomal' translocations
    sel = bedpe[,1] %in% chrs & bedpe[,4] %in% chrs & !(bedpe[,7] %in% BFB_ids)  # Breakage-fusion-bridge RGMTS to be plotted separately
    col =
        ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col,
        ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "-", del_col,
        ifelse( bedpe[sel, 9] == "-" & bedpe[sel, 10] == "+", td_col, tail_tail_col)))

    if (sum(sel) > 0) {
        arc(
            x0  = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            x1  = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
            y   = ifelse(bedpe[sel,9] == bedpe[sel, 10], yrange[2]+.3*yrange_size, yrange[2]+.75*yrange_size),
            yr  = ifelse(bedpe[sel, 10] == "-", 1, -1) * 0.2 * yrange_size,
            col = col,
            lwd = lwd
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75*yrange_size, 0.3*yrange_size),
            y1 = yrange[1],
            col = rgb(t(col2rgb(col)), alpha = 127, max=255),
            lwd = lwd
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
            y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75*yrange_size, 0.3*yrange_size),
            y1 = yrange[1],
            col = rgb(t(col2rgb(col)), alpha = 127, max=255),
            lwd = lwd
        )
    }

    # Then rearrangments where low end %in% chrs and !(high end %in% chrs)
    sel = bedpe[,1] %in% chrs & !(bedpe[,4] %in% chrs)
    if (sum(sel) > 0) {
        # arrows(
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange[1],
            y1 = yrange[2] + yrange_size,
            col = rgb(t(col2rgb("black")), alpha = 127, max=255),
            lwd = lwd  # ,
            # length = arrow_ln
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange[2] + yrange_size,
            x1 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]) + ifelse(bedpe[sel, 9] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            y1 = ifelse(bedpe[sel, 9] == "+", yrange[2] + 1.1 * yrange_size, yrange[2] + 1.3 * yrange_size),
            xpd = NA
        )
        text(
            as.character(bedpe[sel, 4]),
            x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]) + ifelse(bedpe[sel, 9] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            # y = yrange[2] + 1.2 * yrange_size
            y = yrange[2] + ifelse(bedpe[sel, 9] == "+", 1.2, 1.4) * yrange_size,
            xpd = NA
        )
    }

    # Then rearrangments where high end %in% chrs and !(low end %in% chrs)
    sel = !(bedpe[,1] %in% chrs) & bedpe[,4] %in% chrs
    if (sum(sel) > 0) {
        # arrows(
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
            y0 = yrange[1],
            y1 = yrange[2] + yrange_size,
            col = rgb(t(col2rgb("black")), alpha = 127, max=255),
            lwd = lwd  # ,
            # length = arrow_ln
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
            y0 = yrange[2] + yrange_size,
            x1 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]) + ifelse(bedpe[sel, 10] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            y1 = ifelse(bedpe[sel, 10] == "+", yrange[2] + 1.1 * yrange_size, yrange[2] + 1.3 * yrange_size),
            xpd = NA
        )
        text(
            as.character(bedpe[sel, 1]),
            x = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]) + ifelse(bedpe[sel, 10] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            # y = yrange[2] + 1.2 * yrange_size
            y = yrange[2] + ifelse(bedpe[sel, 10] == "+", 1.2, 1.4) * yrange_size,
            xpd = NA
        )
    }


    # Then BFBs, currently only supporting 2 BFB events max
    # Also, in case of BFBs, only plotting 1st end
    if (length(BFB_ids) > 2) {
        stop()
    }
    if (length(BFB_ids) == 2) {
        sel = bedpe[,7] == BFB_ids[2]
        col =
            ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col,
                ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "-", del_col,
                    ifelse( bedpe[sel, 9] == "-" & bedpe[sel, 10] == "+", td_col, 		 tail_tail_col)))
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange[2] + 1.2*yrange_size,
            y1 = yrange[1],
            col = rgb(t(col2rgb(col)), alpha = 127, max=255)
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange[2] + 1.2*yrange_size,
            y1 = yrange[2] + (1.2 + 0.1)*yrange_size,
            col = col
        )

        # The curved arrow
        lines(
            x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(seq(-pi/2, pi/2, length.out=20)) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50,
            y = yrange[2] + (1 + 0.3)*yrange_size +  0.05*yrange_size*sin(seq(-pi/2, pi/2, length.out=20)),
            col = col,
            lwd = 2 * lwd
        )
        x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(-pi/2) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50
        arrows(
            x0 = x0,
            x1 = x0 + ifelse(bedpe[sel,9] == "+", -1, 1) * xrange_size / 50,
            y0 = yrange[2] + (1 + 0.3)*yrange_size +  0.05*yrange_size*sin(-pi/2),
            col = col,
            lwd = 2 * lwd,
            length = arrow_ln
        )
    }
    if (length(BFB_ids) >= 1) {
        abline(h = yrange[2] + c(1, 1.2)*yrange_size)
        if (length(BFB_ids) > 1) abline(h = yrange[2] + 1.4*yrange_size)

        sel = bedpe[,7] == BFB_ids[1]
        col =
            ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col,
                ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "-", del_col,
                    ifelse( bedpe[sel, 9] == "-" & bedpe[sel, 10] == "+", td_col, 		 tail_tail_col)))
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange[2] + 1*yrange_size,
            y1 = yrange[1],
            col = rgb(t(col2rgb(col)), alpha = 127, max=255)
        )
        segments(
            x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            y0 = yrange[2] + 1*yrange_size,
            y1 = yrange[2] + (1 + 0.1)*yrange_size,
            col = col
        )

        # The curved arrow
        lines(
            x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(seq(-pi/2, pi/2, length.out=20)) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50,
            y = yrange[2] + (1 + 0.1)*yrange_size +  0.05*yrange_size*sin(seq(-pi/2, pi/2, length.out=20)),
            col = col,
            lwd = 2 * lwd
        )
        x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(-pi/2) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50
        arrows(
            x0 = x0,
            x1 = x0 + ifelse(bedpe[sel,9] == "+", -1, 1) * xrange_size / 50,
            y0 = yrange[2] + (1 + 0.1)*yrange_size +  0.05*yrange_size*sin(-pi/2),
            col = col,
            lwd = 2 * lwd,
            length = arrow_ln
        )
    }


    # Annotations on top
    if (!is.null(annot) && sum(annot[,1] %in% chrs > 0)) {
        sel = annot[,1] %in% chrs
        segments(
            x0 = chr_cum_lns[as.character(annot[sel, 1])] + annot[sel, 2],
            y0 = 1.1 * yrange_size + yrange[2] + (0:(sum(annot[,1] %in% chrs)-1)) * 0.05 * yrange_size,
            x1 = chr_cum_lns[as.character(annot[sel, 1])] + annot[sel, 3],
            lwd = 2
        )

        text(
            chr_cum_lns[as.character(annot[sel, 1])] + rowMeans(annot[sel, 2:3]),
            1.2 * yrange_size + yrange[2] + (0:(sum(annot[,1] %in% chrs)-1)) * 0.05 * yrange_size,
            labels = annot[sel,4],
            xpd = NA
        )
    }


    #
    # Plot point mutations?
    #
    if (!is.null(muts)) {
        rc = c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
        muts = muts[muts[,1] %in% chrs, ]
        muts[,4] = ifelse(muts[,3] %in% c("A", "G"), rc[muts[,4]], muts[,4])
        muts[,3] = ifelse(muts[,3] %in% c("A", "G"), rc[muts[,3]], muts[,3])
        mut_col = c(
            "C>A" = "orange",
            "C>G" = "brown",
            "C>T" = "red",
            "T>A" = "purple",
            "T>C" = "pink",
            "T>G" = "lightgreen"
        )

        muts = muts[order(muts[,1], muts[,2]), ]
        dist_to_next = c(muts[-1,2] - muts[-nrow(muts), 2])
        dist_to_next[muts[-1, 1] != muts[-nrow(muts), 1]] = NA
        dist_to_next = c(dist_to_next, NA)

        xlim = par("usr")[1:2]
        par(new = T)
        pos = chr_cum_lns[as.character(muts[,1])] + muts[,2]
        plot(
            pos,
            dist_to_next,
            xlim = xlim,
            xaxs = "i",
            axes = F,
            xlab = "",
            ylab = "",
            main = "",
            yaxt = "n",
            xaxt = "n",
            log = "y",
            col = mut_col[paste(muts[,3], muts[,4], sep = ">")],
            pch = 16,
            ylim = c(1, 1e7)
        )
        axis(4)

        legend("topleft", bty = "n", pch = 16, col = mut_col, legend = names(mut_col), ncol = 2)
    }


    # Now plot the B-allele frequency
    par(mar = c(0, 4, 0, 4) + .1)
    plot(
        c(),
        ylim = c(-.6, 1),
        xlim = xlim,
        bty  = "n",
        yaxt = "n",
        xaxt = "n",
        xlab = "",
        ylab = "",
        yaxs = "i",
        xaxs = "i"
    )
    abline(h = seq(0, 1, by = 0.25), lwd = lwd, col = "lightgrey")

    for (c in chrs) {
        if (!is.null(allele_freqs)) {
            # When plotting only one chromosome with zoomed-in image, only plot
            # the data that will be visible in the graph.
            if (length(chrs) == 1 && !all(xlim == c(1, chr_lens[chrs]))) {
                sel = allele_freqs[,1] == c & allele_freqs[,2] > xlim[1] - 1e5 & allele_freqs[,3] < xlim[2] + 1e5
                x = allele_freqs[sel, 3]
                y = allele_freqs[sel, 4]
                points(
                    x = x,
                    y = y,
                    pch = 16,
                    cex = cn_cex
                )
            }
            else {
                sel = allele_freqs[,1] == c
                points(
                    x = chr_cum_lns[c] + allele_freqs[sel, 3],
                    y = allele_freqs[sel, 4],
                    pch = 16,
                    cex = cn_cex
                )
            }
        }
        if (!is.null(allele_freq_segments)) {
            # When plotting only one chromosome with zoomed-in image, only plot
            # the data that will be visible in the graph.
            if (length(chrs) == 1 && !all(xlim == c(1, chr_lens[chrs]))) {
                sel = allele_freq_segments[,1] == c & allele_freq_segments[,3] > xlim[1] & allele_freq_segments[,2] < xlim[2]
                x0 = allele_freq_segments[sel, 2]
                x1 = allele_freq_segments[sel, 3]
                y = allele_freq_segments[sel, 4]
                segments(
                    x0,
                    y,
                    x1,
                    col = ifelse(x0 == x1, "chartreuse2", "blue"),
                    lwd = 2
                )
            }
            else {
                sel = allele_freq_segments[,1] == c
                x0 = chr_cum_lns[c] + allele_freq_segments[sel, 2]
                y0 = allele_freq_segments[sel, 4]
                x1 = chr_cum_lns[c] + allele_freq_segments[sel, 3]
                segments(
                    x0 = x0,
                    y0 = y0,
                    x1 = x1,
                    col = ifelse(x0 == x1, "chartreuse2", "blue"),
                    lwd = 2
                )
            }
        }
    }


    # X axis names and ticks
    par(mgp = par("mgp") + c(0,1,0))
    if (all(xlim == c(1, sum(chr_lens[chrs])))) {
        axis(
            1,
            at = (cumsum(chr_lens[chrs]) + chr_cum_lns)/2,
            labels = ifelse(rep(xaxt, length(chrs)), paste("chr", chrs, " position (Mb)", sep=""), chrs),
            tick = F,
            cex.lab = 1.5
        )
    }
    else if (!is.null(chr_lim)) {
        axis(
            1,
            at = mean(xlim),
            labels = paste("chr", chr_lim, " position (Mb)", sep=""),
            tick = F,
            cex.lab = 1.5
        )
    }
    else {
        axis(
            1,
            at = mean(xlim),
            labels = paste("chr", chrs, " position (Mb)", sep=""),
            tick = F,
            cex.lab = 1.5
        )
    }
    par(mgp = par("mgp") - c(0,1,0))

    if (xaxt) {
        if (length(chrs) > 1) {
            if (all(xlim == c(1, sum(chr_lens[chrs])))) {
                for (c in chrs) {
                    pretty_ticks = pretty(c(1, chr_lens[c]))
                    pretty_ticks = pretty_ticks[which(pretty_ticks < chr_lens[c])]
                    axis(1, at = pretty_ticks + chr_cum_lns[c], labels = pretty_ticks/1e6)
                }
            }
            else {
                if (is.null(chr_lim) || !(chr_lim %in% names(chr_lens))) {
                    stop()
                }

                pretty_ticks = pretty(xlim - chr_cum_lns[chr_lim])

                axis(1, at = pretty_ticks + chr_cum_lns[chr_lim], labels = pretty_ticks/1e6)
            }
        }
        else {
            if (all(xlim == c(1, sum(chr_lens[chrs])))) {
                axis(1, at = axisTicks(usr=c(1, chr_lens[chrs]), log=F), labels = axisTicks(usr=c(1, chr_lens[chrs]), log=F)/1e6)
            }
            else {
                pretty_ticks = pretty(xlim)
                axis(1, at = pretty_ticks, labels = pretty_ticks / 1e6)
            }
        }
    }


    # Finally the ideogram
    if (xlim[2] - xlim[1] < 10e6) {
        print("Ideogram plotting disabled because xlim[2] - xlim[1] < 10e6")

        ideogram = F
    }
    if (ideogram) {
        for (c in chrs) {
            paintCytobands(
                c,
                pos = c(1 + chr_cum_lns[c], -0.1),
                units = "bases",
                width = 0.5,
                length.out = chr_lens[c],
                legend = F,
                xpd = NA
            )
        }
    }
}
