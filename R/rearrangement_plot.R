#' Function for plotting campbellgrams
#'
#' Copyright: 2019 Yilong Li <yilong.li.yl3@gmail.com>


#' Add an arc to the currently active plot.
#'
#' @param x X-axis coordinate of the centre of the arc.
#' @param y Y-axis coordinate of the centre of the arc.
#' @param xr X-axis radius.
#' @param yr Y-axis radius.
#' @param col colour of the arc.
#' @param lwd line width of the arc.
draw_arc = function(x, y, xr, yr, col, lwd) {
    x_points = seq(x - xr, x + xr, length.out = 200)
    y_points = y + yr * sqrt(1 - ( (x_points - x) / xr) ^ 2 )
    lines(x_points, y_points, col = col, lwd = lwd)
}


#' Add horizontal arcs to the current plot.
#'
#' Draw the top or bottom half of a horizontal ellipsoid, i.e. one where both
#' focal points are on the same Y-axis coordinate.
#'
#' @param x0 start points of an arc.
#' @param x1 end points of an arc.
#' @param y Y-axis coordinates of the arc.
#' @param yr Y-axis radiuses. The sign of these values determine whether each
#'   arc is drawn above or below \emph{y}.
#' @param col colour of the arcs.
#' @param lwd line width of the arcs.
draw_arcs = function(x0, x1, y, yr, col, lwd) {
    x = (x0 + x1) / 2  # Center of arc
    xr = x - x0  # x-radius of arc
    mapply_args = append(
        list(FUN = draw_arc),
        as.list(data.frame(x, y, xr, yr, col, lwd))
    )
    do.call(mapply, mapply_args)
}


#' Windowed means by chromosomal coordinate.
#'
#' Compute non-overlapping windowed means by setting window size based on
#' chromosomal coordinates.
#'
#' @param coords coordinates of the data points.
#' @param cns copy number levels for each coordinate.
#' @param min_pos minimum \emph{coordinate} to accept. Coordinates lower than
#'   this are discarded.
#' @param max_pos maximum \emph{coordinate} to accept. Coordinates higher than
#'   this are discarded.
#' @param win_size window size.
#'
#' @return list of centre coordinates of the windows ("window_coords") and
#'   their average values ("window_means").
window_means = function(coords, cns, min_pos, max_pos, win_size) {
    cut_levels = cut(
      coords,
      seq(min_pos, max_pos + win_size, win_size),
      labels = F,
      include_lowest = T,
      right = F
    )

    window_coords = win_size / 2  +  seq(min_pos, max_pos, win_size)
    cut_values = by(cns, cut_levels, function(x) mean(x, na.rm = T))
    window_means = cut_values[
        as.character(1:ceiling( (max_pos - min_pos + 1) / win_size ))
    ]

    return(list(window_coords, window_means))
}


#' Helper function for validating chrs_used.
set_chrs_used = function(chrs_used, ref_genome) {
    chrs_used = as.character(chrs_used)
    if (any(!(chrs_used %in% ref_genome[, 1]))) {
        msg = paste("Chromosome name(s) in chrs_used are not present on the ",
                    "first column of parameter ref_genome.", sep = "")
        stop(msg)
    }
    return(chrs_used)
}


#' Helper function for setting chrs_shown based on provided parameter values.
#'
#' Also performs additional validations.
set_chrs_shown = function(chrs_used, chrs_shown) {
    if (length(chrs_shown) > 1) {
        stop("Parameter 'chrs_shown' must be of length 1.")
    }
    if (is.null(chrs_shown)) {
        return(chrs_used)
    } else {
        return(as.character(chrs_shown))
    }
}


#' Helper function for setting xlim.
#'
#' If only one chromosome is to be shown, then set xlim as the entire chromosome
#' unless xlim is already provided. Otherwise, set xlim as NULL.
set_xlim = function(xlim, chrs_shown, chr_lens) {
    if (length(chrs_shown) == 1) {
        if (is.null(xlim)) {
            xlim = c(1, chr_lens[chrs_shown])
        }
        return(xlim)
    } else {
        if (!is.null(xlim)) {
            message(paste("Parameter xlim ignored, since more than one ",
                          "chromosome is to be shown", sep = "\t"))
        }
        return(NULL)
    }
}


#' Subset cn_bedgraph to actually visible regions.
subset_cn_bedgraph = function(cn_bedgraph, chrs_shown, xlim) {
    if (length(chrs_shown) == 1) {
        idx = (cn_bedgraph[, 1] == chrs_shown
               & cn_bedgraph[, 2] <= xlim[2]
               & cn_bedgraph[, 3] >= xlim[1])
    } else {
        idx = cn_bedgraph[, 1] %in% chrs_shown
    }
    return(cn_bedgraph[idx, ])
}


#' Helper function for computing a copy number axis range for
#' \code{campbellgram()}.
compute_yrange = function(cn_bedgraph) {
    if (is.null(cn_bedgraph)) {
        # No copy number data for computing yrange - just return default values.
        return(c(0, 4))
    } else {
        return(c(0, quantile(cn_bedgraph, p = 0.999, na.rm = T)))
    }
}


#' Helper function for setting CN averaging window size dynamically.
set_cn_win_size = function(xlim) {
    if (xlim[2] - xlim[1] >= 1e+07) {
        return(10000)
    } else if (xlim[2] - xlim[1] >= 1e+06) {
        return(1000)
    } else {
        return(500)
    }
}


ARC_COLS = c(
    td_col = "darkorange4",
    del_col = "darkslateblue",
    tail_tail_col = "cadetblue4",
    head_head_col = "chartreuse2"
)


#' Plot a campbellgram
#'
#' Generate a campbellgram with rearrangement junctions at the top of the panel,
#' and optionally absolute copy number, copy number segmentation, B allele
#' frequency and point mutation rainfall plot at the bottom of the panel.
#'
#' Parameter \emph{chrs} controls which junctions are plotted as arcs and which
#' ones are plotted as "inter-chromosomal" junctions. Parameter \emph{chr_lim}
#' controls which chromosomes are actually plotted. For instance, by setting
#' \emph{chr_lim = "1"} and \emph{chrs = c("1", "2")}, chromosomes 1 will be
#' plotted, any intra-chromosomal junctions in chromosome 1 or inter-chromosomal
#' junctions between chromosome 1 and chromosome 2 will be visualised as arcs,
#' and all other inter-chromosomal junctions are visualised as inter-chromosomal
#' junctions.
#'
#' One way to generate the \emph{ref_genome} table is using
#' \emph{samtools faidx}.
#'
#' @param ref_genome data frame with chromosome name on the first column and
#'   chromosome length on the second column.
#' @param bedpe SV junctions in BEDPE format
#'   (\link{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}).
#'   TODO(yl3): default
#' @param chrs_used rearrangement junctions within and between these chromosomes
#'   are visualised as intra-chromosomal junctions with arcs. Other
#'   rearrangement junctions are visualised as vertical lines. Names for
#'   \emph{chrs_used} must be defined in \emph{ref_genome}.
#' @param chrs_shown either a single chromosome name or NULL: chromosomes to be
#'   displayed in the campbellgram. By default (NULL), all chromosomes in
#'   \emph{chrs_used} are shown. If provided, the one provided chromosome is
#'   displayed instead. Whether a rearrangement junction is displayed as an arc
#'   or a vertical line is controlled using parameter \emph{chrs_used}.
#' @param xlim X-axis limits for chromosomal coordinates for zooming into a sub-
#'   chromosomal locus. This parameter is ignored if \emph{chrs_displ} contains
#'   more than one chromosome.
#' @param cn_bedgraph absolute copy number (for instance as opposed to log R
#'   ratio) in bedGraph format
#'   (\link{https://genome.ucsc.edu/goldenPath/help/bedgraph.html}).
#' @param segments copy number segments to overlay in bedGraph format
#'   (\link{https://genome.ucsc.edu/goldenPath/help/bedgraph.html}).
#' @param yrange Y-axis range for copy number.
#'   Default: \code{c(0, quantile(cn_bedgraph[,4], 0.999))}.
#' @param ideogram: whether to plot an ideogram. By default, an ideogram is
#'   shown if the X-axis spans more than 1 Mb.
#' @param cn_win_size window size for averaging copy number before plotting
#'   them. By default, if more than 1e7 bp is displayed, then window size is
#'   10 kb, otherwise if more than 1e6 bp is displayed, then window size is
#'   1 kb, otherwise window size is set to 500 bp.
#' @param arc_cols named vector with labels "td_col", "del_col", "tail_tail_col"
#'   and "head_head_col" for the colors of the arcs of different rearrangement
#'   junction orientations.
campbellgram = function(ref_genome, bedpe, chrs_used, chrs_shown = NULL,
                        xlim = NULL, cn_bedgraph = NULL, segments = NULL,
                        yrange = NULL, ideogram = NULL, cn_win_size = NULL,
                        arc_cols = ARC_COLS, cn_cex = 0.3, lwd = 0.75,
                        BFB_ids = c(), arrow_ln = 0.15, annot = NULL,
                        main = NULL, muts = NULL, allele_freqs = NULL,
                        allele_freq_segments = NULL, xaxt = T) {
    # Internally, the X-axis is composed of all chromosomes in chrs_used stacked
    # side by side. The actual part of the chromosome shown is controlled by
    # chrs_displ and/or xlim.

    # Create chromosome lengths vector. Ensure chromosome names are characters.
    ref_genome[, 1] = as.character(ref_genome[, 1])
    chr_lens = setNames(ref_genome[, 2], ref_genome[, 1])
    chrs_used = set_chrs_used(chrs_used, ref_genome)
    chrs_shown = set_chrs_shown(chrs_used, chrs_shown)

    # Compute the cumulative chromosomal coordinate for each chromosome.
    chr_coord_offset = c(0, cumsum(chr_lens[chrs_used])[-length(chrs_used)])
    names(chr_coord_offset) = chrs_used

    # Subset cn_bedgraph to coordinates actually to be plotted. Set yrange using
    # subsetted copy number data.
    xlim = set_xlim(xlim, chrs_shown, chr_lens)
    if (!is.null(cn_bedgraph)) {
        plotted_cn = subset_cn_bedgraph(cn_bedgraph, chrs_shown, xlim)
    }
    if (is.null(yrange)) {
        yrange = compute_yrange(plotted_cn)
    }

    # Set CN window size dynamically
    if (is.null(cn_win_size)) {
        cn_win_size = set_cn_win_size(xlim)
    }

    xrange_size = cumsum(chr_lens[chrs_used])  # TODO(yl3): what is this?
    yrange_size = yrange[2] - yrange[1]

    # Create the plot par(mar = c(5, 4, 2, 2) + .1)
    layout(matrix(c(rep(1, 7), 2)))
    par(mar = c(0, 4, 4, 4) + 0.1, oma = c(5 + 0.1, 0, 0, 0))
    
    plot(c(), ylim = c(yrange[1] - 0.5, yrange[2] + 1.4 * yrange_size), xlim = xlim, 
        bty = "n", yaxt = "n", xaxt = "n", xlab = "", ylab = "", yaxs = "i", xaxs = "i", 
        main = main)
    
    
    # Shaded grid
    for (i in yrange[1]:yrange[2]) {
        polygon(c(1, sum(chr_lens[chrs]), sum(chr_lens[chrs]), 1), i - 0.5 + c(0, 
            0, 1, 1), col = rgb(0.1, 0.1, 0.1, ifelse(i%%2 == 0, 0.1, 0.05)), lty = 0)
    }
    
    
    # Line to separate chromosomes
    if (length(chrs) > 1) {
        segments(x0 = cumsum(chr_lens[chrs])[-length(chrs)], y0 = yrange[1] - 0.5, 
            y1 = yrange[2] + 0.5)
    }
    
    
    # Plot CN
    win_size = cn_win_size
    for (c in chrs) {
        if (!is.null(cn_bedgraph)) {
            # sel = cn[,1] == c When plotting only one chromome with zoomed-in image, only
            # plot the data that will be visible in the graph.
            if (length(chrs) == 1 && !all(xlim == c(1, chr_lens[chrs]))) {
                sel = cn[, 1] == c & cn[, 2] > xlim[1] - 1e+05 & cn[, 3] < xlim[2] + 
                  1e+05
                x = win_size/2 + seq(xlim[1] - 1e+05, xlim[2] + 1e+05, win_size)
                y = window_means(rowMeans(cn[sel, 2:3]), cn[sel, 4], xlim[1] - 1e+05, 
                  xlim[2] + 1e+05, win_size)
                points(x = x, y = y, pch = 16, cex = cn_cex)
            } else {
                sel = cn[, 1] == c
                points(x = chr_cum_lns[c] + win_size/2 + seq(1, chr_lens[c], win_size), 
                  y = window_means(rowMeans(cn[sel, 2:3]), cn[sel, 4], 1, chr_lens[c], 
                    win_size), pch = 16, cex = cn_cex)
            }
        }
        
        if (!is.null(segments)) {
            sel = segments[, 1] == c
            
            segments(x0 = chr_cum_lns[c] + segments[sel, 2] + 1, x1 = chr_cum_lns[c] + 
                segments[sel, 3], y0 = segments[sel, 4], lwd = 2, col = "blue")
        }
    }
    axis(2, at = axisTicks(usr = yrange, log = F), las = 2)
    title(ylab = "Copy number")
    
    
    # Plot rearrangements: First dotted lines
    segments(x0 = c(1, 1), x1 = c(sum(chr_lens[chrs]), sum(chr_lens[chrs])), y0 = c(yrange[2] + 
        0.3 * yrange_size, yrange[2] + 0.75 * yrange_size), lty = 3)
    abline(h = yrange[2] + 0.5)
    
    # Then 'intra-chromosomal' translocations
    sel = bedpe[, 1] %in% chrs & bedpe[, 4] %in% chrs & !(bedpe[, 7] %in% BFB_ids)  # Breakage-fusion-bridge RGMTS to be plotted separately
    col = ifelse(bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col, ifelse(bedpe[sel, 
        9] == "+" & bedpe[sel, 10] == "-", del_col, ifelse(bedpe[sel, 9] == "-" & 
        bedpe[sel, 10] == "+", td_col, tail_tail_col)))
    
    if (sum(sel) > 0) {
        arc(x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]), 
            x1 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]), 
            y = ifelse(bedpe[sel, 9] == bedpe[sel, 10], yrange[2] + 0.3 * yrange_size, 
                yrange[2] + 0.75 * yrange_size), yr = ifelse(bedpe[sel, 10] == "-", 
                1, -1) * 0.2 * yrange_size, col = col, lwd = lwd)
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 
            2:3]), y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75 * yrange_size, 
            0.3 * yrange_size), y1 = yrange[1], col = rgb(t(col2rgb(col)), alpha = 127, 
            max = 255), lwd = lwd)
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 
            5:6]), y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75 * yrange_size, 
            0.3 * yrange_size), y1 = yrange[1], col = rgb(t(col2rgb(col)), alpha = 127, 
            max = 255), lwd = lwd)
    }
    
    # Then rearrangments where low end %in% chrs and !(high end %in% chrs)
    sel = bedpe[, 1] %in% chrs & !(bedpe[, 4] %in% chrs)
    if (sum(sel) > 0) {
        # arrows(
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 
            2:3]), y0 = yrange[1], y1 = yrange[2] + yrange_size, col = rgb(t(col2rgb("black")), 
            alpha = 127, max = 255), lwd = lwd)
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 
            2:3]), y0 = yrange[2] + yrange_size, x1 = chr_cum_lns[as.character(bedpe[sel, 
            1])] + rowMeans(bedpe[sel, 2:3]) + ifelse(bedpe[sel, 9] == "+", -1, +1) * 
            (par("usr")[2] - par("usr")[1])/50, y1 = ifelse(bedpe[sel, 9] == "+", 
            yrange[2] + 1.1 * yrange_size, yrange[2] + 1.3 * yrange_size), xpd = NA)
        text(as.character(bedpe[sel, 4]), x = chr_cum_lns[as.character(bedpe[sel, 
            1])] + rowMeans(bedpe[sel, 2:3]) + ifelse(bedpe[sel, 9] == "+", -1, +1) * 
            (par("usr")[2] - par("usr")[1])/50, y = yrange[2] + ifelse(bedpe[sel, 
            9] == "+", 1.2, 1.4) * yrange_size, xpd = NA)
    }
    
    # Then rearrangments where high end %in% chrs and !(low end %in% chrs)
    sel = !(bedpe[, 1] %in% chrs) & bedpe[, 4] %in% chrs
    if (sum(sel) > 0) {
        # arrows(
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 
            5:6]), y0 = yrange[1], y1 = yrange[2] + yrange_size, col = rgb(t(col2rgb("black")), 
            alpha = 127, max = 255), lwd = lwd)
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 
            5:6]), y0 = yrange[2] + yrange_size, x1 = chr_cum_lns[as.character(bedpe[sel, 
            4])] + rowMeans(bedpe[sel, 5:6]) + ifelse(bedpe[sel, 10] == "+", -1, 
            +1) * (par("usr")[2] - par("usr")[1])/50, y1 = ifelse(bedpe[sel, 10] == 
            "+", yrange[2] + 1.1 * yrange_size, yrange[2] + 1.3 * yrange_size), xpd = NA)
        text(as.character(bedpe[sel, 1]), x = chr_cum_lns[as.character(bedpe[sel, 
            4])] + rowMeans(bedpe[sel, 5:6]) + ifelse(bedpe[sel, 10] == "+", -1, 
            +1) * (par("usr")[2] - par("usr")[1])/50, y = yrange[2] + ifelse(bedpe[sel, 
            10] == "+", 1.2, 1.4) * yrange_size, xpd = NA)
    }
    
    
    # Then BFBs, currently only supporting 2 BFB events max Also, in case of BFBs,
    # only plotting 1st end
    if (length(BFB_ids) > 2) {
        stop()
    }
    if (length(BFB_ids) == 2) {
        sel = bedpe[, 7] == BFB_ids[2]
        col = ifelse(bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col, 
            ifelse(bedpe[sel, 9] == "+" & bedpe[sel, 10] == "-", del_col, ifelse(bedpe[sel, 
                9] == "-" & bedpe[sel, 10] == "+", td_col, tail_tail_col)))
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 
            2:3]), y0 = yrange[2] + 1.2 * yrange_size, y1 = yrange[1], col = rgb(t(col2rgb(col)), 
            alpha = 127, max = 255))
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 
            2:3]), y0 = yrange[2] + 1.2 * yrange_size, y1 = yrange[2] + (1.2 + 0.1) * 
            yrange_size, col = col)
        
        # The curved arrow
        lines(x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 
            2:3]) + xrange_size/50 * cos(seq(-pi/2, pi/2, length.out = 20)) * ifelse(bedpe[sel, 
            9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1) * xrange_size/50, 
            y = yrange[2] + (1 + 0.3) * yrange_size + 0.05 * yrange_size * sin(seq(-pi/2, 
                pi/2, length.out = 20)), col = col, lwd = 2 * lwd)
        x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]) + 
            xrange_size/50 * cos(-pi/2) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 
            9] == "+", -1, 1) * xrange_size/50
        arrows(x0 = x0, x1 = x0 + ifelse(bedpe[sel, 9] == "+", -1, 1) * xrange_size/50, 
            y0 = yrange[2] + (1 + 0.3) * yrange_size + 0.05 * yrange_size * sin(-pi/2), 
            col = col, lwd = 2 * lwd, length = arrow_ln)
    }
    if (length(BFB_ids) >= 1) {
        abline(h = yrange[2] + c(1, 1.2) * yrange_size)
        if (length(BFB_ids) > 1) 
            abline(h = yrange[2] + 1.4 * yrange_size)
        
        sel = bedpe[, 7] == BFB_ids[1]
        col = ifelse(bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col, 
            ifelse(bedpe[sel, 9] == "+" & bedpe[sel, 10] == "-", del_col, ifelse(bedpe[sel, 
                9] == "-" & bedpe[sel, 10] == "+", td_col, tail_tail_col)))
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 
            2:3]), y0 = yrange[2] + 1 * yrange_size, y1 = yrange[1], col = rgb(t(col2rgb(col)), 
            alpha = 127, max = 255))
        segments(x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 
            2:3]), y0 = yrange[2] + 1 * yrange_size, y1 = yrange[2] + (1 + 0.1) * 
            yrange_size, col = col)
        
        # The curved arrow
        lines(x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 
            2:3]) + xrange_size/50 * cos(seq(-pi/2, pi/2, length.out = 20)) * ifelse(bedpe[sel, 
            9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1) * xrange_size/50, 
            y = yrange[2] + (1 + 0.1) * yrange_size + 0.05 * yrange_size * sin(seq(-pi/2, 
                pi/2, length.out = 20)), col = col, lwd = 2 * lwd)
        x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]) + 
            xrange_size/50 * cos(-pi/2) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 
            9] == "+", -1, 1) * xrange_size/50
        arrows(x0 = x0, x1 = x0 + ifelse(bedpe[sel, 9] == "+", -1, 1) * xrange_size/50, 
            y0 = yrange[2] + (1 + 0.1) * yrange_size + 0.05 * yrange_size * sin(-pi/2), 
            col = col, lwd = 2 * lwd, length = arrow_ln)
    }
    
    
    # Annotations on top
    if (!is.null(annot) && sum(annot[, 1] %in% chrs > 0)) {
        sel = annot[, 1] %in% chrs
        segments(x0 = chr_cum_lns[as.character(annot[sel, 1])] + annot[sel, 2], y0 = 1.1 * 
            yrange_size + yrange[2] + (0:(sum(annot[, 1] %in% chrs) - 1)) * 0.05 * 
            yrange_size, x1 = chr_cum_lns[as.character(annot[sel, 1])] + annot[sel, 
            3], lwd = 2)
        
        text(chr_cum_lns[as.character(annot[sel, 1])] + rowMeans(annot[sel, 2:3]), 
            1.2 * yrange_size + yrange[2] + (0:(sum(annot[, 1] %in% chrs) - 1)) * 
                0.05 * yrange_size, labels = annot[sel, 4], xpd = NA)
    }
    
    
    # Plot point mutations?
    if (!is.null(muts)) {
        rc = c(A = "T", C = "G", G = "C", T = "A")
        muts = muts[muts[, 1] %in% chrs, ]
        muts[, 4] = ifelse(muts[, 3] %in% c("A", "G"), rc[muts[, 4]], muts[, 4])
        muts[, 3] = ifelse(muts[, 3] %in% c("A", "G"), rc[muts[, 3]], muts[, 3])
        mut_col = c(`C>A` = "orange", `C>G` = "brown", `C>T` = "red", `T>A` = "purple", 
            `T>C` = "pink", `T>G` = "lightgreen")
        
        muts = muts[order(muts[, 1], muts[, 2]), ]
        dist_to_next = c(muts[-1, 2] - muts[-nrow(muts), 2])
        dist_to_next[muts[-1, 1] != muts[-nrow(muts), 1]] = NA
        dist_to_next = c(dist_to_next, NA)
        
        xlim = par("usr")[1:2]
        par(new = T)
        pos = chr_cum_lns[as.character(muts[, 1])] + muts[, 2]
        plot(pos, dist_to_next, xlim = xlim, xaxs = "i", axes = F, xlab = "", ylab = "", 
            main = "", yaxt = "n", xaxt = "n", log = "y", col = mut_col[paste(muts[, 
                3], muts[, 4], sep = ">")], pch = 16, ylim = c(1, 1e+07))
        axis(4)
        
        legend("topleft", bty = "n", pch = 16, col = mut_col, legend = names(mut_col), 
            ncol = 2)
    }
    
    
    # Now plot the B-allele frequency
    par(mar = c(0, 4, 0, 4) + 0.1)
    plot(c(), ylim = c(-0.6, 1), xlim = xlim, bty = "n", yaxt = "n", xaxt = "n", 
        xlab = "", ylab = "", yaxs = "i", xaxs = "i")
    abline(h = seq(0, 1, by = 0.25), lwd = lwd, col = "lightgrey")
    
    for (c in chrs) {
        if (!is.null(allele_freqs)) {
            # When plotting only one chromosome with zoomed-in image, only plot the data that
            # will be visible in the graph.
            if (length(chrs) == 1 && !all(xlim == c(1, chr_lens[chrs]))) {
                sel = allele_freqs[, 1] == c & allele_freqs[, 2] > xlim[1] - 1e+05 & 
                  allele_freqs[, 3] < xlim[2] + 1e+05
                x = allele_freqs[sel, 3]
                y = allele_freqs[sel, 4]
                points(x = x, y = y, pch = 16, cex = cn_cex)
            } else {
                sel = allele_freqs[, 1] == c
                points(x = chr_cum_lns[c] + allele_freqs[sel, 3], y = allele_freqs[sel, 
                  4], pch = 16, cex = cn_cex)
            }
        }
        if (!is.null(allele_freq_segments)) {
            # When plotting only one chromosome with zoomed-in image, only plot the data that
            # will be visible in the graph.
            if (length(chrs) == 1 && !all(xlim == c(1, chr_lens[chrs]))) {
                sel = allele_freq_segments[, 1] == c & allele_freq_segments[, 3] > 
                  xlim[1] & allele_freq_segments[, 2] < xlim[2]
                x0 = allele_freq_segments[sel, 2]
                x1 = allele_freq_segments[sel, 3]
                y = allele_freq_segments[sel, 4]
                segments(x0, y, x1, col = ifelse(x0 == x1, "chartreuse2", "blue"), 
                  lwd = 2)
            } else {
                sel = allele_freq_segments[, 1] == c
                x0 = chr_cum_lns[c] + allele_freq_segments[sel, 2]
                y0 = allele_freq_segments[sel, 4]
                x1 = chr_cum_lns[c] + allele_freq_segments[sel, 3]
                segments(x0 = x0, y0 = y0, x1 = x1, col = ifelse(x0 == x1, "chartreuse2", 
                  "blue"), lwd = 2)
            }
        }
    }
    
    
    # X axis names and ticks
    par(mgp = par("mgp") + c(0, 1, 0))
    if (all(xlim == c(1, sum(chr_lens[chrs])))) {
        axis(1, at = (cumsum(chr_lens[chrs]) + chr_cum_lns)/2, labels = ifelse(rep(xaxt, 
            length(chrs)), paste("chr", chrs, " position (Mb)", sep = ""), chrs), 
            tick = F, cex.lab = 1.5)
    } else if (!is.null(chr_lim)) {
        axis(1, at = mean(xlim), labels = paste("chr", chr_lim, " position (Mb)", 
            sep = ""), tick = F, cex.lab = 1.5)
    } else {
        axis(1, at = mean(xlim), labels = paste("chr", chrs, " position (Mb)", sep = ""), 
            tick = F, cex.lab = 1.5)
    }
    par(mgp = par("mgp") - c(0, 1, 0))
    
    if (xaxt) {
        if (length(chrs) > 1) {
            if (all(xlim == c(1, sum(chr_lens[chrs])))) {
                for (c in chrs) {
                  pretty_ticks = pretty(c(1, chr_lens[c]))
                  pretty_ticks = pretty_ticks[which(pretty_ticks < chr_lens[c])]
                  axis(1, at = pretty_ticks + chr_cum_lns[c], labels = pretty_ticks/1e+06)
                }
            } else {
                if (is.null(chr_lim) || !(chr_lim %in% names(chr_lens))) {
                  stop()
                }
                
                pretty_ticks = pretty(xlim - chr_cum_lns[chr_lim])
                
                axis(1, at = pretty_ticks + chr_cum_lns[chr_lim], labels = pretty_ticks/1e+06)
            }
        } else {
            if (all(xlim == c(1, sum(chr_lens[chrs])))) {
                axis(1, at = axisTicks(usr = c(1, chr_lens[chrs]), log = F), labels = axisTicks(usr = c(1, 
                  chr_lens[chrs]), log = F)/1e+06)
            } else {
                pretty_ticks = pretty(xlim)
                axis(1, at = pretty_ticks, labels = pretty_ticks/1e+06)
            }
        }
    }
    
    
    # Finally the ideogram
    if (xlim[2] - xlim[1] < 1e+07) {
        print("Ideogram plotting disabled because xlim[2] - xlim[1] < 10e6")
        
        ideogram = F
    }
    if (ideogram) {
        for (c in chrs) {
            paintCytobands(c, pos = c(1 + chr_cum_lns[c], -0.1), units = "bases", 
                width = 0.5, length.out = chr_lens[c], legend = F, xpd = NA)
        }
    }
}
