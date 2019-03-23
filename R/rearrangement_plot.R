#' Function for plotting campbellgrams
#'
#' Copyright: 2019 Yilong Li <yilong.li.yl3@gmail.com>


# Colours for arcs
ARC_COLS = c(
    td_col = "darkorange4",
    del_col = "darkslateblue",
    tail_tail_col = "cadetblue4",
    head_head_col = "chartreuse2"
)


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
compute_cn_yrange = function(cn_bedgraph) {
    if (is.null(cn_bedgraph)) {
        # No copy number data for computing cn_yrange - just return default
        # values.
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


#' Helper function for creating a plot outline for a campbellgram.
#'
#' Create an empty canvas to which all the SV and CN etc. data can be added.
#'
#' @param cn_yrange Range for copy number to be displayed.
#' @param xlim X-axis limits in user (not necessarily chromosomal) coordinates.
#' @param main Title for the plot.
#' @param chr_lens named vector of chromosome lengths.
#'
#' @export
make_campbellgram_outline = function(cn_yrange, xlim, main, chr_lens) {
    par(mar = c(5, 4, 2, 2) + .1)

    cn_yrange_size = cn_yrange[2] - cn_yrange[1]
    ideogram_width = cn_yrange_size / 10  # Space allocated for an ideogram.
    SV_region_width = 1.4 * cn_yrange_size  # Space allocated for SV arcs.

    # Add space in the plot for SV junctions and the ideogram
    plot_yrange = c(cn_yrange[1] - ideogram_width,
                    cn_yrange[2] + SV_region_width)
    plot(c(), ylim = plot_yrange, xlim = xlim, bty = "n", yaxt = "n",
         xaxt = "n", xlab = "", ylab = "", yaxs = "i", xaxs = "i", main = main)

    # Add a shaded grid for the CN section.
    xmin = par("usr")[1]
    xmax = par("usr")[2]
    odd_cn_col = rgb(0.1, 0.1, 0.1, 0.05)
    even_cn_col = rgb(0.1, 0.1, 0.1, 0.1)
    for (i in floor(cn_yrange[1]):floor(cn_yrange[2])) {
        polygon(
            x = c(xmin, xmax, xmax, xmin),
            y = i - 0.5 + c(0, 0, 1, 1),
            col = ifelse(i %% 2 == 0, even_cn_col, odd_cn_col),
            lty = 0
        )
    }

    # Add lines to separate chromosomes.
    if (length(chrs_shown) > 1) {
        segments(
            x0 = cumsum(chr_lens[chrs_shown])[-length(chrs_shown)],
            y0 = yrange[1] - 0.5,
            y1 = yrange[2] + 0.5
        )
    }

    # Add dotted lines where the SV arcs emanate from.
    segments(
        x0 = c(1, 1),
        x1 = rep(sum(chr_lens[chrs_shown]), 2),
        y0 = cn_yrange[2] + c(0.3, 0.75) * cn_yrange_size,
        lty = 3
    )
    abline(h = cn_yrange[2] + 0.5)
}


#' Helper function for adding CN data to a campbellgram.
#'
#' @param cn_win_size window size in base pairs for copy number rolling window
#'   averaging.
#' @param chrs_shown chromosomes for which the CN data are to be plotted.
#' @param plotted_cn copy number to be plotted in a bedGraph format data frame.
#' @param cn_yrange Y-axis range for copy number.
#' @param cn_cex marker size (see \code{?par}) for CN dots.
#' @param chr_coord_offset a named vector of X-axis offsets for each chromosome.
#'   Names must include chromosomes in \emph{chrs_shown}.
#'
#' @export
add_cns_to_campbellgram = function(chrs_shown, cn_win_size, plotted_cn,
                                   cn_yrange, cn_cex, chr_coord_offset) {
    is_zoomed = (length(chrs_shown) == 1
                 && !all(xlim == c(1, chr_lens[chrs_shown])))
    if (is_zoomed) {
        # Plotting a single zoomed chromosome.
        idx = (plotted_cn[, 1] == chrom
               & plotted_cn[, 2] >= xlim[1]
               & plotted_cn[, 3] <= xlim[2])
        average_cn_by_window = window_means(
            rowMeans(plotted_cn[idx, 2:3]),
            plotted_cn[idx, 4],
            xlim[1] - cn_win_size,
            xlim[2] + cn_win_size,
            cn_win_size
        )
        x = average_cn_by_window$window_coords
        y = average_cn_by_window$window_means
        y = ifelse(y < cn_yrange | y > cn_yrange, NA, y)
        points(x, y, pch = 16, cex = cn_cex)
    } else {
        # Plotting one or multiple entire chromosomes.
        for (chrom in chrs_shown) {
            idx = plotted_cn[, 1] == chrom
            average_cn_by_window = window_means(
                rowMeans(plotted_cn[idx, 2:3]),
                plotted_cn[idx, 4],
                1,
                chr_lens[chrom],
                cn_win_size
            )
            x = chr_coord_offset[chrom] + average_cn_by_window$window_coords
            y = average_cn_by_window$window_means
            y = ifelse(y < cn_yrange | y > cn_yrange, NA, y)
            points(x, y, pch = 16, cex = cn_cex)
        }
    }

    # Add Y-axis label and Y-axis ticks.
    axis(2, at = axisTicks(usr = cn_yrange, log = F), las = 2)
    mtext("Copy number", side = 2, line = 1, at = mean(cn_yrange))
}


#' Helper function for adding CN segments to a campbellgram.
#'
#' @param segments absolute copy number segments in bedGraph format.
#' @param chrs_shown chromosomes for which the CN data are to be plotted.
#' @param cn_yrange Y-axis range for copy number.
#' @param chr_coord_offset a named vector of X-axis offsets for each chromosome.
#'   Names must include chromosomes in \emph{chrs_shown}.
#' @param lwd line width.
#' @param col line color.
#'
#' @export
add_cn_segs_to_campbellgram = function(segments, chrs_shown, cn_yrange,
                                       chr_coord_offset, lwd,
                                       col = "mediumslateblue") {
    for (chrom in chrs_shown) {
        idx = (segments[, 1] == chrs_shown
               & segments[, 4] >= cn_yrange[1]
               & segments[, 4] <= cn_yrange[2])
        if (any(idx)) {
            segments(
                x0 = chr_coord_offset[chrom] + segments[idx, 2] + 1,
                x1 = chr_coord_offset[chrom] + segments[idx, 3],
                y0 = segments[idx, 4],
                col = col
            )
        }
    }
}


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
#' @param cn_yrange Y-axis range for copy number.
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
                        cn_yrange = NULL, ideogram = NULL, cn_win_size = NULL,
                        arc_cols = ARC_COLS, cn_cex = 0.3, lwd = 0.75,
                        arrow_ln = 0.15, annot = NULL, main = NULL, muts = NULL,
                        allele_freqs = NULL, allele_freq_segments = NULL,
                        xaxt = T) {
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
    } else {
        plotted_cn = NULL
    }
    if (is.null(cn_yrange)) {
        cn_yrange = compute_cn_yrange(plotted_cn)
        cn_yrange_size = cn_yrange[2] - cn_yrange[1]
    }

    # Set CN window size dynamically
    if (is.null(cn_win_size)) {
        cn_win_size = set_cn_win_size(xlim)
    }

    xrange_size = cumsum(chr_lens[chrs_used])  # TODO(yl3): what is this?

    # Create the campbellgram plot outline.
    make_campbellgram_outline(cn_yrange, xlim, main, chr_lens)

    # Plot CN data points and segments.
    if (!is.null(plotted_cn)) {
        add_cns_to_campbellgram(cn_win_size, plotted_cn, cn_yrange, cn_cex,
                                chr_coord_offset)
    }
    if (!is.null(segments)) {
        add_cn_segs_to_campbellgram(segments, chrs_shown, cn_yrange,
                                    chr_coord_offset, lwd)
    }

    # Plot rearrangements: First dotted lines
    # Then 'intra-chromosomal' translocations
    sel = bedpe[, 1] %in% chrs & bedpe[, 4] %in% chrs
    col = ifelse(bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col, ifelse(bedpe[sel,
        9] == "+" & bedpe[sel, 10] == "-", del_col, ifelse(bedpe[sel, 9] == "-" &
        bedpe[sel, 10] == "+", td_col, tail_tail_col)))

    if (sum(sel) > 0) {
        arc(x0 = chr_coord_offset[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
            x1 = chr_coord_offset[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
            y = ifelse(bedpe[sel, 9] == bedpe[sel, 10], yrange[2] + 0.3 * yrange_size,
                yrange[2] + 0.75 * yrange_size), yr = ifelse(bedpe[sel, 10] == "-",
                1, -1) * 0.2 * yrange_size, col = col, lwd = lwd)
        segments(x0 = chr_coord_offset[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel,
            2:3]), y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75 * yrange_size,
            0.3 * yrange_size), y1 = yrange[1], col = rgb(t(col2rgb(col)), alpha = 127,
            max = 255), lwd = lwd)
        segments(x0 = chr_coord_offset[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel,
            5:6]), y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75 * yrange_size,
            0.3 * yrange_size), y1 = yrange[1], col = rgb(t(col2rgb(col)), alpha = 127,
            max = 255), lwd = lwd)
    }

    # Then rearrangments where low end %in% chrs and !(high end %in% chrs)
    sel = bedpe[, 1] %in% chrs & !(bedpe[, 4] %in% chrs)
    if (sum(sel) > 0) {
        # arrows(
        segments(x0 = chr_coord_offset[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel,
            2:3]), y0 = yrange[1], y1 = yrange[2] + yrange_size, col = rgb(t(col2rgb("black")),
            alpha = 127, max = 255), lwd = lwd)
        segments(x0 = chr_coord_offset[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel,
            2:3]), y0 = yrange[2] + yrange_size, x1 = chr_coord_offset[as.character(bedpe[sel,
            1])] + rowMeans(bedpe[sel, 2:3]) + ifelse(bedpe[sel, 9] == "+", -1, +1) *
            (par("usr")[2] - par("usr")[1])/50, y1 = ifelse(bedpe[sel, 9] == "+",
            yrange[2] + 1.1 * yrange_size, yrange[2] + 1.3 * yrange_size), xpd = NA)
        text(as.character(bedpe[sel, 4]), x = chr_coord_offset[as.character(bedpe[sel,
            1])] + rowMeans(bedpe[sel, 2:3]) + ifelse(bedpe[sel, 9] == "+", -1, +1) *
            (par("usr")[2] - par("usr")[1])/50, y = yrange[2] + ifelse(bedpe[sel,
            9] == "+", 1.2, 1.4) * yrange_size, xpd = NA)
    }

    # Then rearrangments where high end %in% chrs and !(low end %in% chrs)
    sel = !(bedpe[, 1] %in% chrs) & bedpe[, 4] %in% chrs
    if (sum(sel) > 0) {
        # arrows(
        segments(x0 = chr_coord_offset[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel,
            5:6]), y0 = yrange[1], y1 = yrange[2] + yrange_size, col = rgb(t(col2rgb("black")),
            alpha = 127, max = 255), lwd = lwd)
        segments(x0 = chr_coord_offset[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel,
            5:6]), y0 = yrange[2] + yrange_size, x1 = chr_coord_offset[as.character(bedpe[sel,
            4])] + rowMeans(bedpe[sel, 5:6]) + ifelse(bedpe[sel, 10] == "+", -1,
            +1) * (par("usr")[2] - par("usr")[1])/50, y1 = ifelse(bedpe[sel, 10] ==
            "+", yrange[2] + 1.1 * yrange_size, yrange[2] + 1.3 * yrange_size), xpd = NA)
        text(as.character(bedpe[sel, 1]), x = chr_coord_offset[as.character(bedpe[sel,
            4])] + rowMeans(bedpe[sel, 5:6]) + ifelse(bedpe[sel, 10] == "+", -1,
            +1) * (par("usr")[2] - par("usr")[1])/50, y = yrange[2] + ifelse(bedpe[sel,
            10] == "+", 1.2, 1.4) * yrange_size, xpd = NA)
    }


    # Annotations on top
    if (!is.null(annot) && sum(annot[, 1] %in% chrs > 0)) {
        sel = annot[, 1] %in% chrs
        segments(x0 = chr_coord_offset[as.character(annot[sel, 1])] + annot[sel, 2], y0 = 1.1 *
            yrange_size + yrange[2] + (0:(sum(annot[, 1] %in% chrs) - 1)) * 0.05 *
            yrange_size, x1 = chr_coord_offset[as.character(annot[sel, 1])] + annot[sel,
            3], lwd = 2)

        text(chr_coord_offset[as.character(annot[sel, 1])] + rowMeans(annot[sel, 2:3]),
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
        pos = chr_coord_offset[as.character(muts[, 1])] + muts[, 2]
        plot(pos, dist_to_next, xlim = xlim, xaxs = "i", axes = F, xlab = "", ylab = "",
            main = "", yaxt = "n", xaxt = "n", log = "y", col = mut_col[paste(muts[,
                3], muts[, 4], sep = ">")], pch = 16, ylim = c(1, 1e+07))
        axis(4)

        legend("topleft", bty = "n", pch = 16, col = mut_col, legend = names(mut_col),
            ncol = 2)
    }


    # X axis names and ticks
    par(mgp = par("mgp") + c(0, 1, 0))
    if (all(xlim == c(1, sum(chr_lens[chrs])))) {
        axis(1, at = (cumsum(chr_lens[chrs]) + chr_coord_offset)/2, labels = ifelse(rep(xaxt,
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
                  axis(1, at = pretty_ticks + chr_coord_offset[c], labels = pretty_ticks/1e+06)
                }
            } else {
                if (is.null(chr_lim) || !(chr_lim %in% names(chr_lens))) {
                  stop()
                }

                pretty_ticks = pretty(xlim - chr_coord_offset[chr_lim])

                axis(1, at = pretty_ticks + chr_coord_offset[chr_lim], labels = pretty_ticks/1e+06)
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
            quantsmooth::paintCytobands(c, pos = c(1 + chr_coord_offset[c], -0.1), units = "bases",
                width = 0.5, length.out = chr_lens[c], legend = F, xpd = NA)
        }
    }
}
