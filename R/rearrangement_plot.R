#' Function for plotting campbellgrams
#'
#' Copyright: 2019 Yilong Li <yilong.li.yl3@gmail.com>


# Colours for arcs
ARC_COLS = c(
    td_col = "darkorange4",
    del_col = "darkslateblue",
    head_head_col = "chartreuse2",
    tail_tail_col = "cadetblue4"
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
        as.list(data.frame(x, y, xr, yr, col, lwd, stringsAsFactors = F))
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


#' @describeIn campbellgram Helper function for setting chrs_shown based on
#'   provided parameter values.
#'
#' Also performs additional validations. Currently either all chrs_used
#' chromosomes are shown, or one of the chromosomes in chrs_used.
set_chrs_shown = function(chrs_used, chrs_shown) {
    if (length(chrs_shown) > 1) {
        stop("Currently parameter 'chrs_shown' must be of length 1.")
    }
    if (is.null(chrs_shown)) {
        return(chrs_used)
    } else {
        if (!(chrs_shown %in% chrs_used)) {
            stop("'chrs_shown' is not in 'chrs_used'.")
        }
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
        return(c(1, sum(chr_lens[chrs_shown])))
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
#' @param chrs_shown chromosomes to be shown.
#'
#' @export
make_campbellgram_outline = function(cn_yrange, xlim, main, chr_lens,
                                     chrs_shown) {
    par(mar = c(5, 4, 4, 4) + .1)
    cn_yrange_size = cn_yrange[2] - cn_yrange[1]

    # Some HARDCODED variables. See also function draw_sv_arc().
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
            y0 = cn_yrange[1] - 0.5,
            y1 = cn_yrange[2] + 0.5
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

    # Add Y-axis label and Y-axis ticks.
    axis(2, at = axisTicks(usr = cn_yrange, log = F), las = 2)
    mtext("Copy number", side = 2, line = 2, at = mean(cn_yrange))
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


#' Helper function for colouring SVs based on breakpoint orientations.
#'
#' @param low_end_dir Breakpoint direction for SV's low end, "+" or "-".
#' @param high_end_dir Breakpoint direction for SV's high end, "+" or "-".
#' @param arc_cols Arc colours - a vector with colour values for keys "td_col",
#'   "del_col", "tail_tail_col" and "head_head_col".
choose_sv_colour = function(low_end_dir, high_end_dir, arc_cols) {
    sapply(
        paste(low_end_dir, high_end_dir, sep = ""),
        function(value) {
            switch(
                value,
                "++" = arc_cols["tail_tail_col"],
                "+-" = arc_cols["del_col"],
                "-+" = arc_cols["td_col"],
                "--" = arc_cols["head_head_col"]
            )
        }
    )
}


#' Add alpha to a vector of colours.
add_alpha = function(col, alpha, max = 255) {
    rgb(t(col2rgb(col)), alpha = alpha, max = max)
}


#' Helper function for drawing different types of SV arcs.
#'
#' Also draws vertical segments for each SV arc.
draw_sv_arc = function(sv_bedpe, chr_coord_offset, arc_cols, cn_yrange, lwd) {
    # Some HARDCODED variables. See also function make_campbellgram_outline().
    inv_sv_pos_multiplier = 0.3
    noninv_sv_pos_multiplier = 0.75
    arc_radius_multiplier = 0.2

    cn_yrange_size = diff(cn_yrange)

    # Extract some values from the SV BEDPE table.
    chrom_l = sv_bedpe[, 1]
    chrom_h = sv_bedpe[, 4]
    bkpt_pos_l = rowMeans(sv_bedpe[, 2:3])
    bkpt_pos_h = rowMeans(sv_bedpe[, 5:6])
    dir_l = sv_bedpe[, 9]
    dir_h = sv_bedpe[, 10]
    sv_col = choose_sv_colour(dir_l, dir_h, arc_cols)
    arc_y_pos = ifelse(dir_l == dir_h,
                       cn_yrange[2] + inv_sv_pos_multiplier * cn_yrange_size,
                       cn_yrange[2] + noninv_sv_pos_multiplier * cn_yrange_size)
    y_radius =
        ifelse(dir_h == "-", 1, -1) * arc_radius_multiplier * cn_yrange_size

    # Draw arcs.
    draw_arcs(
        x0 = chr_coord_offset[chrom_l] + bkpt_pos_l,
        x1 = chr_coord_offset[chrom_h] + bkpt_pos_h,
        y = arc_y_pos,
        yr = y_radius,
        col = sv_col,
        lwd = lwd
    )

    # Draw segments going into the CN region.
    segments(
        x0 = c(
            chr_coord_offset[chrom_l] + bkpt_pos_l,
            chr_coord_offset[chrom_h] + bkpt_pos_h
        ),
        y0 = rep(arc_y_pos, 2),
        y1 = cn_yrange[1],
        col = add_alpha(rep(sv_col, 2), 127),
        lwd = lwd
    )
}


#' Helper function for drawing intra-chromosomal SV arcs.
draw_intrachr_sv_arc = function(sv_bedpe, chrs_used, chr_coord_offset, arc_cols,
                                cn_yrange, lwd) {
    idx = (sv_bedpe[, 1] %in% chrs_used
           & sv_bedpe[, 4] %in% chrs_used)
    if (sum(idx) == 0) {
        return()
    }
    sv_bedpe = sv_bedpe[idx, ]
    draw_sv_arc(sv_bedpe, chr_coord_offset, arc_cols, cn_yrange, lwd)
}


#' Draw inter-chromosomal rearrangements at specified locations.
draw_interchr_sv_arc = function(chr_coord_offset, sv_chr, sv_chr_pos, sv_dir,
                                cn_yrange, lwd) {
    cn_yrange_size = diff(cn_yrange)
    x = chr_coord_offset[sv_chr] + sv_chr_pos
    x_delta = ifelse(sv_dir == "+", -1, +1) * diff(par("usr")[1:2]) / 50
    seg_bot = cn_yrange[1]
    seg_top = cn_yrange[2] + cn_yrange_size
    y_delta = ifelse(sv_dir == "+", 0.1, 0.3) * cn_yrange_size
    text_y_delta = ifelse(sv_dir == "+", 0.2, 0.4) * cn_yrange_size
    transparent_black = rgb(t(col2rgb("black")), alpha = 127, max = 255)
    segments(x0 = x, y0 = seg_bot, y1 = seg_top, col = transparent_black,
             lwd = lwd)
    segments(x0 = x, y0 = seg_top, x1 = x + x_delta, y1 = seg_top + y_delta,
             xpd = T)
    text(x + x_delta, seg_top + text_y_delta, sv_chr)
}


#' Find inter-chromosomal SVs from BEDPE and plot them.
draw_interchr_sv_from_bedpe = function(sv_bedpe, chrs_used,
                                       chr_coord_offset, cn_yrange, lwd) {
    chr_low = sv_bedpe[, 1]
    chr_high = sv_bedpe[, 4]

    # Find the low end inter-chromosomal SVs.
    is_inter_chr_low_end = chr_low %in% chrs_used & !(chr_high %in% chrs_used)
    chrom = chr_low[is_inter_chr_low_end]
    pos = rowMeans(sv_bedpe[is_inter_chr_low_end, 2:3])
    dir = sv_bedpe[is_inter_chr_low_end, 9]

    # Append the additional high end inter-chromosomal SVs.
    is_inter_chr_high_end = chr_high %in% chrs_used & !(chr_low %in% chrs_used)
    chrom = c(chrom, chr_high[is_inter_chr_high_end])
    pos = c(pos, rowMeans(sv_bedpe[is_inter_chr_high_end, 5:6]))
    dir = c(dir, sv_bedpe[is_inter_chr_high_end, 10])

    # Draw the inter-chromosomal lines.
    draw_interchr_sv_arc(chr_coord_offset, chrom, pos, dir, cn_yrange, lwd)
}


#' Add annotations to the top of campbellgram.
#'
#' Draws line segments and adds a text annotation for each BED entry.
#'
#' @param annot an annotation data frame in BED format.
#' @param chrs_used see ?campbellgram.
#' @param chr_coord_offset see ?campbellgram.
#' @param cn_yrange Range for copy number to be displayed.
add_annotations = function(annot, chrs_used, chr_coord_offset, cn_yrange) {
    if (is.null(annot) || sum(annot[, 1] %in% chrs_used) == 0) {
        return()
    }
    yrange_size = cn_yrange[2] - cn_yrange[1]
    chrom = as.character(annot[, 1])
    idx = chrom %in% chrs_used
    x_start = chr_coord_offset[chrom[idx]] + annot[idx, 2]
    x_end = chr_coord_offset[chrom[idx]] + annot[idx, 3]

    # The Y-axis coordinate of the annotation line segments. This depends on a
    # HARDCODED multiplier. See also make_campbellgram_outline()
    OFFSET = 1.1
    segs_y_coords = (cn_yrange[2]  +  (1 + OFFSET) * yrange_size
                     + 0.05 * (0:sum(idx) - 1) * yrange_size)
    text_y_coords = (cn_yrange[2]  +  (1 + 2 * OFFSET) * yrange_size
                     + 0.05 * (0:sum(idx) - 1) * yrange_size)
    segments(x0 = x_start, x1 = x_end, y0 = segs_y_coords, lwd = 2)
    text(x_start/2 + x_end/2, text_y_coords, labels = annot[idx, 4], xpd = T)
}


#' @describeIn campbellgram Plot a rainfall plot of mutation.
#'
#' Overlays a rainfall plot on top of a campbellgram.
#'
#' Every single mutation is displayed on the plot. The mutation distance
#' (Y-axis) of each mutation is computed as the \emph{harmonic mean} of the
#' distances to the previous and next mutation. Note that this differs from the
#' original rainfall plot version, which plotted n - 1 mutations, with the
#' distance computed as each mutation's distance to its next mutation.
#'
#' In the original rainfall plot version, if a cluster had k mutations, only
#' k - 1 mutations would be visualised in the cluster, since the k'th mutation
#' would have a large distance to the next mutation. With a harmonic mean, which
#' emphasises its smallest members, all k mutations would be shown.
add_mutation_rainfall_plot = function(muts, chrs_used, chr_coord_offset,
                                      mut_col = NULL, mut_ylim = NULL) {

    # Pairwise harmonic mean for vectors of same length.
    pairwise_harmonic_mean = function(x, y) {
        ifelse(is.na(x) & is.na(y), NA,
        ifelse(is.na(x), y,
        ifelse(is.na(y), x,
               2 / (1 / x + 1 / y))))
    }

    # Set up some basic variables.
    complement_of = c(A = "T", C = "G", G = "C", T = "A")
    if (is.null(mut_col)) {
        mut_col = c("C>A" = "orange", "C>G" = "brown", "C>T" = "red",
                    "T>A" = "purple", "T>C" = "pink", "T>G" = "lightgreen")
    }

    # Change all original nucleotides to the pyrimidine strand & prepare muts
    # otherwise.
    chrom = muts[, 1]
    idx = chrom %in% chrs_used
    muts = muts[idx, ]
    chrom = chrom[idx]
    is_purine = muts[, 3] %in% c("A", "G")
    mut_type = ifelse(
        is_purine,
        paste(complement_of[muts[, 3]], complement_of[muts[, 4]], sep = ">"),
        paste(muts[, 3], muts[, 4], sep = ">"))
    # Order based on chromosomal order in chrs_used.
    # chrom_order = setNames(1:length(chrs_used), chrs_used)
    # muts = muts[order(chrom_order[muts[, 1]], muts[, 2]), ]

    # Compute the coordinates to be plotted.
    pos = chr_coord_offset[chrom] + muts[, 2]
    intermut_dists = pmax(diff(pos), 1)
    # Mask "mutation distances" of mutations between chromosomes.
    intermut_dists[chrom[-1] != chrom[-length(chrom)]] = NA
    mut_dist_harmonic_means = pairwise_harmonic_mean(
        c(intermut_dists, NA),
        c(NA, intermut_dists))

    # Shrink the vertical height of the plot and add the mutations to the plot.
    y_bottom_ndc = grconvertY(0, "user", "ndc")
    y_top_ndc = grconvertY(par("usr")[4] - 0.05 * diff(par("usr")[3:4]), "user",
                           "ndc")
    par(plt = c(par("plt")[1:2], y_bottom_ndc, y_top_ndc), new = T)
    xlim = par("usr")[1:2]

    plot(pos, mut_dist_harmonic_means, xlim = xlim, xaxs = "i", axes = F,
         xlab = "", ylab = "", main = "", yaxt = "n", xaxt = "n", log = "y",
         col = mut_col[mut_type], pch = 16, ylim = mut_ylim, cex = 0.75)
    axis(4)
    mtext("Inter-mutation distance (bp)", 4, line = 2)
}


#' @describeIn add_xlab_and_ticks Helper function for adding the X axis
#'   (chromosome name) labels.
add_xlab = function(xlim, chrs_used, chrs_shown, chr_lens, chr_coord_offset,
                    plot_xaxt) {
    par(mgp = par("mgp") + c(0, 1, 0))
    if (all(xlim == c(1, sum(chr_lens[chrs_used])))) {
        # xlim encompasses all chromosomes in chrs_used.
        label_pos = (cumsum(chr_lens[chrs_used]) + chr_coord_offset) / 2
        labels = ifelse(rep(plot_xaxt, length(chrs_used)),
                        paste("chr", chrs_used, " position (Mb)", sep = ""),
                        chrs_used)
        axis(1, at = label_pos, labels = labels, tick = F, cex.lab = 1.5)
    } else {
        # We are here when xlim does not cover all chromosomes in chrs_used. In
        # this case only one chromosome is shown either through chrs_shown or
        # chrs_used.
        if (!is.null(chrs_shown)) {
            label = paste("chr", chrs_shown, " position (Mb)", sep = "")
        } else {
            label = paste("chr", chrs_used, " position (Mb)", sep = "")
        }
        axis(1, at = mean(xlim), labels = label, tick = F, cex.lab = 1.5)
    }
    par(mgp = par("mgp") - c(0, 1, 0))
}


#' @descrbeIn add_xlab_and_xticks add_xticks
add_xticks = function(xlim, chr_lens, chrs_used, chrs_shown, chr_coord_offset) {
    if (all(xlim == c(1, sum(chr_lens[chrs_used])))) {
        # xlim encompasses all chromosomes in chrs_used.
        if (length(chrs_used) > 1) {
            for (c in chrs_used) {
                pretty_ticks = pretty(c(1, chr_lens[c]))
                pretty_ticks = pretty_ticks[which(pretty_ticks < chr_lens[c])]
                axis(1, at = pretty_ticks + chr_coord_offset[c],
                     labels = pretty_ticks / 1e6)
            }
        } else {
            tick_pos = axisTicks(usr = c(1, chr_lens[chrs_used]), log = F)
            axis(1, at = tick_pos, labels = tick_pos / 1e6)
        }
    } else {
        # We are here when xlim does not cover all chromosomes in chrs_used. In
        # this case only one chromosome is shown either through chrs_shown or
        # chrs_used.
        if (!is.null(chrs_shown)) {
            tick_pos = pretty(xlim - chr_coord_offset[chrs_shown])
            axis(1, at = tick_pos + chr_coord_offset[chrs_shown],
                 labels = tick_pos / 1e6)
        } else {
            tick_pos = pretty(xlim)
            axis(1, at = tick_pos, labels = tick_pos / 1e6)
        }
    }
}


#' @describeIn campbellgram Add X axis ticks.
#'
#' @param xlim the current X-axis coordinate limits.
#' @param chrs_used see ?campbellgram.
#' @param chrs_shown see ?campbellgram.
#' @param plot_xaxt see ?campbellgram.
add_xlab_and_xticks = function(xlim, chrs_used, chrs_shown, chr_lens,
                               chr_coord_offset, plot_xaxt) {
    add_xlab(xlim, chrs_used, chrs_shown, chr_lens, chr_coord_offset, plot_xaxt)
    if (plot_xaxt) {
        add_xticks(xlim, chr_lens, chrs_used, chrs_shown, chr_coord_offset)
    }
}


#' @describeIn campbellgram Helper function for adding an ideogram to a
#'   campbellgram.
add_ideogram = function(ideogram, xlim, chr_lens, chrs_used, cn_yrange,
                        chr_coord_offset) {
    if (xlim[2] - xlim[1] < 1e7) {
        message("Ideogram plotting disabled because xlim[2] - xlim[1] < 1e7")
        ideogram = F
    }
    if (ideogram) {
        cn_yrange_size = cn_yrange[2] - cn_yrange[1]
        # Currently, ideogram_width is a constant that's set to
        # cn_yrange_size / 10. See make_campbellgram_outline().
        ideogram_top = cn_yrange[1] + cn_yrange_size / 100
        ideogram_width = cn_yrange_size * 9 / 100
        for (chrom in chrs_used) {
            ideogram_xy = c(1 + chr_coord_offset[chrom], -ideogram_top)
            quantsmooth::paintCytobands(
                chrom, pos = ideogram_xy, units = "bases",
                width = ideogram_width, length.out = chr_lens[chrom],
                legend = F, xpd = NA)
        }
    }
}


#' Plot a campbellgram.
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
#' @param cn_cex marker size (see \code{?par}) for CN dots.
#' @param lwd line width for all SV arcs and other line segments.
#' @param muts a 'data.frame' object with the following columns: chromosome,
#'   position, base mutated from, base mutated to. Bases must be in
#'   {A, C, G, T}.
#' @param mut_col a vector of >= 6 colours for C>A, C>G, C>T, T>A, T>C and T>G
#'   mutations. Default: c("orange", "brown", "red", "purple", "pink",
#'   "lightgreen").
#' @param mut_ylim Y-axis limits for the rainfall plot.
#' @param plot_xaxt whether to plot X axis ticks or not.
campbellgram = function(ref_genome, bedpe, chrs_used, chrs_shown = NULL,
                        xlim = NULL, cn_bedgraph = NULL, segments = NULL,
                        cn_yrange = NULL, ideogram = T, cn_win_size = NULL,
                        arc_cols = ARC_COLS, cn_cex = 0.3, lwd = 0.75,
                        annot = NULL, main = NULL, muts = NULL, mut_col = NULL,
                        mut_ylim = NULL, plot_xaxt = T) {
    # Internally, the X-axis is composed of all chromosomes in chrs_used stacked
    # side by side. The actual part of the chromosome shown is controlled by
    # chrs_displ and/or xlim.

    # Create chromosome lengths vector. Ensure chromosome names are characters.
    ref_genome[, 1] = as.character(ref_genome[, 1])
    chr_lens = setNames(ref_genome[, 2], ref_genome[, 1])
    chrs_used = set_chrs_used(chrs_used, ref_genome)
    chrs_shown = set_chrs_shown(chrs_used, chrs_shown)
    bedpe[, 1] = as.character(bedpe[, 1])
    bedpe[, 4] = as.character(bedpe[, 4])
    annot[, 1] = as.character(annot[, 1])
    muts[, 1] = as.character(muts[, 1])

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

    # Create the campbellgram plot outline.
    make_campbellgram_outline(cn_yrange, xlim, main, chr_lens, chrs_shown)

    # Plot CN data points and segments.
    if (!is.null(plotted_cn)) {
        add_cns_to_campbellgram(cn_win_size, plotted_cn, cn_yrange, cn_cex,
                                chr_coord_offset)
    }
    if (!is.null(segments)) {
        add_cn_segs_to_campbellgram(segments, chrs_shown, cn_yrange,
                                    chr_coord_offset, lwd)
    }

    # Draw SVs where both breakpoints are in chrs_used.
    draw_intrachr_sv_arc(bedpe, chrs_used, chr_coord_offset, ARC_COLS,
                         cn_yrange, lwd)

    # Then SVs where one of the chromosomes %in% chrs_used.
    draw_interchr_sv_from_bedpe(bedpe, chrs_used, chr_coord_offset,
                                cn_yrange, lwd)

    # Add gene/genomic location annotations on top.
    add_annotations(annot, chrs_used, chr_coord_offset, cn_yrange)

    # Finally add the ideogram.
    add_ideogram(ideogram, xlim, chr_lens, chrs_used, cn_yrange,
                 chr_coord_offset)

    # Add X axis names and ticks.
    add_xlab_and_xticks(xlim, chrs_used, chrs_shown, chr_lens, chr_coord_offset,
                        plot_xaxt)

    # Plot point mutations? Needs to be done last since this changes the Y-axis!
    if (!is.null(muts)) {
        add_mutation_rainfall_plot(muts, chrs_used, chr_coord_offset, mut_col,
                                   mut_ylim)
    }
}
