## internal functions
###########################
bc_plot_draw_pair <- function(d, log_coord, highlight, count_marks, alpha) {

    sample_name_x <- names(d)[2]
    sample_name_y <- names(d)[3]

    d[, 2:3] <- d[, 2:3] + 1

    g <- ggplot(d) + 
        aes_string(x=sample_name_x, y=sample_name_y) +
        geom_point() + theme_bw() + guides(color="none", size="none", alpha="none") +
        labs( 
            x = paste0("count.", sample_name_x, " + 1"),
            y = paste0("count.", sample_name_y, " + 1")) + aes(alpha=alpha)

    if (!is.null(highlight)) {
        boolColors <- c("TRUE"="red", "FALSE"="black")
        boolScale <- scale_colour_manual(name="myboolean", values=boolColors)
        g <- g + aes(color=barcode_seq %in% highlight) + boolScale
    }

    if (!is.null(count_marks)) {
        count_marks <- rep_len(count_marks, 2)
        count_left_top <- sum(d[, 2] < count_marks[1] & d[, 3] >= count_marks[2])
        count_left_bottom <- sum(d[, 2] < count_marks[1] & d[, 3] < count_marks[2])
        count_right_top <- sum(d[, 2] >= count_marks[1] & d[, 3] >= count_marks[2])
        count_right_bottom <- sum(d[, 2] >= count_marks[1] & d[, 3] < count_marks[2])
        g <- g + 
            geom_hline(yintercept=count_marks[2], alpha=0.3, color='red') + 
            geom_vline(xintercept=count_marks[1], alpha=0.3, color='red') +
            annotate(geom="text", x=count_marks[1], y=count_marks[2],
                label=count_left_bottom, hjust=1+1/nchar(count_left_bottom),
                vjust=1+1) +
            annotate(geom="text", x=count_marks[1], y=count_marks[2],
                label=count_left_top, hjust=1+1/nchar(count_left_top),
                vjust=0-1) +
            annotate(geom="text", x=count_marks[1], y=count_marks[2],
                label=count_right_bottom, hjust=0-1/nchar(count_right_bottom),
                vjust=1+1) +
            annotate(geom="text", x=count_marks[1], y=count_marks[2],
                label=count_right_top, hjust=0-1/nchar(count_right_top),
                vjust=0-1)
    }

    if (log_coord) {
        g <- g + scale_y_log10() + scale_x_log10()
    }
    g
}


bc_plot_draw_single <- function(d, log_coord, highlight, count_marks, alpha) {
    # d: data.frame with two columns, barcode_seq and count
    #
    g <- ggplot(d) +
        geom_density(aes(x=count), fill=NA, color='black')

    if (log_coord)
        g <- g + scale_x_log10()

    if (!is.null(highlight)) {
        boolColors <- c("TRUE"="red", "FALSE"="black")
        boolScale <- scale_colour_manual(name="myboolean", values=boolColors)
        g <- g + aes(color=barcode_seq %in% highlight) + boolScale 
    }

    if (!is.null(count_marks)) {
        count_left <- sum(d$count < count_marks)
        count_right <- sum(d$count >= count_marks)
        g <- g +
            geom_vline(xintercept=count_marks, alpha=0.3, color='red') +
            annotate(geom="text", x=count_marks, y=0, label=count_left,
                hjust=1+1/nchar(count_left)) +
            annotate(geom="text", x=count_marks, y=0, label=count_right,
                hjust=0-1/nchar(count_right))
    }
    g <- g + geom_rug(aes(x=count, y=0), position=position_jitter(height=0),
        alpha=alpha, sides='b') +
        labs(y="density") + guides(color="none", size="none", alpha="none")
    g
}

## bc_plot*
###########################

#' @rdname bc_plot_mutual
#' @exportMethod bc_plot_mutual
setMethod("bc_plot_mutual", c("BarcodeObj"), function(
    barcodeObj,
    count_marks=NULL,
    highlight=NULL,
    log_coord=TRUE,
    alpha=0.7
    ) {

    sample_names <- bc_names(barcodeObj)

    if (is.null(barcodeObj@cleanBc))
        stop("No cleanBc available please run bc_cure* first.")

    if (!is.null(count_marks) & length(count_marks) < length(sample_names)) {
        count_marks <- rep_len(count_marks, length(sample_names))
        message("-message----\nbc_plot_mutual: count_marks is less than sample number, it will be repeatedly used.\n------------")
    }

    if (length(sample_names) < 2) {
        stop("The input BarcodeObj only contains less than 2 samples.")
    } else {
        g_list <- apply(utils::combn(seq_along(sample_names), 2), 2, function(x_i) {
            d <- bc_2dt(bc_subset(barcodeObj, sample = x_i))
            d <- data.table::dcast(d, barcode_seq ~ sample_name, sep = "_", value.var = "count")
            names(d) <- make.names(names(d))
            d[is.na(d)] <- 0
            bc_plot_draw_pair(d, log_coord, highlight, count_marks[x_i], alpha)
        })

    }
    egg::ggarrange(plots=g_list)
})

#' @rdname bc_plot_single
#' @exportMethod bc_plot_single
setMethod("bc_plot_single", c("BarcodeObj"), function(
    barcodeObj, 
    sample_names=bc_names(barcodeObj), 
    count_marks=NULL, 
    highlight=NULL,
    log_coord=TRUE,
    alpha=0.7
    ) {

    if (!is.null(count_marks) & length(count_marks) < length(sample_names)) {
        count_marks <- rep_len(count_marks, length(sample_names))
        message("-message----\nbc_plot_single: count_marks is less than sample number, it will be repeatedly used.\n------------")
    }

    g_list <- lapply(seq_along(sample_names), function(i) {
        d <- bc_2df(barcodeObj[, sample_names[i]])
        g <- bc_plot_draw_single(d, log_coord=log_coord, highlight=highlight,
            count_marks=count_marks[i], alpha=alpha)
        g + labs(title=sample_names[i])
    })

    egg::ggarrange(plots=g_list)
})

#' @rdname bc_plot_pair
#' @exportMethod bc_plot_pair
setMethod("bc_plot_pair", c("BarcodeObj"), function(
    barcodeObj,
    sample_x,
    sample_y,
    count_marks_x=NULL,
    count_marks_y=count_marks_x,
    highlight=NULL,
    log_coord=TRUE,
    alpha=0.7
    ) {
    if (is.null(barcodeObj@cleanBc))
        stop("No cleanBc available please run bc_cure* first.")

    if (length(sample_y) > length(sample_x)) {
        message("-message----\nbc_plot_pair: sample_x is shorter than sample_y, sample_x will be repeatedly used.\n------------")
        sample_x <- rep_len(sample_x, length(sample_y))
    } 

    if (length(sample_x) > length(sample_y)) {
        message("-message----\nbc_plot_pair: sample_y is shorter than sample_x, sample_y will be repeatedly used.\n------------")
        sample_y <- rep_len(sample_y, length(sample_x))
    }

    if (!is.null(count_marks_x)) {
        if (length(count_marks_x) < length(sample_x)) {
            count_marks_x <- rep_len(count_marks_x, length(sample_x))
            message("-message----\nbc_plot_pair: count_marks_x length is less than sample number, it will be repeatedly used.\n------------")
        }

        if (length(count_marks_y) < length(sample_y)) {
            count_marks_y <- rep_len(count_marks_y, length(sample_y))
            message("-message----\nbc_plot_pair: count_marks_y length is less than sample number, it will be repeatedly used.\n------------")
        }
    }

    g_list <- lapply(seq_along(sample_x), function(i) {
        sample_xy <- c(sample_x[i], sample_y[i])
        bc_sub <- barcodeObj[, sample_xy]
        d <- bc_2dt(bc_sub)
        d <- data.table::dcast(d, barcode_seq ~ sample_name, value.var = "count")
        d <- d[, c("barcode_seq", bc_names(bc_sub)), with=FALSE]
        names(d) <- make.names(names(d))
        d[is.na(d)] <- 0
        bc_plot_draw_pair(d, log_coord, highlight, c(count_marks_x[i],
                count_marks_y[i]), alpha)
    })

    egg::ggarrange(plots=g_list)
})



