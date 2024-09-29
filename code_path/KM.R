ggsurvplot2=function (fit, data = NULL, fun = NULL, color = NULL, palette = NULL,
    linetype = 1, conf.int = FALSE, pval = FALSE, pval.method = FALSE,
    test.for.trend = FALSE, surv.median.line = "none", risk.table = FALSE,
    cumevents = FALSE, cumcensor = FALSE, tables.height = 0.25,
    group.by = NULL, facet.by = NULL, add.all = FALSE, combine = FALSE,
    ggtheme = theme_survminer(), tables.theme = ggtheme, ...) {
    if (length(group.by) > 2)
        stop("group.by should be of length 1 or 2.")
    opts_df <- list(fit = fit, fun = fun, color = color, palette = palette,
        linetype = linetype, conf.int = conf.int, ggtheme = ggtheme,
        ...)
    opts <- list(fit = fit, data = data, fun = fun, color = color,
        palette = palette, linetype = linetype, conf.int = conf.int,
        pval = pval, pval.method = pval.method, test.for.trend = test.for.trend,
        surv.median.line = surv.median.line, risk.table = risk.table,
        cumevents = cumevents, cumcensor = cumcensor, tables.height = tables.height,
        ggtheme = ggtheme, tables.theme = tables.theme, ...)
    if (survminer:::.is_list(fit)) {
        if (combine)
            ggsurv <- do.call(ggsurvplot_combine, opts)
        else ggsurv <- do.call(ggsurvplot_list, opts)
    }
    else if (is.data.frame(fit))
        ggsurv <- do.call(ggsurvplot_df, opts_df)
    else if (survminer:::.is_survfit(fit)) {
        if (!is.null(group.by)) {
            opts$group.by <- group.by
            ggsurv <- do.call(ggsurvplot_group_by, opts)
        }
        else if (!is.null(facet.by)) {
            opts$facet.by <- facet.by
            ggsurv <- do.call(ggsurvplot_facet, opts)
        }
        else if (add.all) {
            ggsurv <- do.call(ggsurvplot_add_all, opts)
        }
        else {
            if (is.null(fit$strata)) {
                if (missing(conf.int)) {
                  opts$conf.int = TRUE
                  opts$conf.int.fill = "strata"
                }
            }
            ggsurv <- do.call(ggsurvplot_core2, opts)
        }
    }
return(ggsurv)
}
ggsurvplot_core2= function (fit, data = NULL, fun = NULL, color = NULL, palette = NULL,
    linetype = 1, break.x.by = NULL, break.y.by = NULL, break.time.by = NULL,
    surv.scale = c("default", "percent"), xscale = 1, conf.int = FALSE,
    conf.int.fill = "gray", conf.int.style = "ribbon", conf.int.alpha = 0.3,
    censor = TRUE, censor.shape = "+", censor.size = 4.5, pval = FALSE,
    pval.size = 5, pval.coord = c(NULL, NULL), test.for.trend = FALSE,
    pval.method = FALSE, pval.method.size = pval.size, pval.method.coord = c(NULL,
        NULL), log.rank.weights = c("survdiff", "1", "n", "sqrtN",
        "S1", "S2", "FH_p=1_q=1"), title = NULL, xlab = "Time",
    ylab = "Survival probability", xlim = NULL, ylim = NULL,
    axes.offset = TRUE, legend = c("top", "bottom", "left", "right",
        "none"), legend.title = "Strata", legend.labs = NULL,
    fontsize = 4.5, font.family = "", tables.height = 0.25, tables.y.text = TRUE,
    tables.col = "black", tables.y.text.col = TRUE, risk.table = FALSE,
    risk.table.pos = c("out", "in"), risk.table.title = NULL,
    risk.table.col = tables.col, risk.table.fontsize = fontsize,
    risk.table.y.text = tables.y.text, risk.table.y.text.col = tables.y.text.col,
    risk.table.height = tables.height, surv.plot.height = 0.75,
    ncensor.plot.height = tables.height, cumevents.height = tables.height,
    cumcensor.height = tables.height, ncensor.plot = FALSE, ncensor.plot.title = NULL,
    cumevents = FALSE, cumevents.col = tables.col, cumevents.title = NULL,
    cumevents.y.text = tables.y.text, cumevents.y.text.col = tables.y.text.col,
    cumcensor = FALSE, cumcensor.col = tables.col, cumcensor.title = NULL,
    cumcensor.y.text = tables.y.text, cumcensor.y.text.col = tables.y.text.col,
    surv.median.line = c("none", "hv", "h", "v"), ggtheme = theme_survminer(),
    tables.theme = ggtheme, ...) {
    if (!inherits(fit, "survfit"))
        stop("Can't handle an object of class ", class(fit))
    surv.median.line <- match.arg(surv.median.line)
    stopifnot(log.rank.weights %in% c("survdiff", "1", "n", "sqrtN",
        "S1", "S2", "FH_p=1_q=1"))
    log.rank.weights <- match.arg(log.rank.weights)
    if (ncensor.plot & cumcensor) {
        warning("Both ncensor.plot and cumsensor are TRUE.",
            "In this case, we consider only cumcensor.", call. = FALSE)
        ncensor.plot <- FALSE
    }
    if (cumcensor)
        ncensor.plot.height <- cumcensor.height
    if (is.null(ncensor.plot.title))
        ncensor.plot.title <- "Number of censoring"
    if (is.null(cumcensor.title))
        cumcensor.title <- "Cumulative number of censoring"
    if (is.null(cumevents.title))
        cumevents.title <- "Cumulative number of events"
    risk.table.pos <- match.arg(risk.table.pos)
    risktable <- survminer:::.parse_risk_table_arg(risk.table)
    risk.table <- risktable$display
    risk.table.type <- risktable$type
    extra.params <- list(...)
    .expand <- ggplot2::waiver()
    if (!axes.offset)
        .expand <- c(0, 0)
    data <- survminer:::.get_data(fit, data = data, complain = FALSE)
    d <- surv_summary(fit, data = data)
    if (!is.null(fit$start.time))
        d <- subset(d, d$time >= fit$start.time)
    xmin <- ifelse(survminer:::.is_cloglog(fun), min(c(1, d$time)), 0)
    if (!is.null(fit$start.time))
        xmin <- fit$start.time
    xmax <- survminer:::.get_default_breaks(d$time, .log = survminer:::.is_cloglog(fun)) %>%
        max()
    if (is.null(xlim))
        xlim <- c(xmin, xmax)
    p <- ggsurvplot_df(d, fun = fun, color = color, palette = palette,
        linetype = linetype, break.x.by = break.x.by, break.time.by = break.time.by,
        break.y.by = break.y.by, surv.scale = surv.scale, xscale = xscale,
        conf.int = conf.int, conf.int.fill = conf.int.fill, conf.int.style = conf.int.style,
        conf.int.alpha = conf.int.alpha, censor = censor, censor.shape = censor.shape,
        censor.size = censor.size, title = title, xlab = xlab,
        ylab = ylab, xlim = xlim, ylim = ylim, axes.offset = axes.offset,
        legend = legend, legend.title = legend.title, legend.labs = legend.labs,
        ggtheme = ggtheme, ...)
    pms <- attr(p, "parameters")
    color <- surv.color <- pms$color
    pval <- surv_pvalue(fit, method = log.rank.weights, data = data,
        pval = pval, pval.coord = pval.coord, pval.method.coord = pval.method.coord,
        test.for.trend = test.for.trend)
return(pval)
}