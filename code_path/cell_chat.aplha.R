
netVisual_spatial=function (net, coordinates, meta, sample.use = NULL, color.use = NULL, 
    title.name = NULL, sources.use = NULL, targets.use = NULL, 
    idents.use = NULL, remove.isolate = FALSE, remove.loop = TRUE, 
    top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
    vertex.size.max = NULL, vertex.label.cex = 5, vertex.label.color = "black", 
    edge.weight.max = NULL, edge.width.max = 8, edge.curved = 0.2, 
    alpha.edge = 0.6, arrow.angle = 5, arrow.size = 0.2, alpha.image = 0.15, 
    point.size = 1.5, legend.size = 5) 
{
    cells.level <- rownames(net)
    labels <- meta$labels
    samples <- meta$samples
    if (ncol(coordinates) == 2) {
        colnames(coordinates) <- c("x_cent", "y_cent")
        if (length(unique(samples)) > 1) {
            if (is.null(sample.use)) {
                stop("`sample.use` should be provided for visualizing signaling on each individual sample.")
            }
            else if (sample.use %in% unique(samples)) {
                coordinates = coordinates[samples == sample.use, 
                  ]
                labels = labels[samples == sample.use]
            }
            else {
                stop("Please check the input `sample.use`, which should be the element in `meta$samples`.")
            }
        }
        temp_coordinates = coordinates
        coordinates[, 1] = temp_coordinates[, 2]
        coordinates[, 2] = temp_coordinates[, 1]
    }
    else {
        stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
    }
    num_cluster <- length(cells.level)
    node_coords <- matrix(0, nrow = num_cluster, ncol = 2)
    for (i in c(1:num_cluster)) {
        node_coords[i, 1] <- median(coordinates[as.character(labels) == 
            cells.level[i], 1])
        node_coords[i, 2] <- median(coordinates[as.character(labels) == 
            cells.level[i], 2])
    }
    rownames(node_coords) <- cells.level
    if (is.null(vertex.size.max)) {
        if (length(unique(vertex.weight)) == 1) {
            vertex.size.max <- 5
        }
        else {
            vertex.size.max <- 15
        }
    }
    options(warn = -1)
    thresh <- stats::quantile(net, probs = 1 - top)
    net[net < thresh] <- 0
    if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use))) {
        if (is.null(rownames(net))) {
            stop("The input weighted matrix should have rownames!")
        }
        df.net <- reshape2::melt(net, value.name = "value")
        colnames(df.net)[1:2] <- c("source", "target")
        if (!is.null(sources.use)) {
            if (is.numeric(sources.use)) {
                sources.use <- cells.level[sources.use]
            }
            df.net <- subset(df.net, source %in% sources.use)
        }
        if (!is.null(targets.use)) {
            if (is.numeric(targets.use)) {
                targets.use <- cells.level[targets.use]
            }
            df.net <- subset(df.net, target %in% targets.use)
        }
        if (!is.null(idents.use)) {
            if (is.numeric(idents.use)) {
                idents.use <- cells.level[idents.use]
            }
            df.net <- filter(df.net, (source %in% idents.use) | 
                (target %in% idents.use))
        }
        df.net$source <- factor(df.net$source, levels = cells.level)
        df.net$target <- factor(df.net$target, levels = cells.level)
        df.net$value[is.na(df.net$value)] <- 0
        net <- tapply(df.net[["value"]], list(df.net[["source"]], 
            df.net[["target"]]), sum)
    }
    net[is.na(net)] <- 0
    if (remove.loop) {
        diag(net) <- 0
    }
    if (remove.isolate) {
        idx1 <- which(Matrix::rowSums(net) == 0)
        idx2 <- which(Matrix::colSums(net) == 0)
        idx <- intersect(idx1, idx2)
        net <- net[-idx, ]
        net <- net[, -idx]
        node_coords <- node_coords[-idx, ]
        cells.level <- cells.level[-idx]
    }
    g <- graph_from_adjacency_matrix(net, mode = "directed", 
        weighted = T)
    edgelist <- get.edgelist(g)
    edges <- data.frame(node_coords[edgelist[, 1], ], node_coords[edgelist[, 
        2], ])
    colnames(edges) <- c("X1", "Y1", "X2", "Y2")
    node_coords = data.frame(node_coords)
    node_idents = factor(cells.level, levels = cells.level)
    node_family = data.frame(node_coords, node_idents)
    if (is.null(color.use)) {
        color.use = scPalette(length(igraph::V(g)))
        names(color.use) <- cells.level
    }
    if (is.null(vertex.weight.max)) {
        vertex.weight.max <- max(vertex.weight)
    }
    vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
        5
    if (is.null(edge.weight.max)) {
        edge.weight.max <- max(igraph::E(g)$weight)
    }
    if (weight.scale == TRUE) {
        igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
            edge.width.max
    }
    else {
        igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
    }
    gg <- ggplot(data = node_family, aes(X1, X2)) +  geom_point(aes(X1, 
        X2, colour = node_idents), data = node_family, size = vertex.weight, 
        show.legend = TRUE) + geom_curve(aes(x = X1, 
        y = Y1, xend = X2, yend = Y2), data = edges, size = igraph::E(g)$width, 
        curvature = edge.curved, alpha = 1, arrow = arrow(angle = arrow.angle, 
            type = "closed", length = unit(arrow.size, "inches")), 
        colour = color.use[edgelist[, 1]]) + scale_color_manual(values = color.use) + 
        guides(color = guide_legend(override.aes = list(size = legend.size))) + 
        xlab(NULL) + ylab(NULL) + coord_fixed() + theme(aspect.ratio = 1) + 
        theme(legend.key = element_blank()) + theme(panel.background = element_blank(), 
        axis.ticks = element_blank(), panel.border = element_blank(), 
        axis.text = element_blank(), legend.title = element_blank())
    gg <- gg + geom_point(aes(x_cent, y_cent), data = coordinates, 
        colour = color.use[labels], alpha = 1, size = point.size, 
        show.legend = FALSE)
    gg <- gg + scale_y_reverse()
    if (vertex.label.cex > 0) {
        gg <- gg + ggnetwork::geom_nodetext_repel(aes(label = node_idents), 
            color = "black", size = vertex.label.cex)
    }
    if (!is.null(title.name)) {
        gg <- gg + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5, 
            vjust = 0))
    }
    gg
    return(gg)
}