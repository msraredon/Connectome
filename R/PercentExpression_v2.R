PercentExpression_v2 <- function (object, assays = 'RNA', features = NULL, return.seurat = FALSE,
    add.ident = NULL, slot = "counts", use.scale = FALSE, use.counts = FALSE,
    verbose = TRUE, ...)
{
    if (use.scale) {
        .Deprecated(msg = "'use.scale' is a deprecated argument, please use the 'slot' argument instead")
        slot <- "scale.data"
    }
    if (use.counts) {
        .Deprecated(msg = "'use.counts' is a deprecated argument, please use the 'slot' argument instead")
        if (use.scale) {
            warning("Both 'use.scale' and 'use.counts' were set; using counts",
                call. = FALSE, immediate. = TRUE)
        }
        slot <- "counts"
    }
    fxn.average <- function(x) {sum(x>0)/length(x)}
    object.assays <- Seurat:::FilterObjects(object = object, classes.keep = "Assay")
    # assays <- assays %||% object.assays
    ident.orig <- Idents(object = object)
    orig.levels <- levels(x = Idents(object = object))
    ident.new <- c()
    if (!all(assays %in% object.assays)) {
        assays <- assays[assays %in% object.assays]
        if (length(assays) == 0) {
            stop("None of the requested assays are present in the object")
        }
        else {
            warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
        }
    }
    if (!is.null(x = add.ident)) {
        new.data <- FetchData(object = object, vars = add.ident)
        new.ident <- paste(Idents(object)[rownames(x = new.data)],
            new.data[, 1], sep = "_")
        Idents(object, cells = rownames(new.data)) <- new.ident
    }
    data.return <- list()
    for (i in 1:length(x = assays)) {
        data.use <- GetAssayData(object = object, assay = assays[i],
            slot = slot)
        features.assay <- features
        if (length(x = intersect(x = features, y = rownames(x = data.use))) <
            1) {
            features.assay <- rownames(x = data.use)
        }
        data.all <- data.frame(row.names = features.assay)
        for (j in levels(x = Idents(object))) {
            temp.cells <- WhichCells(object = object, idents = j)
            features.assay <- unique(x = intersect(x = features.assay,
                y = rownames(x = data.use)))
            if (length(x = temp.cells) == 1) {
                data.temp <- (data.use[features.assay, temp.cells])
                if (slot == "data") {
                  data.temp <- expm1(x = data.temp)
                }
            }
            if (length(x = temp.cells) > 1) {
                data.temp <- apply(X = data.use[features.assay,
                  temp.cells, drop = FALSE], MARGIN = 1, FUN = fxn.average)
            }
            data.all <- cbind(data.all, data.temp)
            colnames(x = data.all)[ncol(x = data.all)] <- j
            if (verbose) {
                message(paste("Finished calculating percentage expressing", assays[i],
                  "for cluster", j))
            }
            if (i == 1) {
                ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
            }
        }
        names(x = ident.new) <- levels(x = Idents(object))
        data.return[[i]] <- data.all
        names(x = data.return)[i] <- assays[[i]]
    }
    if (return.seurat) {
        toRet <- CreateSeuratObject(counts = data.return[[1]],
            project = "Average", assay = names(x = data.return)[1],
            ...)
        if (length(x = data.return) > 1) {
            for (i in 2:length(x = data.return)) {
                toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]])
            }
        }
        if (DefaultAssay(object = object) %in% names(x = data.return)) {
            DefaultAssay(object = toRet) <- DefaultAssay(object = object)
        }
        Idents(toRet, cells = colnames(x = toRet)) <- ident.new[colnames(x = toRet)]
        Idents(object = toRet) <- factor(x = Idents(object = toRet),
            levels = as.character(x = orig.levels), ordered = TRUE)
        toRet <- NormalizeData(object = toRet, verbose = verbose)
        toRet <- ScaleData(object = toRet, verbose = verbose)
        return(toRet)
    }
    else {
        return(data.return)
    }
}
