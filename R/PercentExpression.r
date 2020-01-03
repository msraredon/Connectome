PercentExpression <- function (object, genes.use = NULL, return.seurat = FALSE, add.ident = NULL,
                                 use.scale = FALSE, use.raw = TRUE, show.progress = TRUE,
                                 ...)
{
    ident.orig <- object@ident
    orig.levels <- levels(x = object@ident)
    ident.new <- c()
    if (!is.null(x = add.ident)) {
        new.data <- FetchData(object = object, vars.all = add.ident)
        new.ident <- paste(object@ident[rownames(x = new.data)],
                           new.data[, 1], sep = "_")
        object <- SetIdent(object = object, cells.use = rownames(x = new.data),
                           ident.use = new.ident)
    }
    if (return.seurat) {
        assays.use <- c("RNA", names(x = object@assay))
    }
    else {
        assays.use <- "RNA"
    }
    slot.use <- "data"
    fxn.average <- function(x) {sum(expm1(x>0))/length(x)}
    if (use.scale) {
        slot.use <- "scale.data"
        fxn.average <- function(x) {sum(x>0)/length(x)}
    }
    if (use.raw) {
        slot.use <- "raw.data"
        fxn.average <- function(x) {sum(x>0)/length(x)}
    }
    data.return <- list()
    for (i in 1:length(x = assays.use)) {
        data.use <- GetAssayData(object = object, assay.type = assays.use[i],
                                 slot = slot.use)
        genes.assay <- genes.use
        if (length(x = intersect(x = genes.use, y = rownames(x = data.use))) <
            1) {
            genes.assay <- rownames(x = data.use)
        }
        data.all <- data.frame(row.names = genes.assay)
        for (j in levels(x = object@ident)) {
            temp.cells <- WhichCells(object = object, ident = j)
            genes.assay <- unique(x = intersect(x = genes.assay,
                                                y = rownames(x = data.use)))
            if (length(x = temp.cells) == 1) {
                data.temp <- (data.use[genes.assay, temp.cells])
            }
            if (length(x = temp.cells) > 1) {
                data.temp <- apply(X = data.use[genes.assay,
                                                temp.cells], MARGIN = 1, FUN = fxn.average)
            }
            data.all <- cbind(data.all, data.temp)
            colnames(x = data.all)[ncol(x = data.all)] <- j
            if (show.progress) {
                print(paste0("Finished calculating percent expression", assays.use[i],
                             " for cluster ", j))
            }
            if (i == 1) {
                ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
            }
        }
        names(x = ident.new) <- levels(x = object@ident)
        data.return[[i]] <- data.all
        names(x = data.return)[i] <- assays.use[[i]]
    }
    if (return.seurat) {
        toRet <- CreateSeuratObject(raw.data = data.return[[1]],
                                    project = "Average", min.cells = 0, min.genes = 0,
                                    is.expr = 0, ...)
        if (length(x = data.return) > 1) {
            for (i in 2:length(x = data.return)) {
                toRet <- SetAssayData(object = toRet, assay.type = names(x = data.return)[i],
                                      slot = "raw.data", new.data = data.return[[i]])
            }
        }
        toRet <- SetIdent(object = toRet, cells.use = toRet@cell.names,
                          ident.use = ident.new[toRet@cell.names])
        toRet@ident <- factor(x = toRet@ident, levels = as.character(x = orig.levels),
                              ordered = TRUE)
        toRet <- NormalizeData(toRet, display.progress = show.progress)
        toRet <- ScaleData(toRet, display.progress = show.progress)
        return(toRet)
    }
    else {
        return(data.return[[1]])
    }
}
