#' SignalingTSNE
#'
#' In beta with James Garritano.
#'
#' <Description here once developed>
#'
#' @param object A Seurat 3.0 object.
#' @param mechanisms A dataframe of desired signaling mechanisms to inventigate at once. Should be two matched columns, with the first column the ligand and the second column the receptor. #James I think this is the best way to do this for now. Can modify later. Also, if we later want to add ability to modify relative weighting of different mechanisms (unclear to me right now) we could add the option to include a third column which sums to 1 and specifies this
#' @param sending.cells Cells to consider when looking at ligand expression.  Default all cells in object.
#' @param receiving.cells Cells to consdier when looking at receptor expression.  Default all cells in object.
#' @param plot.centroids Whether to plot centroids for combined signaling affinity and each individual mechanism. Default = ?
#' @param plot.trajectories Whether to plot mean trajectories for each receiving cell type, for each mechanism considered. Default = ?
#' @param plot.forces Whether to plot force vectors which explain the position of each single-mechanism centroid. Default = ?
#' @param output.gif Whether to output a GIF of cells moving through their trajectories. Default = ?
#' @param ... Additional parameters passed to any nested functions. # James this may or may not be needed, I tend to add just in case.
#' @export

SignalingTSNE <- function(object,
                          mechanisms,
                          sending.cells = NULL,
                          receiving.cells = NULL,
                          plot.centroids = T,
                          plot.trajectories = T,
                          plot.forces = F,
                          output.gif = F,...){

}
