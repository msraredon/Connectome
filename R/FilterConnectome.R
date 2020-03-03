#' FilterConnectome
#'
#' Filters a connectomic edgelist output from CreateConnectome according to user inputs. Defaults are set to reasonable initial parameters for data clean-up and exploration.
#'
#' @param connectome A connectomic edgelist output from CreateConnectome
#' @param min.pct Minimum fraction of cells within a given cluster expressing the ligand or receptor. Defaults to 0.10, allows NULL.
#' @param min.exp Minimum normalized expression level of ligand and receptor. Defaults to 0, allows NULL.
#' @param min.z Minimum z-score for ligand and receptor. Defaults to 0, allows NULL.
#' @param max.p Maximum p-value for ligand and receptor. Defaults to NULL. Filtration on this column requires prior p-value calculation.
#' @param modes.include String or vector signifying mode(s) of interest in include. Defaults to all modes.
#' @export

FilterConnectome <- function(connectome,
                              min.pct = 0.10,
                              max.p = 0.01,
                              min.exp = 0,
                              min.z = 0,
                              modes.include = NULL,...){

  if (!is.null(min.pct)){
    connectome <- subset(connectome, percent.source > min.pct & percent.target > min.pct)
  }

  if (!is.null(min.exp)){
    connectome <- subset(connectome, ligand.expression > min.exp & recept.expression > min.exp)
  }

  if (!is.null(min.z)){
    connectome <- subset(connectome, ligand.scale > min.z & recept.scale > min.z)
  }

  if (!is.null(max.p)){
    if ('lig.p' %in% colnames(connectome) & 'rec.p' %in% colnames(connectome)){
      connectome <- subset(connectome, lig.p < max.p & rec.p < max.p)
    }else{message(paste("\np-values not available; p-value filtration was not performed"))}
  }

  if (!is.null(modes.include)){
    if (length(modes.include) == 1){
      connectome <- subset(connectome, mode == modes.include)
    }else{
      connectome <- subset(connectome, mode %in% modes.include)
    }
  }

message(paste("\nConnectome filtration completed"))
return(connectome)
}
