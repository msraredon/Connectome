#' FilterConnectome
#'
#' Filters a connectomic edgelist output from CreateConnectome according to user inputs. Defaults are set to reasonable initial parameters for data clean-up and exploration.
#'
#' @param connectome A connectomic edgelist output from CreateConnectome
#' @param min.pct Minimum fraction of cells within a given cluster expressing the ligand or receptor.
#' @param max.p Maximum p-value for ligand and receptor. Filtration on this column requires prior p-value calculation.
#' @param min.DOR Minimum log-normalized Diagnostic Odds Ratio for the ligand or receptor for its cell type within an edge.
#' @param min.exp Minimum normalized expression level of ligand and receptor.
#' @param min.z Minimum z-score for ligand and receptor.
#' @param modes.include String or vector signifying mode(s) of interest in include.
#' @param sources.include Source nodes of interest. Output will be limited to edges coming from these sources.
#' @param targets.include Target nodes of interest. Output will be limited to edges landing on these targets.
#' @param mechanisms.include Ligand - Receptor pairs of interest. The character string should match entries in the 'pair' column of the connectome.
#' @param verbose Whether to output feedback to user
#' @param remove.na Whether to remove edges containing 'NA' (no mapping to original object - only useful if investigating orphan ligands and receptors)

#' @export

FilterConnectome <- function(connectome,
                              min.pct = NULL,
                              max.p = NULL,
                              min.DOR = NULL,
                              min.exp = NULL,
                              min.z = NULL,
                              modes.include = NULL,
                              sources.include = NULL,
                              targets.include = NULL,
                              mechanisms.include = NULL,
                              features = NULL,
                              verbose = T,
                              remove.na = F){

  pre.filter <- nrow(connectome)

  # Percentage values
  if (!is.null(min.pct)){
    connectome <- subset(connectome, percent.source > min.pct & percent.target > min.pct)
  }

  # Expression values
  if (!is.null(min.exp)){
    connectome <- subset(connectome, ligand.expression > min.exp & recept.expression > min.exp)
  }

  # DOR values
  if (!is.null(min.DOR)){
    connectome <- subset(connectome, DOR.source > min.DOR & DOR.target > min.DOR)
  }

  # Scaled values
  if (!is.null(min.z)){
    connectome <- subset(connectome, ligand.scale > min.z & recept.scale > min.z)
  }

  # P-values
  if (!is.null(max.p)){
    if ('p_val_adj.lig' %in% colnames(connectome) & 'p_val_adj.rec' %in% colnames(connectome)){
      connectome <- subset(connectome, p_val_adj.lig < max.p & p_val_adj.rec < max.p)
    }else{message(paste("\np-values not available; p-value filtration was not performed"))}
  }

  # Modes
  if (!is.null(modes.include)){
      connectome <- subset(connectome, mode %in% modes.include)
    }


  # Nodes
  if (!is.null(sources.include)){
      connectome <- subset(connectome, source %in% sources.include)
    }
  if (!is.null(targets.include)){
      connectome <- subset(connectome, target %in% targets.include)
    }

  # Features
  if (!is.null(features)){
    connectome <- subset(connectome,ligand %in% features | receptor %in% features)
  }

  # Mechanisms
  if (!is.null(mechanisms.include)){
    connectome <- subset(connectome, pair %in% mechanisms.include)
  }
  
  # Remove NAs
  if (remove.na){
    connectome <- connectome[!is.na(connectome$percent.source),]
    connectome <- connectome[!is.na(connectome$percent.target),]
  }

post.filter <- nrow(connectome)

if (verbose){
            message(paste("\nPre-filter edges: ",as.character(pre.filter)))
            message(paste("\nPost-filter edges: ",as.character(post.filter)))
            message("\nConnectome filtration completed")
          }

return(connectome)
}
