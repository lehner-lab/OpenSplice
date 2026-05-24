## 06_shared.R
## Shared helper functions for the 06_cis_regulatory_elements analysis scripts.
##
## Sourced by:
##   06.6_transition.R
##   06.8_clustering_plot2.R

library(data.table)

#' Count state transitions per exon per region
#'
#' @param annotated_summary  data.table/data.frame with columns:
#'   exon_id, start, exon_length, type_annotation
#' @return data.table with columns: exon_id, region, n_transition, transition_type
calculate_state_transitions <- function(annotated_summary) {

  setDT(annotated_summary)

  valid_states <- c("O", "S", "E", "N")

  possible_transitions <- na.omit(
    as.vector(outer(valid_states, valid_states,
                    FUN = Vectorize(function(a, b) {
                      if (a != b) paste(a, b, sep = " > ") else NA_character_
                    })))
  )

  all_results <- list()

  for (exon in unique(annotated_summary$exon_id)) {
    sub_dt            <- annotated_summary[exon_id == exon][order(start)]
    intron_down_start <- unique(sub_dt$exon_length) + 71

    sub_dt[, region := fifelse(
      start <= 70, "Intron_up",
      fifelse(start >= intron_down_start, "Intron_down", "Exon")
    )]

    region_results <- sub_dt[, {
      state_shift <- shift(type_annotation, type = "lead")
      transition  <- paste(type_annotation, state_shift, sep = " > ")
      transition  <- transition[
        !is.na(state_shift) & type_annotation != state_shift
      ]
      transition_table <- table(factor(transition, levels = possible_transitions))
      data.table(
        transition_type = names(transition_table),
        n_transition    = as.integer(transition_table)
      )
    }, by = .(exon_id, region)]

    all_results[[exon]] <- region_results
  }

  final_dt <- rbindlist(all_results, fill = TRUE)
  return(final_dt[, .(exon_id, region, n_transition, transition_type)])
}
