#' Get a compiled Stan model by name
#'
#' @description
#' Retrieves and compiles a Stan model from the `tickTack` package by its
#' logical name. Supported models are `"timing_betabinomial"`,
#' `"timing_binomial"`, `"timing"`, `"clustering"`, and `"hierarchical"`.
#'
#' @param model_name A character string specifying the name of the Stan model
#'   to load. Must be one of `"timing_betabinomial"`, `"timing_binomial"`,
#'   `"timing"`, `"clustering"`, or `"hierarchical"`.
#'
#' @return A compiled `CmdStanModel` object as returned by
#'   \code{\link[cmdstanr]{cmdstan_model}}.
#'
#' @importFrom cmdstanr cmdstan_model
#' @importFrom utils capture.output
#'
#' @keywords internal

get_model <- function(model_name) {
  all_paths <- list(
    "base"        = "draft_base_model.stan",
    "constraint"  = "model_1_constraint.stan",
    "pseudotime"  = "model_1_pseudotime.stan",
    "realtime"    = "model_1_realtime.stan",
    "model1"      = "model_1.stan",
    "pseudotime2" = "model_2_pseudotime.stan"
  )

  if (!(model_name %in% names(all_paths))) {
    stop(paste0(
      "model_name '", model_name, "' not recognized. ",
      "Available models: ", paste(names(all_paths), collapse = ", ")
    ))
  }

  model_path <- system.file("stan", all_paths[[model_name]],
                            package = "SCtickTack", mustWork = TRUE)
  tmp <- utils::capture.output(
    suppressMessages(model <- cmdstanr::cmdstan_model(model_path))
  )
  model
}



#' Validate and normalise a sample name prefix
#'
#' @description
#' Checks that `sample_prefix` is either `NULL` or a non-empty character
#' string, and ensures it ends with `"_"` or `" "`. If the last character is
#' neither, a trailing `"_"` is appended automatically.
#'
#' @param sample_prefix Either `NULL` (returns `""`) or a character string
#'   to be used as a prefix for sample names.
#'
#' @return A character string: `""` when `sample_prefix` is `NULL`, otherwise
#'   the validated (and possibly `"_"`-suffixed) prefix.
#'
#' @keywords internal
validate_sample_prefix <- function(sample_prefix) {

  if (is.null(sample_prefix)) {
    return("")
  }

  if (!is.character(sample_prefix)) {
    stop("The parameter \"sample_prefix\" must be either NULL or a string.")
  }

  if (length(sample_prefix) > 0) {
    last_char <- substr(sample_prefix, nchar(sample_prefix),
                        nchar(sample_prefix))
    if (last_char != "_" && last_char != " ") {
      sample_prefix <- paste0(sample_prefix, "_")
    }
  }

  return(sample_prefix)
}


#' Simulate DLP+ sampling
#'
#' @description
#' This function simulates DLP+ sampling by collecting a sample for
#' each of the cells in a tissue rectangle. ProCESS supports multiple
#' DLP+ samples, however, the overall number of collected cells must
#' be less than 2^16.
#'
#' @param tissue_simulation The ProCESS tissue simulation from which the
#'   DLP+ sample must be collected.
#' @param bottom_left The coordinate of the bottom left corner of the
#'   DLP+ sample.
#' @param top_right The coordinate of the top right corner of the
#'   DLP+ sample.
#' @param sample_prefix A prefix for the DLP+ sample (default: `NULL`).
#' @export
#'
#' @examples
#' set.seed(0)
#'
#' sim <- TissueSimulation()
#' sim$add_mutant(name = "A", growth_rates = 0.1, death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(100)
#'
#' # sampling tissue
#' n_w <- n_h <- 10
#' ncells <- 5
#'
#' # adding second mutant
#' sim$add_mutant(name = "B", growth_rates = 0.3, death_rates = 0.0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#' sim$run_up_to_time(300)
#'
#' # find a tissue rectangle containing 5 cells of type A and B at least
#' bbox <- sim$search_sample(c("A" = ncells, "B" = ncells), n_w, n_h)
#'
#' # collect a DLP+ sample from the tissue simulation `sim`
#' DLP.sample(sim, bbox$lower_corner, bbox$upper_corner, sample_prefix="test")
#'
#' @seealso DLP.coverage()
DLP.sample <- function(tissue_simulation, bottom_left, top_right,
                       sample_prefix = NULL) {
  if (class(tissue_simulation)[1] != "Rcpp_TissueSimulation") {
    stop(paste("The parameter \"tissue_simulation\" must be a",
               "ProCESS's TissueSimulation object."))
  }

  sample_prefix <- validate_sample_prefix(sample_prefix)

  if (length(bottom_left) != 2 || bottom_left[1] < 0 || bottom_left[2] < 0) {
    stop(paste("The parameter \"bottom_left\" must be a 2-dimensional vector",
               "containing the coordinates of the bottom left cell of",
               "the DLP sample."))
  }

  if (length(top_right) != 2 || top_right[1] < 0 || top_right[2] < 0) {
    stop(paste("The parameter \"top_right\" must be a 2-dimensional vector",
               "containing the coordinates of the top right cell of",
               "the DLP sample."))
  }

  num_of_samples <- nrow(tissue_simulation$get_samples_info())

  num_of_future_samples <- ((top_right[1]-bottom_left[1]+1)
                            * (top_right[2]-bottom_left[2]+1))

  if (num_of_samples+num_of_future_samples >= 2^(16)) {
    stop(paste("Too many cells in the sample. Currently, ProCESS supports",
               "up to 2^16 cells in DLP samples."))
  }

  for (y in seq(bottom_left[2], top_right[2])) {
    for (x in seq(bottom_left[1], top_right[1])) {
      cell_pos <- c(x, y)
      cells <- tissue_simulation$get_cells(c(x, y), c(x,y))

      if (nrow(cells)>0) { # if the cell is a tumour cell
        cell_sample_name <- paste0(sample_prefix, x, "-", y)
        tissue_simulation$sample_cells(cell_sample_name, cell_pos, cell_pos)
      }
    }
  }
}

#' Compute DLP+ adjusted coverage
#'
#' @description
#' This function computes the coverage to use in `simulate_seq()` to
#' properly simulate a DLP+ sequencing in ProCESS.DLP+ sequencing
#' is simulated in ProCESS by associating each DPL+ sample cell to a
#' sample. However, the coverage parameter in `simulate_seq()` refers to
#' the coverage per sample, i.e., the coverage per cell. Thus, we need
#' to compute the coverage of the ProCESS samples before sequencing a
#' full DPL+ sample. Luckly, the coverage is uniform over cells in DLP+.
#' Hence, we can compute the cell coverage as the DPL+ sample coverage
#' divided by the number of cells in the DLP+ sample.
#'
#' @param phylo_forest The ProCESS phylogenetic forest from which the
#'   DLP+ sample must be collected.
#' @param overall_coverage The aimed coverage for the entire DPL+ sample.
#' @param sample_prefix The prefix of a DLP+ sample (default: `NULL`).
#' @export
#'
#' @examples
#' set.seed(0)
#'
#' sim <- TissueSimulation()
#' sim$add_mutant(name = "A", growth_rates = 0.1, death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(100)
#'
#' # sampling tissue
#' n_w <- n_h <- 10
#' ncells <- 5
#'
#' # adding second mutant
#' sim$add_mutant(name = "B", growth_rates = 0.3, death_rates = 0.0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#' sim$run_up_to_time(300)
#'
#' # find a tissue rectangle containing 5 cells of type A and B at least
#' bbox <- sim$search_sample(c("A" = ncells, "B" = ncells), n_w, n_h)
#'
#' # collect a DLP+ sample from the tissue simulation `sim`
#' DLP.sample(sim, bbox$lower_corner, bbox$upper_corner, sample_prefix="test")
#'
#' forest <- sim$get_sample_forest()
#'
#' rm(sim)
#'
#' # placing mutations
#' m_engine <- MutationEngine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8),
#'                     drivers = list(SNV("22", 46510210, "C", "A", allele = 1),
#'                                    "DGCR8 P26L"))
#' m_engine$add_mutant(mutant_name="B", passenger_rates=c(SNV=5e-9),
#'                     drivers = list(list("DGCR8 A18V", allele = 1)))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#'
#' phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#'
#' rm(forest)
#'
#' # compute adjusted coverage for simulating DLP+ sequencing
#' cov <- DLP.coverage(phylo_forest, 10, sample_prefix = "test")
#' cov
#'
#' # simulating DLP+ sequencing
#' seq_results <- simulate_seq(phylo_forest, coverage = cov, write_SAM = F,
#'                             with_normal_sample = FALSE)
#'
#' head(seq_results$mutations)
#' @seealso DLP.sample()
DLP.coverage <- function(phylo_forest, overall_coverage,
                         sample_prefix = NULL) {

  if (class(phylo_forest)[1] != "Rcpp_PhylogeneticForest") {
    stop(paste("The parameter \"phylo_forest\" must be a",
               "ProCESS's PhylogeneticForest object."))
  }

  if (class(overall_coverage) != "numeric") {
    stop("The parameter \"overall_coverage\" must be a numeric value.")
  }

  sample_prefix <- validate_sample_prefix(sample_prefix)

  sample_info <- phylo_forest$get_samples_info()

  num_of_DLP_samples <- sum(grepl(paste0("^", sample_prefix),
                                  sample_info$name))
  rm(sample_info)

  return(overall_coverage/num_of_DLP_samples)
}

get_segment_info_cell = function(data, chr, cell, new_from, new_to, keep_columns){
  data %>%
    dplyr::filter(cell_id == cell,
                  chr == !!chr,
                  begin <= new_from,
                  end >= new_to) %>%
    select(all_of(keep_columns))
}

DLP.uniform_segments <- function(x){
  out = lapply(x$chr %>% unique(), function(chr) {

    cli::cli_alert_info("Iterating on {.field {chr}}")
    cli::cli_h2("Iterating on {.field {chr}}")

    old_segments = sapply(x$cell_id %>% unique(), function(c) {
      x %>%
        dplyr::filter(cell_id == c) %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(seg_id) %>%
        unique() %>%
        length()
    })

    cli::cli_alert_info("Number of original segments in individual {.cls CNAqc} objects in {.field {chr}}:")
    #cli::cli_ul(paste(names(old_segments), old_segments, sep = " = "))
    cat("\n")

    # Chromosome-specific new breakpoints
    new_breakpoints = c(
      x %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(begin),
      x %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(end)) %>%
      unique() %>%
      sort()

    cli::cli_alert_info("Found {.field {length(new_breakpoints)}} breakpoints")
    #cat("\n")

    #  Separate new breakpoints into segment from and to values
    new_from = new_breakpoints[ !new_breakpoints == dplyr::last(new_breakpoints)] # last element will not be included in the from column
    new_to = new_breakpoints[!new_breakpoints == dplyr::first(new_breakpoints)] # first element will not be included in the to column

    # iterate on the new breakpoints to subset the cna piled up
    lapply(new_breakpoints[-1] %>% seq_along(), function(i) {

      # iterate for each sample
      lapply(x$cell_id %>% unique(), function(c) {

        # define which columns must be kept in the new table
        not_wanted = c("begin", "end", "length", "size", "seg_id", "chr", "cell_id", "n")
        wanted = setdiff(colnames(x), not_wanted)

        # get the information for the sample in the new segment
        tmp = get_segment_info_cell(x,
                                    chr = chr,
                                    cell = c,
                                    new_from = new_from[i],
                                    new_to = new_to[i],
                                    keep_columns = wanted)

        # do some checking on the result
        if (nrow(tmp) == 0) { # there is no information on copy number on the new segment: insert NA as value of all the columns, except from, to, segment_id and sample_id

          tmp_v2 = rep(NA, ncol(tmp))
          names(tmp_v2) = colnames(tmp)

          tmp = tmp_v2 %>%
            tibble::as_tibble_row()
        }

        # create a tibble with the information on the new breakpoints and include the previously retrieved information
        tidyr::tibble(chr = chr,
                      begin = new_from[i],
                      end = new_to[i],
                      cell_id = c) %>%
          dplyr::bind_cols(tmp) %>%
          dplyr::mutate(seg_id = paste(chr, begin, end, sep = ":"))

      }) %>% do.call(bind_rows, .)
    }) %>% do.call(bind_rows, .)
  })  %>% do.call(bind_rows, .)

  out <- out %>%
    filter(end - begin > 1)
  return(out)
}
