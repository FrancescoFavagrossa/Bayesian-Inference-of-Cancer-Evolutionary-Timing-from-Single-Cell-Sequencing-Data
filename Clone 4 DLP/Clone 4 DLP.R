
set.seed(123)

library(ProCESS)
library(tidyverse)
library(cmdstanr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(stringr)

genome_dir <- '/Users/francescofavagrossa/Documents/Projects/Thesis/GRCh38'
wd <- '/Users/francescofavagrossa/Documents/Projects/Thesis/Clone 4 DLP'
setwd(wd)
source('/Users/francescofavagrossa/Documents/Projects/Thesis/Helpers and Utils/Helpers.R')
source('/Users/francescofavagrossa/Documents/Projects/Thesis/Helpers and Utils/utils.R')

show_plot <- function(plot) {
  if (interactive()) print(plot)
  invisible(plot)
}

################################################################################
############################ ProCESS SIMULATION ################################
################################################################################

# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------


mu_SNV <- 5e-8
mu_INDELs <- 0

w_bar_prior_mean <- 0.06
w_bar_prior_sd <- 0.35
sigma_log_w_rate <- 2
clone_abundance_nb_phi <- 100

clone_names <- c("A", "B", "C", "D", "E")
parent <- c(NA, 1, 2, 3, 4)

n_clones <- length(clone_names)

w_bar <- rlnorm(
  n = 1,
  meanlog = log(w_bar_prior_mean),
  sdlog = w_bar_prior_sd
)

sigma_log_w <- abs(rnorm(
  n = 1,
  mean = 0,
  sd = 1 / sigma_log_w_rate
))
# Random part of the per-clone growth spread (a deterministic floor is added to
# each increment below to guarantee the clones are genuinely distinct).
sigma_log_w <- 0.25

w <- numeric(n_clones)
log_w <- numeric(n_clones)

w[1] <- w_bar
log_w[1] <- log(w[1])

for (i in 2:n_clones) {
  # Floor on each increment so consecutive clones along a lineage are genuinely
  # distinct in growth (the clock-growth coupling needs separated w to break the
  # w-t degeneracy; small random increments leave clones near-identical).
  increment <- 0.3 + abs(rnorm(1, mean = 0, sd = sigma_log_w))
  log_w[i] <- log_w[parent[i]] + increment
  w[i] <- exp(log_w[i])
}

data.frame(
  clone = clone_names,
  parent = c(NA, clone_names[parent[-1]]),
  log_w = log_w,
  w = w
)

# ------------------------------------------------------------
# 2. Cancer evolution simulation
# ------------------------------------------------------------

sim <- TissueSimulation(
  name          = "4C Prop",
  save_snapshots = FALSE,
  width  = 5000,
  height = 5000
)

# Volumetric (exponential) growth: every cell can divide, not only the tumour
# border. Required for consistency with the Stan model, which assumes
# exponential clone growth (terminal_clone_mass = exp(w * delta_t * (1 - t)))
# and clonal-mutation accumulation proportional to w * delta_t. The default
# border-growth model grows ~polynomially and biases the inference.
sim$border_growth_model <- FALSE

# With volumetric (exponential) growth the founder keeps dividing for the whole
# run, so intermediates are grown to a small size and only the sampled leaf clone
# is grown large (mirrors Clone 6 DLP). This keeps the tissue and the DLP
# sampling box within ProCESS's 2^16-cell limit and the leaf clone abundant
# enough to sample.
sim$add_mutant(name = "A", growth_rates = w[1], death_rates = 0)
sim$place_cell("A", 2500, 2500)
sim$run_up_to_size(species = "A", 5)

sim$add_mutant(name = "B", growth_rates = w[2], death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("A"), "B")
sim$run_up_to_size(species = "B", 5)

sim$add_mutant(name = "C", growth_rates = w[3], death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("B"), "C")
sim$run_up_to_size(species = "C", 5)

sim$add_mutant(name = "D", growth_rates = w[4], death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("C"), "D")
sim$run_up_to_size(species = "D", 100)

Tf <- sim$get_clock()


# ------------------------------------------------------------
# 3. Tissue Plot
# ------------------------------------------------------------

clones <- data.frame(
  mutant = c("A", "B", "C", "D"),
  color  = c("grey70", "dodgerblue", "orange", "red")
)

color_map <- stats::setNames(clones$color, clones$mutant)
forest_color_map <- c(
  A = "#56B4E9",
  B = "#0072B2",
  C = "#009E73",
  D = "#D55E00"
)

make_forest_species_color_map <- function(forest, clone_color_map) {
  forest$get_nodes() %>%
    mutate(species = paste0(.data$mutant, .data$epistate)) %>%
    distinct(.data$species, .data$mutant) %>%
    mutate(color = clone_color_map[.data$mutant]) %>%
    { stats::setNames(.$color, .$species) }
}

make_colored_forest_edges <- function(forest, forest_plot) {
  layout <- forest_plot$data %>%
    select(.data$name, .data$x, .data$y) %>%
    mutate(name = as.character(.data$name))
  
  forest$get_nodes() %>%
    transmute(
      from = ifelse(is.na(.data$ancestor), "WT", as.character(.data$ancestor)),
      to = as.character(.data$cell_id),
      species = paste0(.data$mutant, .data$epistate)
    ) %>%
    left_join(layout %>% rename(x_from = x, y_from = y), by = c("from" = "name")) %>%
    left_join(layout %>% rename(x_to = x, y_to = y), by = c("to" = "name")) %>%
    filter(!is.na(.data$x_from), !is.na(.data$y_from), !is.na(.data$x_to), !is.na(.data$y_to))
}

cells_plot <- sim$get_cells() %>%
  filter(mutant %in% clones$mutant)

pad <- 2

clone_abundance <- cells_plot %>%
  count(mutant, name = "n_tissue") %>%
  right_join(clones, by = "mutant") %>%
  mutate(
    n_tissue = tidyr::replace_na(n_tissue, 0L),
    fraction = n_tissue / sum(n_tissue)
  ) %>%
  arrange(match(mutant, clones$mutant))

print(clone_abundance)

p <- ggplot(cells_plot, aes(x = position_x, y = position_y, fill = mutant)) +
      geom_tile() +
      scale_fill_manual(values = color_map, drop = FALSE) +
      coord_equal(
        xlim = range(cells_plot$position_x) + c(-pad, pad),
        ylim = range(cells_plot$position_y) + c(-pad, pad)
      ) +
      theme_classic() +
      labs(
        title = "Tissue simulation",
        x = "x position",
        y = "y position",
        fill = "Clone"
      )

ggsave(
  filename = "figures 4 DLP/Tissue_Plot.png",
  plot = p, width = 16, height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 4. Sampling
# ------------------------------------------------------------

# DLP samples are capped at 2^16 cells, and search_sample returns a box of
# exactly n_w x n_h positions, so n_w * n_h must stay below 65536.
n_w <- 250
n_h <- 250

ncells_by_clone <- cells_plot %>%
  count(mutant, name = "n") %>%
  right_join(clones, by = "mutant") %>%
  mutate(n = tidyr::replace_na(n, 0L)) %>%
  arrange(match(mutant, clones$mutant))

sample_targets <- stats::setNames(
  ncells_by_clone$n / 10,
  ncells_by_clone$mutant
)

bbox1 <- sim$search_sample(sample_targets, n_w, n_h)

p_cells <- ggplot(cells_plot, aes(x = position_x, y = position_y, fill = mutant)) +
  geom_tile() +
  annotate(
    "rect",
    xmin = bbox1$lower_corner[1],
    xmax = bbox1$upper_corner[1],
    ymin = bbox1$lower_corner[2],
    ymax = bbox1$upper_corner[2],
    fill = NA,
    color = "black",
    linewidth = 0.9
  ) +
  scale_fill_manual(values = color_map, drop = FALSE) +
  coord_equal(
    xlim = range(cells_plot$position_x) + c(-pad, pad),
    ylim = range(cells_plot$position_y) + c(-pad, pad)
  ) +
  theme_classic() +
  labs(
    title = "Tissue simulation with sampled region",
    x = "x position",
    y = "y position",
    fill = "Clone"
  )

print(p_cells)

ggsave(
  filename = "figures 4 DLP/Tissue_Sampled_Plot.png",
  plot = p_cells, width = 16, height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

DLP.sample(sim, bbox1$lower_corner, bbox1$upper_corner, sample_prefix="DLP Clone 4")

forest <- sim$get_sample_forest()

forest_nodes <- forest$get_nodes()

prop <- forest_nodes %>%
  filter(!is.na(sample)) %>%
  count(mutant, name = "sample_n") %>%
  mutate(
    sample_percent = 100 * sample_n / sum(sample_n)
  )

sim_prop <- sim$get_counts() %>%
  select(mutant, overall) %>%
  mutate(
    simulation_percent = 100 * overall / sum(overall, na.rm = TRUE)
  )

comparison <- sim_prop %>%
  left_join(prop, by = "mutant") %>%
  mutate(
    sample_n = ifelse(is.na(sample_n), 0, sample_n),
    sample_percent = ifelse(is.na(sample_percent), 0, sample_percent),
    simulation_percent = round(simulation_percent, 2),
    sample_percent = round(sample_percent, 2)
  )

comparison

# ------------------------------------------------------------
# 5. Forest Plot
# ------------------------------------------------------------

forest_species_color_map <- make_forest_species_color_map(forest, forest_color_map)
p <- plot_forest(forest, color_map = forest_species_color_map)
colored_forest_edges <- make_colored_forest_edges(forest, p)

p <- p +
  geom_segment(
    data = colored_forest_edges,
    aes(
      x = .data$x_from,
      y = .data$y_from,
      xend = .data$x_to,
      yend = .data$y_to,
      color = .data$species
    ),
    inherit.aes = FALSE,
    linewidth = 0.16,
    alpha = 0.45
  ) +
  theme(legend.position = "none")

show_plot(p)

ggsave(
  filename = "figures 4 DLP/Forest_Plot.png",
  plot = p, width = 16, height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 6. Mutation Engine
# ------------------------------------------------------------

setwd(genome_dir)

m_engine <- MutationEngine(
  setup_code = "GRCh38",
  tumour_type = "COADREAD",
  context_sampling = 20
)

setwd(wd)

CNA_A        <- ProCESS::CNA(type = "A", "12", from = 10000000,  len = 4e7)
CNA_C_hidden <- ProCESS::CNA(type = "A", "5",  from = 107707518, len = 2e7)


cna_definitions <- tribble(
  ~cna_name,               ~mutant, ~target_label, ~type, ~chr, ~begin,     ~len,
  "CNA_A_chr12_40Mb",      "A",     "Clone_A",     "A",   "12", 10000000,   40000000,
  "CNA_C_hidden_chr5_20Mb","C",     "Clone_C",     "A",   "5",  107707518,  20000000
)


m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),
                    drivers = list("APC R1450*",CNA_A))

m_engine$add_mutant(mutant_name = "B",
                    passenger_rates = c(SNV = mu_SNV, CNA = 0, indel = mu_INDELs))

m_engine$add_mutant(mutant_name = "C",
                    passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),
                    drivers = list(CNA_C_hidden))

m_engine$add_mutant(mutant_name = "D",
                    passenger_rates = c(SNV = mu_SNV, CNA = 0, indel = mu_INDELs))

m_engine$add_exposure( time = 0, coefficients = c(SBS1 = 0.15, SBS5 = 0.40, SBS18 = 0.15,
                                                  SBS17b = 0.20, ID1 = 0.40, ID2 = 0.40, ID18 = 0.20, SBS88 = 0.10)
)

phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs = 800, num_of_preneoplatic_indels = 200)

forest_species_color_map <- make_forest_species_color_map(forest, forest_color_map)
p_base <- plot_forest(forest, color_map = forest_species_color_map)
colored_forest_edges <- make_colored_forest_edges(forest, p_base)

p <- (
  p_base +
    geom_segment(
      data = colored_forest_edges,
      aes(
        x = .data$x_from,
        y = .data$y_from,
        xend = .data$x_to,
        yend = .data$y_to,
        color = .data$species
      ),
      inherit.aes = FALSE,
      linewidth = 0.16,
      alpha = 0.45
    )
) %>%
  annotate_forest(phylo_forest, MRCAs = FALSE) +
  ggraph::geom_node_point(
    aes(color = .data$species),
    size = 0.65,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  theme(legend.position = "none")

show_plot(p)

ggsave(
  filename = "figures 4 DLP/Forest_Mutations_Plot.png",
  plot = p,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

clones <- tibble::tribble(
  ~mutant, ~parent_mutant, ~cell_label,
  "A",     NA,             "Clone_A",
  "B",     "A",            "Clone_B",
  "C",     "B",            "Clone_C",
  "D",     "C",            "Clone_D"
)

################################################################################
############################ ANALYSE SIMULATION ################################
################################################################################

# ------------------------------------------------------------
# 1. Sampled nodes and clone IDs
# ------------------------------------------------------------

forest_nodes <- phylo_forest$get_nodes()

sampled_nodes <- forest_nodes %>%
  filter(
    !is.na(.data$sample),
    .data$mutant %in% clones$mutant
  ) %>%
  left_join(
    clones %>% select(mutant, cell_label),
    by = "mutant"
  )

n_cells_per_clone <- sampled_nodes %>%
  count(mutant, name = "n_cells_total")

print(n_cells_per_clone)

sampled_clone_order <- clones$mutant

cell_label_map <- sampled_nodes %>%
  select(cell_label, mutant, cell_id)

clone_ids <- lapply(clones$mutant, function(clone_name) {
  sampled_nodes %>%
    filter(mutant == clone_name) %>%
    pull(cell_id)
})

names(clone_ids) <- clones$mutant

print(lengths(clone_ids))


# ------------------------------------------------------------
# 2. CNA regions present in sampled cells
# ------------------------------------------------------------

cna_row <- phylo_forest$get_cell_allelic_fragmentation() %>%
  mutate(karyo = paste0(.data$major, ":", .data$minor)) %>%
  filter(
    .data$karyo %in% c("2:2", "2:0", "2:1"),
    .data$cell_id %in% cell_label_map$cell_id
  )

shared_events_sample <- cna_row %>%
  group_by(.data$chr, .data$begin, .data$end, .data$karyo) %>%
  summarise(
    cells = list(unique(.data$cell_id)),
    n_cells = n_distinct(.data$cell_id),
    .groups = "drop"
  ) %>%
  filter(.data$n_cells >= 1) %>%
  arrange(.data$chr, .data$begin, .data$end, .data$karyo)

cna_regions <- tibble()

for (i in seq_len(nrow(cna_definitions))) {
  event <- shared_events_sample %>%
    filter(
      .data$chr == cna_definitions$chr[[i]],
      .data$begin == cna_definitions$begin[[i]]
    )
  
  if (nrow(event) != 1) {
    stop(
      "Expected exactly one CNA event for ",
      cna_definitions$cna_name[[i]],
      "; found ",
      nrow(event)
    )
  }
  
  cna_regions <- bind_rows(
    cna_regions,
    tibble(
      cna_name = cna_definitions$cna_name[[i]],
      target_label = cna_definitions$target_label[[i]],
      mutant = cna_definitions$mutant[[i]],
      chr = cna_definitions$chr[[i]],
      begin = event$begin[[1]],
      end = event$end[[1]],
      length = event$end[[1]] - event$begin[[1]] + 1
    )
  )
}

print(cna_regions)


# ------------------------------------------------------------
# 3. Parent structure
# ------------------------------------------------------------

parents <- as.integer(rep(0, nrow(clones)))

for (node in seq_len(nrow(clones))) {
  if (node == 1L) {
    parents[[node]] <- 0L
  } else {
    parent_idx <- match(clones$parent_mutant[[node]], clones$mutant)
    
    if (is.na(parent_idx) || parent_idx >= node) {
      stop("Invalid parent order in clones table.")
    }
    
    parents[[node]] <- parent_idx
  }
}

print(parents)


# ------------------------------------------------------------
# 4. Descendant mutants for MRCA calculation
# ------------------------------------------------------------

node_descendant_mutants <- vector("list", nrow(clones))

for (node in seq_len(nrow(clones))) {
  node_descendant_mutants[[node]] <- clones$mutant[
    get_node_descendants(node, parents)
  ]
}

names(node_descendant_mutants) <- clones$mutant

print(node_descendant_mutants)


# ------------------------------------------------------------
# 5. Events table: MRCA + CNA
# ------------------------------------------------------------

event_rows <- list()
event_idx <- 1L

for (node in seq_len(nrow(clones))) {
  mutant_name <- clones$mutant[[node]]
  label <- clones$cell_label[[node]]
  
  ids_for_mrca <- unlist(clone_ids[node_descendant_mutants[[node]]])
  
  node_time <- get_mrca_time(
    phylo_forest,
    ids_for_mrca,
    paste0("MRCA ", label)
  )
  
  # MRCA event
  event_rows[[event_idx]] <- tibble(
    node = as.integer(node),
    event_offset = 1L,
    event_name = paste0("MRCA_", mutant_name),
    event_type = "mrca",
    cell_label = label,
    cna_name = NA_character_,
    true_time = node_time
  )
  
  event_idx <- event_idx + 1L
  
  # CNA events associated with this clone
  node_cnas <- cna_regions %>%
    filter(.data$target_label == label)
  
  if (nrow(node_cnas) > 0) {
    for (j in seq_len(nrow(node_cnas))) {
      event_rows[[event_idx]] <- tibble(
        node = as.integer(node),
        event_offset = as.integer(j + 1L),
        event_name = node_cnas$cna_name[[j]],
        event_type = "cna",
        cell_label = label,
        cna_name = node_cnas$cna_name[[j]],
        true_time = node_time
      )
      
      event_idx <- event_idx + 1L
    }
  }
}

events <- bind_rows(event_rows) %>%
  arrange(node, event_offset) %>%
  mutate(
    true_pseudotime = .data$true_time / Tf,
    event_index = row_number()
  ) %>%
  distinct(event_name, event_type, node, .keep_all = TRUE) %>%
  arrange(node, event_offset) %>%
  mutate(event_index = row_number())

print(events)


# ------------------------------------------------------------
# 6. True event pseudotimes
# ------------------------------------------------------------

true_tevents <- events$true_pseudotime
names(true_tevents) <- events$event_name

print(true_tevents)

# ------------------------------------------------------------
# 7. Common normal regions
# ------------------------------------------------------------

sampled_clone_order <- clones$cell_label

common_normal_regions <- make_common_normal_regions(
  cna_regions = cna_regions,
  sampled_clone_order = sampled_clone_order
)

total_normal_length <- common_normal_regions %>%
  mutate(region_length = end - begin + 1L) %>%
  group_by(cell_label) %>%
  summarise(
    total_11_length = sum(region_length),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 8. Mutations
# ------------------------------------------------------------

mutation_row <- phylo_forest$get_sampled_cell_mutations()

mutation_key_cols <- c("chr", "from", "ref", "alt")

mut_clean <- build_mutation_clean(
  mutation_row,
  cell_label_map,
  sampled_clone_order
)

mutation_presence <- build_mutation_presence(mut_clean)

# ------------------------------------------------------------
# 9. Descendant cell labels for downstream alpha/beta
# ------------------------------------------------------------

node_descendant_labels <- vector("list", nrow(clones))

for (node in seq_len(nrow(clones))) {
  node_descendant_labels[[node]] <- clones$cell_label[
    get_node_descendants(node, parents)
  ]
}

names(node_descendant_labels) <- clones$cell_label

print(node_descendant_labels)

# ------------------------------------------------------------
# 10. Initialise downstream containers
# ------------------------------------------------------------

alpha_beta_table <- tibble()

event_lengths <- numeric(nrow(events))
event_karyotypes <- integer(nrow(events))

mutations <- c()
coeff <- c()

idx_first_mutation <- integer(nrow(events))
mutation_labels <- c()

for (i in seq_len(nrow(events))) {
  label <- events$cell_label[[i]]
  node <- events$node[[i]]
  parent <- parents[[node]]
  required_labels <- node_descendant_labels[[node]]
  compatible_labels <- if (parent == 0L) character(0) else clones$cell_label[[parent]]
  incompatible_labels <- setdiff(setdiff(sampled_clone_order, required_labels), compatible_labels)
  idx_first_mutation[[i]] <- length(mutations) + 1L
  
  if (events$event_type[[i]] == "mrca") {
    m_coeff <- 2.0
    # l_tot,i = regions that were 1:1 during THIS clone's own branch = the
    # clone's own normal regions. Descendant/sibling CNAs are later events and
    # must NOT shrink it; intersecting over the subtree (required_labels) lets
    # brother clones erode the ancestor's normal length and under-counts m.
    branch_normal_regions <- get_diploid_regions(
      common_normal_regions = common_normal_regions,
      region_labels = label
    )
    L_11 <- sum(branch_normal_regions$end - branch_normal_regions$begin + 1L)

    if (length(L_11) != 1 || is.na(L_11) || L_11 <= 0) {
      stop("Invalid 1:1 length for ", label)
    }

    event_lengths[[i]] <- L_11
    event_karyotypes[[i]] <- 0L
    m_count <- 0L

    if (node != 1L) {
      m_table <- count_m(
        mutation_presence = mutation_presence,
        common_normal_regions = common_normal_regions,
        clone_label = required_labels,
        parent_label = clones$cell_label[[parent]],
        region_label = label,
        n_cells_per_clone = n_cells_per_clone
      )
      m_count <- nrow(m_table)
    }

    mutations <- c(mutations, as.integer(m_count))
    coeff <- c(coeff, m_coeff)
    mutation_labels <- c(mutation_labels, paste0("m_", label))
  } else {
    alpha_coeff <- 1.0
    beta_coeff <- 3.0
    region <- cna_regions %>% filter(.data$cna_name == events$cna_name[[i]])
    event_lengths[[i]] <- region$length[[1]]
    event_karyotypes[[i]] <- 1L

    alpha_beta <- count_alpha_beta(
      mut_clean = mut_clean,
      cna_regions = cna_regions,
      cna_name = events$cna_name[[i]],
      target_label = label,
      n_cells_per_clone = n_cells_per_clone,
      required_labels = required_labels,
      excluded_labels = incompatible_labels,
      compatible_labels = compatible_labels,
      redistribute_non_amp = TRUE,
      mutation_presence = mutation_presence
    )

    alpha_beta_table <- bind_rows(alpha_beta_table, alpha_beta)
    mutations <- c(
      mutations,
      as.integer(round(alpha_beta$N_alpha_total[[1]])),
      as.integer(round(alpha_beta$N_beta_total[[1]]))
    )
    coeff <- c(coeff, alpha_coeff, beta_coeff)
    mutation_labels <- c(
      mutation_labels,
      paste0("alpha_", events$cna_name[[i]]),
      paste0("beta_", events$cna_name[[i]])
    )
  }
}

print(alpha_beta_table)

mutation_count_summary <- tibble(
  equation = mutation_labels,
  count = as.integer(mutations)
)
print(mutation_count_summary)

number_of_events <- as.integer(vapply(seq_len(nrow(clones)), function(node) {
  sum(events$node == node)
}, integer(1)))

idx_first_event <- as.integer(vapply(seq_len(nrow(clones)), function(node) {
  min(which(events$node == node))
}, integer(1)))

# ------------------------------------------------------------
# 11. Stan
# ------------------------------------------------------------

stan_data <- list(
  Tf_real = Tf,
  t0_real = 0,
  mu = mu_SNV,
  phi = clone_abundance_nb_phi,
  
  Nnodes = nrow(clones),
  parents = parents,
  number_of_events = number_of_events,
  idx_first_event = idx_first_event,
  
  Nevents = nrow(events),
  karyotypes = as.integer(event_karyotypes),
  lengths = as.numeric(event_lengths),
  idx_first_mutation = as.integer(idx_first_mutation),
  
  Nequations = length(mutations),
  coeff = as.numeric(coeff),
  mutations = as.integer(mutations),
  
  Ncells = as.integer(sum(n_cells_per_clone$n_cells_total)),
  Ncells_per_clone = as.integer(n_cells_per_clone$n_cells_total),
  w_bar_prior_mean = w_bar_prior_mean,
  w_bar_prior_sd = w_bar_prior_sd,
  sigma_log_w_rate = sigma_log_w_rate,
  sigma_nu_rate = 2
)

stan_file <- file.path(
  "/Users/francescofavagrossa/Documents/Projects/Thesis/stan_models/n growth fixed nb clock coupled.stan"
)

model <- cmdstan_model(stan_file)
fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 4000,
  iter_sampling = 6000,
  adapt_delta = 0.99,
  max_treedepth = 14,
  seed = 42
)

print(fit$summary(c("tevents", "w", "w_node_out", "delta_log_w_out", "sigma_log_w")))

true_tevents

# ------------------------------------------------------------
# 12. Posterior checks 
# ------------------------------------------------------------

true_events <- tibble(
  index = seq_along(true_tevents),
  label = names(true_tevents),
  true_value = as.numeric(true_tevents)
)

posterior_tevents <- fit$draws("tevents", format = "draws_df") %>%
  as_tibble() %>%
  select(starts_with("tevents[")) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(
    cols = starts_with("tevents["),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    index = as.integer(str_extract(parameter, "\\d+")),
    label = true_events$label[index]
  ) %>%
  left_join(true_events, by = c("index", "label")) %>%
  mutate(label = factor(label, levels = true_events$label))

p_tevents <- make_posterior_ci_plot(
  posterior_df = posterior_tevents,
  plot_title = "Posterior event times",
  x_label = "Pseudotime",
  x_limits = c(0, 1)
)

show_plot(p_tevents)

ggsave(
  filename = "figures 4 DLP/Posterior events time.png",
  plot = p_tevents,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

true_growth <- tibble(
  index = seq_len(nrow(clones)),
  label = clones$cell_label,
  true_value = c(
    w[1],
    w[2],
    w[3],
    w[4]
  )
)

posterior_growth <- fit$draws("w_node_out", format = "draws_df") %>%
  as_tibble() %>%
  select(starts_with("w_node_out[")) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(
    cols = starts_with("w_node_out["),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    index = as.integer(str_extract(parameter, "\\d+")),
    label = true_growth$label[index]
  ) %>%
  left_join(true_growth, by = c("index", "label")) %>%
  mutate(label = factor(label, levels = true_growth$label))

p_growth <- make_posterior_ci_plot(
  posterior_df = posterior_growth,
  plot_title = "Posterior growth rates",
  x_label = "Growth rate"
)

show_plot(p_growth)

ggsave(
  filename = "figures 4 DLP/Posterior growth rates.png",
  plot = p_growth,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 13. Save bias results (posterior vs true: event times + growth rates)
# ------------------------------------------------------------

source('/Users/francescofavagrossa/Documents/Projects/Thesis/Helpers and Utils/bias_table.R')
save_bias(fit, true_events, true_growth, "Clone 4 DLP")
res <- combine_bias()
print(res$wide, width = Inf)
