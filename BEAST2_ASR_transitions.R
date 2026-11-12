library(ape)
library(treeio)
library(tidytree)
library(phytools)
library(dplyr)

tree_obj <- read.beast("final.tree")
tree_phylo <- tree_obj@phylo
meta <- read.csv("IC2_uniform_year_balanced.csv")
rownames(meta) <- meta$ID

latest_year <- 2025.0

states <- meta[tree_phylo$tip.label, "Type"]
names(states) <- tree_phylo$tip.label

fit_er <- fitMk(tree_phylo, states, model = "ER")
fit_sym <- fitMk(tree_phylo, states, model = "SYM")
fit_ard <- fitMk(tree_phylo, states, model = "ARD")

aic_df <- data.frame(
  Model = c("ER", "SYM", "ARD"),
  logLik = c(logLik(fit_er), logLik(fit_sym), logLik(fit_ard)),
  AIC = c(AIC(fit_er), AIC(fit_sym), AIC(fit_ard)),
  stringsAsFactors = FALSE
) %>%
  mutate(delta_AIC = round(AIC - min(AIC), 3)) %>%
  arrange(AIC)

best_model <- aic_df$Model[1]

fit <- rerootingMethod(tree_phylo, states, model = best_model)
node_states <- colnames(fit$marginal.anc)[apply(fit$marginal.anc, 1, which.max)]
all_states <- c(states, node_states)

tree_df <- as_tibble(tree_obj)
n_tips <- length(tree_phylo$tip.label)

height_col <- if ("height_mean" %in% colnames(tree_df)) "height_mean" else "height"

node_years <- setNames(
  round(latest_year - tree_df[[height_col]], 2),
  tree_df$node
)

edges <- tree_phylo$edge
events <- list()

for (i in 1:nrow(edges)) {
  parent_node <- edges[i, 1]
  child_node <- edges[i, 2]

  from_state <- all_states[parent_node]
  to_state <- all_states[child_node]

  if (from_state != to_state) {
    child_label <- if (child_node <= n_tips) {
      tree_phylo$tip.label[child_node]
    } else {
      paste0("Node_", child_node)
    }

    parent_yr <- node_years[as.character(parent_node)]
    child_yr <- node_years[as.character(child_node)]
    mid_yr <- round((parent_yr + child_yr) / 2, 2)

    events[[length(events) + 1]] <- data.frame(
      Transition_ID = paste0("Trans_", length(events) + 1),
      Model_Used = best_model,
      Parent_Node = parent_node,
      Child_Node = child_node,
      Child_Label_or_Sample = child_label,
      Source_Niche = from_state,
      Sink_Niche = to_state,
      Direction = paste(from_state, "->", to_state),
      Parent_Year = parent_yr,
      Child_Year = child_yr,
      Transition_Year_Mid = mid_yr,
      Branch_Duration_Years = round(child_yr - parent_yr, 2),
      stringsAsFactors = FALSE
    )
  }
}

transition_table <- bind_rows(events) %>%
  arrange(Transition_Year_Mid)

write.csv(transition_table, "IC2_Cross_Niche_Transition_Events_Detailed.csv", row.names = FALSE)
write.csv(aic_df, "IC2_Transition_Model_AIC_Comparison.csv", row.names = FALSE)