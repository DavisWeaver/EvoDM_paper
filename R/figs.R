#Always load the python stuff in a clean session\
reticulate::use_python("C:/Program Files/Python37/python.exe", required = TRUE) #python path may vary

reticulate::source_python("load.py")

#packages
library(here)
library(stringr)
library(magrittr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(visNetwork)
library(scales)
library(ggraph)
library(igraph)
library(tidygraph)
library(dplyr)

#Helper functions
source(here("R", "prep_figs.R"))
source(here("R", "clean_evol_out.R"))
#set theme
theme_set(theme_classic())


#Do we want a different intro figure
#I actually think this belongs in a different plot
# g3 <- ggplot(df_big_wide, aes(x = index_main, y = fitness_ma, color = condition)) +
#   geom_line() + 
#   labs(x = "Evolutionary time", 
#        y = "Population fitness \n (150 step moving average)", 
#        tag = "B") +
#   theme(text = element_text(size = 16))

##############Mira Performance Plot####################
#All code to generate mira performance plot
mira_perf_plot <- function(out, df_combos_long, df_combos, reset_every = 20, 
                           replicates = 75, quad_mut = FALSE, total_res = FALSE) {
  
  mem_df <- out[[7]] %>% 
    filter(replicate <= replicates,
           quadmut == quad_mut,
           totalres == totalres)
  df_big <- out[[3]]
  df <- out[[2]] %>% 
    filter(reset == reset_every, replicate %in% 1:replicates,
           quadmut == quad_mut,
           totalres == total_res)
  
  
  mem_short <- mem_df %>% group_by(original_ep) %>% 
    summarise(
      upper_bound = quantile(fitness, probs = c(0.95)), 
      lower_bound = quantile(fitness, probs = c(0.05)),
      fitness = mean(fitness)
    )
  
  g1 <- ggplot(mem_short, aes(x = original_ep, y = fitness)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.3) + 
    theme(text = element_text(size = 16)) + 
    labs(x = "original episode", y = "population fitness", tag = "A")
  
  
  g2 <- ggplot(df, aes(x = replicate, y = benefit)) + 
    geom_boxplot() + 
    xlab("Replicate") + 
    ylab("Fitness Improvement \n (random - evoDM)") +
    labs(tag = "B") + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size = 16), 
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  g4 <- ggplot(df, aes(x = replicate, y = distance_optimal)) + 
    geom_boxplot() + 
    xlab("Replicate") + 
    ylab("Fitness decrement \n (evoDM - optimal)") +
    labs(tag = "D") + 
    theme(axis.text.x = element_blank(),
          text = element_text(size = 16),
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  #now histogram 
  df_hist <- df %>% 
    pivot_longer(cols = c(evodm, naive, optimal_policy), 
                 names_to =  "condition", 
                 values_to = "population_fitness") 
  
  g3 <- ggplot(df_hist, aes(x = population_fitness, fill = condition)) + 
    geom_density(alpha = 0.6) + 
    geom_density(data = df_combos, aes(x=average_fitness, 
                                       fill = "two-drug \n combinations"), alpha = 0.6) + 
    scale_fill_viridis_d(labels = c("evodm", "random", "optimal", "two-drug \n combinations")) + 
    scale_x_continuous(expand = c(0,0), limits = c(0,3.3)) + 
    scale_y_continuous(expand = c(0,0)) + 
    labs(y = "density", x = "population fitness", tags = "C") + 
    theme(text = element_text(size = 16),
          legend.position = "top", 
          legend.title = element_blank())
  
  g <- plot_grid(g1, g2, g3, g4, nrow = 2, rel_heights = c(1.25,1), align = "vh", axis = "blr")
  
  return(g)
  
}

mira_perf_sens <- function(out) {
  mem_df <- out[[7]] %>% 
    mutate(
      condition = case_when(
        quadmut ~ "quadmut",
        totalres ~ "totalres",
        TRUE ~ "base_case")
      ) %>% 
    filter(original_ep < 500)
  
  mem_short <- mem_df %>% group_by(original_ep, condition) %>% 
    summarise(
      upper_bound = quantile(fitness, probs = c(0.95)), 
      lower_bound = quantile(fitness, probs = c(0.05)),
      fitness = mean(fitness)
    )
  
  g1 <- ggplot(mem_short, aes(x = original_ep, y = fitness, color = condition)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.3) + 
    theme(text = element_text(size = 16)) + 
    labs(x = "original episode", y = "population fitness", tag = "A")
  
  return(g1)
  
}

###############Mira Policy Plot################# 

mira_policy_plot <- function(out) {
  df_big_wide <- out[[4]]
  
  policies <- out[[1]] %>% filter(reset == 20) %>% 
    mutate(replicate = fct_reorder(replicate, benefit, .desc=TRUE), 
           replicate_num = as.numeric(replicate))
  
  joint_prob = out[[6]]
  
  out <- prep_opt_policy(df_big_wide)
  opt_props <- out[[1]]
  opt_props_sum <- out[[2]]
  opt_props_sum$replicate_num = -1
  
  
  g1a <- ggplot(policies, aes(x = replicate_num, y = factor(drug), fill = prob_selection)) + 
    geom_tile() + 
    scale_fill_gradient(low = "navy", high = "yellow") + 
    scale_x_continuous(expand = c(0,0), limits = c(-5,100)) + 
    scale_y_discrete(expand = c(0,0)) + 
    labs(y = "drug", x = "Replicate", tag = "A") + 
    theme(axis.text.x = element_blank(), 
          text = element_text(size = 16),
          legend.title = element_blank(),
          legend.position = "top",
          legend.margin = margin(), 
          legend.key.width = unit(2.5, units = "cm")) 
  
  g_opt <- ggplot(opt_props_sum, aes(x = replicate_num, y = factor(drug), fill = prob_selection)) + 
    geom_tile() + 
    scale_fill_gradient(low = "navy", high = "yellow") + 
    theme(axis.line = element_blank(),
          axis.text = element_blank(), 
          legend.position = "none",
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          margin()) + 
    scale_y_discrete(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) + theme_nothing()
  
  g1 <- g1a + 
    annotation_custom(ggplotGrob(g_opt), xmin = -5, xmax = 0, ymin = 0.5, ymax = 15.5)
  
  ###Policy Network
  g <- compute_combo_graph(joint_probs = joint_prob)
  
  g2 <- ggraph(g) + 
    geom_node_point(size = 20, fill = "navy", colour = "navy") + 
    geom_edge_fan(aes(edge_width = weight, 
                      label = paste0("prob = ", round(weight, digits = 3))), strength = 2,
                  arrow = arrow(ends = "last", length = unit(1, "cm"), angle = 15, type = "closed"), 
                  show.legend = FALSE,
                  angle_calc = 'along',
                  label_dodge = unit(-8, 'mm'),
                  label_size = 6,
                  end_cap = circle(10, 'mm')) +
    geom_node_text(aes(label = label), size = 10, nudge_x = -0.23, nudge_y = 0.015) +
    labs(tag = "B") + 
    theme(text = element_text(size = 16))
  
  #now take a look at action by state heatmap
  policies_wide =  out[[8]]
  g3 <- ggplot(policies_wide, aes(x = drug, y = state, fill = prob_selection)) + 
    geom_tile()
    

  plot_grid(g1,g2)
  #Plot network figure
 
}

#########Two-drug combos plot######## 
plot_2drug <- function(df_combos_long, df_combos_big) {
  
  #more plots of the combos
  #another plot about the two-drug combos
  g1 <- ggplot(df_combos_long, aes(x = combo, y = average_fitness)) + 
    #geom_boxplot() + geom_quasirandom() + 
    #geom_violin() + 
    geom_boxplot() + 
    geom_rect(aes(ymin = 0, ymax = 0.25, fill = correl, 
                  xmin = as.numeric(combo) -0.5, xmax = as.numeric(combo) +0.5)) + 
    theme(text = element_text(size = 16), 
          axis.text.x = element_blank()) + 
    scale_y_continuous(expand = c(0,0)) + 
    labs(x = "two-drug combination", 
         y = "mean population fitness", 
         fill = "landscape \n correlation") + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")
  
  g2 <- ggplot(df_combos_big, aes(x = factor(starting_state), y = average_fitness)) + 
    geom_boxplot() +geom_quasirandom(alpha = 0.4) +
    labs(x = "starting genotype", 
         y = "mean population fitness",) + 
    theme(text = element_text(size = 16)) 
  
  g <- plot_grid(g1,g2)
  return(g)
}

################hyperparameter tuning######################

###gamma, batch size, lr, update_target_every all need to be represented here

plot_hp_tune <- function() {
  
  mems_lr <- prep_lr() #%>% 
  #filter(lr > 1e-5)
  
  g1 <- ggplot(mems_lr, aes(x = original_ep, y = fitness, color = lr)) + 
    geom_line() + 
    theme(text = element_text(size = 16)) #+ 
  facet_wrap(~lr)
  
  g2 <-ggplot(mems_lr, aes(x = factor(lr), y = fitness)) + 
    geom_violin()
  
  mems_batchsize <- prep_batch_size()
  
  g3 <- ggplot(mems_batchsize, aes(x = original_ep, y = fitness, color = batch_size)) + 
    geom_line() + 
    theme(text = element_text(size = 16))
  
  g4 <- ggplot(mems_batchsize, aes(x = factor(batch_size), y = fitness)) + 
    geom_violin()
  
  mems_updatetarget <- prep_update_target()
  g5 <- ggplot(mems_batchsize, aes(x = original_ep, y = fitness, color = update_target)) + 
    geom_line() + 
    theme(text = element_text(size = 16))
  g <- plot_grid(g1,g3)
  return(g)
}

#### Procedurally generated landscapes figure#####
#### MDP parameter sweep######
plot_random_mdp <- function(){
  load(file = here("data", "results", "normalized_mdp_sweep_random.Rda"))
  df_prob = out[[1]]
  df_sum = out[[2]]
  df_sum <- df_sum %>% group_by(sigma, num_drugs, N) %>% 
    summarise(benefit = mean(benefit), sd = sd(benefit))
  
  g1 <- ggplot(df_sum, aes(x = sigma, y = num_drugs, fill = benefit)) + 
    geom_tile() + 
    facet_wrap(~N, labeller = labeller(.cols = label_both)) + 
    scale_fill_viridis_c() + 
    labs(x = "sigma (epistasis coefficient)",
         y = "number of available drugs", 
         fill = "fitness benefit (mdp - random)") + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    theme(legend.position = "top", 
          legend.key.width = unit(2.5, "cm"),
          text = element_text(size = 26)) + 
    guides(fill = guide_colorbar(title.position = "top",
                                 title.hjust = unit(-0, "cm")))
  
  return(g1)
  
}

####MDP Mira stuff####
plot_mira_mdp <- function(){
  load(here("data", "results", "not_normalized_mdp_sweep.Rda"))
  df_prob <-out[[1]]
  df_sum <- out[[2]]
  
  g1 <- ggplot(df_sum, aes(x = gamma, y = fitness)) + geom_line() + 
    labs(y = "Average Fitness") + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    theme(text = element_text(size= 16)) 
  
  g2 <- ggplot(df_prob, aes(x = gamma, y = prob_selection, fill = factor(drug))) + 
    geom_area() + 
    labs(y = "probability of selection", fill = "drug") + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    theme(text = element_text(size= 16),
          legend.position = "top") 
  return(plot_grid(g1, g2))
}

##Load/ clean data
##Heatmap describing repeatability of results
out_hm <- prep_mira_heatmap(total_eps = 500)

out_combos <- prep_twodrug_sweep()
df_combos_long <- out_combos[[1]]
df_combos <- out_combos[[3]]
df_combos_big <- out_combos[[2]]

###Actually plot everything

#Performance Plot
#base case
g <- mira_perf_plot(out = out_hm, df_combos_long = df_combos_long, 
                    df_combos = df_combos, reset = 20, replicates = 100)
png(filename = here("figs", "mira_landscapes20.png"), 
    width = 1000, height = 600)
g
dev.off()

#quadmut starting position
g <- mira_perf_plot(out = out_hm, df_combos_long = df_combos_long, 
                    df_combos = df_combos, reset = 20, replicates = 100, 
                    quad_mut = TRUE)
png(filename = here("figs", "mira_landscapes20_quadmut.png"), 
    width = 1000, height = 600)
g
dev.off()


#optimize over total resistance instead of single drug resistance
g <- mira_perf_plot(out = out_hm, df_combos_long = df_combos_long, 
                    df_combos = df_combos, reset = 20, replicates = 100, 
                    quad_mut = FALSE, total_res = TRUE)
png(filename = here("figs", "mira_landscapes20_totalres.png"), 
    width = 1000, height = 600)
g
dev.off()
#policy plot


g <- mira_policy_plot(out = out_hm)
png(filename = here("figs", "mira_policy.png"), 
    width = 1000, height = 600)
g
dev.off()


#random landscapes figure

##Supps##

#Sensitivity analyses around training on mira landscapes
g <- mira_perf_sens(out=out_hm)
png(filename = here("figs", "mira_perf_sens.png"), 
    width = 1000, height = 600)
g
dev.off()


#two-drug sweep
g <- plot_2drug(df_combos_long = df_combos_long, df_combos_big = df_combos_big)
png(filename = here("figs", "two_drug_sweep.png"), 
    width = 1000, height = 600)
g
dev.off()

#MDP Random Sweep 4
png(filename = here("figs", "random_mdp_sweep.png"), width = 800, height = 800)
plot_random_mdp()
dev.off()

#MDP Mira Sweep
png(filename = here("figs", "mira_mdp_sweep.png"), 
    width = 800, height = 400)
plot_mira_mdp()
dev.off()

#hyper-parameter tuning
png(filename = here("figs", "hp_tune_plot.png"), 
    width = 1000, height = 600)
plot_hp_tune()
dev.off()
