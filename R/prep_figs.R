
#Figure that looks at policies for a given agent
prep_performance_overview <- function(N=5, landscapes = "random", files = NULL, 
                                      lr_range= NULL) {
  
  
  add_filestring_vars <- function(file, df) {
    
    df$landscapes <- str_extract(file, "\\d+(?=N)") #grab the landscape set from the filename
    df$train_input <- str_extract(file, "(?<=_)[:alpha:]{2,3}(?=[\\d\\.])")
    df$replicate <- str_extract(file, "(?<=fit)\\d+")
    df$totalres <- str_detect(file, "totalres")
    df$quadmut <- str_detect(file, "quadmut")
    
    if(landscapes == "lr") {
      lr_num <- as.numeric(str_extract(file, "(?<=lr)\\d+"))
      df$lr <- lr_range[lr_num]
      
    }
    
    if(landscapes == "updatetarget") {
      df$update_target <- as.numeric(str_extract(file, "(?<=updatetarget)\\d+"))
    }
    
    if(landscapes == "batchsize") {
      df$batch_size <- as.numeric(str_extract(file, "(?<=batchsize)\\d+"))
    }
    if(landscapes == "gamma") {
      df$gamma <- as.numeric(str_extract(file, "(?<=gamma)\\d+\\.\\d+"))
    }
    
    if(landscapes == "mira") {
      df$reset <- str_extract(files[i], "(?<=reset)\\d+")
      df$replicate_mod <- str_extract(file, "(?<=mira)\\d+")
    }
    return(df)
  }
  
  #check if the user provided a vector of filenames
  if(is.null(files)) {
    files <- list.files(here("data", "results"), pattern = landscapes)
  }
  #omit the hyperparameter sweeps if we are looking to prep the mira heatmap figure
  if(landscapes == "mira") {
    files <- files[!grepl("lr\\d", files)] #totalres has an lr in it hilariously
    files <- files[!grepl("updatetarget", files)]
    files <- files[!grepl("batchsize", files)]
    files <- files[!grepl("gamma", files)]
  }
  df_list <- list()
  policies <- list()
  mems <- list()
  for(i in 1:length(files)) {
    #clean  the master_memory dataframes
    load(here("data", "results", files[i]))
    
    test <- try(out_list[[1]], silent = TRUE)
    if(inherits(test, 'try-error')) {
      next
      df = cleaned_out[[1]] %>%
        tidyr::pivot_wider(names_from = drug_regime, values_from = moving_prob)
      df <- add_filestring_vars(file = files[i], df)
    } else {
      cleaned_out <- out_list[[1]]
      df <- out_list[[2]]
      df <- add_filestring_vars(file = files[i], df) 
      df$condition = "evodm"
      
      
      #these are the other conditions
      df2 <- cleaned_out[[1]] 
      df2 <- add_filestring_vars(file = files[i], df2) %>% 
        filter(condition != "evodm") %>% 
        tidyr::pivot_wider(names_from = drug_regime, values_from = moving_prob)
      
      df <- bind_rows(df, df2)
    }
    
    #pivot so we don't have num_drugs * the number of obs we want
    #df <- tidyr::pivot_wider(df, names_from = drug_regime, values_from = moving_prob)
    df_list[[i]] <- df
    
    #join the policy dfs.
    policy <- cleaned_out[[2]]
    policy <- add_filestring_vars(file = files[i], df = policy)
    #now handle the optimal policies
    opt_policy <- cleaned_out[[3]] %>% 
      as.data.frame() %>% 
      mutate(state = 1:nrow(cleaned_out[[3]])) %>% 
      pivot_longer(cols = -state, 
                   names_to = "action", 
                   names_prefix = "V", 
                   values_to = "optimal_action") %>%
      mutate(action = as.numeric(action))
    
    policy <- left_join(policy, opt_policy)
    
    policies[[i]] <- policy
    
    mems_i <- out_list[[3]][[2]]
    mems_i <- add_filestring_vars(file = files[i], mems_i)
    mems[[i]] <- mems_i
  }
  
  df <- bind_rows(df_list)
  policies <- bind_rows(policies)
  mems <- bind_rows(mems)
  
  #more work to do with df
  df <- df %>% 
    filter((train_input != "sv" | train_input == "sv" & condition == "evodm")) %>% 
    mutate(condition = ifelse(train_input == "sv", "evodm_sv", condition))
  return(list(df, policies, opt_policy, mems))
}

prep_mira_heatmap <- function(total_eps = 500) {
  
  #drug index table
  drug_index <- data.frame(action = 1:15, 
                           drug = c("AMP", "AM", "CEC", "CTX", "ZOX", "CXM", 
                                    "CRO", "AMC", "CAZ", "CTT", "SAM", "CPR",
                                    "CPD", "TZP", "FEP"))
  out <- prep_performance_overview(N=4, landscapes = "mira")
  
  #first prep the performance of all 100 runs for the lollipop
  df <- out[[1]] %>% filter(!is.na(replicate), 
                            condition != "evodm" | ep_number > total_eps) %>%
    mutate(episode = ifelse(ep_number>total_eps, ep_number-total_eps, ep_number)) %>%
    mutate(replicate = as.numeric(replicate),
           replicate_mod = as.numeric(replicate_mod),
           replicate = replicate + (replicate_mod-1))
  #compute joint probability of a given drug being selected for all replicates
  
  
  joint_probs <- compute_joint_probability(df)
  
  #finish cleaning performance data
  df <- df  %>% 
    group_by(replicate, condition, reset, episode, quadmut, totalres) %>% 
    summarise(fitness = mean(average_fitness)) %>% 
    pivot_wider(names_from = "condition", values_from = "fitness") %>% 
    mutate(benefit = naive - evodm,
           distance_optimal = optimal_policy - evodm,
           optimal_benefit = naive - optimal_policy
    ) %>% 
    ungroup() 
  
  #need replicate to still be numeric at this stage
  df_small <- df %>% 
    group_by(replicate, reset) %>% 
    summarise(benefit = mean(benefit))
  
  #re-order to be in descending order of benefit 
  df <- df %>%
    mutate(
      replicate = factor(replicate),
      replicate = fct_reorder(replicate, benefit, .desc=TRUE)
    ) 
  
  #now grab data on just one replicate for a snapshot
  #use the middle replicate
  rep_target <- levels(df$replicate)[round(length(levels(df$replicate))/2)]
  
  df_full <- out[[1]] %>% 
    filter(replicate == rep_target)
  
  df_long <- df_full %>% 
    pivot_longer(cols = matches("\\d"), names_to = "action", 
                 values_to = "moving_prob") %>% 
    rename(drug_regime = drug) %>% 
    mutate(action = as.numeric(action)) %>%
    left_join(drug_index) %>% 
    filter(condition == "evodm")
  
  
  ##Process the policies data - many of the same steps as for the performance data
  policies <- out[[2]] %>% 
    filter(!is.na(replicate)) %>% 
    mutate(replicate = as.numeric(replicate),
           replicate_mod = as.numeric(replicate_mod),
           replicate = replicate + (replicate_mod-1)) %>%
    filter(episode == total_eps, !quadmut, !totalres) 
  
  #
  policies_state <- policies %>% 
    group_by(state, action) %>% 
    summarise(prob_selection = mean(prob_selection)) %>% 
    left_join(drug_index)
    
    
  #finish summarising the policy df
  policies <- policies %>% 
    group_by(action, replicate, reset) %>% 
    summarise(prob_selection = mean(prob_selection)) %>% 
    ungroup() %>% 
    left_join(drug_index) %>%
    left_join(df_small) %>% 
    mutate(replicate = factor(replicate),
           replicate = fct_reorder(replicate, benefit, .desc = TRUE))
  
  
  #add drug index to optimal policy
  drug_index <- drug_index %>% mutate(action = action-1)
  opt_policy <- out[[3]] %>% 
    rename(time_step = action,
           action = optimal_action) %>% 
    left_join(drug_index)
  
  mems <- out[[4]] %>%
    mutate(replicate= as.numeric(replicate))
  
  
  return(list(policies, df, df_long, df_full, opt_policy, joint_probs, mems,policies_state))
  #now compute 
}
#function to compute the joint probability of a drug being selected in the 500 episode 
compute_joint_probability <- function(df) {
  #this is going to be different for every reset
  #nutso looping structure here
  
  reset_vec <- as.numeric(unique(df$reset))
  action_vec <- unique(df$drug)
  #list for outer loop
  out_list <- list()
  for(i in 1:length(reset_vec)) {
    df_i <- df %>% filter(reset == reset_vec[i], condition == "evodm")
    df_i <- df_i %>% group_by(ep_number, replicate) %>% 
      summarise(drug = paste(drug, collapse = ",", sep = ","))
    
    out_inner <- list()
    for(j in 1:length(action_vec)) {
      
      str_pattern = paste0("(?<=", action_vec[j], "),\\d+")
      out = str_extract_all(df_i$drug, pattern = str_pattern)
      #setup counter
      #oh my god another loop help me 
      event_counts = replicate(15, 0)
      drug2 = 1:15
      for(z in 1:length(out)) {
        events <- as.numeric(gsub(",", "", out[[z]]))
        for(y in 1:length(events)) {
          event_counts[events[y]] =  event_counts[events[y]] + 1 #count all the times action b followed action a
        }
      }
      
      #get number of total events so we can compute the joint prob
      num_events = length(unique(df_i$replicate)) * length(unique(df_i$ep_number)) * (reset_vec[i]-2)
      df_j <- data.frame(drug1 = action_vec[j], drug2 = drug2, event_counts = event_counts, 
                         joint_prob = event_counts/num_events)
      out_inner[[j]] <- df_j
      
      
    }
    #clean and store results of inner loop
    df_i <- bind_rows(out_inner)
    df_i$reset <- reset_vec[[i]]
    out_list[[i]] <- df_i
  }
  #clean and process results of two loops
  df <- bind_rows(out_list)
  return(df)
}

compute_combo_graph <- function(joint_probs, reset = 20){
  
  drug_index <- data.frame(name = 1:15, 
                           label = c("AMP", "AM", "CEC", "CTX", "ZOX", "CXM", 
                                     "CRO", "AMC", "CAZ", "CTT", "SAM", "CPR",
                                     "CPD", "TZP", "FEP"),
                           label_long = c("Ampicillin", "Amoxicillin", "Cefaclor", "Cefotaxime",
                                          "Ceftizoxime", "Cefuroxime", "Ceftriaxone",
                                          "Amoxicillin + Clavulanic acid", "Ceftazidime",
                                          "Cefotetan", "Ampicillin + Sulbactam", 
                                          "Cefprozil", "Cefpodoxime", 
                                          "Pipercillin + Tazobactam", "Cefepime")
  )
  
  #Filter joint_probs to only have the reset distance under study. 
  df_i <- filter(joint_probs, reset == reset)
  
  #put together an adjacency matrix
  edge_mat <- df_i %>% select(drug1, drug2, joint_prob) %>% 
    mutate(joint_prob = ifelse(joint_prob < 0.01, 0, joint_prob)) %>%
    pivot_wider(names_from = drug2, values_from = joint_prob) %>% 
    arrange(-desc(drug1)) %>% 
    select(-drug1) %>%
    as.matrix()
  
  g <- igraph::graph_from_adjacency_matrix(edge_mat, weighted=  TRUE, mode = "directed")
  E(g)$width <- E(g)$weight*100
  
  
  g <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    mutate(degree = centrality_degree(), 
           name = as.integer(name)) %>% 
    left_join(drug_index) %>%
    filter(degree > 0)
  
  
  
  return(g)
}

prep_lr <- function(lr_range = 10^seq(1,3,by= 0.2)/1e7) {
  #  #drug index table
  
  out <- prep_performance_overview(N=4, landscapes = "lr", lr_range = lr_range)
  
  #first prep the performance of all 70 runs for the lollipop
  mems <- out[[4]] %>% 
    group_by(lr, original_ep) %>% 
    summarise(fitness = mean(fitness))
  
  return(mems)
}

prep_update_target <- function() {
  out <- prep_performance_overview(N=4, landscapes = "updatetarget", lr_range = NULL)
  
  #first prep the performance of all 70 runs for the lollipop
  mems <- out[[4]] %>% 
    group_by(update_target, original_ep) %>% 
    summarise(fitness = mean(fitness))
  return(mems)
}

prep_batch_size <- function() {
  out <- prep_performance_overview(N=4, landscapes = "batchsize", lr_range = NULL)
  
  #first prep the performance of all 70 runs for the lollipop
  mems <- out[[4]] %>% 
    group_by(batch_size, original_ep) %>% 
    summarise(fitness = mean(fitness))
  return(mems)
}

prep_gamma <- function() {
  
  out <- prep_performance_overview(N=4, landscapes = "gamma", lr_range = NULL)
  
  df <- out[[1]] %>% filter(!is.na(replicate), ep_number %in% 450:500) %>%
    group_by(replicate, gamma, condition) %>% 
    summarise(fitness = mean(average_fitness)) %>% 
    pivot_wider(names_from = "condition", values_from = "fitness") %>% 
    mutate(benefit = naive - evodm,
           distance_optimal = optimal_policy - evodm, 
           replicate = factor(replicate)) %>% 
    ungroup() %>%
    mutate(replicate = fct_reorder(replicate, benefit, .desc=TRUE)) 
  return(df)
  
}

prep_opt_policy <- function(df) {
  drug_index <- data.frame(action = 1:15, 
                           drug = c("AMP", "AM", "CEC", "CTX", "ZOX", "CXM", 
                                    "CRO", "AMC", "CAZ", "CTT", "SAM", "CPR",
                                    "CPD", "TZP", "FEP"))
  df_opt <- df %>% filter(condition == "optimal_policy") 
  
  #Put together figure showing probability of a drug being selected over the course of an episode (for the optimal)
  opt_freq_table <- (table(df_opt$drug, df_opt$action_number))
  opt_props <- as.data.frame(opt_freq_table/Matrix::colSums(opt_freq_table)) %>%
    rename(action = Var1, time_step = Var2, prob_selection = Freq) %>% 
    mutate(action = as.numeric(as.character(action)),
           time_step = as.numeric(time_step),
           prob_selection = as.numeric(prob_selection)) %>% # have to add an as.character for some reason
    left_join(drug_index)
  
  opt_props_sum <- opt_props %>% 
    group_by(action, drug) %>% 
    summarise(prob_selection = mean(prob_selection)) 
  
  opt_props_sum <- left_join(drug_index, opt_props_sum) %>% 
    mutate(prob_selection = ifelse(is.na(prob_selection), 0, prob_selection))
  
  return(list(opt_props, opt_props_sum))
}

prep_twodrug_sweep <- function(update = FALSE) {
  if(file.exists("data/results/twodrugsweep.Rda") & !update) {
    load("data/results/twodrugsweep.Rda")
    return(out)
  }
  drug_index1 <- data.frame(drug1 = 0:14, 
                            drug_name1 = c("AMP", "AM", "CEC", "CTX", "ZOX", "CXM", 
                                           "CRO", "AMC", "CAZ", "CTT", "SAM", "CPR",
                                           "CPD", "TZP", "FEP"))
  drug_index2 <- data.frame(drug2 = 0:14, 
                            drug_name2 = c("AMP", "AM", "CEC", "CTX", "ZOX", "CXM", 
                                           "CRO", "AMC", "CAZ", "CTT", "SAM", "CPR",
                                           "CPD", "TZP", "FEP"))
  #run function to sweep through all two drug policies and compute performance over 100 episodes
  out <- policy_sweep(episodes = as.integer(100))
  #loop through output list
  out2 <- list()
  for(i in 1:length(out)) {
    #clean mem for every policy
    mem_i <- clean_mem(out[[i]][[1]])
    
    #add which drugs made up the alternating two drug policy
    combo_i <- unlist(out[[i]][[2]])
    mem_i$drug1 <- combo_i[1]
    mem_i$drug2 <- combo_i[2]
    mem_i$starting_state <- out[[i]][[3]]
    
    #done
    out2[[i]] <- mem_i
  }
  
  #get landscape correlations
  df <- compute_ls_cor()
  
  df_long <- bind_rows(out2) %>% 
    left_join(df) %>% 
    left_join(drug_index1) %>% left_join(drug_index2) %>% 
    group_by(drug_name1, drug_name2, ep_number, starting_state) %>% 
    summarise(average_fitness = mean(average_fitness), correl = mean(correl)) %>%
    ungroup() %>%
    mutate(combo =  paste(drug_name1, drug_name2, sep = ",")) 
  
  df_hist <- df_long %>% group_by(drug_name1, drug_name2, combo, starting_state) %>% 
    summarise(average_fitness = mean(average_fitness), 
              correl = mean(correl))
  
  #order fct
  df_long <- df_long %>% 
    mutate(combo = fct_reorder(combo, average_fitness)) %>%
    filter(starting_state == 0)
  
  df_hist_wt <- df_hist %>% 
    filter(starting_state == 0)
  
  out <- list(df_long, df_hist, df_hist_wt)
  save(out, file = "data/results/twodrugsweep.Rda")
  return(out)
  
}


## Function to compute pairwise correlation of each mira landscape
#out: df with pairwise mira correlations

compute_ls_cor <- function(){
  drugs <- define_mira_landscapes() 
  
  #get ready to loop through
  drug1_vec <- vector(mode = "numeric", length = length(drugs)^2)
  drug2_vec <- vector(mode = "numeric", length = length(drugs)^2)
  cor_vec <- drug1_vec
  
  counter = 0
  for(i in 1:length(drugs)) {
    for(j in 1:length(drugs)) {
      counter = counter+1
      cor_vec[counter] <- cor(drugs[[i]], drugs[[j]], method = "pearson")
      drug1_vec[counter] <- i
      drug2_vec[counter] <- j
    }
  }
  df <- data.frame(drug1 = drug1_vec-1, drug2 = drug2_vec-1, correl = cor_vec)
  return(df)
  
  
}

prep_network <- function(joint_prob) {
  drug_index <- data.frame(id = 1:15, 
                           label = c("AMP", "AM", "CEC", "CTX", "ZOX", "CXM", 
                                     "CRO", "AMC", "CAZ", "CTT", "SAM", "CPR",
                                     "CPD", "TZP", "FEP"),
                           label_long = c("Ampicillin", "Amoxicillin", "Cefaclor", "Cefotaxime",
                                          "Ceftizoxime", "Cefuroxime", "Ceftriaxone",
                                          "Amoxicillin + Clavulanic acid", "Ceftazidime",
                                          "Cefotetan", "Ampicillin + Sulbactam", 
                                          "Cefprozil", "Cefpodoxime", 
                                          "Pipercillin + Tazobactam", "Cefepime")
  )
  edges <- joint_prob %>% 
    rename(from = drug1, to = drug2, value = event_counts) %>% filter(value > 2000) %>% 
    mutate(value = value/10)
  nodes <- data.frame(id = unique(c(edges$from, edges$to))) %>% 
    left_join(drug_index)
  return(list(nodes, edges))
}


