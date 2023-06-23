##################################
#Functions for coalition formation
##################################

#Note that for the following functions, D2 simply refers to the fact that these rules were first specified in the second Deliverable for the grant that initially funded this project (see ReadMe)

#Nodes that experience excess demand equal to at least 20% of their current gross demand will consider coalition membership or altering their membership
#In addition, if a node is lacking access to any suppliers, they will automatically seek coalition membership (note, this condition should never actually occur)
D2_add_rule <- function(model_state){
  model_state$sink_nodes %>%
    lapply(get_node_suppliers, model_state = model_state) %>% 
    sapply(length) -> tmp; 
  res1 <- model_state$sink_nodes[tmp==0]
  model_state$node_params %>% 
    filter(timeStep == (model_state$current_step - 1), currentGrossDemand > 0,
           ExcessDemand / currentGrossDemand >= 0.2) %>% pull("nodeID") -> res2
  
  return(c(res1,res2) %>% unique() %>% sort())
}

#Nodes that experience excess demand equal to at least 20% of their current gross demand will reevaluate their current coalition membership
#To limit unrealistic fluidity, nodes must remain part of a given coalition for at least 3 time steps before they will consider leaving/switching
D2_leave_rule <- function(model_state){
  model_state$node_params %>% 
    filter(timeStep %in% (model_state$current_step - 1), currentGrossDemand > 0,
           coalitionTime >= 3,
           ExcessDemand / currentGrossDemand >= 0.2) %>%
    pull("nodeID")
}

#Nodes that are interested in joining a coalition identify a candidate set of sink nodes with which they may want to join
#These sink nodes are: (1) selected from the subset of sink nodes that are connected to the focal node's suppliers
#and (2) interested in forming a coalition themselves
r2_potential_01 <- function(seeker_nodes,a_set, candidate_nodes=seeker_nodes ){
  as.list(seeker_nodes) %>% 
    purrr::lmap(function(xx) 
      c(xx, r2_neighbor(xx[[1]], a_set) %>% 
          purrr::keep(~.x %in% candidate_nodes) %>%
          purrr::map(~ sort(c(xx[[1]],.x)) ))) %>% unique()
}

#Used to identify a node's 1-cofaces (i.e., its direct connections) and the 1-cofaces of its 1-cofaces
#In other words, a node can identify all nodes within two steps from itself in the network
r2_neighbor <- function(nodeID, a_set) {
  cofaceList(a_set,nodeID) %>% purrr::map(~ cofaceList(a_set,.x)) %>% 
    unlist() %>% unique() %>% sort() %>% 
    purrr::keep(~  ! .x %in% nodeID)
}

#For each possible coalition involving a node's identified candidates, the node estimates the strength of that coalition
#Coalition options are ranked according to their strength and then reduced via select_coalitions
simple_sinkCoalition_formation_stage1 <- function(model_state, seeking_nodes){
  model_state %>% 
    with(r2_potential_01(seeking_nodes ,assocSkel)) -> coal_cand
  coal_strength = sapply(coal_cand, 
                         simple_coalition_strength_V2, 
                         model_state=model_state)
  coal_cand[order(coal_strength, decreasing = T)]  %>%
    select_coalitions()
}

#This function produces an estimate of each coalition's ability to procure sufficient inventory to meet its demand
#It takes into account a coalition's total demand, total liquid funds, the amount of inventory held by its suppliers and its suppliers' price points
#In addition, "coalitions" made up of a single node will downgrade their strength estimate according to its baseline reluctance to join a coalition along with recent evidence of its performance when acting alone
#Nodes that repeatedly fail to acquire sufficient supplies to meet their demand will be increasingly likely to downgrade their own ability when acting alone and thus prioritize joining/forming a coalition
simple_coalition_strength_V2 <- function(model_state, coalition){
  if(length(coalition) ==0 )return(0)
  if(length(coalition) == 1 && model_state$current_step > 5) {
    if("reluctance" %in% names(model_state$node_params)){
      nodeDataTemp <- model_state$node_params[model_state$node_params$timeStep < model_state$current_step & model_state$node_params$timeStep > model_state$current_step - 6 & 
                                                model_state$node_params$nodeID == coalition & model_state$node_params$coalitionTime == 0, ]
      if(nrow(nodeDataTemp) >0) {
        estimate_adjust = 1 - mean(nodeDataTemp$ExcessDemand/nodeDataTemp$currentGrossDemand)
        reluctance = 1 - ((1 + exp(-(sum(nodeDataTemp$ExcessDemand/nodeDataTemp$currentGrossDemand >= 0.2))))^-(sus_node(model_state$node_params, timeStep = model_state$current_step, nodeID = coalition)$reluctance))
        estimate_adjust = estimate_adjust * reluctance
        
      } else{
        estimate_adjust = 1
      }
    } else {
      estimate_adjust = 1
    }
  } else{
    estimate_adjust = 1
  }
  #Calculate buy price for coalition
  buy_price = 
    (sum(coalition %>% purrr::map_dbl(sus_node_budget, model_state = model_state)) - 
       sum(simple_coalition_cost(model_state, coalition, coalition)))/
    sum((coalition %>% purrr::map_dbl(sus_node_demand, model_state = model_state )))
  
  total_demand = sum(coalition %>% purrr::map_dbl(sus_node_demand, model_state = model_state))
  
  #Determine potential sellers 
  sellers <- coalition %>% 
    purrr::map(get_node_suppliers, 
               model_state = model_state ) %>% 
    unlist() %>% unique() %>% sort()
  if(length(sellers)==0) return(0)  # Coalition power is 0 if no sellers exist
  sus_node(nP = model_state$node_params,
           timeStep = model_state$current_step,
           nodeID = sellers )%>% group_by(nodeID) %>%
    group_modify(~ trace_potential_supplies(model_state, .y$nodeID) %>% select(-nodeID) ) %>%
    ungroup() %>% arrange(src_ID,pricePerUnit) %>% group_by(src_ID) %>% slice(1) %>%
    ungroup() -> tmp
  tmp %>% mutate(
    sell_power=pmin(availableSupplies, total_demand * (buy_price/pricePerUnit))) %>%
    pull("sell_power") %>% sum()/length(sellers) -> tmp
  tmp <- (tmp / total_demand) * estimate_adjust
}

#Nodes are allocated to coalitions, starting with the top-ranked potential coalition
#A node cannot be allocated to multiple coalitions
#Thus, as coalitions are formed, nodes that have not yet been allocated to a coalition will have to settle for the highest-ranked option of unallocated nodes remaining to them
select_coalitions <- function(candidate_list) {
  res <- list()
  det <- c()
  for(node_set in candidate_list){
    if(!any(node_set %in% det)){
      res <- c(res, list(node_set))
      det <- c(det, node_set) %>% unique()
    }
  }
  res
}

#Nodes that are dissatisfied with their current coalition are split from their current coalition
#Note, that they may opt to rejoin if their current coalition is still estimated as the best option by both parties
simple_sinkCoalition_reduction <- function(model_state, seeking_nodes){
  current_sink_coalitions(model_state = model_state) %>%
    simplex_splitting(seeking_nodes)
}

#Function that splits the dissatisfied coalition member from the rest of the coalition
simplex_splitting <- function(cand_sets, exc_vert) {
  exc_vert <- exc_vert %>% purrr::keep(~ .x %in% as.node_set(cand_sets) )
  res <- 
    cand_sets %>% lapply(function(xx) xx %>% purrr::keep(~ ! .x %in% exc_vert) ) %>%
    purrr::keep(~ length(.x) > 0)
  c( res, as.list(exc_vert)) 
}

#A second round of proposed coalition ranking occurs, based on the coalition-strength algorithm defined above
#This round involves all valid potential coalitions identified by the coalition_expansion function (see below)
#The ranked coalitions are then reduced via select_coalitions (see above).
#The outcome of this process finalizes the coalitions that are used on the current time step
simple_sinkCoalition_formation_stage2 <- function(model_state, old_coal, prop_coal){
  
  coalition_expansion(model_state$assocSkel,
                      c(old_coal, prop_coal),
                      prop_coal ) -> st2_cand
  coal_strength_2 = sapply(st2_cand, 
                           simple_coalition_strength_V2, 
                           model_state=model_state)
  st2_cand[order(coal_strength_2, decreasing = T)]  %>% select_coalitions()
}

#Identifies potential coalitions made up of every (valid) combination of proposed coalitions (from stage 1) and current coalitions
#For current coalitions that are considering splitting, all possible outcomes are considered when identifying potential coalitions
#For example, if node 33 is considering splitting from 44 and 55, then coalitions: "33", "44_55", and "33_44_55" are included as possible candidates
coalition_expansion <- function(a_set, cand_list_A, cand_list_B =list()){
  if(is.list(cand_list_A)) cand_list_A <- as.char.coal_set(cand_list_A)
  if(length(cand_list_B) == 0)  return(cand_list_A %>% unique() %>% as.coal_set.char)
  if(is.list(cand_list_B)) cand_list_B <- as.char.coal_set(cand_list_B)
  tmp <- function(cand_set_A, cand_set_B) {
    cand_set_A = as.coal_set.char(cand_set_A)[[1]]
    cand_set_B = as.coal_set.char(cand_set_B)[[1]]
    if(cand_set_B %>%
       purrr::keep(~any(cand_set_A %in% r2_neighbor(.x,a_set)) ) %>% length() >0 )
      return(as.char.coal_set(c(cand_set_A,cand_set_B)))
    NA_character_
  }
  print(cand_list_A)
  print(cand_list_B)
  exp_sets <- expand.grid(aa=cand_list_A,bb=cand_list_B)
  exp_sets %>% with(mapply(tmp,cand_set_A=aa,cand_set_B=bb )) %>%
    purrr::keep(~!is.na(.x)) %>% c(cand_list_A,cand_list_B) %>% unique %>%
    as.coal_set.char
}