
######################################################
#Functions for setting up and tracking the model state
######################################################

shockSimulation_v1.3 <- 
  function(n_source, 
           n_transfer, 
           n_sink, 
           edge.src_out,
           edge.transfer,
           edge.sink_in,
           disruptionType,
           disruptionMagnitude,
           supplyDisruption,
           meanSourcePrice,
           meanTransferMarkup,
           meanReplenishment,
           meanGD, 
           inventoryType, 
           purchase_Strategy,
           coalitions = "Prohibited",
           coalitionReluctance = "Uniform",
           t_final) {
    
    #Provide unique numeric identifiers for nodes (src = source nodes, mid = transfer nodes, snk = sink nodes)
    srcID = 1:n_source
    midID = (1:n_transfer) + as.integer(n_source)
    snkID = (1:n_sink) + as.integer(n_source+n_transfer)
    
    #Create the initial dyadic supply network
    sc_mat =  create_supply_network(N = n_source + n_transfer + n_sink, 
                                    sourceNodes = srcID, 
                                    transferNodes = midID, 
                                    sinkNodes = snkID, 
                                    sourceEdgeMax = edge.src_out, 
                                    transferToSourceEdgeMax = edge.transfer[1], 
                                    transferToTransferEdgeMax = edge.transfer[2], 
                                    transferToSinkEdgeMax = edge.transfer[3], 
                                    sinkEdgeMax = edge.sink_in)
    
    
    if(supplyDisruption == "supplyReduction") {
      sourceNodes_Post <- srcID
      transferNodes_Post <- midID
      sinkNodes_Post <- snkID
      scMatrix_Post <- sc_mat
    } else{    
      
      #Create the updated list of nodes, supply network, and simplicial set 1-skeleton that will be used post-disruption
      postDisruptionNodeLists <-
        update_node_list_post_disruption(sourceNodeList = srcID, 
                                         transferNodeList = midID, 
                                         sinkNodeList = snkID, 
                                         disruptionType = disruptionType, 
                                         disruptionMagnitude)
      
      sourceNodes_Post <- postDisruptionNodeLists$sourceNodes_Post
      transferNodes_Post <- postDisruptionNodeLists$transferNodes_Post
      sinkNodes_Post <- postDisruptionNodeLists$sinkNodes_Post
      
      
      scMatrix_Post <- 
        create_disrupted_supply_network(adjMatrix = sc_mat, 
                                        sourceNodeList_Pre = srcID, 
                                        sourceNodeList_Post = sourceNodes_Post, 
                                        transferNodeList_Pre = midID, 
                                        transferNodeList_Post = transferNodes_Post, 
                                        sinkNodeList_Pre = snkID, 
                                        sinkNodeList_Post = sinkNodes_Post, 
                                        disruptionType = disruptionType,
                                        N =  n_source + n_transfer + n_sink,
                                        sinkEdgeMax = edge.sink_in)
      
      oneSkelList_Pre <- create_supply_simplicial_list(supplyMatrix = sc_mat)
      oneSkelList_Post <- create_supply_simplicial_list(supplyMatrix = scMatrix_Post[[2]])
      
      oneSkel_Post <- 
        create_disrupted_simplicial_skeleton(simplices = oneSkelList_Post,
                                             adjMatrix = scMatrix_Post[[1]], 
                                             sourceNodeList_Pre = srcID, 
                                             sourceNodeList_Post = sourceNodes_Post, 
                                             transferNodeList_Pre = midID, 
                                             transferNodeList_Post = transferNodes_Post, 
                                             sinkNodeList_Pre = snkID, 
                                             sinkNodeList_Post = sinkNodes_Post, 
                                             disruptionType = disruptionType)
    }
    
    #Set suppliers' (source and transfer nodes) minimum price per unit
    sourcePrice <- 
      round(rnorm(n_source, meanSourcePrice, 0.1 * meanSourcePrice), digits = 2)
    
    transferMarkUp <-
      1 + round(rnorm(n_transfer,  meanTransferMarkup,  0.1 * meanTransferMarkup), digits = 2)
    
    transferPrice <- round(meanSourcePrice * transferMarkUp, digits = 2)
    
    pricePerUnit <-  c(sourcePrice, transferPrice, rep(0, n_sink))
    
    #Set the capacity for sink nodes to replenish their stock
    #This represents manufacturing capability or access to upstream suppliers not included in the simulation
    meanReplenish <- 
      c(round(rnorm(n_source,  meanReplenishment, 0.2 * meanReplenishment)),  
        rep(0, n_transfer + n_sink))
    
    #Generate the initial estimate for each source node of their manufacturing capability on the firs time step
    anticReplenish <- anticipatedReplenishment(adjMatrix = sc_mat, 
                                               meanReplenishmentValues = meanReplenish)
    
    #Generate the starting inventory for suppliers when durable goods are used
    if(inventoryType == "Durable"){
      currentSupply <-  c(sample(1:50, n_source, replace = TRUE), 
                          sample(1:50, n_transfer, replace = TRUE), 
                          rep(0, n_sink))
    } else{
      currentSupply <- 0
    }
    
    #Specify the type of goods that a supplier is selling
    if(inventoryType == "Durable") {
      inventory <- c(rep("Durable", n_source), rep("NA", n_transfer + n_sink))
    }
    if(inventoryType == "Perishable") {
      inventory <- c(rep("Perishable", n_source), rep("NA", n_transfer + n_sink))
    }
    

    
    timeSinceSwitch = rep(0, n_source + n_transfer + n_sink)
    
    #Specify the typical gross demand experienced by sink nodes and the actual gross demand for the first time step
    meanGrossDemand <- c(rep(0, n_source+n_transfer), 
                         round(rnorm(n_sink, meanGD, 0.2 * meanGD)))
    
    currentGrossDemand <- currentGD(adjMatrix = sc_mat, meanGD = meanGrossDemand)
    
    #Generate the initial expected budget needed per step by a sink node (proportional to their typical gross demand and source nodes' mean pricing)
    budgetPerStep <- round(meanGrossDemand * (1.08 * meanSourcePrice), digits = 2)
    
    #Parameter that governs a node's reluctance to join coalitions by limiting the effects of supply shortfalls it experiences on its evaluation of its own performance
    if(coalitionReluctance == "Variable"){
      reluctance = c(rep(0, n_source + n_transfer), 
                     rtruncnorm(n = n_sink, a = 0.01, b = 100, mean = 50, sd = 20))
    }
    
    if(coalitionReluctance == "Uniform"){
      reluctance = c(rep(0, n_source + n_transfer), 
                     rep(50, n_sink))
    }
    
    #Specify the strategy used by sink nodes during purchasing
    if(purchase_Strategy == "pricePriority") {
      purchaseStrategy <- rep("pricePriority", n_source + n_transfer + n_sink)
    }
    if(purchase_Strategy == "random") {
      purchaseStrategy <- c(rep("pricePriority", n_source + n_transfer), rep("random", n_sink))
    }
    if(purchase_Strategy == "betHedging") {
      purchaseStrategy <- c(rep("pricePriority", n_source + n_transfer), rep("betHedging", n_sink))
    }
    if(purchase_Strategy == "mixed") {
      purchaseStrategy <- c(rep("pricePriority", n_source + n_transfer), 
                            sample(c("pricePriority","betHedging"), n_sink, replace = TRUE))
    }
    
    #Create the node-level parameters dataframe
    #Many variables are defined above
    #Vol stands for volume (e.g., OrderVolPlaced = volume of units that a node has ordered)
    #SupplyReserved are units held by a supplier that have been held in reserve to fulfill received orders
    #ExcessDemand refers to units of unmet demand at the end of a round
    #Nodes are assigned a lump sum of liquid funds at the start of the simulation to purchase goods
    #This liquid fund pool depletes over time as nodes purchase inventory
    #UsableLiquidFunds refers to the amount a node sets aside for purchasing on the current time step
    nP <-
      data.frame(
        "Coalitions" = coalitions,
        "nodeID" =  c(srcID,midID,snkID),
        entityID = c(srcID,midID,snkID) %>% as.list() %>% as.char.coal_set(),
        "nodeType" =  c(rep("Source", n_source),
                        rep("Transfer", n_transfer),
                        rep("Sink",n_sink)),
        "timeStep" = 1,
        "currentSupply" = currentSupply,
        "pricePerUnit" = pricePerUnit,
        "meanReplenish" = meanReplenish,
        "potentialSupply" = 0,
        "currentGrossDemand" = currentGrossDemand,
        "budgetPerStep" = budgetPerStep,
        "anticReplenish" = anticReplenish,
        "meanGrossDemand" = meanGrossDemand,
        "OrderVolReceived" = 0,
        "OrderVolPlaced" = 0,
        "DeliveryVolReceived" = 0,
        "DeliveryVolMade" = 0,
        "SupplyReserved" = 0,
        "PaymentMade" = 0,
        "PaymentReceived" = 0,
        "ExcessDemand" = 0,
        "supplyGenerated" = 0,
        "desiredBuffer" = 0,
        "initialLiquidFunds" = budgetPerStep * t_final,
        "liquidFunds" = budgetPerStep * t_final,
        "usableLiquidFunds" = budgetPerStep,
        "inventoryType" = inventory,
        "purchaseStrategy" = purchaseStrategy,
        "timeSinceSwitch" = timeSinceSwitch,
        "failureResponse" = 0,
        "supplyLost" = 0,
        "reluctance" = reluctance,
        stringsAsFactors = F)
    
    #Note, the node-level parameters for t = 0 are needed for some functions (i.e., buffer inventory calculation at t = 1), but are ultimately discarded pre-analysis
    nP0 <- nP
    nP0$timeStep <- 0
    nP <- rbind(nP0, nP)
    row.names(nP) <- with(nP, sprintf("t%03d_n%02d",timeStep,nodeID))
    
    #Create the model state
    model_state = modelState_v1.3(nP, sc_mat, inventoryType = inventoryType)
    shock_params <-list(  )
    shock_params$sourceNodes_Post <- sourceNodes_Post
    shock_params$transferNodes_Post <- transferNodes_Post
    shock_params$sinkNodes_Post <- sinkNodes_Post
    shock_params$disruptionType <- disruptionType
    shock_params$meanGD = meanGD
    if(supplyDisruption == "supplyReduction") {
      shock_params$oneSkel_Post <- model_state$oneSkel
      shock_params$scMatrix_Post <- sc_mat
    } else {
    shock_params$oneSkel_Post <- oneSkel_Post
    shock_params$scMatrix_Post <- scMatrix_Post[[1]]
    }
    shock_params$t_final <- t_final
    return(list(model_state=model_state, shock_params=shock_params))
  }

modelState_v1.3 <- function(nodeParams, scMatrix, inventoryType){
  #Initialize a data frame to store order records
  orderRecord <-
    data.frame(
      "TimeStep" = integer(),
      "Customer" = character(),
      "Supplier" = character(),
      "Volume" = integer(),
      "pricePerUnit" = numeric())
  
  srcID = nodeParams %>% filter(nodeType == "Source") %>% pull(nodeID) %>% unique()
  midID = nodeParams %>% filter(nodeType == "Transfer") %>% pull(nodeID) %>% unique()
  snkID = nodeParams %>% filter(nodeType == "Sink") %>% pull(nodeID) %>% unique()
  
  #Create the initial supply chain simplicial 1-skeleton
  oneSkel_Pre <- scMatrix %>% create_supply_simplicial_list %>% create_supply_simplicial_skeleton
  st_coal <- nodeParams %>% pull(nodeID) %>% unique() %>% as.list() %>% simplex_tree()
  
  networkSize =  nodeParams %>% pull(nodeID) %>% max() 
  
  #Specify matrix of values indicating the costs of working as a coalition for each unique node pair
  c_cost = matrix(1e6,nrow=networkSize,ncol=networkSize)
  diag(c_cost) <- 0
  c_cost_d = matrix(0,nrow=networkSize,ncol=networkSize)
  
  #Initialize the model_state object that tracks nodes' parameters, orders, the structure of the supply network, etc.
  model_state<-list(node_budget = liquid_budget, 
                    source_nodes = srcID, 
                    transfer_nodes = midID, 
                    sink_nodes = snkID,
                    current_step = max(nodeParams$timeStep),
                    oneSkel = oneSkel_Pre, 
                    oneSkel_Initial = simplextree::clone(oneSkel_Pre),
                    assocSkel = simplextree::clone(oneSkel_Pre), 
                    coalSkel = st_coal,
                    scMatrix = scMatrix,
                    node_params = nodeParams, 
                    orderRecord = orderRecord,
                    coalition_cost = c_cost,
                    coalition_min_cost = c_cost,
                    coalition_cost_change = c_cost_d)
  
  #Specify order and distribution functions
  if(inventoryType == "Perishable") {
    model_state$node_demand <- estimated_simple_demand
    model_state$place_orders <- place_orders_with_budgets
  } else{
    model_state$node_demand <- estimated_demand_for_durables
    model_state$place_orders <- place_orders_with_budgets
  }
  model_state$distribute_stock <- simple_distribute_stock
  model_state$distribute_stock_coalition <- fail_distribute_stock_coalition
  model_state$make_delivery  <- make_delivery_with_budgets
  
  #Decision rules for coalition-seeking are specified in the function: sinkCoalitionModelState
  model_state$add_rule <- function(model_state) return(integer(0))
  model_state$leave_rule <- function(model_state) return(integer(0))
  
  return(model_state)
}

update_shock_params <- function(model_state, shock_params, edge.sink_in) {
  
  scMatrix_Post <- 
    create_disrupted_supply_network(adjMatrix = model_state$scMatrix, 
                                    sourceNodeList_Pre = model_state$source_nodes, 
                                    sourceNodeList_Post = shock_params$sourceNodes_Post, 
                                    transferNodeList_Pre = model_state$transfer_nodes, 
                                    transferNodeList_Post = shock_params$transferNodes_Post, 
                                    sinkNodeList_Pre = model_state$sink_nodes, 
                                    sinkNodeList_Post = shock_params$sinkNodes_Post, 
                                    disruptionType = shock_params$disruptionType,
                                    N =  length(model_state$source_nodes) + length(model_state$transfer_nodes) + 
                                      length(model_state$sink_nodes),
                                    sinkEdgeMax = edge.sink_in)
  
  oneSkelList_Post <- create_supply_simplicial_list(supplyMatrix = scMatrix_Post[[2]])
  
  oneSkel_Post <- 
    create_disrupted_simplicial_skeleton(simplices = oneSkelList_Post,
                                         adjMatrix = scMatrix_Post[[1]], 
                                         sourceNodeList_Pre = model_state$source_nodes, 
                                         sourceNodeList_Post = shock_params$sourceNodes_Post, 
                                         transferNodeList_Pre = model_state$transfer_nodes, 
                                         transferNodeList_Post = shock_params$transferNodes_Post, 
                                         sinkNodeList_Pre = model_state$sink_nodes, 
                                         sinkNodeList_Post = shock_params$sinkNodes_Post, 
                                         disruptionType = shock_params$disruptionType)
  shock_params$oneSkel_Post <- oneSkel_Post
  shock_params$scMatrix_Post <- scMatrix_Post[[1]]
  return(shock_params)
}

sinkCoalitionModelState <-
  function(model_state,
           coal_rule.add=D2_add_rule, 
           coal_rule.leave=D2_leave_rule,
           coal_cost.init=0, 
           coal_cost.min=0, 
           coal_cost.delta=0, 
           coal_distribute_rule = "Subsidized", 
           coalitions = "Prohibited") {
    
    mS = model_state
    networkSize = with(mS, max(source_nodes,transfer_nodes,sink_nodes))
    
    #Create the simplex tree (st) for use in determining potential coalition membership
    st_coal = diag(1,networkSize,networkSize) %>% create_supply_simplicial_list %>%
      create_supply_simplicial_skeleton
    
    #Create matrix that keeps track of dyad-specific costs of coalition formation
    c_cost = matrix(nrow = networkSize, ncol = networkSize)
    c_cost[mS$sink_nodes, mS$sink_nodes] = coal_cost.init
    diag(c_cost[mS$sink_nodes, mS$sink_nodes]) <- 0
    c_cost_min = matrix(nrow = networkSize, ncol = networkSize)
    c_cost_min[mS$sink_nodes, mS$sink_nodes] = coal_cost.min
    diag(c_cost_min[mS$sink_nodes, mS$sink_nodes]) <- 0
    c_cost_d = matrix(nrow = networkSize, ncol = networkSize)
    c_cost_d[mS$sink_nodes, mS$sink_nodes] = coal_cost.delta
    diag(c_cost_d[mS$sink_nodes, mS$sink_nodes]) <- 0
    
    model_state<-list(node_budget = mS$node_budget, 
                      node_demand = mS$node_demand,
                      source_nodes = mS$source_nodes,
                      transfer_nodes = mS$transfer_nodes,
                      sink_nodes = mS$sink_nodes,
                      current_step = mS$current_step,
                      oneSkel = simplextree::clone(mS$oneSkel),
                      assocSkel = simplextree::clone(mS$oneSkel),
                      coalSkel = st_coal,
                      scMatrix = mS$scMatrix,
                      node_params = mS$node_params,
                      orderRecord = mS$orderRecord,
                      coalition_cost = c_cost,
                      coalition_min_cost = c_cost_min,
                      coalition_cost_change = c_cost_d,
                      coalition_new_node=c(coal_cost.init, 
                                           coal_cost.min, 
                                           coal_cost.delta))
    
    # Adding on coalition rules, parameters and functions
    model_state$node_params <- 
      ensure_param_def(model_state$node_params, "coalition", 
                       sprintf("n%02d", model_state$node_params$nodeID)) 
    model_state$node_params <- ensure_param(model_state$node_params,"coalitionTime")
    
    if (coalitions == "Allowed") {
      model_state$add_rule <- coal_rule.add
    }
    if (coalitions == "Prohibited"){
      model_state$add_rule <- function(model_state) return(integer(0))
    }
    model_state$leave_rule <- coal_rule.leave
    
    model_state$update_coalition_cost <- simple_update_coalition_cost
    model_state$update_coalition_demand <- simple_update_coalition_demand
    model_state$update_coalition_budget <- update_coalition_liquidFunds
    
    if(coal_distribute_rule == "Subsidized") {
      model_state$distribute_stock_coalition <- distribute_stock_coalition_with_subsidies
    }
    if(coal_distribute_rule == "Unsubsidized") {
      model_state$distribute_stock_coalition <- distribute_stock_coalition_without_subsidies
    }
    
    #Copy over any additional parameters that have not been overridden
    model_state <- c(model_state, mS[!names(mS) %in% names(model_state)])
    return(model_state)
  }

update_modelState_wShock <- function(model_state, shock_params, inventoryType, disruptionMagnitude, supplyDisruption){
  nodeParams  = model_state$node_params
  future_modelState = c(model_state,list())
  orderRecord = model_state$orderRecord
  
  # Do shock modifications
  if("t_final" %in% names(shock_params)){
    remaining_steps = shock_params$t_final - model_state$current_step
  } else {
    remaining_steps = -1
  }
  
  futureNodeParams <- 
    update_node_params_post_disruption(nodeParams = model_state$node_params, 
                                       timeStep = model_state$current_step, 
                                       sourceNodeList_Pre = model_state$source_nodes, 
                                       sourceNodeList_Post = shock_params$sourceNodes_Post, 
                                       transferNodeList_Pre = model_state$transfer_nodes, 
                                       transferNodeList_Post = shock_params$transferNodes_Post, 
                                       sinkNodeList_Pre = model_state$sink_nodes, 
                                       sinkNodeList_Post = shock_params$sinkNodes_Post, 
                                       meanGD = shock_params$meanGD, 
                                       disruptionType = shock_params$disruptionType,
                                       supplyDisruption = supplyDisruption,
                                       disruptionMagnitude = disruptionMagnitude,
                                       budgeted_timesteps = remaining_steps)
  
  futureNodeParams$timeSinceSwitch <- futureNodeParams$timeSinceSwitch + 1 
  
  if(failure_Response == "none") {} else {
    failureResponseVect <- rep(0, nrow(futureNodeParams))
    for(n in 1:nrow(futureNodeParams)) {
      if(futureNodeParams[n,"currentGrossDemand"] > 0) {
        if(futureNodeParams[n, "ExcessDemand"]/futureNodeParams[n, "currentGrossDemand"] > 0.2) {
          failureResponseVect[n] <- 1
        }
      }
    }
    futureNodeParams$failureResponse <- failureResponseVect
  }
  
  #Set anticipated replenishment and gross demand for next time step (as during simulation set up)
  futureNodeParams$anticReplenish <-
    rnorm_multi_uint(futureNodeParams$meanReplenish, 0.1)
  
  futureNodeParams$currentGrossDemand <-
    rnorm_multi_uint(futureNodeParams$meanGrossDemand, 0.1)
  
  if(inventoryType == "Perishable"){
    futureNodeParams$supplyLost <- futureNodeParams$currentSupply
    futureNodeParams$currentSupply <- 0
  }
  
  futureNodeParams$liquidFunds[which(futureNodeParams$nodeType == "Sink")] <- futureNodeParams$liquidFunds[which(futureNodeParams$nodeType == "Sink")] -  
    futureNodeParams$PaymentMade[which(futureNodeParams$nodeType == "Sink")]
  futureNodeParams$usableLiquidFunds <- round(futureNodeParams$liquidFunds / (t_final - futureNodeParams$timeStep + 1), digits = 2)
  
  #All other parameters are set to 0
  futureNodeParams$OrderVolPlaced <- 0
  futureNodeParams$OrderVolReceived <- 0
  futureNodeParams$DeliveryVolReceived <- 0
  futureNodeParams$DeliveryVolMade <- 0
  futureNodeParams$SupplyReserved <- 0
  futureNodeParams$PaymentMade <- 0
  futureNodeParams$PaymentReceived <- 0
  futureNodeParams$ExcessDemand <- 0
  futureNodeParams$supplyGenerated <- 0
  futureNodeParams$desiredBuffer <- 0
  
  #Update timeStep to loop_t + 1
  futureNodeParams$timeStep <- rep(model_state$current_step + 1, dim(futureNodeParams)[1])
  row.names(futureNodeParams) <- with(futureNodeParams, sprintf("t%03d_n%02d",timeStep,nodeID))
  
  #Append node parameters for time loop_t + 1 to node-level data frame
  nodeParams <- rbind(nodeParams, futureNodeParams)
  
  cat("Pre-Shock Node-level Parameters:\n")
  print(nodeParams[nodeParams$timeStep == t_shock, ])
  
  cat("Shock adjusted Node-level Parameters:\n")
  print(nodeParams[nodeParams$timeStep == t_shock+1, ])
  
  #Now update future model_state
  future_modelState$source_nodes <- shock_params$sourceNodes_Post
  future_modelState$transfer_nodes <- shock_params$transferNodes_Post
  future_modelState$sink_nodes <- shock_params$sinkNodes_Post
  future_modelState$current_step = model_state$current_step+1
  future_modelState$node_params <- nodeParams
  
  #Update simplicial sets
  {
    future_modelState$oneSkel=shock_params$oneSkel_Post
    future_modelState$coalSkel=clone(future_modelState$coalSkel)
    future_modelState$assocSkel=clone(future_modelState$assocSkel)
    if (shock_params$disruptionType == "Supply") {
      nodeRemovals <- setdiff(model_state$source_nodes, future_modelState$source_nodes)
      remove(future_modelState$coalSkel, nodeRemovals %>% as.list())
    } else if (shock_params$disruptionType == "Transport") {
      nodeRemovals <- setdiff(model_state$transfer_nodes, future_modelState$transfer_nodes)
      remove(future_modelState$coalSkel, nodeRemovals %>% as.list())
    } else if (shock_params$disruptionType == "Demand") {
      sinkAdditions <- setdiff(future_modelState$sink_nodes, model_state$sink_nodes)
      insert(future_modelState$coalSkel, sinkAdditions %>% as.list())
      insert(future_modelState$assocSkel, maximal(shock_params$oneSkel_Post) %>% as.list())
      cat("New coalition cost update code for new nodes inserted")
      networkSize_old = with(model_state, max(source_nodes,transfer_nodes,sink_nodes))
      networkSize = with(future_modelState, max(source_nodes,transfer_nodes,sink_nodes))
      if("coalition_new_node" %in% names(model_state)){
        coal_cost.init=model_state$coalition_new_node[1]
        coal_cost.min=model_state$coalition_new_node[2]
        coal_cost.delta=model_state$coalition_new_node[3]
      } else {
        coal_cost.init=1e6
        coal_cost.min=1e6
        coal_cost.delta=0
      }
      
      #Create matrix that keeps track of dyad-specific costs of coalition formation
      c_cost = matrix(nrow = networkSize, ncol = networkSize)
      if("coalition_cost" %in% names(model_state))
        c_cost[1:networkSize_old,1:networkSize_old] = model_state$coalition_cost
      c_cost[sinkAdditions, future_modelState$sink_nodes] = coal_cost.init
      c_cost[future_modelState$sink_nodes, sinkAdditions] = coal_cost.init
      diag(c_cost[sinkAdditions, sinkAdditions])<-0
      
      c_cost_min = matrix(nrow=networkSize, ncol=networkSize)
      if("coalition_min_cost" %in% names(model_state))
        c_cost_min[1:networkSize_old, 1:networkSize_old] = model_state$coalition_min_cost
      c_cost_min[sinkAdditions, future_modelState$sink_nodes] = coal_cost.min
      c_cost_min[future_modelState$sink_nodes, sinkAdditions] = coal_cost.min
      diag(c_cost_min[sinkAdditions, sinkAdditions])<-0
      
      c_cost_d = matrix(nrow=networkSize, ncol=networkSize)
      if("coalition_cost_change" %in% names(model_state))
        c_cost_d[1:networkSize_old, 1:networkSize_old] = model_state$coalition_cost_change
      c_cost_d[sinkAdditions, future_modelState$sink_nodes] = coal_cost.delta
      c_cost_d[future_modelState$sink_nodes, sinkAdditions] = coal_cost.delta
      diag(c_cost_d[sinkAdditions, sinkAdditions])<-0
      
      future_modelState$coalition_cost = c_cost
      future_modelState$coalition_min_cost = c_cost_min
      future_modelState$coalition_cost_change = c_cost_d
    }
  }
  future_modelState$scMatrix=shock_params$scMatrix_Post
  return(future_modelState)
}

update_node_params_post_disruption <- function(nodeParams,
                                               timeStep,
                                               sourceNodeList_Pre,
                                               sourceNodeList_Post,
                                               transferNodeList_Pre,
                                               transferNodeList_Post,
                                               sinkNodeList_Pre,
                                               sinkNodeList_Post,
                                               meanGD,
                                               disruptionType = c("Supply",
                                                                  "Transport",
                                                                  "Demand"),
                                               supplyDisruption,
                                               disruptionMagnitude,
                                               budgeted_timesteps=-1) {
  #Copy over current node-level data
  nodeParamsTemp <-
    nodeParams[which(nodeParams$timeStep == timeStep),]
  
  if(supplyDisruption == "supplyReduction" & disruptionType == "Supply") {
  nodeParamsTemp$meanReplenish <- round(nodeParamsTemp$meanReplenish * (1 - disruptionMagnitude))
  }
  
  else{
  #Add or remove nodes depending on the disruption type
  if (disruptionType == "Supply") {
    sourceRemovals <- setdiff(sourceNodeList_Pre, sourceNodeList_Post)
    nodeParamsTemp <-
      nodeParamsTemp[!nodeParamsTemp$nodeID %in% sourceRemovals,]
  } else
    (
      if (disruptionType == "Transport") {
        transferRemovals <-
          setdiff(transferNodeList_Pre, transferNodeList_Post)
        nodeParamsTemp <-
          nodeParamsTemp[!nodeParamsTemp$nodeID %in% transferRemovals,]
      } else
        (
          if (disruptionType == "Demand") {
            if(budgeted_timesteps < 0  && "initialLiquidFunds" %in% names(nodeParamsTemp))
              stop("For demand disruption, budgeted_timesteps must be set to a value >= 0.")
            sinkAdditions <- setdiff(sinkNodeList_Post, sinkNodeList_Pre)
            mGD <-
              round(rnorm(length(sinkAdditions), meanGD, 0.2 * meanGD))
            budget <-
              round(mGD * (1.08 * mean(nodeParamsTemp$pricePerUnit[nodeParamsTemp$nodeType == "Source"])), digits = 2)
            #Specify the strategy used by sink nodes during purchasing
            if(purchase_Strategy == "pricePriority") {
              purchaseStrategy_add <- rep("pricePriority", length(sinkAdditions))
            }
            if(purchase_Strategy == "random") {
              purchaseStrategy_add <- rep("random", length(sinkAdditions))
            }
            if(purchase_Strategy == "betHedging") {
              purchaseStrategy_add <- rep("betHedging", length(sinkAdditions))
            }
            if(purchase_Strategy == "mixed") {
              purchaseStrategy_add <- sample(c("pricePriority","betHedging"), length(sinkAdditions), replace = TRUE)
            }
            str(nodeParamsTemp)
            nodeParamAdditions <-
              data.frame(
                "Coalitions" = unique(nodeParamsTemp$Coalitions)[1],
                "nodeID" = sinkAdditions,
                "entityID" = sinkAdditions %>% as.list() %>% as.char.coal_set(),
                "nodeType" = rep("Sink", length(sinkAdditions)),
                "timeStep" = rep(timeStep, length(sinkAdditions)),
                "currentSupply" = 0,
                "pricePerUnit" = 0.00,
                "meanReplenish" = 0,
                "potentialSupply" = 0,
                "currentGrossDemand" = 0,
                "budgetPerStep" = budget,
                "anticReplenish" = 0,
                "meanGrossDemand" = mGD,
                "OrderVolReceived" = 0,
                "OrderVolPlaced" = 0,
                "DeliveryVolReceived" = 0,
                "DeliveryVolMade" = 0,
                "SupplyReserved" = 0,
                "PaymentMade" = 0,
                "PaymentReceived" = 0,
                "ExcessDemand" = 0,
                "supplyGenerated" = 0,
                "desiredBuffer" = 0)
            if("initialLiquidFunds" %in% names(nodeParamsTemp))
              nodeParamAdditions$initialLiquidFunds = budget * (budgeted_timesteps)
            if("liquidFunds" %in% names(nodeParamsTemp))
              nodeParamAdditions$liquidFunds = budget * (budgeted_timesteps)
            if("usableLiquidFunds" %in% names(nodeParamsTemp))
              nodeParamAdditions$usableLiquidFunds = budget
            if("inventoryType" %in% names(nodeParamsTemp))
              nodeParamAdditions$inventoryType = "NA"
            if("purchaseStrategy" %in% names(nodeParamsTemp))
              nodeParamAdditions$purchaseStrategy = purchaseStrategy_add
            if("timeSinceSwitch" %in% names(nodeParamsTemp))
              nodeParamAdditions$timeSinceSwitch = 0
            if("failureResponse" %in% names(nodeParamsTemp))
              nodeParamAdditions$failureResponse = 0
            if("supplyLost" %in% names(nodeParamsTemp))
              nodeParamAdditions$supplyLost = 0
            if("reluctance" %in% names(nodeParamsTemp))
              nodeParamAdditions$reluctance = 0
            if("coalition"  %in% names(nodeParamsTemp))
              nodeParamAdditions$coalition = sprintf("n%02d", sinkAdditions)
            if("coalitionTime"  %in% names(nodeParamsTemp))
              nodeParamAdditions$coalitionTime = 0
            #if("Coalitions" %in% names(nodeParamsTemp))
            #  nodeParamAdditions$Coalitions = unique(nodeParamsTemp$Coalitions)[1]
            
            str(nodeParamAdditions)
            nodeParamsTemp <- rbind(nodeParamsTemp, nodeParamAdditions)
            str(nodeParamsTemp)
          }
        )
    )
  }
  return(nodeParamsTemp)
}

###########################################################
#Functions for creating supply networks and simplicial sets
###########################################################

create_supply_network <-
  function(N,
           sourceNodes,
           transferNodes,
           sinkNodes,
           sourceEdgeMax,
           transferToSourceEdgeMax,
           transferToTransferEdgeMax,
           transferToSinkEdgeMax,
           sinkEdgeMax) {
    
    #Set dimensions (number of rows and columns) of an empty adjacency matrix
    scMatrix <- matrix(0, nrow = N, ncol = N)
    
    #Generate edges between nodes, representing pre-exisiting supplier-customer relationships
    #Each source nodes produces >= 1 outgoing edges to randomly selected transfer and sink nodes
    for (i in min(sourceNodes):max(sourceNodes)) {
      nEdges <- sample(1:sourceEdgeMax, 1)
      customerNodes <- c(transferNodes, sinkNodes)
      customers <- sample(customerNodes, nEdges, replace = FALSE)
      for (c in customers) {
        scMatrix[i, c] <- 1
      }
    }
    
    #Each transfer node produces >= 1 incoming edges from sources, >= 1 outgoing edges to sinks,
    #and >= 0 horizontal edges to fellow transfer nodes
    for (i in min(transferNodes):max(transferNodes)) {
      nEdgesSource <- sample(1:transferToSourceEdgeMax, 1)
      nEdgesSink <- sample(1:transferToSinkEdgeMax, 1)
      nEdgesTransfer <- sample(0:transferToTransferEdgeMax, 1)
      sources <- sample(sourceNodes, nEdgesSource, replace = FALSE)
      sinks <- sample(sinkNodes, nEdgesSink, replace = FALSE)
      transfer <-
        sample(transferNodes, nEdgesTransfer, replace = FALSE)
      for (s in sources) {
        scMatrix[s, i] <- 1
      }
      for (c in sinks) {
        scMatrix[i, c] <- 1
      }
      for (t in transfer) {
        scMatrix[i, t] <- 1
      }
    }
    
    #Each sink node generates >= 1 edges to suppliers (either source or transfer nodes)
    for (i in min(sinkNodes):max(sinkNodes)) {
      nEdges <- sample(1:sinkEdgeMax, 1)
      supplierNodes <- c(transferNodes, sourceNodes)
      suppliers <- sample(supplierNodes, nEdges, replace = FALSE)
      for (s in suppliers) {
        scMatrix[s, i] <- 1
      }
    }
    return(scMatrix)
  }

#From the populated adjacency matrix, generate the list of 1-simplices comprising the pre-exisiting supply network
create_supply_simplicial_list <- function(supplyMatrix) {
  oneSkelList <- list()
  for (r in 1:nrow(supplyMatrix)) {
    for (c in 1:nrow(supplyMatrix)) {
      if (supplyMatrix[r, c] >= 1) {
        listIndex <- length(oneSkelList) + 1
        oneSkelList[[listIndex]] <- c(r, c)
      }
    }
  }
  return(oneSkelList)
}

#From the list of 1-simplices, create a 1-dimensional simplicial complex (or 1-skeleton)
#Note, that this is equivalent to a traditional network, but allows for dynamic formation and dissolution of n-dimensional coalitions
create_supply_simplicial_skeleton <- function(simplicialList) {
  oneSkel <- simplex_tree()
  oneSkel %>% insert(simplicialList)
  return(oneSkel)
}

#Specify the list of nodes that will be included in the simulation post-disruption
update_node_list_post_disruption <- function(sourceNodeList,
                                             transferNodeList,
                                             sinkNodeList,
                                             disruptionType = c("Supply",
                                                                "Demand",
                                                                "Transport"),
                                             disruptionMagnitude) {
  if (disruptionType == "Supply") {
    sourceRemovals <- sample(sourceNodeList,
                             round(disruptionMagnitude * length(sourceNodeList)),
                             replace = FALSE)
    return(
      list(
        sourceNodes_Post = sourceNodeList[!sourceNodeList %in% sourceRemovals],
        transferNodes_Post = transferNodeList,
        sinkNodes_Post = sinkNodeList
      )
    )
  } else
    (
      if (disruptionType == "Transport") {
        transferRemovals <- sample(transferNodeList,
                                   round(disruptionMagnitude * length(transferNodeList)),
                                   replace = FALSE)
        return(
          list(
            sourceNodes_Post = sourceNodeList,
            transferNodes_Post = transferNodeList[!transferNodeList %in% transferRemovals],
            sinkNodes_Post = sinkNodeList
          )
        )
      } else
        (
          if (disruptionType == "Demand") {
            demandAdditions <- c(sinkNodeList,
                                 (max(sinkNodeList) + 1):(max(sinkNodeList) +
                                                            round(
                                                              disruptionMagnitude * length(sinkNodeList)
                                                            )))
            return(
              list(
                sourceNodes_Post = sourceNodeList,
                transferNodes_Post = transferNodeList,
                sinkNodes_Post = demandAdditions
              )
            )
          }
        )
    )
}

#Generated the post-disruption supply network
create_disrupted_supply_network <- function(adjMatrix,
                                            sourceNodeList_Pre,
                                            sourceNodeList_Post,
                                            transferNodeList_Pre,
                                            transferNodeList_Post,
                                            sinkNodeList_Pre,
                                            sinkNodeList_Post,
                                            disruptionType = c("Supply",
                                                               "Demand",
                                                               "Transport"),
                                            N,
                                            sinkEdgeMax) {
  
  #For supply and transport disruptions, if a sink node is fully disconnected from suppliers following disruption,
  #it gets to form a "free" connection with a new supplier
  if (disruptionType == "Supply") {
    sourceRemovals <- setdiff(sourceNodeList_Pre, sourceNodeList_Post)
    if(length(sourceNodeList_Pre) == length(sourceNodeList_Post)){
      scMatrixTemp <- adjMatrix
    } 
    else{
    scMatrixTemp <- adjMatrix[-sourceRemovals,-sourceRemovals]
    }
    supplierNodes <-
      c(transferNodeList_Post, sourceNodeList_Post)
    for(i in (length(supplierNodes)+1):(length(supplierNodes)+length(sinkNodeList_Post))) {
      if(sum(scMatrixTemp[,i]) == 0){
        newSupplier <- sample(supplierNodes, 1)
        newSupplierAdj <- sum(newSupplier > sourceRemovals)
        scMatrixTemp[(newSupplier-newSupplierAdj),i] <- 1
        adjMatrix[newSupplier, (i+length(sourceRemovals))]
      }
    }
    scMatrixTempIncomingSourceLinks <- scMatrixTemp[1:length(sourceNodeList_Post), ]
    for(i in (length(sourceNodeList_Post) +1):(length(sourceNodeList_Post)+length(transferNodeList_Post))) {
      if(sum(scMatrixTempIncomingSourceLinks[,i]) == 0) {
        newSource <- sample(sourceNodeList_Post, 1)
        newSourceAdj <- sum(newSource > sourceRemovals)
        scMatrixTemp[(newSource-newSourceAdj), i] <- 1
        adjMatrix[newSource, (i+length(sourceRemovals))]
      }
    }
    return(list(scMatrixTemp, adjMatrix))
  } else
    (
      if (disruptionType == "Transport") {
        transferRemovals <-
          setdiff(transferNodeList_Pre, transferNodeList_Post)
        scMatrixTemp <-
          adjMatrix[-transferRemovals,-transferRemovals]
        supplierNodes <-
          c(transferNodeList_Post, sourceNodeList_Post)
        for(i in (length(supplierNodes)+1):(length(supplierNodes)+length(sinkNodeList_Post))) {
          if(sum(scMatrixTemp[,i]) == 0){
            newSupplier <- sample(supplierNodes, 1)
            newSupplierAdj <- sum(newSupplier > transferRemovals)
            scMatrixTemp[(newSupplier-newSupplierAdj),i] <- 1
            adjMatrix[newSupplier, (i+length(transferRemovals))]
          }
        }
        customerNodes <-
          c(transferNodeList_Post, sinkNodeList_Post)
        for(i in 1:length(sourceNodeList_Post)) {
          if(sum(scMatrixTemp[i,]) == 0) {
            newCustomer <- sample(customerNodes, 1)
            newCustomerAdj <- sum(newCustomer > transferRemovals)
            scMatrixTemp[i, newCustomer-newCustomerAdj] <- 1
            adjMatrix[i, newCustomer] <- 1
          }
        }
        return(list(scMatrixTemp,adjMatrix))
      } else
        (
          if (disruptionType == "Demand") {
            sinkAdditions <- setdiff(sinkNodeList_Post, sinkNodeList_Pre)
            scMatrixTemp <-
              matrix(0, nrow = (N + length(sinkAdditions)), ncol = (N + length(sinkAdditions)))
            scMatrixTemp[1:nrow(adjMatrix), 1:ncol(adjMatrix)] <-
              adjMatrix
            #Each new sink node generates >= 1 edges to suppliers (either source or transfer nodes)
            for (i in min(sinkAdditions):max(sinkAdditions)) {
              nEdges <- sample(1:sinkEdgeMax, 1)
              supplierNodes <-
                c(transferNodeList_Post, sourceNodeList_Post)
              suppliers <-
                sample(supplierNodes, nEdges, replace = FALSE)
              
              for (s in suppliers) {
                scMatrixTemp[s, i] <- 1
              }
            }
            return(list(scMatrixTemp,scMatrixTemp))
          }
        )
    )
}

#Create the simplicial skeleton to be used in the simulation post-disruption
create_disrupted_simplicial_skeleton <- function(simplices,
                                                 adjMatrix,
                                                 sourceNodeList_Pre,
                                                 sourceNodeList_Post,
                                                 transferNodeList_Pre,
                                                 transferNodeList_Post,
                                                 sinkNodeList_Pre,
                                                 sinkNodeList_Post,
                                                 disruptionType = c("Supply",
                                                                    "Demand",
                                                                    "Transport")) {
  if (disruptionType == "Supply") {
    sourceRemovals <- setdiff(sourceNodeList_Pre, sourceNodeList_Post)
    sourceRemovalsList <- sourceRemovals %>% as.list()
    stTemp <- simplex_tree()
    stTemp %>% insert(simplices)
    stTemp %>% remove(sourceRemovalsList)
    return(stTemp)
  } else
    (if (disruptionType == "Transport") {
      transferRemovals <-
        setdiff(transferNodeList_Pre, transferNodeList_Post)
      transferRemovalsList <- transferRemovals %>% as.list()
      stTemp <- simplex_tree()
      stTemp %>% insert(simplices)
      stTemp %>% remove(transferRemovalsList)
      return(stTemp)
    } else
      (
        if (disruptionType == "Demand") {
          sinkAdditions <- setdiff(sinkNodeList_Post, sinkNodeList_Pre)
          additionsList <- list()
          for (r in 1:nrow(adjMatrix)) {
            for (c in sinkAdditions) {
              if (adjMatrix[r, c] >= 1) {
                listIndex <- length(additionsList) + 1
                additionsList[[listIndex]] <- c(r, c)
              }
            }
          }
          stTemp <- simplex_tree()
          stTemp %>% insert(simplices)
          stTemp %>% insert(additionsList)
          return(stTemp)
        })
    )
}

anticipatedReplenishment <-
  function(adjMatrix, meanReplenishmentValues) {
    
    anticReplenish <- rep(0, nrow(adjMatrix))
    for (i in 1:nrow(adjMatrix)) {
      if (meanReplenishmentValues[i] > 0) {
        anticReplenish[i] <-
          max(0, round(rnorm(
            1,
            meanReplenishmentValues[i],
            (0.1 * meanReplenishmentValues[i])
          )))
      } else{
        anticReplenish[i] <- 0
      }
    }
    return(anticReplenish)
  }

currentGD <- function(adjMatrix, meanGD) {
  
  currentGrossDemand <- rep(0, nrow(adjMatrix))
  for (i in 1:nrow(adjMatrix)) {
    if (meanGD[i] > 0) {
      currentGrossDemand[i] <-
        max(0, round(rnorm(1, meanGD[i], (0.1 * meanGD[i]))))
    } else{
      currentGrossDemand[i] <- 0
    }
  }
  return(currentGrossDemand)
}

##################################################
#Functions for determining node demand and budgets
##################################################

simple_budget <- function(model_state, nodeID, ...) {
  if(is.numeric(nodeID)) nodeID =as.list(nodeID)
  with(sus_coalition(model_state, entityID = nodeID),
       budgetPerStep-PaymentMade)
}

liquid_budget <- function(model_state, nodeID, ...) {
  if(is.numeric(nodeID)) nodeID =as.list(nodeID)
  with(sus_coalition(model_state, entityID = nodeID),
       usableLiquidFunds-PaymentMade)
}

#Demand for perishable goods is simply a node's current gross demand minus any inventory it already holds
simple_demand <- function(model_state, nodeID, ...) {
  if(is.numeric(nodeID)) nodeID =as.list(nodeID)
  with(sus_coalition(model_state, entityID = nodeID),
       currentGrossDemand-currentSupply)
}

#Estimated demand for perishable goods is simply a node's estimated current gross demand minus any inventory it already holds
estimated_simple_demand <- function(model_state, nodeID, ...) {
  if(is.numeric(nodeID)) nodeID =as.list(nodeID)
  with(sus_coalition(model_state, entityID = nodeID),
       max(0, round(rnorm(1, currentGrossDemand, 0.1*currentGrossDemand), digits = 0) - currentSupply))
}

#Demand for durable goods allows for the possibility of stockpiling extra units
#This is done by inflating a node's actual demand by 1.1x
demand_for_durables <- function(model_state, nodeID, ...) {
  if(is.numeric(nodeID)) nodeID =as.list(nodeID)
  with(sus_node(model_state$node_params, timeStep = model_state$current_step, 
                nodeID = as.node_set(nodeID)), sum(round((currentGrossDemand - currentSupply) * 1)))
}

#Demand for durable goods allows for the possibility of stockpiling extra units
#This is done by inflating a node's actual demand by 1.1x
estimated_demand_for_durables <- function(model_state, nodeID, ...) {
  if(is.numeric(nodeID)) nodeID =as.list(nodeID)
  with(sus_node(model_state$node_params, timeStep = model_state$current_step, 
                nodeID = as.node_set(nodeID)), 
       sum(round((round(rnorm(1, currentGrossDemand, 0.1*currentGrossDemand), digits = 0) - currentSupply) * 1)))
}

##############################################
#Functions for placing and distributing orders
##############################################

get_node_suppliers <- function(model_state, nodeID, ...){
  with(model_state,
       cofaceList(st = oneSkel, nodeID = nodeID) %>%
         purrr::keep(~.x %in% c(transfer_nodes, source_nodes)) %>%
         purrr::keep(~ ! .x %in% nodeID))
  
}

#Place orders using liquid funds
place_orders_with_budgets <- function(model_state, entityID, total_volume, supplierParams){
  
  cat("Stub for placing orders for:",entityID,"\n" )
  cat("n SRC:", nrow(supplierParams), "demand:", total_volume, "\n")
  
  if(is.na(total_volume)){
    if(!is.null(model_state$coalition_params))
      print(with(model_state, coalition_params[coalition_params$timeStep==current_step, ] ))
    print(with(model_state, node_params[node_params$timeStep==current_step, ]))
    model_state$coalSkel %>% maximal() %>% as.list() %>% as.char.coal_set() %>% print()
    stop("Cannot place order with total_volume == NA")
  }
  
  #Tracks how many funds a coalition/node has already set aside towards purchases
  fundsAllocated <- 0
  
  if(sus_coalition(model_state, entityID)$purchaseStrategy == "betHedging") {
    
    orderVolPlaced_PrevRound <- 0
    #Work through potential suppliers so long as: 
    #(i) orders placed have not yet met current demand, 
    #(ii) there are suppliers remaining, 
    #(iii) there are funds remaining in usableLiquidFunds
    
    while (sus_coalition(model_state, entityID)$OrderVolPlaced < total_volume &&
           nrow(supplierParams) > 0 && 
           (sus_coalition(model_state, entityID)$usableLiquidFunds - fundsAllocated) > min(supplierParams[,"pricePerUnit"])) {
      
      src_node = supplierParams[1, "entityID"]
      cat("src_nod:", src_node, "ordered:", sus_coalition(model_state, entityID)$OrderVolPlaced, "desired", total_volume , "\n")
      
      #Check if the next supplier with the lowest price per unit has supplies to purchase (either currently held or able to be ordered from a secondary supplier)
      if (with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved + potentialSupply) > 0) {
        cat("A")
        
        #Target order volume is equal either to remaining demand or the rest of the supplier's stock, whichever is lower
        targetOrderVol <-
          min(ceiling((total_volume - sus_coalition(model_state, entityID)$OrderVolPlaced)/2),
              with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved + potentialSupply))
        if (targetOrderVol * sus_coalition(model_state, entityID = src_node)$pricePerUnit <= (sus_coalition(model_state, entityID)$usableLiquidFunds - fundsAllocated)) {
          orderVol <- targetOrderVol
        } else{
          orderVol <- targetOrderVol - 
            ceiling(((targetOrderVol * sus_coalition(model_state, entityID = src_node)$pricePerUnit) -
                       (sus_coalition(model_state, entityID)$usableLiquidFunds - fundsAllocated)) /
                      sus_coalition(model_state, entityID = src_node)$pricePerUnit)
        }
        
        #Record order details in the node-level parameters
        model_state <- 
          sus_add_order(model_state, entityID, src_node, orderVol,
                        sus_coalition(model_state, entityID = src_node)$pricePerUnit)
        
        #Calculate funds a node/coalition has allocated towards orders already
        fundsAllocated <- fundsAllocated + (orderVol * sus_coalition(model_state, entityID = src_node)$pricePerUnit)
        
        #If the order volume exceeds a supplier's current and anticipated stock, the remainder must be drawn from its secondary suppliers
        #In practice, this only applies to transfer nodes
        if (orderVol > with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved)) {
          
          cat("B")
          
          #Determine the amount a supplier needs to still obtain to fulfill an order
          supplierDemand <-
            orderVol - with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved) 
          
          #Update SupplyReserved for supplier
          sus_coalition(model_state, entityID=src_node)$SupplyReserved <-
            orderVol - supplierDemand + sus_coalition(model_state, entityID = src_node)$SupplyReserved
          
          #Identify the suppliers of the original supplier
          src_suppliers <- locate_suppliers(model_state, src_node)
          src_suppliers <-
            src_suppliers[!src_suppliers %in% c(src_node, entityID)]
          
          if (length(src_suppliers) > 0) {
            #Update secondary supplier information and order them according to price
            src_supplierParams <- sus_coalition(model_state, src_suppliers)
            src_supplierParams <- src_supplierParams[order(src_supplierParams$pricePerUnit), ]
            
            #Work through secondary suppliers so long as orders placed have not yet met supplierDemand and secondary suppliers remain
            while (supplierDemand > 0 &&
                   nrow(src_supplierParams) > 0) {
              #Check if the secondary supplier with the lowest price has supplies to be purchased
              if (with(src_supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved) > 0) {
                
                altsrc_node <- src_supplierParams[1, "entityID"]
                
                #If supplies are available, place order equal either to remaining demand or the available suppliers, whichever is lower
                supplierOrderVol <-
                  min(supplierDemand, with(src_supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved))
                
                #Update supplier demand
                supplierDemand <- supplierDemand - supplierOrderVol
                model_state <- 
                  sus_add_order(model_state, src_node, altsrc_node, supplierOrderVol,
                                sus_coalition(model_state, entityID=altsrc_node)$pricePerUnit)
                
                #Update supply reserved for secondary supplier
                sus_coalition(model_state, entityID=altsrc_node)$SupplyReserved <-
                  sus_coalition(model_state, entityID=altsrc_node)$SupplyReserved 
              }
              #Remove current secondary supplier from future consideration for the current primary supplier
              src_supplierParams <- src_supplierParams[-1,]
            }
          }
        } else{
          #Update SupplyReserved for supplier
          sus_coalition(model_state, entityID=src_node)$SupplyReserved <-
            orderVol + sus_coalition(model_state, entityID=src_node)$SupplyReserved 
        }
      }
      supplierParams <- supplierParams[-1,]
      cat("Z\n")
    }
    
    while(sus_coalition(model_state, entityID)$OrderVolPlaced < total_volume &&
          sus_coalition(model_state, entityID)$OrderVolPlaced > orderVolPlaced_PrevRound) {
      
      orderVolPlaced_PrevRound <- sus_coalition(model_state, entityID)$OrderVolPlaced
      
      #Determine possible suppliers, drawn from both source and transfer nodes
      suppliers <- locate_suppliers(model_state, entityID)
      #cat("Stub for placing orders for:",i,"\n" )
      print(suppliers)
      
      #Update transfer nodes' supplier information (e.g. current stock levels)
      sus_node(model_state$node_params,timeStep = model_state$current_step )$potentialSupply <-
        update_potential_supplies(adjMatrix = model_state$scMatrix,
                                  transferNodeList = model_state$transfer_nodes,
                                  sourceNodeList = model_state$source_nodes,
                                  st = model_state$oneSkel,
                                  timeStep = model_state$current_step,
                                  nodeParams = model_state$node_params)
      
      #Create list of potential suppliers and order them according to price per unit
      supplierParams <- sus_coalition(model_state, suppliers)
      supplierParams <- supplierParams[order(supplierParams$pricePerUnit),]
      
      #After ordering suppliers, if "distributorPriority" is active, place distributors first
      if(sus_coalition(model_state, entityID)$failureResponse == 1 && 
         failure_Response == "distributorPriority") {
        supplierParams <- supplierParams[order(supplierParams$nodeType, decreasing = TRUE),]
      }
      
      while (sus_coalition(model_state, entityID)$OrderVolPlaced < total_volume &&
             nrow(supplierParams) > 0 && 
             (sus_coalition(model_state, entityID)$usableLiquidFunds - fundsAllocated) > min(supplierParams[,"pricePerUnit"])) {
        
        src_node = supplierParams[1, "entityID"]
        cat("src_nod:", src_node, "ordered:", sus_coalition(model_state, entityID)$OrderVolPlaced, "desired", total_volume , "\n")
        
        #Check if the next supplier with the lowest price per unit has supplies to purchase (either currently held or able to be ordered from a secondary supplier)
        if (with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved + potentialSupply) > 0) {
          cat("A")
          
          #Target order volume is equal either to remaining demand or the rest of the supplier's stock, whichever is lower
          targetOrderVol <-
            min(ceiling((total_volume - sus_coalition(model_state, entityID)$OrderVolPlaced)/2),
                with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved + potentialSupply))
          if (targetOrderVol * sus_coalition(model_state, entityID = src_node)$pricePerUnit <= (sus_coalition(model_state, entityID)$usableLiquidFunds - fundsAllocated)) {
            orderVol <- targetOrderVol
          } else{
            orderVol <- targetOrderVol - 
              ceiling(((targetOrderVol * sus_coalition(model_state, entityID = src_node)$pricePerUnit) -
                         (sus_coalition(model_state, entityID)$usableLiquidFunds - fundsAllocated)) /
                        sus_coalition(model_state, entityID = src_node)$pricePerUnit)
          }
          
          #Record order details in the node-level parameters
          model_state <- 
            sus_add_order(model_state, entityID, src_node, orderVol,
                          sus_coalition(model_state, entityID = src_node)$pricePerUnit)
          
          #Calculate funds a node/coalition has allocated towards orders already
          fundsAllocated <- fundsAllocated + (orderVol * sus_coalition(model_state, entityID = src_node)$pricePerUnit)
          
          #If the order volume exceeds a supplier's current and anticipated stock, the remainder must be drawn from its secondary suppliers
          #In practice, this only applies to transfer nodes
          if (orderVol > with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved)) {
            
            cat("B")
            
            #Determine the amount a supplier needs to still obtain to fulfill an order
            supplierDemand <-
              orderVol - with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved) 
            
            #Update SupplyReserved for supplier
            sus_coalition(model_state, entityID=src_node)$SupplyReserved <-
              orderVol - supplierDemand + sus_coalition(model_state, entityID = src_node)$SupplyReserved
            
            #Identify the suppliers of the original supplier
            src_suppliers <- locate_suppliers(model_state, src_node)
            src_suppliers <-
              src_suppliers[!src_suppliers %in% c(src_node, entityID)]
            
            if (length(src_suppliers) > 0) {
              #Update secondary supplier information and order them according to price
              src_supplierParams <- sus_coalition(model_state, src_suppliers)
              src_supplierParams <- src_supplierParams[order(src_supplierParams$pricePerUnit), ]
              
              #Work through secondary suppliers so long as orders placed have not yet met supplierDemand and secondary suppliers remain
              while (supplierDemand > 0 &&
                     nrow(src_supplierParams) > 0) {
                #Check if the secondary supplier with the lowest price has supplies to be purchased
                if (with(src_supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved) > 0) {
                  
                  altsrc_node <- src_supplierParams[1, "entityID"]
                  
                  #If supplies are available, place order equal either to remaining demand or the available suppliers, whichever is lower
                  supplierOrderVol <-
                    min(supplierDemand, with(src_supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved))
                  
                  #Update supplier demand
                  supplierDemand <- supplierDemand - supplierOrderVol
                  model_state <- 
                    sus_add_order(model_state, src_node, altsrc_node, supplierOrderVol,
                                  sus_coalition(model_state, entityID=altsrc_node)$pricePerUnit)
                  
                  #Update supply reserved for secondary supplier
                  sus_coalition(model_state, entityID=altsrc_node)$SupplyReserved <-
                    sus_coalition(model_state, entityID=altsrc_node)$SupplyReserved 
                }
                #Remove current secondary supplier from future consideration for the current primary supplier
                src_supplierParams <- src_supplierParams[-1,]
              }
            }
          } else{
            #Update SupplyReserved for supplier
            sus_coalition(model_state, entityID=src_node)$SupplyReserved <-
              orderVol + sus_coalition(model_state, entityID=src_node)$SupplyReserved 
          }
        }
        supplierParams <- supplierParams[-1,]
        cat("Z\n")
      }
      
    }
    
  } else {
    
    #Work through potential suppliers so long as: 
    #(i) orders placed have not yet met current demand, 
    #(ii) there are suppliers remaining, 
    #(iii) there are funds remaining in usableLiquidFunds
    
    while (sus_coalition(model_state, entityID)$OrderVolPlaced < total_volume &&
           nrow(supplierParams) > 0 && 
           (sus_coalition(model_state, entityID)$usableLiquidFunds - fundsAllocated) > min(supplierParams[,"pricePerUnit"])) {
      
      src_node = supplierParams[1, "entityID"]
      cat("src_nod:", src_node, "ordered:", sus_coalition(model_state, entityID)$OrderVolPlaced, "desired", total_volume , "\n")
      
      #Check if the supplier with the lowest price per unit has supplies to purchase (either currently held or able to be ordered from a secondary supplier)
      if (with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved + potentialSupply) > 0) {
        cat("A")
        
        #Target order volume is equal either to remaining demand or the rest of the supplier's stock, whichever is lower
        targetOrderVol <-
          min(total_volume - sus_coalition(model_state, entityID)$OrderVolPlaced,
              with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved + potentialSupply))
        if (targetOrderVol * sus_coalition(model_state, entityID = src_node)$pricePerUnit <= (sus_coalition(model_state, entityID)$usableLiquidFunds - fundsAllocated)) {
          orderVol <- targetOrderVol
        } else{
          orderVol <- targetOrderVol - 
            ceiling(((targetOrderVol * sus_coalition(model_state, entityID = src_node)$pricePerUnit) -
                       (sus_coalition(model_state, entityID)$usableLiquidFunds - fundsAllocated)) /
                      sus_coalition(model_state, entityID = src_node)$pricePerUnit)
        }
        
        #Record order details in the node-level parameters
        model_state <- 
          sus_add_order(model_state, entityID, src_node, orderVol,
                        sus_coalition(model_state, entityID = src_node)$pricePerUnit)
        
        #Calculate funds a node/coalition has allocated towards orders already
        fundsAllocated <- fundsAllocated + (orderVol * sus_coalition(model_state, entityID = src_node)$pricePerUnit)
        
        #If the order volume exceeds a supplier's current and anticipated stock, the remainder must be drawn from its secondary suppliers
        #In practice, this only applies to transfer nodes
        if (orderVol > with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved)) {
          
          cat("B")
          
          #Determine the amount a supplier needs to still obtain to fulfill an order
          supplierDemand <-
            orderVol - with(supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved) 
          
          #Update SupplyReserved for supplier
          sus_coalition(model_state, entityID=src_node)$SupplyReserved <-
            orderVol - supplierDemand + sus_coalition(model_state, entityID = src_node)$SupplyReserved
          
          #Identify the suppliers of the original supplier
          src_suppliers <- locate_suppliers(model_state, src_node)
          src_suppliers <-
            src_suppliers[!src_suppliers %in% c(src_node, entityID)]
          
          if (length(src_suppliers) > 0) {
            #Update secondary supplier information and order them according to price
            src_supplierParams <- sus_coalition(model_state, src_suppliers)
            src_supplierParams <- src_supplierParams[order(src_supplierParams$pricePerUnit), ]
            
            #Work through secondary suppliers so long as orders placed have not yet met supplierDemand and secondary suppliers remain
            while (supplierDemand > 0 &&
                   nrow(src_supplierParams) > 0) {
              #Check if the secondary supplier with the lowest price has supplies to be purchased
              if (with(src_supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved) > 0) {
                
                altsrc_node <- src_supplierParams[1, "entityID"]
                
                #If supplies are available, place order equal either to remaining demand or the available suppliers, whichever is lower
                supplierOrderVol <-
                  min(supplierDemand, with(src_supplierParams[1,,drop=F], currentSupply + anticReplenish - SupplyReserved))
                
                #Update supplier demand
                supplierDemand <- supplierDemand - supplierOrderVol
                model_state <- 
                  sus_add_order(model_state, src_node, altsrc_node, supplierOrderVol,
                                sus_coalition(model_state, entityID=altsrc_node)$pricePerUnit)
                
                #Update supply reserved for secondary supplier
                sus_coalition(model_state, entityID=altsrc_node)$SupplyReserved <-
                  sus_coalition(model_state, entityID=altsrc_node)$SupplyReserved 
              }
              #Remove current secondary supplier from future consideration for the current primary supplier
              src_supplierParams <- src_supplierParams[-1,]
            }
          }
        } else{
          #Update SupplyReserved for supplier
          sus_coalition(model_state, entityID=src_node)$SupplyReserved <-
            orderVol + sus_coalition(model_state, entityID=src_node)$SupplyReserved 
        }
      }
      supplierParams <- supplierParams[-1,]
      cat("Z\n")
    }
  }
  return(model_state)
}

#Estimate the amount of buffer stock that a transfer or source node seeks to maintain
estimate_target_buffer <-
  function(memory,
           timeStep,
           adjMatrix,
           bufferPercent,
           nodeParams) {
    
    nodeParamsTemp <- nodeParams
    
    if (memory > 0) {
      if (timeStep > memory) {
        maxOrderVolReceived <- rep(0, nrow(adjMatrix))
        meanOrderVolReceived <- rep(0, nrow(adjMatrix))
        timeRange <- c(timeStep:c(timeStep - memory))
        
        #Modified to work with disruptions
        #Desired buffers are calculated for all nodes at once
        #This is done twice to allow source nodes to account for new orders placed by transfer nodes seeking to meet their buffer
        for (i in 1:nrow(adjMatrix)) {
          orderVolReceivedTemp <- c()
          for (r in timeRange) {
            orderVolReceivedTemp <-
              c(orderVolReceivedTemp,
                nodeParamsTemp$OrderVolReceived[nodeParamsTemp$timeStep == r][i])
          }
          maxOrderVolReceived[i] <- max(orderVolReceivedTemp)
          meanOrderVolReceived[i] <- mean(orderVolReceivedTemp)
        }
        
        bufferSupply <-
          round(maxOrderVolReceived - meanOrderVolReceived, digits = 0)
        return(bufferSupply)
      } else{
        bufferSupply <- rep(0, nrow(adjMatrix))
        for (i in 1:nrow(adjMatrix)) {
          bufferSupply[i] <-
            round(bufferPercent * nodeParamsTemp$OrderVolReceived[nodeParamsTemp$timeStep == timeStep][i],
                  digits = 0)
        }
        return(bufferSupply)
      }
    } else{
      bufferSupply <- rep(0, nrow(adjMatrix))
      for (i in 1:nrow(adjMatrix)) {
        bufferSupply[i] <-
          round(bufferPercent * nodeParamsTemp$OrderVolReceived[nodeParamsTemp$timeStep == timeStep][i],
                digits = 0)
      }
      return(bufferSupply)
    }
  }

current_sink_buyers <- function(model_state){
  maximal(model_state$coalSkel) %>% as.list() %>%
    purrr::map(~purrr::keep(.x,function(xx) xx %in% model_state$sink_nodes)) %>%
    purrr::keep(~length(.x)>0)
}

#For each transfer node, update its record of the total supplies held by its upstream suppliers or fellow transfer nodes to which it is connected
#These potential supplies can be ordered by the transfer node to fulfill orders that it itself receives
#Potential supplies include both supplies that are currently held by a supplier, as well as its anticipated amount of resupply, minus any supplies already reserved for previous orders
update_potential_supplies <-
  function(adjMatrix,
           transferNodeList,
           sourceNodeList,
           st,
           timeStep,
           nodeParams) {
    
    potentialSupply <- rep(0, nrow(adjMatrix))
    for (n in 1:nrow(adjMatrix)) {
      if (nodeParams$nodeType[nodeParams$timeStep == timeStep][n] == "Transfer") {
        suppliersT <-
          cofaceList(st = st, nodeID = nodeParams$nodeID[nodeParams$timeStep == timeStep][n])
        suppliersT <-
          suppliersT[suppliersT %in% c(transferNodeList, sourceNodeList) &
                       !suppliersT == nodeParams$nodeID[nodeParams$timeStep == timeStep][n]]
        potentialSupply[n] <-
          sum((nodeParams$currentSupply[nodeParams$timeStep == timeStep &
                                          nodeParams$nodeID %in% suppliersT] +
                 nodeParams$anticReplenish[nodeParams$timeStep == timeStep &
                                             nodeParams$nodeID %in% suppliersT]) -
                nodeParams$SupplyReserved[nodeParams$timeStep == timeStep &
                                            nodeParams$nodeID %in% suppliersT]
          )
      } else {
        potentialSupply[n] <- 0
      }
    }
    return(potentialSupply)
  }

simple_distribute_stock <- function(model_state, supplierID, records){
  
  #Largest orders are prioritized
  records <- records %>% arrange(-Volume)
  
  cat("Stub for fulfilling orders by:",supplierID,"\n" )
  cat("n orders:",nrow(records),
      "units:",sum(records$Volume),
      "value:",sum(records$Volume*records$pricePerUnit), "\n" )
  
  #While supply exists and orders remain to be filled, distribute supplies
  while (sus_coalition(model_state, supplierID)$currentSupply > 0 &&
         nrow(records) > 0) {
    #Either deliver the full order amount if supplies permit or otherwise deliver whatever stock remains
    deliverVolume <-
      min(sus_coalition(model_state, supplierID)$currentSupply, records$Volume[1])
    model_state <- 
      model_state$make_delivery(model_state, 
                                from=supplierID,
                                to=records$Customer[1],
                                unit_volume= deliverVolume,
                                total_price=round(deliverVolume*records$pricePerUnit[1],digits = 2))
    
    #Remove transaction from supplier's list
    records <- records[-1,]
  }
  return(model_state)
}

fail_distribute_stock_coalition <- function(model_state, supplierID, records){
  stop("No coalitions should be distributing stock!")
}

#Handles final payment by sink nodes for deliveries received
make_delivery_with_budgets <- function(model_state, from, to, unit_volume, total_price){
  
  cat("Stub for transaction:",from,"=>",to,"\n" )
  sus_coalition(model_state, entityID = from)$currentSupply <-
    sus_coalition(model_state, entityID = from)$currentSupply - unit_volume
  sus_coalition(model_state, entityID = from)$DeliveryVolMade <-
    sus_coalition(model_state, entityID = from)$DeliveryVolMade + unit_volume
  
  sus_coalition(model_state, entityID = to)$currentSupply <-
    sus_coalition(model_state, entityID = to)$currentSupply + unit_volume
  sus_coalition(model_state, entityID = to)$DeliveryVolReceived <-
    sus_coalition(model_state, entityID = to)$DeliveryVolReceived + unit_volume
  
  sus_coalition(model_state, entityID = from)$PaymentReceived <-
    sus_coalition(model_state, entityID = from)$PaymentReceived + total_price
  sus_coalition(model_state, entityID = to)$PaymentMade <-
    sus_coalition(model_state, entityID = to)$PaymentMade + total_price
  
  if(sus_coalition(model_state, entityID = to)$nodeType == "Sink") {
  sus_coalition(model_state, entityID = to)$usableLiquidFunds <- 
    sus_coalition(model_state, entityID = to)$usableLiquidFunds -  total_price
  }
  
  stopifnot(sus_coalition(model_state, entityID=from)$currentSupply >=0)
  return(model_state)
}

#When subsidized distribution is active, coalition members with excess liquid funds will effectively purchase units
#For members that don't have sufficient funds to purchase their portion of the joint order
distribute_stock_coalition_with_subsidies <- function(model_state, supplierID, records){
  cat("Stub for fulfilling orders by:",supplierID,"\n" )
  cat("n orders:",nrow(records),
      "units:",sum(records$Volume),
      "value:",sum(records$Volume*records$pricePerUnit), "\n" )
  
  records$VolumePart <- 
    proportion_vol_fully(avail=sus_coalition(model_state,supplierID)$DeliveryVolReceived,
                         parts = records[, "Volume"],
                         total = sum(records[, "Volume"]))
  
  records$liquidFunds <- sus_node(model_state$node_params, 
                                  timeStep = model_state$current_step, 
                                  nodeID = as.coal_set.char(supplierID)[[1]])$usableLiquidFunds
  records$SubsidyNeeded <-
    (records$pricePerUnit * records$VolumePart)-records$liquidFunds
  records <- records[order(records$SubsidyNeeded, decreasing = TRUE),]
  totalSubsidiesNeeded <- sum(records$SubsidyNeeded[which(records$SubsidyNeeded > 0)])
  totalSubsidiesAvailable <- sum(abs(records$SubsidyNeeded[which(records$SubsidyNeeded <= 0)]))
  
  for(r in 1:nrow(records))
    print(records)
  #While supply exists and orders remain to be filled, distribute supplies
  while (sus_coalition(model_state, supplierID)$currentSupply > 0 &&
         nrow(records) > 0) {
    
    #Either deliver the full order amount if supplies permit or otherwise deliver whatever stock remains
    deliverVolume <-
      min(sus_coalition(model_state, supplierID)$currentSupply, records$VolumePart[1])
    remainingVolume <-
      min(sus_coalition(model_state, supplierID)$currentSupply-deliverVolume, sum(records$VolumePart[-1]))
    remainingPayment <-
      with(sus_coalition(model_state,supplierID), PaymentMade-PaymentReceived)
    if(records$SubsidyNeeded[1] > 0){
      memberPayment <- records$liquidFunds[1]
    } else{
      memberPayment <- (deliverVolume * records$pricePerUnit[1]) + 
        round(totalSubsidiesNeeded * (abs(records$SubsidyNeeded[1])/totalSubsidiesAvailable), digits = 2)
    }
    model_state <- 
      model_state$make_delivery(model_state, 
                                from=supplierID,
                                to=records$Customer[1],
                                unit_volume= deliverVolume,
                                total_price= memberPayment)
    
    #Remove transaction from supplier's list
    records <- records[-1,]
  }
  print(sus_coalition(model_state, supplierID))
  print(sus_coalition(model_state, supplierID%>%as.coal_set.char()%>%
                        unlist()%>%as.list()%>%as.char.coal_set()))
  return(model_state)
}

#When unsubsidized distribution is active, coalition members will purchase from the joint order in proportion to their liquid funds
#This means that a node with more funds may purchase more supplies than it can effectively use, but will not subsidize the purchases of its fellow members
distribute_stock_coalition_without_subsidies <- function(model_state, supplierID, records){
  cat("Stub for fulfilling orders by:",supplierID,"\n" )
  cat("n orders:",nrow(records),
      "units:",sum(records$Volume),
      "value:",sum(records$Volume*records$pricePerUnit), "\n" )
  
  records$VolumePart <- 
    proportion_vol_fully(avail=sus_coalition(model_state,supplierID)$DeliveryVolReceived,
                         parts = records[, "Volume"],
                         total = sum(records[, "Volume"]))
  records$liquidFunds <- sus_node(model_state$node_params, 
                                  timeStep = model_state$current_step, 
                                  nodeID = as.coal_set.char(supplierID)[[1]])$usableLiquidFunds
  records$SubsidyNeeded <-
    (records$pricePerUnit * records$VolumePart)-records$liquidFunds
  records <- records[order(records$SubsidyNeeded, decreasing = TRUE),]
  totalSubsidiesNeeded <- sum(records$SubsidyNeeded[which(records$SubsidyNeeded > 0)])
  totalSubsidiesAvailable <- sum(abs(records$SubsidyNeeded[which(records$SubsidyNeeded <= 0)]))
  records$VolumeByFunds <- 
    round((records$liquidFunds/sum(records$liquidFunds))*sum(records$VolumePart))
  for(r in 1:nrow(records))
    print(records)
  #While supply exists and orders remain to be filled, distribute supplies
  while (sus_coalition(model_state, supplierID)$currentSupply > 0 &&
         nrow(records) > 0) {
    #Either deliver the full order amount if supplies permit or otherwise deliver whatever stock remains
    deliverVolume <-
      min(sus_coalition(model_state, supplierID)$currentSupply, records$VolumeByFunds[1])
    
    remainingVolume <-
      min(sus_coalition(model_state, supplierID)$currentSupply-deliverVolume, sum(records$VolumeByFunds[-1]))
    remainingPayment <-
      with(sus_coalition(model_state,supplierID), PaymentMade-PaymentReceived)
    
    model_state <- 
      model_state$make_delivery(model_state, 
                                from=supplierID,
                                to=records$Customer[1],
                                unit_volume= deliverVolume,
                                total_price= remainingPayment-round(remainingVolume*records$pricePerUnit[1],digits = 2))
    
    #Remove transaction from supplier's list
    records <- records[-1,]
  }
  print(sus_coalition(model_state, supplierID))
  print(sus_coalition(model_state, supplierID%>%as.coal_set.char()%>%
                        unlist()%>%as.list()%>%as.char.coal_set()))
  return(model_state)
}

locate_suppliers <- function(model_state, entityID, include_transfers=T ){
  if( length(entityID) != 1) stop("Multiple entity supplier search not implemented.")
  if(is.character(entityID)) entityID = as.coal_set.char(entityID)
  all_sellers <-current_sellers(model_state, include_transfers)
  entityID[[1]] %>% purrr::map(~ cofaceList(model_state$oneSkel,.x)) %>% 
    unlist() %>% unique() %>% sort() %>% 
    purrr::keep(~  ! .x %in% entityID[[1]]) -> seller_nodes
  all_sellers %>% purrr::keep(~any(.x%in%seller_nodes) ) %>% names()
}

current_sellers <- function(model_state, include_transfers=T) {
  if(!include_transfers){
    return(maximal(model_state$coalSkel) %>% as.list() %>%
             purrr::map(~purrr::keep(.x,function(xx) (xx %in% model_state$source_nodes))) %>%
             purrr::keep(~length(.x)>0) %>% as.char.coal_set() %>% as.coal_set.char())
  }
  maximal(model_state$coalSkel) %>% as.list() %>%
    purrr::map(~purrr::keep(.x,function(xx) !(xx %in% model_state$sink_nodes))) %>%
    purrr::keep(~length(.x)>0) %>% as.char.coal_set() %>% as.coal_set.char()
}


#Once all supplies have been distributed on a given time step, sink nodes consume available supplies
#Supplies are consumed either until current gross demand has been met or until supplies have been exhausted
consume_supplies <- function(nodeID, timeStep, nodeParams) {
  
  currentNodeParams <- nodeParams[nodeParams$timeStep == timeStep, ]
  
  consumeAmount <-
    min(currentNodeParams$currentGrossDemand[currentNodeParams$nodeID == nodeID],
        currentNodeParams$currentSupply[currentNodeParams$nodeID == nodeID])
  
  #Update customer's excess demand (i.e. gross demand - consumed amount)
  currentNodeParams$ExcessDemand[currentNodeParams$nodeID == nodeID] <-
    currentNodeParams$currentGrossDemand[currentNodeParams$nodeID == nodeID] - consumeAmount
  
  #Update customer's supply
  currentNodeParams$currentSupply[currentNodeParams$nodeID == nodeID] <-
    currentNodeParams$currentSupply[currentNodeParams$nodeID == nodeID] - consumeAmount
  
  return(currentNodeParams)
}

consume_supplies_multi <- function(sN, tS, nP, shuffle=TRUE, updateCurrent=TRUE){
  if(shuffle) sN <- c(sample(sN, length(sN), replace = FALSE))
  for (i in sN) {
    consumptionOutput <- consume_supplies(nodeID = i, timeStep = tS, nodeParams = nP)
    
    #Update node-level parameters
    nP[nP$timeStep == tS, ] <- as.data.frame(consumptionOutput)
  }
  if(updateCurrent)
    return(nP[nP$timeStep == tS, ])
  return(nP)
}

trace_potential_supplies <- function(model_state, nodeID){
  if(nodeID %in% model_state$transfer_nodes ) {
    sellers <- c(nodeID, get_node_suppliers(model_state = model_state,
                                            nodeID = nodeID))
  } else {
    sellers <- nodeID
  }
  data.frame(nodeID=nodeID, 
             src_ID=sellers, 
             pricePerUnit= sus_node(model_state$node_params,
                                    model_state$current_step,
                                    nodeID)$pricePerUnit) %>%
    inner_join(available_supplies(model_state,sellers), by ="src_ID", copy=T) -> res
  res  
}

available_supplies <- function(model_state, nodeID){
  sus_node(model_state$node_params,model_state$current_step,nodeID) %>%
    mutate(availableSupplies= currentSupply+ anticReplenish  -SupplyReserved ) %>%
    select(src_ID=nodeID,availableSupplies)
}

######################################################
#Utility functions for manipulating the node dataframe
######################################################

sus_node <- function(nP, timeStep, nodeID) {
  if(missing(timeStep)) return(nP[nP$nodeID %in% nodeID,,drop=FALSE])
  if(missing(nodeID)) return(nP[nP$timeStep %in% timeStep,,drop=FALSE])
  return(nP[nP$timeStep %in% timeStep & nP$nodeID %in% nodeID ,,drop=FALSE])
}

sus_coalition <- function(model_state, entityID, delta_step=0) {
  timeStep= model_state$current_step + delta_step
  if(missing(entityID)) stop("missing(entityID) is not implemented.")
  if(is.numeric(entityID)) {
    warning("entityID should be in coalition format. Treating as a list of vertices.")
    entityID=as.coal_set.char(as.char.coal_set(as.list(entityID)))
  } else if(is.character(entityID)){
    entityID=as.coal_set.char(entityID)
  } else {
    entityID = entityID %>% as.char.coal_set() %>% as.coal_set.char()
  }
  
  xx <- names(entityID)
  if(any(grepl("COAL_COAL_",xx))){
    print(entityID)
    print(xx)
    stop("entityID format Violation at selector")
  }
  
  if(all( (entityID %>% sapply(length) )==1) ) {
    #No coalitions
    nP = sus_node(model_state$node_params,timeStep = timeStep, nodeID = unlist(entityID))
    if(any(grepl("COAL_COAL_",row.names(nP)))){
      print(nP)
      stop("entityID format Violation at getter node results")
    }
    return(nP)
  }
  
  if(all( (entityID %>% sapply(length) ) > 1) ) {
    #All coalitions
    #Check if coalition_params exists
    if(is.null(model_state$coalition_params)){
      model_state$node_params[numeric(0),,drop=FALSE] %>% 
        add_row(entityID= xx, timeStep= model_state$current_step + delta_step ) %>% 
        mutate(rowname=paste0(entityID,"__",timeStep) ) %>% tibble::column_to_rownames() -> res
      if(any(grepl("COAL_COAL_",row.names(res)))){
        print(res)
        stop("entityID format Violation at getter coalition init")
      }
      return(res)
    }
    
    cP_all=model_state$coalition_params
    cP = cP_all[ cP_all$timeStep %in% timeStep, ,drop=FALSE]
    
    #Check if entities are available 
    if(!all(xx %in% cP$entityID)){
      if(length(timeStep)>1 && length(xx)>1) warning("Indexing a matrix with missing elements is unreliable")
      add_xx=xx[!xx %in% cP$entityID]
      cP_add = cP[numeric(0),,drop=FALSE]  %>% 
        add_row(entityID= add_xx, timeStep= model_state$current_step + delta_step ) %>% 
        mutate(rowname=paste0(entityID,"__",timeStep) ) %>% tibble::column_to_rownames()
      cP=rbind(cP,cP_add)
    }
    
    if(any(grepl("COAL_COAL_",row.names(cP)))){
      print(cP)
      print(cP[cP$entityID %in% xx,,drop=FALSE])
      stop("entityID format Violation at getter coalition select")
    }
    
    return(cP[cP$entityID %in% xx,,drop=FALSE])
  }
  
  # mixed case scenario, split up as needed
  key_noC <- xx[ (entityID %>% sapply(length) )==1 ]
  key_C <- xx[ (entityID %>% sapply(length) )>1 ]
  return(rbind(sus_coalition(model_state = model_state, entityID = key_noC, delta_step = delta_step),
               sus_coalition(model_state = model_state, entityID = key_C, delta_step = delta_step)))
}

as.char.coal_set <- function(x, ...) {
  if(is.list(x)) {
    if(length(x)==0) return(character(0))
    return(sapply(x,function(xx) 
      sprintf("COAL_%s", paste0(xx %>% unique %>% sort,collapse = "_")) ))
  }
  if(is.character(x)) return(x)
  if(is.factor(x)) return(as.character(x))
  sprintf("COAL_%s", paste0(x %>% unique %>% sort,collapse = "_"))
}

as.coal_set.char <- function(x, ...) {
  if(is.factor(x)) x<-as.character(x)
  nn <- substring(x,6) %>% strsplit("_") %>% lapply(as.numeric)
  names(nn) <- x
  nn
}

ensure_param <- function(nP_start, param_name){
  if(! param_name %in% names(nP_start) )
    nP_start[[param_name]] <- numeric(nrow(nP_start))
  return(nP_start)
}

ensure_param_def <- function(nP_start, param_name, def_value){
  if(! param_name %in% names(nP_start) )
    nP_start[[param_name]] <- rep(def_value,length.out=nrow(nP_start))
  return(nP_start)
}

as.node_set <- function(x) {
  if(is.character(x)){
    x <- as.coal_set.char(x)
  }
  x %>% unlist %>% unique %>% sort
}

sus_add_order <- function(model_state, customerID, supplierID, order_volume, unit_price){
  #Update order volumes placed by the customer
  sus_coalition(model_state, customerID)$OrderVolPlaced <-
    order_volume + sus_coalition(model_state, customerID)$OrderVolPlaced
  
  #Update order volumes received by the supplier
  sus_coalition(model_state, entityID=supplierID)$OrderVolReceived <-
    order_volume + sus_coalition(model_state, entityID=supplierID)$OrderVolReceived
  
  #Update order records
  model_state$orderRecord <-model_state$orderRecord %>%
    add_row(
      TimeStep = model_state$current_step,
      Customer = customerID,
      Supplier = supplierID,
      Volume = order_volume,
      pricePerUnit=unit_price
    )
  model_state
}

proportion_vol_fully <- function(avail, parts, total)  
  rev(diff(c(0, floor(min(avail, total) * cumsum(rev(parts)) / total) )))

rnorm_multi_uint <- function(mean_vec, sd_prop) {
  resL <- length(mean_vec)
  if(resL ==0 ) return(numeric(0))
  sd_vec = sd_prop * mean_vec
  
  res <- rep(0, resL)
  for (i in 1:resL) {
    if (mean_vec[i] > 0) {
      res[i] <-
        max(0, round(rnorm(1, mean_vec[i], sd_vec[i])))
    } else{
      res[i] <- 0
    }
  }
  return(res)
}

`sus_node<-` <- function(nP, timeStep, nodeID, value) {
  if(missing(timeStep))  {
    nP[nP$nodeID %in% nodeID,] <- value
    return(nP)
  }
  if(missing(nodeID))  {
    nP[nP$timeStep %in% timeStep,] <- value
    return(nP)
  }
  nP[nP$timeStep %in% timeStep & nP$nodeID %in% nodeID , ] <- value
  nP
}

sus_node_demand <- function(model_state, nodeID, ...){
  cat("Demand lookup",
      nodeID,
      model_state$node_demand(model_state, nodeID, ...),
      "\n")
  model_state$node_demand(model_state, nodeID, ...)
}

`sus_coalition<-` <- function(model_state, entityID, delta_step=0, value) {
  timeStep= model_state$current_step + delta_step
  if(missing(entityID)) stop("missing(entityID) is not implemented.")
  if(is.numeric(entityID)) {
    warning("entityID should be in coalition format. Treating as a list of vertices.")
    entityID=as.coal_set.char(as.char.coal_set(as.list(entityID)))
  } else if(is.character(entityID)){
    entityID=as.coal_set.char(entityID)
  } else {
    entityID = entityID %>% as.char.coal_set() %>% as.coal_set.char()
  }
  xx <- names(entityID)
  if(!all(value$timeStep %in% timeStep)) stop("Time step changes are not implemented.")
  if(!all(value$entityID %in% xx)) stop("entityID changes are not implemented.")
  
  key_noC <- xx[ (entityID %>% sapply(length) )==1 ]
  key_C <- xx[ (entityID %>% sapply(length) )>1 ]
  if(length(key_noC)>0){
    val_noC <- value[value$entityID %in% key_noC,, drop=F] # %>%tibble::rownames_to_column() 
    model_state$node_params[rownames(val_noC), ] <- val_noC 
  }
  if(length(key_C)>0){
    if(is.null(model_state$coalition_params)){
      model_state$coalition_params <- value[value$entityID %in% key_C,, drop=F] 
    } else {
      val_C <-  value[value$entityID %in% key_C,, drop=F] 
      model_state$coalition_params[rownames(val_C), ] <- val_C 
    }
    
  }
  if(any(grepl("COAL_COAL_",model_state$node_params$entityID))){
    print(value)
    print(xx)
    print(key_noC)
    print(key_C)
    stop("entityID format Violation at setter")
  }
  if(!is.null(model_state$coalition_params))
    if(any(grepl("COAL_COAL_",model_state$coalition_params$entityID))){
      print(value)
      print(xx)
      print(key_noC)
      print(key_C)
      stop("entityID format Violation at setter!")
    }
  return(model_state)
}

sus_node_budget <- function(model_state, nodeID, ...){
  model_state$node_budget(model_state, nodeID, ...)
}

####################################
#Functions for coalition maintenance
####################################

simple_update_coalition_cost <- function(model_state, nodeID, coalition){
  cat(sprintf("Coalition cost stub: n%02d\t%s\n", nodeID,as.char.coal_set(coalition)))
  
  #Represent cost as a payment made
  sus_node(model_state$node_params,
           timeStep = model_state$current_step, 
           nodeID = nodeID)$PaymentMade =
    sus_node(model_state$node_params,
             timeStep = model_state$current_step, 
             nodeID = nodeID)$PaymentMade +
    simple_coalition_cost(model_state, nodeID, coalition)
  
  #Now update the coalition cost matrix for this node:
  for(other_node in coalition){
    model_state$coalition_cost[nodeID,other_node] <-
      max(model_state$coalition_min_cost[nodeID,other_node],
          model_state$coalition_cost[nodeID,other_node]-
            model_state$coalition_cost_change[nodeID,other_node])
  }
  return(model_state)
}

simple_update_coalition_demand <- function(model_state, coalition){
  cat(sprintf("Coalition demand stub: %s\n", as.char.coal_set(coalition)))
  sus_coalition(model_state, 
                entityID=as.char.coal_set(coalition))$currentGrossDemand <-
    with(sus_node(model_state$node_params,
                  timeStep = model_state$current_step, 
                  nodeID = coalition),
         sum(currentGrossDemand - currentSupply))
  
  sus_coalition(model_state, 
                entityID = as.char.coal_set(coalition))$currentSupply <- 0
  return(model_state)
}

update_coalition_liquidFunds <- function(model_state, coalition, t_final){
  cat(sprintf("Coalition budget stub: %s\n", as.char.coal_set(coalition)))
  sus_coalition(model_state, 
                entityID=as.char.coal_set(coalition))$usableLiquidFunds <-
    round(with(sus_node(model_state$node_params,
                        timeStep = model_state$current_step, 
                        nodeID = coalition),
               sum(usableLiquidFunds)), digits = 2)
  return(model_state)
}

#Carries out the actual coalition formation process, based on the strength algorithm and joining rules
update_sinkCoalitions_simple <- function(model_state, add_nodes, remove_nodes, t_final){
  current_coalitions <- current_sink_coalitions(model_state)
  coalition_nodes <- as.node_set(current_coalitions)
  
  #Identify nodes that are considering leaving/modifying coalition membership
  remove_nodes <- remove_nodes %>%
    purrr::keep(~.x %in% coalition_nodes) %>%
    purrr::keep(~.x %in% model_state$sink_nodes)
  #Identify nodes that are considering joining a coalition and that are not currently part of any coalition
  add_nodes <- add_nodes %>% purrr::keep(~!.x %in% coalition_nodes) %>%
    purrr::keep(~.x %in% model_state$sink_nodes)
  if(length(add_nodes) >0 || length(remove_nodes)>0){
    #Carry out the coalition formation process
    st1_coal <- simple_sinkCoalition_formation_stage1(model_state,add_nodes)
    red_coal <- simple_sinkCoalition_reduction(model_state,remove_nodes)
    st2_coal <- simple_sinkCoalition_formation_stage2(model_state,
                                                      c(current_coalitions, red_coal),
                                                      st1_coal)
    #Remove nodes (and cofaces) reconsidering the coalition
    remove(model_state$coalSkel,remove_nodes %>% as.list())
    #but keep the vertices
    insert(model_state$coalSkel,remove_nodes %>% as.list())
    #Now add in all selected coalitions
    insert(model_state$coalSkel,st2_coal)
    insert(model_state$assocSkel,st2_coal)
  } else {
    st2_coal <- current_coalitions
  }
  
  for(c_i in st2_coal){
    c_i_label = as.char.coal_set(c_i)
    if(length(c_i) > 1 )  {
      sus_node(model_state$node_params, model_state$current_step, c_i)$coalition <- c_i_label
      #Set default coalition time
      sus_node(model_state$node_params,model_state$current_step, c_i)$coalitionTime =1
    } else {
      sus_node(model_state$node_params,model_state$current_step, c_i)$coalition <- sprintf("n%02d", c_i)
      #Set default coalition time
      sus_node(model_state$node_params,model_state$current_step, c_i)$coalitionTime = 0
    }
    if(sus_node(model_state$node_params,model_state$current_step - 1, c_i[1])$coalition == c_i_label &&
       length(c_i) > 1)
      sus_node(model_state$node_params,model_state$current_step, c_i)$coalitionTime<-
        sus_node(model_state$node_params,model_state$current_step - 1, c_i)$coalitionTime +1
    #Initialize basic coalition paramters
    if(length(c_i) > 1)  {
      sus_coalition(model_state, entityID = c_i_label)$nodeType = "Sink"
      sus_coalition(model_state, 
                    entityID = c_i_label)[c("OrderVolReceived",
                                            "OrderVolPlaced",
                                            "DeliveryVolReceived",
                                            "DeliveryVolMade",
                                            "PaymentMade",
                                            "PaymentReceived")] <- 0
      
      for(c_i_node in c_i) {
        model_state <- model_state$update_coalition_cost(model_state, c_i_node, c_i)
      }
      model_state <- model_state$update_coalition_demand(model_state, c_i)
      model_state <- model_state$update_coalition_budget(model_state, c_i, t_final = t_final)
      sus_coalition(model_state, entityID = c_i_label)$coalition = c_i_label
      sus_coalition(model_state, entityID = c_i_label)$coalitionTime =
        sus_node(model_state$node_params,model_state$current_step,c_i[1])$coalitionTime
    }
  }
  return(model_state)
}

simple_coalition_cost <- function(model_state, nodeID, coalition){
  return(model_state$coalition_cost[nodeID, coalition, drop=F] %>% rowSums())
}

###################################
#Useful coalition-related functions
###################################

#Generate a list of all nodes to which a given node is directly connected (i.e. it's 1-faces)
cofaceList <- function(st, nodeID) {
  oneCofaces <- as.list(cofaces(st, nodeID))
  stTemp <- simplex_tree()
  stTemp %>% insert(oneCofaces)
  return(as.vector(stTemp$vertices))
}

#Identify the coalitions currently active
current_sink_coalitions <- function(model_state){
  maximal(model_state$coalSkel) %>% as.list() %>%
    purrr::map(~purrr::keep(.x,function(xx) xx %in% model_state$sink_nodes)) %>%
    purrr::keep(~length(.x)>1)
}

###############################################
#Functions to perform the simulation iterations
###############################################

simulationLoop_v1.3 <- 
  function(t_start, 
           t_end, 
           model_state,
           memory = 0,
           desiredBufferPercentage,
           iteration_handler = simple_iteration_handler, 
           inventoryType, 
           purchase_Strategy, 
           failure_Response){
    
    if(any(grepl("COAL_COAL_", model_state$node_params$entityID)))
      stop("entityID format Violation at START!")
    if(!is.null(model_state$coalition_params))
      if(any(grepl("COAL_COAL_", model_state$coalition_params$entityID)))
        stop("entityID format Violation at START!")
    
    for(loop_t in t_start:t_end){
      {
        
        #Sink nodes adjust purchasing strategy if their current strategy is not working for them
        if(purchase_Strategy == "mixed") {
          model_state = adjust_purchase_strategy(model_state)
        }
        
        if(failure_Response == "addSupplier") {
          model_state = add_supplier(model_state)
        }
        
        #First, determine if coalitions may form or dissolve this time step
        #This occurs if at least one node is considering joining/forming a coalition, leaving a coalition, or if there is already at least one coalition already formed
        add_rule_nodes <- model_state$add_rule(model_state)
        leave_rule_nodes <- model_state$leave_rule(model_state)
        if(length(add_rule_nodes) > 0 || length(leave_rule_nodes) > 0 ||
           length(current_sink_coalitions(model_state)) > 0 ){
          model_state <- iteration_handler$update_coalitions(model_state, add_rule_nodes, leave_rule_nodes, t_final)
        }
        
        if(any(grepl("COAL_COAL_", model_state$node_params$entityID)))
          stop("entityID format Violation at update coalitions!")
        if(!is.null(model_state$coalition_params))
          if(any(grepl("COAL_COAL_", model_state$coalition_params$entityID)))
            stop("entityID format Violation at update coalitions!")
        
        #Place initial orders for sink nodes and coalitions
        model_state = iteration_handler$place_sink_orders(model_state, failure_Response)
        
        if(any(grepl("COAL_COAL_", model_state$node_params$entityID)))
          stop("entityID format Violation at sink orders! nodes")
        if(!is.null(model_state$coalition_params))
          if(any(grepl("COAL_COAL_",model_state$coalition_params$entityID)))
            stop("entityID format Violation at sink orders! coalitions")
        
        #All initial orders have been placed
        #Now have transfer and source nodes place buffer orders
        model_state = iteration_handler$place_buffer_orders(model_state, memory, desiredBufferPercentage, inventoryType)
        
        if(any(grepl("COAL_COAL_", model_state$node_params$entityID)))
          stop("entityID format Violation at buffer orders!")
        if(!is.null(model_state$coalition_params))
          if(any(grepl("COAL_COAL_",model_state$coalition_params$entityID)))
            stop("entityID format Violation at buffer orders!")
        
        #Next, generate source nodes' actual supply based on their anticipated supply
        for (i in current_sellers(model_state, include_transfers = F) %>% names()) {
          #Only source nodes (i.e. nodes that replenish supplies) generate supply
          if (sus_coalition(model_state, i)$meanReplenish > 0) {
            # Calculate total desired replenishment volume 
            # (equal to a node's desired buffer plus current order volumes received minus any supplies it currently holds)
            desiredReplenish <- 
              max(0, with(sus_coalition(model_state, i),desiredBuffer + OrderVolReceived - currentSupply))
            
            # Generate supply by drawing from a normal distribution with mean equal to a node's anticipated replenishment this turn and standard deviation equal to 5% of this value
            potentialGenerated <- 
              with(sus_coalition(model_state, i),
                   max(0, round(rnorm(1,anticReplenish,(0.05 *anticReplenish)))))
            
            #The actual supply generated is either the total desired replenishment volume or the volume generated (whichever is less)
            supplyGenerated <- min(desiredReplenish, potentialGenerated)
          } else{
            supplyGenerated <- 0
          }
          sus_coalition(model_state, i)$currentSupply = sus_coalition(model_state, i)$currentSupply + supplyGenerated
          sus_coalition(model_state, i)$supplyGenerated = sus_coalition(model_state, i)$supplyGenerated + supplyGenerated
        }
        
        if(any(grepl("COAL_COAL_", model_state$node_params$entityID)))
          stop("entityID format Violation at generate!")
        if(!is.null(model_state$coalition_params))
          if(any(grepl("COAL_COAL_", model_state$coalition_params$entityID)))
            stop("entityID format Violation at generate!")
        
        #Supply has been generated
        #Distribute deliveries; source nodes act first, followed by transfer nodes
        model_state <- iteration_handler$distribute_stock(model_state, inventoryType)
        
        if(any(grepl("COAL_COAL_", model_state$node_params$entityID)))
          stop("entityID format Violation at distribute!")
        if(!is.null(model_state$coalition_params))
          if(any(grepl("COAL_COAL_", model_state$coalition_params$entityID)))
            stop("entityID format Violation at distribute!")
        
        #Consume units of inventory up to current demand or until inventory is exhausted
        model_state <- iteration_handler$consume_supplies(model_state)
        
        if(any(grepl("COAL_COAL_", model_state$node_params$entityID)))
          stop("entityID format Violation at consume!")
        if(!is.null(model_state$coalition_params))
          if(any(grepl("COAL_COAL_", model_state$coalition_params$entityID)))
            stop("entityID format Violation at consume!")
      }    
      
      # Generate values for the next time stamp
      {
        #If the end of the simulation has not yet been reached, dispatch to iteration handler's next_modelState
        if (loop_t < t_end) {
          model_state <- iteration_handler$next_modelState(model_state, inventoryType, t_final, failure_Response)
          if(any(grepl("COAL_COAL_", model_state$node_params$entityID)))
            stop("entityID format Violation at next_modelState!")
          if(!is.null(model_state$coalition_params))
            if(any(grepl("COAL_COAL_", model_state$coalition_params$entityID)))
              stop("entityID format Violation at next_modelState!")
        } else{
          break
        }
      }
    }
    return(model_state)
  }

#Sets up the node-level parameters for the next time step
next_modelState_simple <- function(model_state, inventoryType, t_final, failure_Response){
  
  nodeParams = model_state$node_params
  loop_t = model_state$current_step
  cat("Setting up model_state for iteration:",loop_t+1,"\n")
  #First, copy over the current node-level data
  futureNodeParams <- nodeParams[which(nodeParams$timeStep == loop_t),]
  cat("step ")
  #Update timeStep to loop_t + 1
  futureNodeParams$timeStep <- loop_t + 1
  
  #Update nodes' time since switch
  futureNodeParams$timeSinceSwitch <- futureNodeParams$timeSinceSwitch + 1
  
  if(failure_Response == "none") {} else {
    failureResponseVect <- rep(0, nrow(futureNodeParams))
    for(n in 1:nrow(futureNodeParams)) {
      if(futureNodeParams[n,"currentGrossDemand"] > 0) {
        if(futureNodeParams[n, "ExcessDemand"]/futureNodeParams[n, "currentGrossDemand"] > 0.2) {
          failureResponseVect[n] <- 1
        }
      }
    }
    futureNodeParams$failureResponse <- failureResponseVect
  }
  
  #Generate anticipated replenishment for source nodes (as during simulation set up)
  cat("antic ")
  futureNodeParams$anticReplenish <-
    rnorm_multi_uint(futureNodeParams$meanReplenish, 0.1)
  
  #Set gross demand for next time step (as during simulation set up)
  cat("cGD ")
  futureNodeParams$currentGrossDemand <-
    rnorm_multi_uint(futureNodeParams$meanGrossDemand, 0.1)
  
  #If the inventory type is perishable, any unused supplies are lost
  if(inventoryType == "Perishable"){
    futureNodeParams$supplyLost <- futureNodeParams$currentSupply
    futureNodeParams$currentSupply <- 0
  }
  
  #Determine usable liquid funds for next time step
  if("liquidFunds" %in% names(futureNodeParams)){
    futureNodeParams$liquidFunds[which(futureNodeParams$nodeType == "Sink")] <- futureNodeParams$liquidFunds[which(futureNodeParams$nodeType == "Sink")] - 
      futureNodeParams$PaymentMade[which(futureNodeParams$nodeType == "Sink")]
    futureNodeParams$usableLiquidFunds <- round(futureNodeParams$liquidFunds / (t_final - futureNodeParams$timeStep + 1), digits = 2)
  }
  cat("000 ")
  
  #All other parameters are set to 0
  futureNodeParams$OrderVolPlaced <- 0
  futureNodeParams$OrderVolReceived <- 0
  futureNodeParams$DeliveryVolReceived <- 0
  futureNodeParams$DeliveryVolMade <- 0
  futureNodeParams$SupplyReserved <- 0
  futureNodeParams$PaymentMade <- 0
  futureNodeParams$PaymentReceived <- 0
  futureNodeParams$supplyGenerated <- 0
  futureNodeParams$ExcessDemand <- 0
  futureNodeParams$desiredBuffer <- 0
  cat("rowID ")
  row.names(futureNodeParams) <- with(futureNodeParams, sprintf("t%03d_n%02d",timeStep,nodeID))
  #Append node parameters for time loop_t + 1 to node-level data frame
  model_state$node_params <- rbind(nodeParams, futureNodeParams)
  model_state$current_step=loop_t+1
  
  cat("fin\n")
  return(model_state)
}

simple_iteration_handler <- list()
simple_iteration_handler$next_modelState <- next_modelState_simple
simple_iteration_handler$update_coalitions <- update_sinkCoalitions_simple
simple_iteration_handler$place_sink_orders  <- function(model_state, failure_Response){
  cat("Stub for placing orders at step:",model_state$current_step,"\n" )
  current_sink_buyers(model_state) %>% as.char.coal_set() -> sN
  
  #Randomize order in which sink nodes/coalitions place orders
  sN <- sN %>%sample()
  for (i in sN) { 
    #Calculate net demand for node i on this time step (current gross demand - current supply)
    demand <- model_state$node_demand(model_state,i)
    
    if(sus_node(model_state$node_params, timeStep = model_state$current_step, nodeID = as.coal_set.char(i))$failureResponse == 1 && 
       failure_Response == "inflateOrder") {
      demand <- round(demand * 1.2, digits = 0)
    }
    
    #Determine possible suppliers, drawn from both source and transfer nodes
    suppliers <- locate_suppliers(model_state,i)
    if (length(suppliers) > 0) {
      cat("Stub for placing orders for:",i,"\n" )
      print(suppliers)
      
      #Update transfer nodes' supplier information (e.g. current stock levels)
      sus_node(model_state$node_params,timeStep = model_state$current_step )$potentialSupply <-
        update_potential_supplies(adjMatrix = model_state$scMatrix,
                                  transferNodeList = model_state$transfer_nodes,
                                  sourceNodeList = model_state$source_nodes,
                                  st = model_state$oneSkel,
                                  timeStep = model_state$current_step,
                                  nodeParams = model_state$node_params)
      
      #Create list of potential suppliers and order them according to price per unit
      supplierParams <- sus_coalition(model_state, suppliers)
      if(sus_node(model_state$node_params, timeStep = model_state$current_step, nodeID = as.coal_set.char(i))$purchaseStrategy == "random") {
        if (dim(supplierParams)[1] > 1) {
          rows <- sample(nrow(supplierParams))
          supplierParams <- supplierParams[rows, ] 
        }
      } else{
        supplierParams <- supplierParams[order(supplierParams$pricePerUnit),]
      }
      
      #After ordering suppliers, if "distributorPriority" is active, place distributors first
      if(sus_node(model_state$node_params, timeStep = model_state$current_step, nodeID = as.coal_set.char(i))$failureResponse == 1 && 
         failure_Response == "distributorPriority") {
        supplierParams <- supplierParams[order(supplierParams$nodeType, decreasing = TRUE),]
      }
      
      #Focal node places orders with suppliers, prioritizing lowest prices per unit
      #Suppliers that do not hold sufficient stock to meet an order can place orders in turn with upstream suppliers
      print(demand)
      print(supplierParams)
      if(is.null(model_state$place_orders)) 
        stop("Not implemented: model_state$place_orders")
      model_state <- model_state$place_orders(model_state,
                                              entityID=i,
                                              total_volume=demand,
                                              supplierParams = supplierParams)
    }
  }
  return(model_state);
}

simple_iteration_handler$place_buffer_orders  <-
  function(model_state, memory, desiredBufferPercentage, inventoryType){
    
    cat("Stub for buffer orders at step:",model_state$current_step,"\n" )
    if(inventoryType == "Durable") {
      
      #If memory = 0 (default), nodes' desired buffer supply is equal to a percentage of the 
      #order volume they received this turn (desiredBufferPercentage)
      #Otherwise, target buffer is calculated as (maximum orders received - average orders received) over a number of recent time steps equal to memory
      
      bufferSupply <- estimate_target_buffer(memory = memory, 
                                             timeStep = model_state$current_step, 
                                             adjMatrix = model_state$scMatrix, 
                                             bufferPercent = desiredBufferPercentage, 
                                             nodeParams = model_state$node_params)
      sus_node(model_state$node_params,timeStep = model_state$current_step)$desiredBuffer <- bufferSupply
      
      #Randomly determine the order in which transfer nodes place orders to maintain their desired buffer amount
      check_nodes <- model_state$transfer_nodes %>% as.list() %>% 
        as.char.coal_set() %>% sample()
      
      for (i in check_nodes) {
        
        # Calculate the total order volume needed to maintain desired buffer supply and fulfill all received order volumes
        bufferOrder <- 
          with(sus_coalition(model_state, i),
               desiredBuffer + OrderVolReceived - currentSupply)
        
        #Determine possible suppliers, drawn only from source nodes
        suppliers <- locate_suppliers(model_state, i, include_transfers = F)
        
        if (length(suppliers) > 0) {
          supplierParams <- sus_coalition(model_state, suppliers)
          supplierParams <- supplierParams[order(supplierParams$pricePerUnit),]
          
          # Focal node places orders with suppliers, prioritizing lowest prices per unit
          # Suppliers that do not hold sufficient stock to meet an order can place orders in turn with upstream suppliers
          print(model_state$node_demand(model_state,i))
          print(supplierParams)
          if(is.null(model_state$place_orders)) 
            stop("Not implemented: model_state$place_orders")
          model_state <- 
            model_state$place_orders(model_state, 
                                     entityID=i, 
                                     total_volume=bufferOrder,
                                     supplierParams = supplierParams)
        }
      }
      
      #After having received additional orders from transfer nodes, source nodes now calculate their desired buffer supply
      #The process here is the same as that used by the transfer nodes above
      bufferSupply <- estimate_target_buffer(memory = memory, 
                                             timeStep = model_state$current_step, 
                                             adjMatrix = model_state$scMatrix, 
                                             bufferPercent = desiredBufferPercentage, 
                                             nodeParams = model_state$node_params)
      sus_node(model_state$node_params,timeStep = model_state$current_step)$desiredBuffer <- bufferSupply
    }
    
    #If inventory type is Perishable, there is no need to build a buffer supply, as unused stock is lost at end of turn
    if(inventoryType == "Perishable") {
      sus_node(model_state$node_params, timeStep = model_state$current_step)$desiredBuffer <- 0
    }
    
    return(model_state) 
  }

simple_iteration_handler$distribute_stock  <- function(model_state, inventoryType){
  str(model_state)
  
  #Nodes now begin distributing stock according to the orders received  with source nodes distributing first
  src_nodes <- current_sellers(model_state, include_transfers = F) %>% names()
  mid_nodes <- current_sellers(model_state, include_transfers = T) %>% names() %>%
    purrr::keep(~ ! .x %in% src_nodes)
  
  for (i in c(sample(src_nodes),sample(mid_nodes))) {
    
    #A supplier first aggregates all order volumes received from the same customer this turn at the same price!
    recordAggregated <-
      model_state$orderRecord %>%
      filter(Supplier == i, TimeStep == model_state$current_step) %>%
      group_by(Supplier,TimeStep,Customer,pricePerUnit) %>%
      summarise(Volume=sum(Volume)) %>% ungroup()
    
    #Check if there are orders to be fulfilled by the supplier
    if (nrow(recordAggregated) > 0) {
      model_state <- model_state$distribute_stock(model_state, supplierID=i, records=recordAggregated)
    }
  }
  
  #Next coalitions distribute received supplies across coalition members
  if(length(current_sink_coalitions(model_state)) > 0){
    # treat coalition members as making an order for the actual demand at the price required to
    
    sink_coalitions <- current_sink_coalitions(model_state) 
    for(coal_i in sink_coalitions){
      if(inventoryType == "Perishable"){
        sus_coalition(model_state, 
                      entityID=as.char.coal_set(coal_i))$OrderVolReceived <-
          with(sus_node(model_state$node_params,
                        timeStep = model_state$current_step, 
                        nodeID = coal_i),
               sum(currentGrossDemand-currentSupply))
        
        sus_node(model_state$node_params,
                 timeStep = model_state$current_step, 
                 nodeID = coal_i)$OrderVolPlaced <- 
          with(sus_node(model_state$node_params,
                        timeStep = model_state$current_step, 
                        nodeID = coal_i),
               currentGrossDemand-currentSupply)
      } else{
        if(inventoryType == "Durable"){
          sus_coalition(model_state, 
                        entityID=as.char.coal_set(coal_i))$OrderVolReceived <-
            with(sus_node(model_state$node_params,
                          timeStep = model_state$current_step, 
                          nodeID = coal_i),
                 sum(round((currentGrossDemand-currentSupply)*1.1)))
          
          sus_node(model_state$node_params,
                   timeStep = model_state$current_step, 
                   nodeID = coal_i)$OrderVolPlaced <- 
            with(sus_node(model_state$node_params,
                          timeStep = model_state$current_step, 
                          nodeID = coal_i),
                 round((currentGrossDemand-currentSupply)*1.1))
        }
      }
      
      unit_price <- 
        with(sus_coalition(model_state, entityID = as.char.coal_set(coal_i)),
             PaymentMade / currentSupply)
      sus_coalition(model_state, entityID = as.char.coal_set(coal_i))$pricePerUnit <- unit_price
      
      model_state$orderRecord <-
        model_state$orderRecord %>%
        add_row(
          TimeStep = model_state$current_step,
          Customer = paste0("COAL_",sus_node(model_state$node_params,
                                             timeStep = model_state$current_step, 
                                             nodeID = coal_i)$nodeID),
          Supplier = as.char.coal_set(coal_i),
          Volume = sus_node(model_state$node_params,
                            timeStep = model_state$current_step, 
                            nodeID = coal_i)$OrderVolPlaced,
          pricePerUnit=unit_price)
      #Assumes coalitions refuse delivery of units members will not want
      if(with(sus_coalition(model_state, 
                            entityID=as.char.coal_set(coal_i)),OrderVolReceived < currentSupply) ){
        print(sus_coalition(model_state, entityID=as.char.coal_set(coal_i)))
        print(sus_node(model_state$node_params,
                       timeStep = model_state$current_step, 
                       nodeID = coal_i))
        stop("Model error: Coalition accepted more units than members will accept!")
      }
    }
    with(model_state,
         print(coalition_params[coalition_params$timeStep==current_step,]) )
    sink_coalitions <- as.char.coal_set(sink_coalitions)
    print(model_state$orderRecord %>% 
            filter(TimeStep==model_state$current_step,
                   Supplier %in% sink_coalitions))
    
    for (i in sink_coalitions) {
      
      #Aggregate all order volumes received from the same customer this turn
      recordAggregated <- model_state$orderRecord %>% 
        filter(TimeStep==model_state$current_step,
               Supplier %in% i)
      
      #Check if there are orders to be filled
      if (nrow(recordAggregated) > 0) {
        
        if(is.null(model_state$distribute_stock_coalition))
          stop("Not implemented: model_state$distribute_stock_coalition")
        model_state <- 
          model_state$distribute_stock_coalition(model_state, supplierID=i, records=recordAggregated)
        if(with(sus_coalition(model_state,i), round(PaymentReceived) != round(PaymentMade)))
          warning("Funds in != Funds out for ", i," at ",model_state$current_step)
        if(with(sus_coalition(model_state,i),currentSupply<0))
          stop("Non-existant units delivered by ", i," at ",model_state$current_step)
        if(with(sus_coalition(model_state,i),currentSupply>0))
          warning("Extra units trashed by ", i," at ",model_state$current_step)
      }
    }
  }
  return(model_state);
}

simple_iteration_handler$consume_supplies <- function(model_state){
  #All sink nodes consume until supplies are depleted or they've reached their desired amount
  sus_node(model_state$node_params, timeStep = model_state$current_step) <- 
    consume_supplies_multi(sN = model_state$sink_nodes, tS = model_state$current_step, nP = model_state$node_params)
  return(model_state)
}

#################################################
#Function for obtaining summary stability metrics
#################################################

getStabilityMetricsV1 <- function(nP, simID=0, inventoryType){
  if(! "simID" %in% names(nP))
    nP$simID=simID
  
  orderStatusSummary <-
    nP %>% group_by(simID, timeStep) %>%
    filter(OrderVolPlaced > 0) %>%
    mutate(orderFilled = (OrderVolPlaced == DeliveryVolReceived)) %>%
    summarise(PropOrdersFilled=mean(orderFilled),
              PropOrdersPartFilled=mean(DeliveryVolReceived>0 ),
              NetworkFillRate=sum(DeliveryVolReceived) / sum(OrderVolPlaced))
  
  if(inventoryType == "Durable"){
    sourceStockSummary <-
      nP %>% group_by(simID, timeStep) %>%
      filter(nodeType %in% c("Source")) %>%
      summarise(meanExcessStockSource=mean(currentSupply - desiredBuffer),
                percentExcessStockSource=mean((currentSupply - desiredBuffer)/(currentSupply - desiredBuffer + DeliveryVolMade), na.rm = TRUE),
                sdExcessStockSource=sd(currentSupply - desiredBuffer), 
                meanDesiredBufferSource=mean(desiredBuffer))
    transferStockSummary <-
      nP %>% group_by(simID, timeStep) %>%
      filter(nodeType %in% c("Transfer" )) %>%
      summarise(meanExcessStockTransfer=mean(currentSupply - desiredBuffer),
                percentExcessStockTransfer=mean((currentSupply - desiredBuffer)/(currentSupply - desiredBuffer + DeliveryVolMade), na.rm = TRUE),
                sdExcessStockTransfer=sd(currentSupply - desiredBuffer), 
                meanDesiredBufferTransfer=mean(desiredBuffer))
  } else{
    sourceStockSummary <-
      nP %>% group_by(simID, timeStep) %>%
      filter(nodeType %in% c("Source")) %>%
      summarise(meanExcessStockSource=mean(supplyLost),
                percentExcessStockSource=mean(supplyLost/(supplyLost + DeliveryVolMade), na.rm = TRUE),
                sdExcessStockSource=sd(supplyLost), 
                meanDesiredBufferSource=mean(desiredBuffer))
    transferStockSummary <-
      nP %>% group_by(simID, timeStep) %>%
      filter(nodeType %in% c("Transfer" )) %>%
      summarise(meanExcessStockTransfer=mean(supplyLost),
                percentExcessStockTransfer=mean(supplyLost/(supplyLost + DeliveryVolMade), na.rm = TRUE),
                sdExcessStockTransfer=sd(supplyLost), 
                meanDesiredBufferTransfer=mean(desiredBuffer))
  }
  
  demandSummary <-
    nP %>% group_by(simID, timeStep) %>%
    filter(nodeType %in% c("Sink" )) %>%
    mutate(sinkOrderFilled = (OrderVolPlaced == DeliveryVolReceived)) %>%
    summarise(
      #meanExcessDemand=mean(ExcessDemand),
              percentShortfall=mean(ExcessDemand / currentGrossDemand, na.rm = TRUE),
              sdPercentShortfall=sd(ExcessDemand / currentGrossDemand, na.rm = TRUE),
              #sinkPropOrdersFilled = mean(sinkOrderFilled), 
              #sinkPropOrdersPartFilled = mean(DeliveryVolReceived > 0),
              meanFillRate = mean(DeliveryVolReceived/OrderVolPlaced, na.rm = TRUE),
              sdFillRate = sd(DeliveryVolReceived/OrderVolPlaced, na.rm = TRUE),
              #sdExcessDemand=sd(ExcessDemand)
              )
  
  #Summary statistics specifically for coalition members on each time step
  demandSummary_CoalMembs <-
    nP %>% group_by(simID, timeStep) %>%
    filter(nodeType %in% c("Sink")) %>%
    filter(coalitionTime > 0) %>%
    mutate(sinkOrderFilled_CM = (OrderVolPlaced == DeliveryVolReceived)) %>%
    summarise(meanExcessDemand_CM = mean(ExcessDemand),
              percentShortfall_CM = mean(ExcessDemand / currentGrossDemand, na.rm = TRUE), 
              sinkPropOrdersFilled_CM = mean(sinkOrderFilled_CM),
              sinkPropOrdersPartFilled_CM = mean(DeliveryVolReceived > 0), 
              meanFillRate_CM = mean(DeliveryVolReceived/OrderVolPlaced, na.rm = TRUE),
              sdExcessDemand_CM = sd(ExcessDemand))
  
  #Summary statistics specifically for lone firms on each time step
  demandSummary_Solo <-
    nP %>% group_by(simID, timeStep) %>%
    filter(nodeType %in% c("Sink")) %>%
    filter(coalitionTime == 0) %>%
    mutate(sinkOrderFilled_S = (OrderVolPlaced == DeliveryVolReceived)) %>%
    summarise(meanExcessDemand_S = mean(ExcessDemand),
              percentShortfall_S = mean(ExcessDemand / currentGrossDemand, na.rm = TRUE), 
              sinkPropOrdersFilled_S = mean(sinkOrderFilled_S),
              sinkPropOrdersPartFilled_S = mean(DeliveryVolReceived > 0), 
              meanFillRate_S = mean(DeliveryVolReceived/OrderVolPlaced, na.rm = TRUE),
              sdExcessDemand_S = sd(ExcessDemand))
  
  budgetSummary <-
    nP %>% group_by(simID, timeStep) %>% filter(budgetPerStep > 0) %>%
    summarise(percentLiquidRemaining=mean(liquidFunds / initialLiquidFunds))
  
  priceSummary <-
    nP %>% group_by(simID, timeStep) %>% filter(budgetPerStep > 0) %>% filter((DeliveryVolReceived - DeliveryVolMade) > 0) %>% 
    summarise(meanPricePerUnit=mean((PaymentMade - PaymentReceived)/(DeliveryVolReceived - DeliveryVolMade)))
  
  stratSummary <- 
    nP %>% group_by(simID, timeStep) %>% filter(nodeType %in% c("Sink")) %>%
    summarise(percentPricePriority = (sum(purchaseStrategy == "pricePriority")/sum(purchaseStrategy == "pricePriority" | 
                                                                                     purchaseStrategy == "betHedging" | 
                                                                                     purchaseStrategy == "random")))
  
  #orderStatusSummary %>% 
  demandSummary %>%
    #full_join(sourceStockSummary, by=c("simID", "timeStep")) %>%
    #full_join(transferStockSummary, by=c("simID", "timeStep")) %>%
    #full_join(demandSummary, by=c("simID", "timeStep")) %>%
    #full_join(demandSummary_CoalMembs, by=c("simID", "timeStep")) %>%
    #full_join(demandSummary_Solo, by=c("simID", "timeStep")) %>%
    full_join(budgetSummary, by=c("simID", "timeStep")) %>%
    full_join(priceSummary, by=c("simID", "timeStep")) %>%
    full_join(stratSummary, by=c("simID", "timeStep")) %>%
    ungroup() %>% rename(TimeStep=timeStep)
}

#Should adjust this so that it takes the average of excess demand / current gross rather than just the most recent time step
adjust_purchase_strategy <- function(model_state) {
  model_state$node_params %>%
    filter(timeStep == (model_state$current_step - 1), currentGrossDemand > 0, 
           ExcessDemand / currentGrossDemand >= 0.2, 
           timeSinceSwitch >= 5) %>% pull("nodeID") -> res1
  for(node_i in res1){
    if(sus_node(model_state$node_params, 
                timeStep = model_state$current_step, 
                nodeID = node_i)$purchaseStrategy == "pricePriority") {
      sus_node(model_state$node_params, timeStep = model_state$current_step, 
               nodeID = node_i)$purchaseStrategy <- "betHedging"
      sus_node(model_state$node_params, timeStep = model_state$current_step, 
               nodeID = node_i)$timeSinceSwitch <- 0
    } else {
      if(sus_node(model_state$node_params, 
                  timeStep = model_state$current_step, 
                  nodeID = node_i)$purchaseStrategy == "betHedging") {
        sus_node(model_state$node_params, timeStep = model_state$current_step, 
                 nodeID = node_i)$purchaseStrategy <- "pricePriority"
        sus_node(model_state$node_params, timeStep = model_state$current_step, 
                 nodeID = node_i)$timeSinceSwitch <- 0
      }
    }
  }
  return(model_state)
}

add_supplier <- function(model_state) {
  model_state$node_params %>% 
    filter(timeStep == (model_state$current_step), currentGrossDemand > 0, 
           failureResponse == 1) %>% pull("nodeID") -> res1
  for(node_i in res1){
    if(runif(1, min = 0, max = 1) <= 0.05){
      
      #Determine possible suppliers, drawn from both source and transfer nodes
      currentSuppliers <- cofaceList(st = model_state$oneSkel, nodeID = node_i)
      currentSuppliers <- currentSuppliers[currentSuppliers %in% c(model_state$source_nodes, model_state$transfer_nodes) & 
                                             currentSuppliers != node_i]
      
      potentSuppliers <- c()
      
      for(c in currentSuppliers) {
        potentSuppliersTemp <- cofaceList(st = model_state$oneSkel_Initial, nodeID = c)
        potentSuppliersTemp <- potentSuppliersTemp[potentSuppliersTemp %in% c(model_state$source_nodes, model_state$transfer_nodes)]
        potentSuppliers <- c(potentSuppliers, potentSuppliersTemp)
      }
      
      potentSuppliers <- unique(potentSuppliers[!potentSuppliers == node_i])
      potentSuppliers <- potentSuppliers[!(potentSuppliers %in% currentSuppliers)]
      
      if(length(potentSuppliers) > 0) {
        if(length(potentSuppliers) == 1) {
          model_state$oneSkel %>% insert(list(c(potentSuppliers[1], node_i)))
        } else{
          if(length(potentSuppliers) > 1){
            newSupplier <- sample(potentSuppliers, 1, replace = FALSE)
            model_state$oneSkel %>% insert(list(c(newSupplier, node_i)))
            supplierAdjustment <- sum(c(model_state$source_nodes, model_state$transfer_nodes, model_state$sink_nodes) < newSupplier) + 1
            customerAdjustment <- sum(c(model_state$source_nodes, model_state$transfer_nodes, model_state$sink_nodes) < node_i) + 1
            model_state$scMatrix[supplierAdjustment, customerAdjustment] <- 1
          }
        }
      }
    }
  }
  return(model_state)
}