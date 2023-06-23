
###########################################################################
#Load libraries
library(tibble)
library(simplextree)
library(dplyr)
library(truncnorm)

#Set working directory
setwd()

#Load custom functions needed for simulation
source("functions_Deliverable2.R")
source("coalitionFormationFunctions_Deliverable3_final.R")

#Create folders in which to store simulation results
#sim_details includes each node's parameters at each time step for each simulation
#sim-summaries includes summary statistics for each time step for each simulation run
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details="Sim-details_purchaseStrat-mixed_responseStrat-addSupplier"
sim_summaries="Sim-summaries_purchaseStrat-mixed_responseStrat-addSupplier"
if(!file.exists(sim_details)) dir.create(sim_details)
if(!file.exists(sim_summaries)) dir.create(sim_summaries)
###########################################################################

#Set random seed to reproduce identical results when running the simulation
set.seed(347569)

#Set when the shock occurs and the final time step of the simulation
t_shock= 50
t_final = 100

#Set number of simulations to run
numberSims = 50
seedVect <- round(rnorm(numberSims, mean = 500000, sd = 200000))

#Set number of each node type within the initial network
n_source_nodes = 20
n_transfer_nodes = 30
n_sink_nodes = 50

#Specify disruption type: "Supply", "Transport", or "Demand"
disruptionType = "Supply"

#When disruption is to supply, disruption can either be "nodeRemoval" or "supplyReduction"
supplyDisruption = "supplyReduction"

#Specify type of good: "Durable" or "Perishable"
inventoryType = "Durable"

#Specify purchase strategy that sink nodes follow: "pricePriority", "betHedging", "random", "mixed"
purchase_Strategy = "mixed"

#Can be "none", "inflateOrder", "distributorPriority", or "addSupplier"
failure_Response = "addSupplier"

#Create empty data frame for simulation summary data
simDataAll <- data.frame(
  "simID" = integer(),
  "TimeStep" = integer(),
  "percentShortfall" = numeric(),
  "meanFillRate" = numeric(),
  "percentLiquidRemaining" = numeric(),
  "meanPricePerUnit" = numeric()
)

#Run requested number of simulations (numberSims)
for (s in 1:numberSims) {
  set.seed(seedVect[s])
  
  #Create the initial model state and post-disruption modifications (e.g., to the network)
  #Can also set here: how many edges nodes tend to form on average to other node types
  #the disruption magnitude (e.g., 0.3 for transfer disruption means 30% of the transfer nodes are removed)
  #the mean price per unit for source nodes, how much transfer nodes mark these prices up,
  #how much inventory source nodes tend to be able to generate each turn, and
  #sink nodes' mean gross demand
  shock_simulation <- 
    shockSimulation_v1.3(n_source = n_source_nodes,
                         n_transfer = n_transfer_nodes,
                         n_sink = n_sink_nodes, 
                         edge.src_out =2,
                         edge.transfer=c(2,1,2),
                         edge.sink_in=2,
                         disruptionType = disruptionType,
                         disruptionMagnitude = 0.3,
                         supplyDisruption = supplyDisruption,
                         meanSourcePrice = 15,
                         meanTransferMarkup = 0.3,
                         meanReplenishment = 150,
                         meanGD = 50, 
                         purchase_Strategy = purchase_Strategy,
                         inventoryType = inventoryType,
                         t_final = t_final)

  #Create the model state object to be used in the simulation
  #If applicable, here is where the rules for joining and leaving coalitions and the costs of coalition formation/maintenence are specified
  sinkCoalition_state <-  sinkCoalitionModelState(model_state = shock_simulation$model_state, 
                                                  coal_rule.add = D2_add_rule, 
                                                  coal_rule.leave = D2_leave_rule,
                                                  coal_cost.init = 0, 
                                                  coal_cost.min = 0, 
                                                  coal_cost.delta = 0)
  print(sinkCoalition_state$node_params %>% head())
  str(sinkCoalition_state)
  
  #Carry out the pre-disruption supply chain simulation
  #Can set desired buffer percentage here for source and transfer nodes
  pre_shock_output <- 
    simulationLoop_v1.3(t_start = 1, 
                        t_end = t_shock,
                        model_state = sinkCoalition_state,
                        desiredBufferPercentage = 0.15,
                        iteration_handler = simple_iteration_handler, 
                        inventoryType = inventoryType, 
                        purchase_Strategy = purchase_Strategy, 
                        failure_Response = failure_Response)
  
  if(!is.null(pre_shock_output$coalition_params)) {
    print(pre_shock_output$coalition_params %>% filter(PaymentMade!=PaymentReceived))
    print(pre_shock_output$coalition_params %>% 
            group_by(entityID) %>% slice(which.max(coalitionTime)) %>% ungroup() %>% as.data.frame())
  }
  
  if(failure_Response == "addSupplier") {
    updatedShockParams <- update_shock_params(model_state = pre_shock_output, 
                                              shock_params = shock_simulation$shock_params, 
                                              edge.sink_in = 2)
    
    #Carry out the specified disruption
    post_shock_state <- update_modelState_wShock(model_state = pre_shock_output,
                                                 shock_params = updatedShockParams,
                                                 inventoryType = inventoryType, 
                                                 disruptionMagnitude = 0.3, 
                                                 supplyDisruption)
    
  } else {
  
  #Carry out the specified disruption
  post_shock_state <- update_modelState_wShock(model_state = pre_shock_output,
                                               shock_params = shock_simulation$shock_params,
                                               inventoryType = inventoryType,
                                               disruptionMagnitude = 0.3,
                                               supplyDisruption)
  }
  
  #Carry out the post-disruption simulation
  post_shock_output <- 
    simulationLoop_v1.3(t_start = t_shock + 1,
                        t_end = t_final,
                        model_state = post_shock_state,
                        desiredBufferPercentage = 0.15,
                        iteration_handler = simple_iteration_handler, 
                        inventoryType = inventoryType, 
                        purchase_Strategy = purchase_Strategy, 
                        failure_Response = failure_Response)
  
  #Extract node-level parameters and order records for complete simulation
  nodeParams  = post_shock_output$node_params
  orderRecord = post_shock_output$orderRecord
  
  #Calculate stability metrics for simulation
  #The rows used to initialize the data (i.e., at time step 0) are first removed
  simData_old <- getStabilityMetricsV1(nodeParams %>% filter(timeStep != 0), s, inventoryType = inventoryType)
  
  cat("FinalNode-level Parameters:\n")
  print(nodeParams[nodeParams$timeStep == t_final, ])
  if(!is.null(post_shock_output$coalition_params)) {
    print(post_shock_output$coalition_params %>% filter(PaymentMade!=PaymentReceived))
    print(post_shock_output$coalition_params %>% 
            group_by(entityID) %>% slice(which.max(coalitionTime)) %>% ungroup() %>% as.data.frame())
  }

  #Add current simulation data to data frame storing the results from the previously run simulations
  simDataAll <- rbind(simDataAll, simData_old)
  
  #Writes a csv with the node-level parameters for each time step in the simulation
  write.csv(nodeParams, file=file.path(sim_details,sprintf( "nPv2_%s_%02i.csv",run_ID,s)))
  
  #Writes a csv including information on the coalitions active on each time step in the simulation (if applicable)
  if(!is.null(post_shock_output$coalition_params)) {
    if("nodeID" %in% names(post_shock_output$coalition_params)){
      kID= which(names(post_shock_output$coalition_params) %in% "nodeID")
    } else { kID=integer(0)}
    write.csv(post_shock_output$coalition_params[,-kID], 
              file=file.path(sim_details,sprintf( "cP_%s_%02i.csv",run_ID,s)))
  }
  
  #Writes a csv including all the purchase orders made during the simulation
  write.csv(orderRecord, file=file.path(sim_details,sprintf( "oR_%s_%02i.csv",run_ID,s)))
}

#Writes a csv including the summary data for each simulation that was run
simDataAll %>% mutate_all(round,digits=4) %>%
  write.csv(file=file.path(sim_summaries,sprintf("simV1_%s.csv",run_ID)),
            row.names = F)

