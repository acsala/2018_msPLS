get_parallel_cv <- function(X, 
                                Y, 
                                lambdas,
                                non_zeros, 
                                label,
                                penalization,
                                max_iterations,
                                tolerance){
  #non-parallel x-fold Cross validation for each alpha, for example alpha = 0, 0.1, ... , 0.9, 1.0
  #
  #input:     X             n*p matrix          - independent variables 
  #           Y             n*q matrix          - dependent variables 
  #           lambda        integer             - factor for Ridge penalty 
  #           labels         vector with n_subset unique elements with length of n - for crossvalidation
  # 
  #output:    abs_cor       vector of integers  - returns the summed abs correlation
  #           stime         int/time            - return algorithms running time
  
  nr_subsets      <-    length(unique(label))
  abs_cors        <-    c()
  iterations_m    <-    c()
  
  # TESTER VALUES
  lambdas <- c(0.1,1)
  non_zeros <- c(20,40)
  tolerance <- 0.001
  max_iterations <- 100
  penalization <- "ust"
  data_sets <- Data
  nr_subsets <- 10
  specs <- s_satpls$model$specs
  
  #get working clusters initialized, otherwise parallel will use only one cpu***************
  nr_cores      <-    detectCores()
  cl            <-    makeCluster(nr_cores)
  registerDoParallel(cl)
  getDoParWorkers()
  print("Number of detected cores:")
  print(nr_cores)
  
  
  #measure time
  stime <- system.time({
    
    #Get the alphas with parallel computing
    x <- foreach(l = lambdas, .combine = c,
                 .export= c("enet",
                            "splspm"),
                 .packages= c()
                 ) %dopar% {
    
    # for (l in 1:length(lambdas)){
      
      sub_abs_cor       <- c()
      sub_results       <- c()
      sub_iterations_m  <- c()
      
      for (nz in 1:length(non_zeros)){
        
        for (i in 1:nr_subsets){
          
          X.train   <- X[label!=i,]
          X.test    <- X[label==i,]
          
          dim(X.train);dim(X.test)
          
          Y.train   <- Y[label!=i,]
          Y.test    <- Y[label==i,]
          
          dim(Y.train);dim(Y.test)
          
          
          data_sets <- cbind(X.train,Y.train)
          
          # Only 2 DATASETS####
          EXPL_X = c(0,0)
          RESP_Y = c(1,0)
          path_matrix = rbind(EXPL_X, RESP_Y)
          
          # blocks of outer model
          blocks = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2])
          
          modes = c("B","A")
          
          
          sub_results[[i]] <- splspm(data_sets,
                                     path_matrix,
                                     blocks,
                                     modes,
                                     scheme="path",
                                     penalization = penalization,
                                     nonzero = non_zeros[nz],
                                     lambda = lambdas[l],
                                     maxiter = max_iterations,
                                     tol = tolerance,
                                     warning_non_convergence = FALSE)
          
          # sub_results[[i]] <- splspm(Data = data_sets, 
          #                    path_matrix = path_matrix, 
          #                    blocks = blocks, 
          #                    modes = modes, 
          #                    scaling = NULL,
          #                    scheme = "path", 
          #                    scaled = TRUE, 
          #                    tol = tolerance, 
          #                    maxiter = max_iterations,
          #                    plscomp = NULL, 
          #                    boot.val = FALSE, 
          #                    br = NULL,
          #                    dataset = TRUE,
          #                    penalization = penalization,
          #                    nonzero = non_zeros[nz],
          #                    lambda = lambdas[l],
          #                    cross_validate = FALSE,
          #                    nr_subsets = 10,
          #                    warning_non_convergence = FALSE)
          
          # sub_results[[i]] = get_weights_p(X = data_sets, 
          #                                  path_matrix = path_matrix, 
          #                                  blocks = blocks, 
          #                                  specs = specs, 
          #                                  penalization = penalization,
          #                                  nonzero = non_zeros[nz],
          #                                  lambda = lambdas[l],
          #                                  warning_non_convergence = FALSE)
          
          ALPHA <- t(X.train)%*% sub_results[[i]]$scores[,1]%*%
            solve(t( sub_results[[i]]$scores[,1])%*%
                    sub_results[[i]]$scores[,1])
          
          XI.test = scale(X.test) %*% ALPHA
          
          #devide with dim(Y.train)[2]
          sub_abs_cor[[i]] <- sum((abs(cor(XI.test,Y.test))))/dim(Y.train)[2]
          
          #sub_iterations_m[[i]] <- sub_results[[i]]$Nr_iterations
          
          
        }#end of subset for loop
        
        abs_cors        <- cbind(abs_cors, sub_abs_cor)
        #iterations_m    <- cbind(iterations_m, sub_iterations_m)
        
      }#end of non_zeros loop
      
      abs_cors
      
    }#end of lambda for loop 
  })[3]#end of measure time
  

  
  #stop cluster
  stopCluster(cl)
  
  
  
}




# #Get the alphas with parallel computing
# x <- foreach(l = lambdas, .combine = c,
#              .export= c(
#                         # "enet",
#                         "splspm",
#                         # "tester",
#                         # "turner",
#                         # "diagram",
#                         # "shape",
#                         # "amap",
#                         # "elasticnet",
#                         # "mvtnorm",
#                         # "auxiliar",
#                         # "check_arguments",
#                         #"check_specifications",
#                         "generate_data",
#                         "get_alpha",
#                         "get_ave",
#                         "get_boots",
#                         "get_cross_validated_penalty_parameters",
#                         "get_dummies",
#                         "get_effects",
#                         # "get_enet",
#                         "get_generals",
#                         "get_gof",
#                         "get_GQI",
#                         "get_inner_summary",
#                         "get_locals_test",
#                         "get_manifests",
#                         "get_metric",
#                         "get_nom_scale",
#                         "get_non_parallel_cv",
#                         "get_num_scale",
#                         "get_ord_scale",
#                         "get_parallel_cv",
#                         "get_path_scheme",
#                         "get_paths",
#                         "get_PLSR",
#                         "get_PLSR_NA",
#                         "get_plsr1",
#                         "get_PLSRdoubleQ",
#                         "get_rank",
#                         "get_rho",
#                         # "get_scaled_data",
#                         # "get_scores",
#                         # "get_split_sets",
#                         # "get_treated_data",
#                         # "get_unidim",
#                         # "get_ust",
#                         # "get_weights_nonmetric",
#                         "get_weights_p",
#                         "innerplot",
#                         "it.reb",
#                         "local.models",
#                         "outerplot",
#                         # "plot.plspm",
#                         # "plspm.fit",
#                         # "plspm.groups",
#                         # "plspm-package",
#                         # "print.rebus",
#                         # "quantiplot",
#                         # "rebus.pls",
#                         # "rebus.test",
#                         # "res.clus",
#                         # "rescale",
#                         # "russett-data",
#                         # "splspm",
#                         # "sPLSPM-package",
#                         # "summary_plspm",
#                         # "test_dataset",
#                         # "test_factors",
#                         # "test_manifest_scaling",
#                         "test_null_weights"
#                         # "unidimensionality"
#                         )) %dopar% {