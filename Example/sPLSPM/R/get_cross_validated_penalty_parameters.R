#prepare for cross validation

get_cross_validated_penalty_parameters <- function(data_sets,
                                                   path_matrix,
                                                   blocks,
                                                   specs,
                                                   penalization,
                                                   ridge_penalty,
                                                   nonzero,
                                                   tolerance,
                                                   max_iterations,
                                                   nr_subsets = nr_subsets
                                                   ){

  best_ridge_penalties  <- c()
  best_nr_of_nonzeros   <- c()
  cv_results_list       <- c()
  aggr_cv_results_list  <- c()
  s_path_matrix <- path_matrix

  s_path_matrix[upper.tri(s_path_matrix)] <- rep(0,length(s_path_matrix[upper.tri(s_path_matrix)]))

  # TESTER VALUES
  # ridge_penalty <- lambda
  # nonzero <- nonzero
  # tolerance <-tol
  # max_iterations <- maxiter
  # penalization <- "enet"
  # data_sets <- Data
  # nr_subsets <- 10
  # specs <- specs
  # specs$modes <- c("A","A","A","B")

  #*********************

  #crossvaldate datasets if they are connected
  for (nr_sets in 1:length(specs$modes)){

    #nr of connections current dataset has
    current_set_path <- path_matrix[,nr_sets]

    aggregate_cv_results <- c()

    #check if dataset itself mode B
    if (specs$modes[nr_sets] == "B"){

      #if dataset is itself mode B and first dataset
      if (sum(current_set_path)>=1 || nr_sets == length(specs$modes)){

        cv_results <- c()

        for (set_connected_checker in 1:length(current_set_path)){

          # sum(path_matrix[1:set_connected_checker,nr_sets]) = 1 # check if it was already used for finding the lambdas
          if(current_set_path[set_connected_checker]==1  ||
             (nr_sets == length(specs$modes) &&
              set_connected_checker == length(current_set_path))){

            # if (sum(path_matrix[1:set_connected_checker,nr_sets])>1){
            #     print("dataset itself mode B and connected (nr set / connected to)")
            #     print(nr_sets)
            #     print(set_connected_checker)
            #     print(specs$modes[c(nr_sets,set_connected_checker)])
            #     print("dataset connected through others check - double work")
            # }


            if (nr_sets>1 &&
                path_matrix[set_connected_checker,nr_sets]==1 &&
                path_matrix[nr_sets,set_connected_checker] == 1 &&
                nr_sets > set_connected_checker){

                print("dataset itself mode B and connected (nr set / connected to)")
                print(nr_sets)
                print(set_connected_checker)
                print(specs$modes[c(nr_sets,set_connected_checker)])
                print("no analysis since it's been already done - double work")

            }else{


                #save current dataset in varaiable
                matrix_X <- data_sets[,blocks[[nr_sets]]]

                dim(matrix_X)

                matrix_Y <- data_sets[,blocks[[set_connected_checker]]]

                dim(matrix_Y)

                print("dataset itself mode B and connected (nr set / connected to)")
                print(nr_sets)
                print(set_connected_checker)
                print(specs$modes[c(nr_sets,set_connected_checker)])

                shuffled <-  get_split_sets(X = matrix_X, Y = matrix_Y,
                                            nr_subsets = nr_subsets)

                X.sampled     <-   shuffled$X.sampled
                Y.sampled     <-   shuffled$Y.sampled
                label         <-   shuffled$labels

                modes2 = specs$modes[c(nr_sets,set_connected_checker)]

                # X = X.sampled
                # Y = Y.sampled
                # lambdas = ridge_penalty
                # non_zeros = nonzero
                # label = label
                # penalization = penalization
                # max_iterations = max_iterations
                # tolerance = tolerance
                # modes2 = modes2

                cv_results[[set_connected_checker]] <- get_non_parallel_cv(X = X.sampled,
                                                  Y = Y.sampled,
                                                  lambdas = ridge_penalty,
                                                  non_zeros = nonzero,
                                                  label = label,
                                                  penalization = penalization,
                                                  max_iterations = max_iterations,
                                                  tolerance = tolerance,
                                                  modes2 = modes2)

                cv_results[[set_connected_checker]]$abs_cors
                cv_results[[set_connected_checker]]$mean_abs_cors

                print(cv_results[[set_connected_checker]])

            }

          }

        }

        s_current_set_path <- s_path_matrix[,nr_sets]

        row_data<- sapply(cv_results[which(s_current_set_path==1)], '[[', 2)

        dim_of_grid_penalization <- length(ridge_penalty) * length(nonzero)

        aggregate_cv_results <- matrix(data =
                                         rep(0,dim_of_grid_penalization * 3),
                                       ncol = 3)

        aggregate_cv_results[,1] <- row_data[1:dim_of_grid_penalization,1]
        aggregate_cv_results[,2] <- row_data[(dim_of_grid_penalization+1):(2*dim_of_grid_penalization),1]

        for (i in 1:dim(row_data)[2]){

          aggregate_cv_results[,3] <- aggregate_cv_results[,3] +
            row_data[((2*dim_of_grid_penalization)+1):dim(row_data)[1],i]

        }

        rownames(aggregate_cv_results)   <- NULL
        colnames(aggregate_cv_results)   <- c("Ridge Penalty",
                                       "Number of Nonzeros", "mean_abs_cors")


        a = aggregate_cv_results[,3]

        best_values     <- aggregate_cv_results[which.max(a),]

        best_ridge   <- best_values[1]
        best_nonzero   <- best_values[2]

        print("aggregated CV results")
        print(aggregate_cv_results)

      }
    }else{
      #if connection to B dataset
      best_ridge    <- 0
      best_nonzero  <- 0
      cv_results    <- c()
    }

    best_ridge_penalties[[nr_sets]]   <- best_ridge
    best_nr_of_nonzeros[[nr_sets]]    <- best_nonzero
    cv_results_list[[nr_sets]]        <- cv_results
    aggr_cv_results_list[[nr_sets]]   <- aggregate_cv_results


  }


  cv_results_list
  best_ridge_penalties
  best_nr_of_nonzeros
  aggr_cv_results_list

  #**********************
  result <- list(cv_results = cv_results_list,
                 best_ridge = best_ridge_penalties,
                 best_nonzero = best_nr_of_nonzeros,
                 aggr_cv_results_list = aggr_cv_results_list
                 )


  result

}




#     #baba
#
#
#     for (set_connected_checker in 1:length(current_set_path)){
#
#       #check if dataset is connected to any other dataset or this
#       #is the first dataset
#       if(current_set_path[set_connected_checker]==1){
#
#         print("dataset itself mode B")
#         print(nr_sets)
#         print(set_connected_checker)
#
#         matrix_Y <- data_sets[,blocks[[set_connected_checker]]]
#
#         shuffled <-  get_split_sets(X = matrix_X, Y = matrix_Y,
#                                     nr_subsets = nr_subsets)
#
#         X.sampled     <-   shuffled$X.sampled
#         Y.sampled     <-   shuffled$Y.sampled
#         label         <-   shuffled$labels
#
#         cv_results <- get_non_parallel_cv(X = X.sampled,
#                                         Y = Y.sampled,
#                                         lambdas = ridge_penalty,
#                                         non_zeros = nonzero,
#                                         label = label,
#                                         penalization = penalization,
#                                         max_iterations = max_iterations,
#                                         tolerance = tolerance)
#
#         cv_results$abs_cors
#         cv_results$mean_abs_cors
#
#         a = cv_results$mean_abs_cors[,3]
#
#         best_values     <- cv_results$mean_abs_cors[which.max(a),]
#
#         best_ridge   <- best_values[1]
#         best_nonzero   <- best_values[2]
#
#       }
#     }
#
#   #if dataset is not mode B
#   } else{
#
#     #check if there is connection to one B dataset
#     if ("B" %in% specs$modes[current_set_path]){
#
#       for (set_connected_checker in 1:length(current_set_path)){
#
#         #check if dataset is connected to mode B dataset
#         if(current_set_path[set_connected_checker]==1 &&
#            specs$modes[set_connected_checker] == "B"){
#
#           #save current dataset that is not mode B but connected to other mode B
#           matrix_X <- data_sets[,blocks[[nr_sets]]]
#
#           #save connected mode B dataset
#           matrix_Y <- data_sets[,blocks[[set_connected_checker]]]
#
#           print("dataset itself not mode B but connected to mode B")
#           print(nr_sets)
#
#           shuffled <-  get_split_sets(X = matrix_X, Y = matrix_Y,
#                                       nr_subsets = nr_subsets)
#
#           X.sampled     <-   shuffled$X.sampled
#           Y.sampled     <-   shuffled$Y.sampled
#           label         <-   shuffled$labels
#
#           cv_results <- get_non_parallel_cv(X = X.sampled,
#                                             Y = Y.sampled,
#                                             lambdas = ridge_penalty,
#                                             non_zeros = nonzero,
#                                             label = label,
#                                             penalization = penalization,
#                                             max_iterations = max_iterations,
#                                             tolerance = tolerance)
#
#           cv_results$abs_cors
#           cv_results$mean_abs_cors
#
#           a = cv_results$mean_abs_cors[,3]
#
#           best_values     <- cv_results$mean_abs_cors[which.max(a),]
#
#           best_ridge   <- best_values[1]
#           best_nonzero   <- best_values[2]
#
#         }
#       }
#
#     }
#     else{
#       #if connection to B dataset
#       best_ridge    <- 0
#       best_nonzero  <- 0
#       cv_results    <- 0
#     }
#
#   }
#
#   best_ridge_penalties[[nr_sets]]   <- best_ridge
#   best_nr_of_nonzeros[[nr_sets]]    <- best_nonzero
#   cv_results_list[[nr_sets]]        <- cv_results
#
# }
