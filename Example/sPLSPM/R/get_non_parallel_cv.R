get_non_parallel_cv <- function(X,
                                Y,
                                lambdas,
                                non_zeros,
                                label,
                                penalization,
                                max_iterations,
                                tolerance,
                                modes2){
  #non-parallel x-fold Cross validation for each alpha, for example alpha = 0, 0.1, ... , 0.9, 1.0
  #
  #input:     X             n*p matrix          - independent variables
  #           Y             n*q matrix          - dependent variables
  #           lambda        integer             - factor for Ridge penalty
  #           labels         vector with n_subset unique elements with length of n - for crossvalidation
  #
  #output:    abs_cor       vector of integers  - returns the summed abs correlation
  #           stime         int/time            - return algorithms running time

  #magic variables for testing
  # lambdas <- lambda
  # non_zeros <- nonzero
  # X = X.sampled
  # Y = Y.sampled
  # tolerance = tol
  # max_iterations = maxiter
  #********************


  nr_subsets      <-    length(unique(label))
  nmodes          <-    modes2
  abs_cors        <-    c()
  iterations_m    <-    c()

  length(lambdas)
  length(non_zeros)

  ## RDA objective function
  if(nmodes[1]!="B" || nmodes[2]!="B"){

    #measure time
    stime <- system.time({
      for (l in 1:length(lambdas)){

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


            ndata_sets <- cbind(X.train,Y.train)

            # Only 2 DATASETS####
            nEXPL_X = c(0,0)
            nRESP_Y = c(1,0)
            npath_matrix = rbind(nEXPL_X, nRESP_Y)

            # blocks of outer model
            nblocks = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2])

            # print out warning message only at the last iteration
            # if (i == nr_subsets){
            #   warning_non_convergence = TRUE
            # } else {
            #   warning_non_convergence = FALSE
            # }


            sub_results[[i]] <- splspm(ndata_sets,
                                       npath_matrix,
                                       nblocks,
                                       nmodes,
                                       scheme="path",
                                       penalization = penalization,
                                       nonzero = non_zeros[nz],
                                       lambda = lambdas[l],
                                       maxiter = max_iterations,
                                       tol = tolerance,
                                       warning_non_convergence = FALSE)

            #calculate alpha by univariate regression
            ALPHA <- t(X.train)%*% sub_results[[i]]$scores[,1]%*%
              solve(t( sub_results[[i]]$scores[,1])%*%
                      sub_results[[i]]$scores[,1])

            XI.test = scale(X.test) %*% ALPHA

            sub_abs_cor[[i]] <- sum((abs(cor(XI.test,Y.test))))/dim(Y.train)[2]

            #sub_iterations_m[[i]] <- sub_results[[i]]$Nr_iterations


          }#end of subset for loop

          abs_cors        <- cbind(abs_cors, sub_abs_cor)
          #iterations_m    <- cbind(iterations_m, sub_iterations_m)

        }#end of non_zeros loop

      }#end of lambda for loop
    })[3]#end of measure time

  }#END OF RDA style
  else{
  #START OF CCA STYLE
    #measure time
    stime <- system.time({
      for (l in 1:length(lambdas)){

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


            ndata_sets <- cbind(X.train,Y.train)

            # Only 2 DATASETS####
            nEXPL_X = c(0,1)
            nRESP_Y = c(1,0)
            npath_matrix = rbind(nEXPL_X, nRESP_Y)

            # blocks of outer model
            nblocks = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2])

            # print out warning message only at the last iteration
            # if (i == nr_subsets){
            #   warning_non_convergence = TRUE
            # } else {
            #   warning_non_convergence = FALSE
            # }


            sub_results[[i]] <- splspm(ndata_sets,
                                       npath_matrix,
                                       nblocks,
                                       nmodes,
                                       scheme="path",
                                       penalization = penalization,
                                       nonzero = non_zeros[nz],
                                       lambda = lambdas[l],
                                       maxiter = max_iterations,
                                       tol = tolerance,
                                       warning_non_convergence = FALSE)

            ALPHA <- t(X.train)%*% sub_results[[i]]$scores[,1]%*%
              solve(t( sub_results[[i]]$scores[,1])%*%
                      sub_results[[i]]$scores[,1])

            XI.test = scale(X.test) %*% ALPHA


              BETA <- t(Y.train)%*% sub_results[[i]]$scores[,2]%*%
                solve(t( sub_results[[i]]$scores[,2])%*%
                        sub_results[[i]]$scores[,2])

              ETA.test = scale(Y.test) %*% BETA

              ##CCA objective function
              #devide with dim(Y.train)[2]
              sub_abs_cor[[i]] <- sum((abs(cor(XI.test,ETA.test))))
              #/dim(XI.test)[1]

            #sub_iterations_m[[i]] <- sub_results[[i]]$Nr_iterations


          }#end of subset for loop

          abs_cors        <- cbind(abs_cors, sub_abs_cor)
          #iterations_m    <- cbind(iterations_m, sub_iterations_m)

        }#end of non_zeros loop

      }#end of lambda for loop
    })[3]#end of measure time

  }


  #Figure out lambdas and non-zeros columns in results
  labels_non_zeros  <- rep(non_zeros, dim(abs_cors)[2]/length(non_zeros))
  labels_non_zeros

  labels_lambdas    <- rep(lambdas, each=dim(abs_cors)[2]/length(lambdas))
  labels_lambdas
  length(labels_lambdas)

  all_abs_cors  <- rbind(labels_lambdas, labels_non_zeros, abs_cors)


  mean_abs_cors <- c()

  for (i in 1:length(labels_lambdas)){

    sub_result    <-  (c(labels_lambdas[i], labels_non_zeros[i], mean(abs_cors[,i], na.rm = TRUE)))
    mean_abs_cors <- rbind(mean_abs_cors, sub_result)

  }
  rownames(mean_abs_cors)   <- NULL
  colnames(mean_abs_cors)   <- c("Ridge Penalty",
                                 "Number of Nonzeros", "mean_abs_cors")
  mean_abs_cors


  #plot(mean_abs_cors[,1],mean_abs_cors[,3], pch=19, col=mean_abs_cors[,2])
  #text(mean_abs_cors[,1],mean_abs_cors[,3], labels=mean_abs_cors[,2], cex= 1.1, pos=2, pch=19, col=mean_abs_cors[,2])


  # #*********************************
  #
  # plot2 <-
  #   ggplot(data=data.frame(mean_abs_cors),
  #          aes(x = factor(Lambda), y = mean_abs_cors,
  #              group = factor(Number.of.Nonzeros),
  #              shape = factor(Number.of.Nonzeros),
  #              color = factor(Number.of.Nonzeros)))+
  #   geom_line() +
  #   geom_point() +
  #   scale_x_discrete("Lambda") +
  #   scale_y_continuous("Mean absolute correlation") +
  #   facet_grid(.~Number.of.Nonzeros )
  #
  # #*********************************

  print("Elapsed time")
  print(stime)
  #Return section**********************
  result        <-    list(abs_cors = abs_cors,
                           mean_abs_cors = mean_abs_cors,
                           stime = stime,
                           nmodes = nmodes
                           # plot2 = plot2,
                           #iterations_m = iterations_m

  )

  result
}
