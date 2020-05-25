#        *         *         *         *         *         *         *         #
# Metadata                                                                  ####
#
# Title:    sparese RDA
# meta:     sparese RDA
#
# code by:  Koos Zwinderman, Attila Csala
# place:    Amsterdam
# date:     2016-11
#        *         *         *         *         *         *         *         #

#' Generate two data sets with some highly correlated variables
#'
#' Generate two data sets with some highly correlated variables.
#' @param number_of_ksi The number of latent variables of the predictive data
#' set.
#' @param number_of_patients The number of observations (rows) in the data sets.
#' @param number_of_Xs_associated_with_ksis Number of variables of the
#' predictive data set that is associated with the predicted data set.
#' @param number_of_not_associated_Xs Number of variables of the predictive data
#' set that is not associated with the predicted data set.
#' @param mean_of_the_regression_weights_of_the_associated_Xs Mean of the
#' regression weights of the predictive varaibles that are associated with the
#' predicted variables.
#' @param sd_of_the_regression_weights_of_the_associated_Xs Standard deviation
#' of the regression weights of the predictive varaibles that are associated
#' with the predicted variables.
#' @param Xnoise_min The lower bound of the unifrom distribution that is used to
#' sample the values for the regression weights of the predictive varaibles that
#' are not associated with the predicted variables.
#' @param Xnoise_max The upper bound of the unifrom distribution that is used to
#' sample the values for the regression weights of the predictive varaibles that
#' are not associated with the predicted variables.
#' @keywords generate data
#' @import mvtnorm
#' @export
#' 
#' @examples
#' dataXY <- generate_data(number_of_ksi = 1,
#'                            number_of_patients = 500,
#'                            number_of_Xs_associated_with_ksis = c(5),
#'                            number_of_not_associated_Xs = 250,
#'                            mean_of_the_regression_weights_of_the_associated_Xs =
#'                              c(0.9),
#'                            sd_of_the_regression_weights_of_the_associated_Xs =
#'                              c(0.05),
#'                            Xnoise_min = -0.3, Xnoise_max = 0.3,
#'                            
#'                            number_of_Ys_associated_with_ksis = c(5),
#'                            number_of_not_associated_Ys = 350,
#'                            mean_of_the_regression_weights_of_the_associated_Ys =
#'                              c(0.9),
#'                            sd_of_the_regression_weights_of_the_associated_Ys =
#'                              c(0.05),
#'                            Ynoise_min = -0.3, Ynoise_max = 0.3)
#'
#' # seperate predictor and predicted sets
#' X <- dataXY$X
#' Y <- dataXY$Y
#' 
#' dim(X);dim(Y)
#' 
#' @import mvtnorm

#******************************************************************************#
#                        Generate Data a'la Koos                       #########
#******************************************************************************#

generate_data <- function(number_of_ksi = 1,
                          number_of_patients = 50,
                          number_of_Xs_associated_with_ksis = c(5),
                          number_of_not_associated_Xs = 250,
                          mean_of_the_regression_weights_of_the_associated_Xs =
                              c(0.7),
                          sd_of_the_regression_weights_of_the_associated_Xs =
                              c(0.05),
                          Xnoise_min = -0.3, Xnoise_max = 0.3,

                          number_of_Ys_associated_with_ksis = c(5),
                          number_of_not_associated_Ys = 350,
                          mean_of_the_regression_weights_of_the_associated_Ys =
                              c(0.7),
                          sd_of_the_regression_weights_of_the_associated_Ys =
                              c(0.05),
                          Ynoise_min = -0.3, Ynoise_max = 0.3){
  #  make the RDA components
  #  number_of_ksi = 3
  #  number_of_patients = 500

  #*Generate latent vector ski from lultivariATE NORMAL distribution with
  #    mean 0, and covariance matrix = 1
  ksi=rmvnorm(number_of_patients,mean=rep(0,number_of_ksi),
              sigma=diag(number_of_ksi))

  # make X data
  #  number_of_Xs_associated_with_ksis = c(5,5,10)
  #  number_of_not_associated_Xs = 120
  #  mean_of_the_regression_weights_of_the_associated_Xs = c(0.6,0.7,0.8)
  #  sd_of_the_regression_weights_of_the_associated_Xs = c(0.05,0.05,0.05)

  #*make empty matrix with nrow and ncol
  x=matrix(NA,
           nrow=number_of_patients,
           ncol=(sum(number_of_Xs_associated_with_ksis)
                 +number_of_not_associated_Xs))
  columncount=0

  #i loops trhough number of ksi's
  for (i in 1:length(number_of_Xs_associated_with_ksis)) {

    #j loops throguh number of associated X's
    #8 since its starts from 1:, it will make the
    #  first j elements associated with ksi
    for (j in 1:number_of_Xs_associated_with_ksis[i]) {

      #calculates a single regression weight for the associated X
      regressionweight =
          rnorm(1,
                mean_of_the_regression_weights_of_the_associated_Xs[i],
                sd_of_the_regression_weights_of_the_associated_Xs[i])

      columncount=columncount+1

      #sd cannot be lower than 0
      #value = rnorm(number_of_patients,mean=regressionweight*ksi[,i],
      #sd=sqrt(max(0,(1-regressionweight^2))+0.001))
      #print(value)
      x[,columncount] = rnorm(number_of_patients,mean=regressionweight*ksi[,i],
                              sd=sqrt(max(0,(1-regressionweight^2))+0.001))

    }

  }

  #fill in number of not associated columns
  for (j in 1:number_of_not_associated_Xs) {

    columncount=columncount+1
    #x[,columncount] = rnorm(number_of_patients,mean=0,sd=1)
    x[,columncount] = runif(number_of_patients,min=Xnoise_min,max=Xnoise_max)
  }

  #cor(x)[1:12,1:12]


  # make Y data
  #   number_of_Ys_associated_with_ksis = c(5,5,10)
  #   number_of_not_associated_Ys = 120
  #
  #   mean_of_the_regression_weights_of_the_associated_Ys = c(0.6,0.7,0.8)
  #   sd_of_the_regression_weights_of_the_associated_Ys = c(0.05,0.05,0.05)
  y=matrix(NA,
           nrow=number_of_patients,ncol=(sum(number_of_Ys_associated_with_ksis)
                                         +number_of_not_associated_Ys))

  columncount=0

  for (i in 1:length(number_of_Ys_associated_with_ksis)) {

    for (j in 1:number_of_Ys_associated_with_ksis[i]) {

      regressionweight =
          rnorm(1,
                mean_of_the_regression_weights_of_the_associated_Ys[i],
                sd_of_the_regression_weights_of_the_associated_Ys[i])
      columncount=columncount+1
      y[,columncount] =
          rnorm(number_of_patients,mean=regressionweight*ksi[,i],
                sd=max(0,(1-regressionweight^2))+0.001)

    }
  }


  for (j in 1:number_of_not_associated_Ys) {
    columncount=columncount+1
    #y[,columncount] = rnorm(number_of_patients,mean=0,sd=1)
    y[,columncount] = runif(number_of_patients,min=Ynoise_min,max=Ynoise_max)
  }

  X <- x
  Y <- y
  data_info <- data.frame(number_of_ksi,
                          number_of_patients,
                          number_of_Xs_associated_with_ksis,
                          number_of_not_associated_Xs,
                          mean_of_the_regression_weights_of_the_associated_Xs,
                          sd_of_the_regression_weights_of_the_associated_Xs,
                          Xnoise_min, Xnoise_max,
                          number_of_Ys_associated_with_ksis,
                          number_of_not_associated_Ys,
                          mean_of_the_regression_weights_of_the_associated_Ys,
                          sd_of_the_regression_weights_of_the_associated_Ys,
                          Ynoise_min, Ynoise_max)

  result <- list(
    X,Y,data_info
  )

  names(result) <- c(
    "X","Y","data_info"
  )

  result

}
