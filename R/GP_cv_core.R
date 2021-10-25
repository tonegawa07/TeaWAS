#' Prediction for data split by rsample::vfold_cv "GP's core"
#'
#' @importFrom dplyr %>%
#'
#' @param fold_split data split by rsample::vfold_cv
#' @param phenotype_name The name of the phenotype to be predicted
#' @param algorithm Algorithm used for prediction
#'
GP_cv_core <- function(fold_split, phenotype_name, algorithm) {
  # training data & test data
  train <- fold_split %>% rsample::analysis()
  test <- fold_split %>% rsample::assessment()

  # train
  train_y <-
    train %>%
    dplyr::select(phenotype_name) %>%
    as.matrix() %>%
    as.numeric()

  train_x <-
    train %>%
    dplyr::select(-phenotype_name) %>%
    as.matrix()

  # test
  test_x <- test %>% dplyr::select(-phenotype_name) %>% as.matrix()

  # rrBLUP ---------------------------------------------------------------------
  if (algorithm == "GBLUP(RR)" || algorithm == "GBLUP(GAUSS)") {
    if (!requireNamespace("rrBLUP", quietly = TRUE)) {
      stop('R package "rrBLUP" is required. Please install it.')
    }
    res <- switch(algorithm,
                  "GBLUP(RR)"    = rrBLUP::kinship.BLUP(y = train_y,
                                                        G.train = train_x,
                                                        G.pred = test_x,
                                                        K.method = "RR"),
                  "GBLUP(GAUSS)" = rrBLUP::kinship.BLUP(y = train_y,
                                                        G.train = train_x,
                                                        G.pred = test_x,
                                                        K.method = "GAUSS"),
                  stop('Only can use "GBLUP(RR)" and "GBLUP(GAUSS)"')
    )
    y_pred <- res$g.pred + rep(res$beta, length(res$g.pred))

  # glmnet ---------------------------------------------------------------------
  } else if (algorithm == "Ridge" || algorithm == "Lasso" || algorithm == "ElasticNet") {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop('R package "glmnet" is required. Please install it.')
    }
    mdl <- switch(algorithm,
                  "Ridge"      = glmnet::cv.glmnet(x = train_x,
                                                   y = train_y,
                                                   alpha = 0,
                                                   standardize = F),
                  "Lasso"      = glmnet::cv.glmnet(x = train_x,
                                                   y = train_y,
                                                   alpha = 1,
                                                   standardize = F),
                  "ElasticNet" = glmnet::cv.glmnet(x = train_x,
                                                   y = train_y,
                                                   alpha = 0.5,
                                                   standardize = F),
                  stop('Only can use "Ridge", "Lasso" and "ElasticNet"')
    )
    y_pred <- stats::predict(mdl, newx = test_x, s = "lambda.min")

  # randomForest ---------------------------------------------------------------
  } else if (algorithm == "RandomForest") {
    if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop('R package "randomForest" is required. Please install it.')
    }
    mdl <- randomForest::randomForest(x = train_x, y = train_y)
    y_pred <- stats::predict(mdl, newdata = test_x)

  }

  output <- data.frame(pred = y_pred)
  colnames(output) <- "pred"
  return(output)
}
