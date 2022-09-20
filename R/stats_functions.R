#' Calculate and Print ROC Plots
#'
#' Calculates a ROC plot from case and control group values. Plots the ROC curve
#' and its sensitivity confidence interval. Returns the ROC curve points and
#' AUC value.
#' @param ctrl_group a vector of values taken by the control samples
#' @param case_group a vector of values taken by the case samples
#' @return A list containing the ROC curve points, the ROC AUC, and the Delong
#' confidence interval for the AUC.
#' @export
calc_roc <- function(ctrl_group, case_group) {
  # calculate the ROC curve
  roc_data <- pROC::roc(controls = ctrl_group,
                  cases = case_group,
                  percent=TRUE, ci = TRUE, print.auc=TRUE)
  # calculate the 95% confidence intervals by bootstrapping
  ciAUC <- pROC::ci.auc(roc_data, parallel = TRUE)
  ciSens <- pROC::ci.se(roc_data, specificities = roc_data$specificities, parallel = TRUE)

  # draw the ROC curve and print the AUC value
  pROC::plot.roc(roc_data, col="#008600", print.auc = TRUE, print.auc.y = 48)

  # plot the confidence intervals as a blue blob
  pROC::plot.roc(ciSens, type = "shape", col="#00860022", no.roc=TRUE)

  # write out the ROC curve and confidence intervals to a CSV file
  curves <- tibble::tibble(threshold = roc_data$thresholds,
                           sensitivity = roc_data$sensitivities,
                   specificity = roc_data$specificities, cisens2.5 = ciSens[,1],
                   cisens50 = ciSens[,2], cisens97.5 = ciSens[,3])
  output <- list(curves = curves, auc = roc_data$auc, ci = roc_data$ci)
  return(output)
}

#' Calculate Optimal ROC Thresholds
#'
#' Calculates the 'best' threshold for a given ROC curve by choosing the
#' threshold with highest (sensitivity+specificity), using sensitivity to break
#' ties.
#' @param cases a vector of values taken by the case samples
#' @param controls a vector of values taken by the control samples
#' @return The optimal threshold value.
#' @importFrom magrittr %>%
#' @export
pick_roc_threshold <- function(cases, controls) {
  if(!requireNamespace("pROC", quietly = TRUE)) {
    stop("Function pick_roc_threshold requires package pROC to work. Please install it.",
         call. = FALSE)
  }
  hs_roc <- pROC::roc(controls = controls, cases = cases, quiet = TRUE,
                ci = FALSE, percent = TRUE, print.auc = FALSE)

  # find the optimum threshold, where sensitivity + specificity is maximized
  df <- tibble::tibble(threshold = hs_roc$thresholds,
                       sensitivity = hs_roc$sensitivities,
                       specificity = hs_roc$specificities) %>%
    dplyr::mutate(sum = sensitivity + specificity) %>%
    dplyr::filter(is.finite(threshold))

  # we pick the 'best' threshold, erring on the side of sensitivity if there are ties
  hs_thresh <- df %>%
    dplyr::top_n(n = 1, sum) %>%
    dplyr::top_n(n = 1, sensitivity) %>%
    .$threshold

  return(hs_thresh)
}

#' Multinomial Deviance Function Stolen from Glmnet
#'
#' The glmnet package doesn't give us the individual sample deviance values, so
#' I've ripped out the relevant function and put it here. This function
#' calculates the multinomial deviance for a glmnet model on a given set of
#' samples. This is particularly useful for calculating model performance on
#' hold-out sets.
#'
#' @param mdl a previously fitted glmnet model
#' @param x a vector of values for which to make predictions
#' @param y the reference category assignments for the samples in x
#' @return a matrix containing the deviance values for each lambda value in the
#'   original model. Matrix has one row per sample and one column per lambda.
#' @export
multinomial_deviance <- function(mdl, x, y) {
  if(class(mdl)[1] != "multnet") {
    stop(simpleError("Tried to call multinomial_deviance with non-multinomial logit model!"))
  }

  prob_min = 1e-05
  prob_max = 1 - prob_min

  lambda = mdl$lambda
  nlambda=length(lambda)

  # number of categories in y
  nc = dim(y)
  # convert y into matrix with 1 col per category
  # single "1" in each row indicates sample category
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
  } else nc = nc[2]

  # 3D matrix: row per sample, col per category, for each lambda value
  predmat = array(NA, c(nrow(y), nc, nlambda))

  # use trained model to predict samples
  # row per sample, col per category, for each lambda value
  preds = predict(object = mdl, newx = x, s = lambda, type = "response")

  # make sure the predictions have a compatible number of lambda values
  nlami = min(dim(preds)[3], nlambda)
  # save predicted probabilities for later
  predmat[ , , seq(nlami)] = preds[ , , seq(nlami)]
  # pad with last prediction if insufficient lambda values
  if(nlami < nlambda) {
    predmat[ , , seq(from=nlami,to=nlambda)] = preds[ , , nlami]
  }

  ywt = apply(y, 1, sum)
  # scale response values if they are part of two groups somehow
  y = y/ywt
  # weights = weights * ywt
  N = nrow(y) - apply(is.na(predmat[, 1, ]), 2, sum)
  # row per sample, 1 in column = group, repeated for all lambda values
  bigY = array(y, dim(predmat))
  # force all values in predmat to be between prob_min and prob_max
  predmat = pmin(pmax(predmat, prob_min), prob_max)

  # log probabilities for each sample/lambda, only for samples' actual groups
  lp = bigY * log(predmat)
  # ly ends up all zeroes unless something was in multiple groups
  ly = bigY * log(bigY)
  ly[bigY == 0] = 0
  # raw deviance, row per sample, column per lambda
  cvraw = apply(2 * (ly - lp), c(1, 3), sum)

  # calculate mean deviance across all samples, one entry per lambda
  # cvm = apply(cvraw, 2, mean, na.rm = TRUE)
  # # calculate SD of deviance across all samples, one entry per lambda
  # cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean,
  #                   na.rm = TRUE)/(N - 1))
  return(cvraw)#list(cvraw = cvraw, cvm = cvm, cvsd = cvsd))
}

#' Binomial Deviance Function Stolen from Glmnet
#'
#' The glmnet package doesn't give us the individual sample deviance values, so
#' I've ripped out the relevant function and put it here. This function
#' calculates the binomial deviance for a glmnet model on a given set of
#' samples. This is particularly useful for calculating model performance on
#' hold-out sets.
#'
#' @param mdl a previously fitted glmnet model
#' @param x a vector of values for which to make predictions
#' @param y the reference category assignments for the samples in x
#' @return a matrix containing the deviance values for each lambda value in the
#'   original model. Matrix has one row per sample and one column per lambda.
#' @export
binomial_deviance <- function (mdl, x, y)
{
  if(class(mdl)[1] != "lognet") {
    stop(simpleError("Tried to call multinomial_deviance with non-multinomial logit model!"))
  }

  prob_min = 1e-05
  prob_max = 1 - prob_min

  # get the lambda sequence used in the model
  lambda = mdl$lambda
  nlambda=length(lambda)

  nc = dim(y)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
  }

  N = nrow(y)

  # prediction matrix - 1 row/sample, 1 col/lambda value
  predmat = matrix(NA, nrow(y), nlambda)
  # nlams = double(nfolds)

  # use trained model to predict samples
  # row per sample, col per category, for each lambda value
  preds = glmnet::predict.glmnet(mdl, x, s = lambda, type = "response")

  # put predictions into the prediction matrix
  nlami = min(ncol(preds),nlambda)
  predmat[, seq(nlami)] = preds[,seq(nlami)]
  if(nlami < nlambda) {
    predmat[which,seq(from=nlami,to=nlambda)]=preds[,nlami]
  }

  # calculate deviance
  ywt = apply(y, 1, sum)
  y = y/ywt
  N = nrow(y) - apply(is.na(predmat), 2, sum)

  predmat = pmin(pmax(predmat, prob_min), prob_max)
  lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
  ly = log(y)
  ly[y == 0] = 0
  ly = drop((y * ly) %*% c(1, 1))
  cvraw = 2 * (ly - lp)

  # return raw deviance, 1 row/sample, 1 col/lambda value
  return(cvraw)

  # cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
  # cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
  #                   w = weights, na.rm = TRUE)/(N - 1))
  # out = list(cvm = cvm, cvsd = cvsd, type.measure=type.measure)
  # if (keep)
  #   out$fit.preval = predmat
  # out
}

#' Absolute Error of glmnet Models
#'
#' This function calculates the absolute error a glmnet model on a given
#' set of samples. This is particularly useful for calculating model performance
#' on hold-out sets. Note that this is NOT producing the mean absolute error!
#'
#' @param mdl a previously fitted glmnet model of any type except "cox"
#' @param x a vector of samples for which to make predictions
#' @param y the actual response values of the samples in x
#' @return a matrix containing the absolute error for each lambda value in
#'   the original model. Matrix has one row per sample and one column per lambda.
#' @export
absolute_error <- function (mdl, x, y)
{
  if(class(mdl)[1] == "coxnet") {
    stop(simpleError("Tried to call mean_absolute_error with a Cox model!"))
  }

  prob_min = 1e-05
  prob_max = 1 - prob_min

  # get the lambda sequence used in the model
  lambda = mdl$lambda
  nlambda=length(lambda)

  # nc = dim(y)
  # if (is.null(nc)) {
  #   y = as.factor(y)
  #   ntab = table(y)
  #   nc = as.integer(length(ntab))
  #   y = diag(nc)[as.numeric(y), ]
  # }

  # N = nrow(y)

  # use trained model to predict samples
  # 1 row/sample, 1 col/lambda value
  predmat = glmnet::predict.glmnet(mdl, x, s = lambda, type = "response")

  # calculate absolute error
  raw_error <- abs(sweep(predmat, MARGIN = 1, STATS = as.array(y), FUN = "-"))
  return(raw_error)
}

#' Calculate Misclassification for a glmnet Model
#'
#' This function calculates the proportion of a given set of samples that are
#' mis-classified by a glmnet model. This is particularly useful for calculating
#' model performance on hold-out sets.
#'
#' @param mdl a previously fitted glmnet model
#' @param x a vector of values for which to make predictions
#' @param y the reference category assignments for the samples in x
#' @return a matrix containing the misclassification error for each lambda value
#'   in the original model. Matrix has one row per sample and one column per
#'   lambda.
#' @export
misclassification_error <- function (mdl, x, y)
{
  if(class(mdl)[2] != "glmnet") {
    stop(simpleError("Tried to call misclassification_error with model not of class glmnet!"))
  }
  if(!(class(mdl)[1] %in% c("lognet", "multnet"))) {
    stop(simpleError("misclassification_error requires a binomial or multinomial model!"))
  }

  # get the lambda sequence used in the model
  lambda = mdl$lambda

  # use trained model to predict samples
  # row per sample, col per lambda value
  preds = stats::predict(mdl, newx = x, s = lambda, type = "class")

  # return misclassification error rates for each lambda value
  mc_error <- colSums(preds != y)/nrow(preds)
  return(mc_error)
}

#' Repeated Cross-Validation of Sparse Logistic Models
#'
#' This function performs repeated k-fold cross-validation of sparse logistic
#' regression models. Unlike the vanilla cv.glmnet method, this version returns
#' the raw deviance values and selected features for all the models. This allows
#' proper selection of the optimal lambda value and the most common features.
#' Note that this method currently only optimizes on deviance (binomial models),
#' multinomial deviance (multinomial models), or absolute error (gaussian
#' models).
#'
#' @param x a matrix with one row per sample and one column per feature
#' @param y a vector containing the true class assignments for all samples in x
#' @param family the model family to use. See [glmnet::glmnet()] for details.
#'   Currently only supports binomial, multinomial, and gaussian.
#' @param eta the proportion of ridge regression to 'mix' into the LASSO
#'   regression. A small amount (~0.05) can help improve model stability.
#' @param nfolds the number of cross-validation folds to use
#' @param nreps the number of times to repeat cross-validation
#'
#' @return A list with the following fields:
#'   \describe{
#'     \item{full_fit}{a glmnet model fit on all the data}
#'     \item{lambda_min}{the lambda value where deviance is minimized when
#'     averaged across all cross-validation folds.}
#'     \item{lambda_1se}{the lambda value giving the simplest model with mean
#'     deviance within 1 standard error of the lambda_min model.}
#'     \item{cv_folds}{the actual cross-validation training folds used.}
#'     \item{cv_fits}{the glmnet models fit on each training fold.}
#'     \item{lambda_perf}{a tibble containing mean and variation in model
#'     deviance when tested on the hold-out sets for each lambda value.}
#'   }
#'
#' @seealso [glmnet::cv.glmnet()] for the original function on which this one
#'   was based.
#'
#' @importFrom magrittr %>%
#' @export
multi.cv.glmnet <- function(x, y, family = "binomial", eta = 0, nfolds = 10,
                            nreps = 10) {
  # select the proper deviance function
  if(family == "binomial") {
    dev_func <- binomial_deviance
  } else if(family == "multinomial") {
    dev_func <- multinomial_deviance
  } else if(family == "gaussian") {
    dev_func <- absolute_error
  } else {
    stop("Unrecognized glmnet family: ", family, call. = FALSE)
  }
  # do a full fit to get the lambda sequence
  init_fit <- glmnet::glmnet(x = x, y = y, alpha = 1 - eta, family = family)
  lambda_seq <- init_fit$lambda

  # split the data into multiple k-fold cross-validation samples
  cv_folds <- caret::createMultiFolds(y = y, k = nfolds, times = nreps)

  # matrix to store the fits so we can extract coefficients later
  cv_fits <- vector(mode = "list", length = length(cv_folds))
  # matrix to store deviances. row per sample, column per lambda, per holdout set
  cv_deviances <- array(data = NA, dim = c(length(y),
                                           length(lambda_seq),
                                           length(cv_folds)))

  # for each training/test set combo
  #TODO: convert to a foreach/dopar loop?
  for(i in seq_along(cv_folds)) {
    fold <- cv_folds[[i]]
    # separate out the train/test sets
    ho_idx <- seq_along(y)[-1*fold]
    xt <- x[fold, ]
    xh <- x[-1*fold, ]
    yt <- y[fold]
    yh <- y[-1*fold]
    # fit the model
    cv_fits[[i]] <- glmnet::glmnet(x = xt, y = yt, alpha = 1 - eta,
                                   family = family, lambda = lambda_seq)
    # calculate deviance on hold-out set
    cv_deviances[ho_idx, , i] <- dev_func(mdl = cv_fits[[i]], x = xh, y = yh)
  }

  # calculate means and standard deviations for each lambda, combining all samples and holdout sets
  cv_means <- apply(cv_deviances, MARGIN = c(2), FUN = mean, na.rm = TRUE)
  cv_sds <- apply(cv_deviances, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)

  cv_lasso_results <- tibble::tibble(Lambda = lambda_seq, Mean = cv_means, SD = cv_sds) %>%
    dplyr::mutate(SEM = SD/sqrt(length(cv_folds))) %>%
    dplyr::mutate(Upper = Mean + SEM, Lower = Mean - SEM)

  # which lambda value gives the lowest average deviance on hold-out sets?
  lambda_min <- cv_lasso_results %>%
    dplyr::slice_min(order_by = Mean, n = 1, with_ties = FALSE) %>%
    .$Lambda

  # which lambda value gives the simplest model within 1-SEM of the minimum?
  dev_thresh <- cv_lasso_results %>%
    dplyr::slice_min(order_by = Mean, n = 1, with_ties = FALSE) %>%
    .$Upper
  lambda_1se <- cv_lasso_results %>%
    dplyr::filter(Lambda >= lambda_min) %>%
    dplyr::filter(Mean <= dev_thresh) %>%
    dplyr::slice_max(order_by = Mean, n = 1, with_ties = FALSE) %>%
    .$Lambda

  output <- list(full_fit = init_fit, lambda_min = lambda_min,
                 lambda_1se = lambda_1se, cv_folds = cv_folds,
                 cv_fits = cv_fits, lambda_perf = cv_lasso_results)
  class(output) <- "multi.cv.glmnet"
  return(output)
}

#' Extract Coefficients from a "multi.cv.glmnet" Object
#'
#' This function extracts the model coefficients at a single lambda value from
#' each of the models in a multi.cv.glmnet object. This is particularly useful
#' for evaluating model stability and selecting robust features for a final
#' model. If the models only consistently select a small subset of features, it
#' is often useful to fit a non-sparse logistic model using just those features.
#'
#' @param object a multi.cv.glmnet object from which to extract coefficients
#' @param s the desired lambda value to use. Special values "lambda.min" and
#'   "lambda.1se" automatically select the corresponding values from the
#'   multi.cv.glmnet object.
#' @param exact if the desired lambda value was not part of the original
#'   sequence, should the models be refit to exactly calculate the coefficients?
#'   See [glmnet::coef.glmnet()] for details.
#' @param ... not used
#'
#' @return a tibble containing every non-zero feature coefficient from every
#'   cross-validation model
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @export
coef.multi.cv.glmnet <- function(object, s = "lambda.1se", exact = FALSE, ...) {
  # check object class
  if(!inherits(object, "multi.cv.glmnet")) {
    stop("This method requires an object of class multi.cv.glmnet!", call. = FALSE)
  }
  # set up desired lambda value
  if(s == "lambda.1se") {
    s <- object$lambda_1se
  } else if(s == "lambda.min") {
    s <- object$lambda_min
  } else if(!is.numeric(s)) {
    stop("Lambda value must be numeric, not ", typeof(s), "!")
  } else if(length(s) != 1) {
    stop("Please supply only a single lambda value!")
  }

  # we have to treat multinomial models differently
  # here's the standard procedure
  if(class(object$full_fit)[1] != "multnet") {
    cv_fits <- foreach::foreach(fit = object$cv_fits, .combine = "rbind",
                                .packages = c("glmnet")) %dopar% {
        t(stats::coef(fit, s = s, exact = exact)[,"s1"])
      }

    # summarize the results
    cv_coefficients <- cv_fits %>%
      tibble::as_tibble() %>%
      tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Feature",
                   values_to = "Coefficient") %>%
      dplyr::filter(Coefficient != 0) %>%
      dplyr::filter(Feature != "(Intercept)")

    return(cv_coefficients)
  } else {
    # special multinomial procedure
    cv_fits <- foreach::foreach(fit = object$cv_fits, .combine = "rbind", .packages = c("glmnet")) %dopar% {
      # pull coefficients and make non-sparse
      df <- stats::coef(fit, s = s, exact = exact) %>%
        lapply(FUN = as.matrix)

      # grab feature names
      rn <- rownames(df[[1]])

      # convert coefficients to a wide matrix
      df <- t(as.matrix(as_tibble(df)))
      # put in the feature names
      colnames(df) <- rn
      # return the matrix
      df
    }

    # summarize the results
    cv_coefficients <- cv_fits %>%
      tibble::as_tibble(rownames = "Level") %>%
      tidyr::pivot_longer(cols = -Level, names_to = "Feature",
                          values_to = "Coefficient") %>%
      dplyr::filter(Coefficient != 0) %>%
      dplyr::filter(Feature != "(Intercept)")

    return(cv_coefficients)
  }
}
