#######################################################################
#' @name phi
#' @title Phi coefficient of correlation 
#' @description The phi coefficient is  a correlation coefficient applied 
#' to dichotomous data. Given a two x two table of counts\cr
#' | a | b | R1 |\cr
#' | c | d | R1 |\cr
#' |---|---|----|\cr
#' |C1 | C2| n  |\cr
#' or a vector c(a,b,c,d) of frequencies.
#' @usage phi(x)
#' @param x a 1 x 4 vector or a matrix 2 x 2 of frequencies.
#' @details The coefficient phi is calculated from \deqn{(ad - bc)/\sqrt{p_qp_2q_1q_2}}
#' where \deqn{p_i} and \deqn{q_i} are the ratios of the dichotomous variables. 
#' @return the value of the phi coefficient correlation.
#' @references Warrens, Matthijs (2008), On Association Coefficients for 2x2 Tables and Properties That Do Not Depend on the Marginal Distributions. Psychometrika, 73, 777-789.
#' @references Yule, G.U. (1912). On the methods of measuring the association between two attributes. Journal of the Royal Statistical Society, 75, 579-652.
#' @export
phi <- function (x) 
{
  stopifnot(prod(dim(x)) == 4 || length(x) == 4)
  if (is.vector(x)) 
    t <- matrix(x, 2)
  rsum <- rowSums(x)
  csum <- colSums(x)
  total <- sum(rsum)
  rsum <- rsum/total
  csum <- csum/total
  v <- prod(rsum, csum)
  phi <- (x[1, 1]/total - csum[1] * rsum[1])/sqrt(v)
  names(phi) <- NULL
  return(round(phi, 2))
}


#######################################################################
#' @name alpha_curve
#' @title Cronbach-Mesbah Curve
#' @description Assess the unidimensionality of a set of items from alpha coefficient.
#' @usage alpha_curve(data)
#' @param data a Dataframe that holds the test response data
#' @details Constructs the curve by iteration:\cr
#' 1. Compute alpha using all items.\cr
#' 2. Remove one item from the scale. The removed item is the one which leaves the scale
#' with its maximum alpha value. If we remove a bad item, the alpha  coefficient  will
#' increase. Whereas if we remove a good item, the alpha must decrease.\cr
#' 3. Repeat until only  two items remain.\cr
#' This function was extracted from CMC_1.0 package.
#' @return The number of items used to calculate the coefficient. 
#' @return The maximum value of the alpha coefficient calculated at each step. 
#' @return The item removed at each step. 
#' @return Tue Cronbach-Mesbah curve plot.
#' @references Cameletti, M. & Caviezel, V. (2010). Checking the unidimensionality in R
#' using the Cronbach-Mesbah curve.
#' @references Mesbah, M. (2010). Statistical quality of life. In "Method and Applications of Statistics in the 
#' Life and Health Sciences", N. BalakrishnanEd., Wiley, pp. 839-864. 
#' @export
alpha_curve=function (data) 
{
  n = nrow(data)
  k = ncol(data)
  max.vec = c()
  which.max.vec = c()
  complete.alpha =alpha_cronbach(data)$alpha
  j = 1
  elements.to.remove.1 = seq(1:k)
  alpha.1 = c()
  for (i in 1:length(elements.to.remove.1)) {
    data.reduced = data[, -elements.to.remove.1[i]]
    alpha.1[i] = as.numeric(alpha_cronbach(data.reduced)$alpha)
  }
  max.vec[j] = max(alpha.1)
  which.max.vec[j] = which.max(alpha.1)
  for (j in 2:(k - 2)) {
    elements.to.remove = matrix(0, nrow = (k - j + 1), ncol = (j - 
                                                                 1))
    for (r in 1:length(which.max.vec)) {
      elements.to.remove[, r] = matrix(rep(which.max.vec[r], 
                                           k - j + 1), ncol = 1)
    }
    elements.to.remove = cbind(elements.to.remove, seq(1:k)[-which.max.vec[1:(j - 
                                                                                1)]])
    alpha = c()
    for (i in 1:nrow(elements.to.remove)) {
      data.reduced = data[, -elements.to.remove[i, ]]
      alpha[i] = alpha_cronbach(data.reduced)$alpha
    }
    max.vec[j] = max(alpha)
    which.max.vec[j] = elements.to.remove[which.max(alpha), 
                                          j]
  }
  output = data.frame(N.Item = seq(2, k), `Alpha Max` = c(rev(max.vec), 
                                                          complete.alpha), `Removed Item` = c(colnames(data)[rev(which.max.vec)], 
                                                                                              "--"))
  plot(output[, 1], output[, 2], t = "b", xlab = "N.Item", 
       ylab = "Alpha Max")
  text(seq(2, k), output[, 2], c(colnames(data)[rev(which.max.vec)], 
                                 ""), pos = 3, cex = 0.6)
  return(output)
}


#######################################################################
#' @name alpha_cronbach
#' @title Cronbach's alpha
#' @description Cronbach's alpha measures how correlated are the 
#' items in a test. Taken from ltm::cronbach.alpha 
#' @usage alpha_cronbach(data,standardized,CI,probs,B,na.rm)
#' @param data a matrix or a Dataframe that holds the test response data
#' @param standardized Standarize data
#' @param CI boolean
#' @param probs Vector of probabilities
#' @param B Vector Size
#' @param na.rm Remove Na values
#' @details the coefficient is calculated \deqn{\alpha = (n/n-1)*(1 - (\sum V_i/V_t))}
#' where \deqn{V_t} is the variance of test scores and \deqn{V_i} is the variance of item scores.
#' It is desirable that the items are closely interrelated (coefficient near 1).
#' This function was extracted from multilevel_2.5 package.
#' @return Cronbach's alpha for the test and the number of individuals of test.
#' @references Cronbach L. J. (1951) Coefficient Alpha and the internal structure of tests. Psychometrika, 16,297-334
#' @export
alpha_cronbach <- function (data, standardized = FALSE, CI = FALSE, probs = c(0.025, 0.975), B = 1000, na.rm = FALSE) {
    if (!inherits(data, "matrix") && !inherits(data, "data.frame"))
      stop("'data' must be either a data.frame or a matrix.\n")
    n <- nrow(data)
    p <- ncol(data)
    nam <- deparse(substitute(data))
    if (p < 2)
      stop("'data' should have more than two columns.\n")
    data <- data.matrix(data)
    if (!na.rm && any(is.na(data))) {
      stop("missing values in 'data'.\n")
    } else {
      alpha <- if (!standardized) {
        VarTot <- var(rowSums(data[complete.cases(data), ]))
        VarInd <- sum(apply(data, 2, sd, na.rm = TRUE)^2)
        (p / (p - 1)) * (1 - (VarInd / VarTot))
      } else {
        mat <- cor(data, use = "complete.obs")
        ave.rho <- mean(mat[upper.tri(mat)])
        (p * ave.rho) / (1 + (p - 1) * ave.rho)
      }        
    }
    out <- list(alpha = alpha, n = n, p = p, standardized = standardized, name = nam)
    if (CI) {
      T.boot <- numeric(B)
      for (i in 1:B) {
        data.boot <- data[sample(1:n, replace = TRUE), ]
        T.boot[i] <- if (!standardized) {
          VarTot <- var(rowSums(data.boot[complete.cases(data.boot), ]))
          VarInd <- sum(apply(data.boot, 2, sd, na.rm = na.rm)^2)
          (p / (p - 1)) * (1 - (VarInd / VarTot))
        } else {
          mat <- if (na.rm) cor(data.boot, use = "complete.obs") else cor(data.boot)
          ave.rho <- mean(mat[upper.tri(mat)])
          (p * ave.rho) / (1 + (p - 1) * ave.rho)
        }
      }
      out$ci <- quantile(T.boot, probs = probs)
      out$probs <- probs
      out$B <- B
    }
    class(out) <- "cronbachAlpha"
    out
  }

#######################################################################
#' @name biserial_cor
#' @title Biserial Correlation
#' @description Point-Biserial correlation coefficient is a correlation coefficient
#' used when one variable is continuous and the other variable is dichotomous.
#' Taken from ltm::biserial_cor
#' @usage biserial_cor(x,y, use = c("all.obs", "complete.obs"), level = 1)
#' @param x a numeric vector representing the continuous variable. 
#' @param y a numeric vector representing the dichotomous variable.
#' @param use is a option for the use of missing values. 
#' @param level which level of y to use.
#' @details It is calculated by applying the Pearson correlation coefficient to the case
#' where one of the variables has dichotomous nature.\cr
#' It is calculated as \deqn{r_{xy} = (\bar{x}_p - \bar{x}_q / S_x)*\sqrt{pq}}
#' Where p  is the proportion of subjects with one of the two possible values of the 
#' variable Y, q is the proportion of subjects with the other possible value, 
#' \deqn{\bar{x}_p} and \deqn{\bar{x}_q} is the average X subjects whose proportion is p and q respectively,
#' and \deqn{S_x} is the standard deviation of all subjects X.
#' This function was adapted from ltm_1.0 package.
#' @examples
#' \dontrun{
#' data <- simulate_dichotomous(size.cluster = c(10),sample.size=1000)
#' biserial_cor(rowSums(data$data), data$data[,1])
#' }
#' @return The value of the point-biserial correlation.
#' @references U.Olsson, F.Drasgow, and N.Dorans (1982). The polyserial correlation coefficient. Psychometrika, 47:337-347.
#' @references Cox. N.R. (1974). Estimation of the Correlation between a Continuous and a Discrete Variable. Biometrics, 30:171-178.
#' @export
biserial_cor <- function (x, y, use = c("all.obs", "complete.obs"), level = 1) 
{
  if (!is.numeric(x)) 
    stop("'x' must be a numeric variable.\n")
  y <- as.factor(y)
  if (length(levs <- levels(y)) > 2) 
    stop("'y' must be a dichotomous variable.\n")
  if (length(x) != length(y)) 
    stop("'x' and 'y' do not have the same length")
  use <- match.arg(use)
  if (use == "complete.obs") {
    cc.ind <- complete.cases(x, y)
    x <- x[cc.ind]
    y <- y[cc.ind]
  }
  ind <- y == levs[level]
  diff.mu <- mean(x[ind]) - mean(x[!ind])
  prob <- mean(ind)
  diff.mu * sqrt(prob * (1 - prob))/sd(x)
}

#######################################################################
#' @name gutt
#' @title Guttman's Lambda
#' @description Six Lower limits of reliability coefficients are presented. 
#' @usage gutt(data)
#' @param data a matrix or a Dataframe that holds the test response data
#' @details Let \deqn{S_j^2} the variances over persons of the n items in the test, and
#' \eqn{S_t^2} the variance over persons of the sum of the items.
#' The firt estimate \eqn{\lambda_1} can be computed from \eqn{L_1 = 1 - (sum{s_j^2}/S_t^2)}
#' Let \eqn{C_2} the sum of squares of the covariances between items, therefore is
#' the sum of \eqn{n(n-1)/2} terms. The bound \eqn{lambda_2} is computed by \eqn{L_2 = L_1 + (\sqrt{n/n-1 C_2}/S_t^2)} 
#' The third lower bound \eqn{\lambda_3} is a modification of \eqn{\lambda_1}, it is computed
#' from the \eqn{L_3 = n/(n-1) L_1}
#' Fourth lower bound \eqn{\lambda_4} has been interpreted as the greatest split half reliability,
#' and requires that the test be scored as twohalves. It is calculated from 
#' \deqn{L_4 = 2(1 - (s_a^2 + s_b^2)/s_t^2)} where \eqn{S_a^2} and \eqn{S_b^2} are the respectives variances
#' of the two parts for the single trial. 
#' For the fifth lower bound \eqn{lambda_5}, let \eqn{C_{2j}} be the sum of the squares of the
#' covariances of item j with the remaining \eqn{n-1} items, and let \eqn{\bar{C}_2} be the largest of
#' the \deqn{C_{2j}}. Then the coefficient can be computed from \eqn{L_5 = L_1 + (2\sqrt{bar{C}_2})/S_t^2}
#' The final bound is based on multiple correlation, let \eqn{e_j^2} be the variance of the errors
#' of estimate of item j from its linear multiple regression on the remaining n-1 items. Then
#' \eqn{lambda_6} can be computed from \eqn{L_6 = 1 - (\sum{e_j^2})/S_t^2} 
#' @return The six coefficients Guttman for the test.
#' @references Guttman, L. (1945). A basis for analyzing test-retest reliability. Psychometrika, 10(4), 255-282.
#' @export
gutt <- function(data){
  test = data
  r <- cov(test)
  n <- ncol(r)
  n.obs <- nrow(r)
  st <- sum(r)
  
  L1 <- 1 - (sum(apply(test, 2, var))/st)
  L2 <- L1 + (sqrt((sum(r^2)/2)*(n/(n-1)))/st)
  L3 <- n/(n-1)*L1
  
  ### 
  
  xy <- matrix(ncol = 2, nrow = n/2)
  r.split <- data.frame(r)
  r.split <- r
  r.split[upper.tri(r.split, diag = TRUE)] <- -999999
  
  for (i in 1:(n/2)) {
    x.m <- which(r.split == max(r.split), arr.ind = TRUE)[1,]
    xy[i, 1] <- x.m[1]
    xy[i, 2] <- x.m[2]
    r.split[(x.m[1]), ] <- -999999
    r.split[, (x.m[1])] <- -999999
    r.split[, (x.m[2])] <- -999999
    r.split[(x.m[2]), ] <- -999999
  }
  A <- xy[, 1]
  B <- xy[, 2]
  lf <- which(1:n %in% c(A, B) == FALSE)
  if (length(c(A, B)) != n) {
    B <- c(B, lf)
  }
  
  
  t1t <- matrix(rep(NA, n), nrow = 1)
  t1t[A] <- 1
  t1t[B] <- 0
  t2 <- (t(t1t) - 1) * -1
  onerow <- matrix(rep(1, n), nrow = 1)
  onevector <- t(onerow)
  
  L4 <- (4 * (t1t %*% r %*% t2))/(onerow %*% r) %*% onevector
  L4 <- as.numeric(L4)
  ## 
  
  c2j <- (r^2)            
  c2 <- NULL
  
  for(j in 1:ncol(test)){
    c2[j] <- sum(c2j[j,-j])
  }
  
  L5 <- L1 + ((2*sqrt(max(c2)))/st)
  
  ##
  
  smc <- 1 - 1/diag(solve(r))
  t1 <- matrix(rep(1, n), ncol = 1)
  t1t <- t(t1)
  
  L6 <- 1 - sum(1 - smc)/(sum(r))
  
  ### 
  
  result <- list(lambda1 = L1, lambda2 = L2, lambda3 = L3, lambda4 = L4, lambda5 = L5, 
                 lambda6 = L6)
  return(result)
  
}

#######################################################################
#' @name Yule
#' @title Yule coefficient of correlation 
#' @description The Yule coefficientis a correlation coefficient applied 
#' to dichotomous data. Given a two x two table of counts\cr
#' | a | b | R1 |\cr
#' | c | d | R1 |\cr
#' |---|---|----|\cr
#' |C1 | C2| n  |\cr
#' or a vector c(a,b,c,d) of frequencies. Taken from psych::Yule.
#' @usage Yule(x, Y = FALSE)
#' @param x a 1 x 4 vector or a matrix 2 x 2 of frequencies.
#' @param Y if Y is true return Yule's Y coefficient of colligation.
#' @details The coefficient of Yule is calculated from (ad - bc)/(ad + bc).
#' This is the number of pairs in agreement (ad) - the number in disagreement (bc) 
#' over the total number of paired observations. 
#' @return the value of the Yule Q coefficient.
#' @references Yule, G.U. (1912). On the methods of measuring the association between two attributes. Journal of the Royal Statistical Society, 75, 579-652.
#' @references Warrens, Matthijs (2008), On Association Coefficients for 2x2 Tables and Properties That Do Not Depend on the Marginal Distributions. Psychometrika, 73, 777-789.
#' @export
Yule <- function (x, Y = FALSE) 
{
  stopifnot(prod(dim(x)) == 4 || length(x) == 4)
  if (is.vector(x)) {
    x <- matrix(x, 2)
  }
  a <- x[1, 1]
  b <- x[1, 2]
  c <- x[2, 1]
  d <- x[2, 2]
  if (Y) {
    Yule <- (sqrt(a * d) - sqrt(b * c))/(sqrt(a * d) + sqrt(b * c))
  }
  else {
    Yule <- (a * d - b * c)/(a * d + b * c)
  }
  return(Yule)
}



#######################################################################
#' @name polyserial_cor
#' @title Polyserial correlation
#' @description Polyserial correlation coefficient is a correlation coefficient
#' used when one variable is continuous and the other variable is dichotomous.
#' Taken and adapted from psych::polyserial
#' @usage polyserial_cor(x,y)
#' @param x a numeric vector representing the continuous variable. 
#' @param y a numeric vector representing the dichotomous variable.
#' @details The coefficient is calculated from \deqn{\rho = r_{xy} * \sqrt{(n - 1)/n} * s_y/\sum{\phi(\tau)}} 
#' where \deqn{r_{xy}} is the coefficient of correlation of Pearson coefficient, S_y is the 
#' standard deviation of Y, and \deqn{\phi(\tau)} are the ordinates of the normal curve at 
#' the normal equivalent of the cut point boundaries between the item responses.
#' This function was adapted from ltm_1.0 package.
#' @examples
#' \dontrun{
#'  x <- rnorm(100)
#'  y <- sample(1:5,100,replace=TRUE)
#'  cor(x, y) 
#'  polyserial_cor(x, y) 
#' }
#' @return the value of the polyserial correlation.
#' @references U.Olsson, F.Drasgow, and N.Dorans (1982). The polyserial correlation coefficient. Psychometrika, 47:337-347.
#' @export
polyserial_cor <- function (x, y) 
{
  min.item <- min(y, na.rm = TRUE) 
  max.item <- max(y, na.rm = TRUE) 
  n.var <- 1  
  n.cases <- length(y) 
  
  dummy <- matrix(rep(min.item:max.item, n.var), ncol = n.var, nrow=length(y))
  colnames(dummy) <- names(y) 
  xdum <- cbind(y, dummy) 
  frequency <- apply(xdum, 2, table) 
  frequency <- t(frequency - 1) 
  responses <- rowSums(frequency) 
  frequency <- frequency/responses 
  frequency <- t(apply(frequency, 1, cumsum)) 
  len <- dim(frequency)[2] 
  tau <- dnorm(qnorm(frequency[, -len, drop = FALSE])) 
  stau <- rowSums(tau) 
  rxy <- cor(x, y, use = "pairwise") 
  sdy <- sd(y) 
  rps <- t(rxy) * sqrt((n.cases - 1)/n.cases) * sdy/stau[1] 
  rps[rps > 1] <- 1
  rps[rps < -1] <- -1
  return(rps)
}

#######################################################################
#' @name an_parallel
#' @title Parallel Analysis
#' @description Performs Horn's parallel analysis for a principal component. 
#' @usage an_parallel(x = NA, iterations = 0, centile = 0, 
#'  seed = 0,mat = NA, n = NA)
#' @param x a matrix or a dataframe that holds the test response data.
#' @param iterations how many iterations are performed during the analysis. Amount of random datasets produced.
#' @param centile a number between 1 and 99 indicating the centile used in estimating bias.
#' @param seed seed used for the random data generation.
#' @param mat a supplied correlation mtrix to use instead of x. If used, then n must be also supplied.
#' @param n the number of observations. Required when a correlation matrix is supplied 
#' with the mat option.
#' @details Implementation of Horn's (1965) tecnique for evaluating the components retained
#' in a principle component analysis (PCA). This procedure is an adaptation of the
#' paran function in the Paran Package.
#' @return Retained Components a scalar integer representing the number of retained components.
#' @return Adjusted Eigenvalues a vector of the estimated adjusted eigenvalues.
#' @return Unadjusted Eigenvalues a vector of the eigenvalues of the observed data from
#' an unrotated principal component analysis.
#' @return Bias a vector of the estimated bias of the unadjusted eigenvalues 
#' @examples
#' \dontrun{
#'  data <- simulate_dichotomous(size.cluster = c(10),sample.size=1000)
#'  an_parallel(data$data, iterations = 100, centile = 99, seed = 12)
#' }
#' @references John L. Horn (1965). A rationale and test for the number of factors 
#' in factor analysis. Psychometrika, Volume 30, Number 2, Page 179.
#' @references Dinno A. 2009. Exploring the Sensitivity of Horn's Parallel Analysis to the 
#' Distributional Form of Simulated Data. Multivariate Behavioral Research. 44(3): 362-388
#' @export

an_parallel <- function(x = NA, iterations = 0, centile = 0, seed = 0, 
                        mat = NA, n = NA) 
{
  if (!is.na(mat[[1]][1]) & !is.na(x[[1]][1])) {
    stop("\nYou must supply either x or mat but not both.")
  }
  if (is.na(mat[[1]][1]) & !is.na(x[[1]][1])) {
    P <- length(as.matrix(x)[1, ]) #numero de variables en una matriz de datos
  }
  if (!is.na(mat[[1]][1]) & is.na(x[[1]][1])) {
    P <- length(mat[1, ]) #numero de variables en una matriz de correlacion
  }
  if (!is.na(mat[[1]][1])) {
    if (length(mat[1, ]) != length(mat[, 1])) {
      stop("\nThe matrix provided with the mat argument is not a correlation matrix.\n")
    }
    if (is.element(FALSE, (diag(mat) == rep(1, P)))) {
      stop("\nThe matrix provided with the mat argument is not a correlation matrix.\nParallel analysis is not compatible with the eigendecomposition of a covariance matrix.")
    }
  }
  
  if (!is.na(x[[1]][1])) {
    N <- length(as.matrix(x)[, 1]) 
    if (!is.na(n)) {
      warning("\nThe n argument is only for use with the mat argument. Ignoring n argument.")
    }
    
    R <- cor(x) 
    cat("\nUsing eigendecomposition of correlation matrix.")
  }
  
  if (!is.na(mat[[1]][1])) {
    if (!is.na(mat[[1]][1]) & is.na(n)) {
      stop("\nYou must also provide the sample size when using the matrix argument.")
    }
    
    
    N <- n 
    R <- mat 
    cat("\nUsing eigendecomposition of provided correlation matrix.")
  }
  
  centile <- round(centile) 
  
  if (centile > 99 | centile < 0) {
    stop("\nYou must specify a centile value between 1 and 99.\n(Specifying centile 0 will use the mean.)")
  }
  
  eigenvalues <- eigen(R, only.values = TRUE, EISPACK = FALSE)[[1]]
  
  Ev <- eigenvalues
  model <- "component"
  models <- "components"
  Model <- "Component"
  Models <- "Components"
  
  if (iterations < 1) {
    iterations <- 30 * P
  }
  
  SimEvs <- matrix(NA, iterations, P)
  
  for (k in 1:iterations) {
    
    Sim <- matrix(NA, N, P)
    
    if (seed != 0) {
      set.seed(seed * k)
    }
    
    Sim <- matrix(rnorm(N * P), N, P)
    
    eigenvalues <- eigen(cor(Sim), only.values = TRUE, EISPACK = FALSE)[[1]]
    Evs <- eigenvalues
    SimEvs[k, ] <- Evs
  }
  
  cat("\n\nResults of Horn's Parallel Analysis for ", model, 
      " retention\n", sep = "")
  if (iterations == 1) {
    if (centile == 0) {
      cat("1 iteration, using the mean estimate", "\n", 
          sep = "")
    }
    if (centile != 0) {
      cat("1 iteration, using the ", centile, " centile estimate", 
          "\n", sep = "")
    }
  }
  if (iterations > 1) {
    if (centile == 0) {
      cat(iterations, " iterations, using the mean estimate", 
          "\n", sep = "")
    }
    if (centile != 0 & centile != 50) {
      cat(iterations, " iterations, using the ", centile, 
          " centile estimate", "\n", sep = "")
    }
    if (centile == 50) {
      cat(iterations, " iterations, using the ", centile, 
          " centile (median) estimate", "\n", sep = "")
    }
  }
  cat("\n--------------------------------------------------", 
      "\n")
  cat(Model, "  Adjusted    Unadjusted    Estimated", "\n")
  cat("            Eigenvalue  Eigenvalue    Bias", "\n")
  cat("--------------------------------------------------", 
      "\n")
  
  RndEv = c(1:P) * NA
  
  if (centile > 0) {
    for (p in 1:P) {
      RndEv[[p]] <- quantile(SimEvs[, p], probs = centile/100)[[1]]
    }
  }
  
  if (centile == 0) {
    for (p in 1:P) {
      RndEv[[p]] <- mean(SimEvs[, p])
    }
  }
  
  if (Ev[[1]] < 1 | RndEv[[1]] < 1) {
    cat("No components passed.", "\n")
    cat("--------------------------------------------------", 
        "\n")
    stop
  }
  
  Bias <- rep(0, P)
  AdjEv <- rep(1, P)
  
  for (p in 1:P) {
    Bias[p] <- RndEv[p] - 1
    AdjEv[p] <- Ev[[p]] - Bias[p]
  }
  
  y <- NA
  
  for (x in 1:P) {
    y <- x
    
    if (AdjEv[x] <= 1) {
      y <- x - 1
      retained <- y
      break
    }
    
  }
  
  y <- P
  
  for (x in 1:y) {
    if (AdjEv[x] >= 0) {
      AdjSpace = " "
    }
    if (AdjEv[x] < 0) {
      AdjSpace = ""
    }
    if (Ev[[x]] >= 0) {
      EvSpace = " "
    }
    if (Ev[[x]] < 0) {
      EvSpace = ""
    }
    if (Bias[x] >= 0) {
      BiasSpace = " "
    }
    if (Bias[x] < 0) {
      BiasSpace = ""
    }
    if (x > 9) {
      xPad = ""
    }
    if (x <= 9) {
      xPad = " "
    }
    AdjFPad = "   "
    if (round(AdjEv[x]) >= 10) {
      AdjFPad = "  "
    }
    if (round(AdjEv[x]) >= 100) {
      AdjFPad <- " "
    }
    SN <- 8
    if (abs(AdjEv[x]) >= 10) {
      SN <- 9
    }
    if (abs(AdjEv[x]) >= 100) {
      SN >= 10
    }
    if (AdjEv[x] < 0) {
      SN <- SN + 1
    }
    EvFPad = "   "
    if (round(Ev[[x]]) >= 10) {
      EvFPad = "  "
    }
    if (round(Ev[[x]]) >= 100) {
      EvFPad = " "
    }
    EvSN <- 8
    if (abs(Ev[[x]]) >= 10) {
      EvSN <- 9
    }
    if (abs(Ev[[x]]) >= 100) {
      EvSN <- 10
    }
    if (abs(Ev[[x]]) >= 5e-07) {
      EvZPad <- ""
    }
    if (abs(Ev[[x]]) < 5e-07) {
      Ev[[x]] <- 0
      EvZPad <- ".000000"
    }
    BiasSN <- 8
    if (Bias[x] >= 10) {
      BiasSN <- 9
    }
    if (Bias[x] >= 100) {
      BiasSN >= 10
    }
    cat(x, xPad, "      ", AdjFPad, AdjSpace, strtrim(AdjEv[x], 
                                                      SN), EvFPad, EvSpace, strtrim(Ev[[x]], EvSN), 
        EvZPad, "     ", BiasSpace, strtrim(Bias[x], 
                                            BiasSN), "\n", sep = "")
  }
  cat("--------------------------------------------------", 
      "\n")
  cat("\nAdjusted eigenvalues > 1 indicate dimensions to retain.\n(", 
      retained, " ", models, " retained)\n\n", sep = "")
  col = c("black", "red", "blue")
  
  AdjEvCol = col[1]
  EvCol = col[2]
  RndEvCol = col[3]
  AdjEvLty = 1
  EvLty = 1
  RndEvLty = 1
  
  par(yaxs = "i", xaxs = "i", lab = c(P, ceiling(max(AdjEv[1], Ev[1], RndEv[1])), 2))
  plot.default(c(1:P), RndEv, type = "o", main = "Parallel Analysis", 
               xlab = "Components", ylab = "Eigenvalues", pch = 20, 
               col = RndEvCol, lty = 1, lwd = 1, 
               xlim = c(0.5,P + 0.5), ylim = c(min(AdjEv, Ev, RndEv) - 0.5, 
                                               ceiling(max(AdjEv[[1]], Ev[[1]], RndEv[[1]]))))
  abline(h = 1, col = "grey", lwd = 0.5)
  points(c(1:P), AdjEv, type = "o", col = AdjEvCol, lty = 1, 
         pch = 21, bg = "white", lwd = 1)
  points(c(1:P), Ev, type = "o", col = EvCol, lty = 1, 
         pch = 20, lwd = 1)
  if (retained >= 1) {
    points(c(1:retained), AdjEv[1:retained], type = "p", 
           pch = 19, col = AdjEvCol, lty = 1, lwd = 1)
  }
  legend("topright", legend = c("Adjusted Ev (retained)", 
                                "Adjusted Ev (unretained)", "Unadjusted Ev", 
                                "Random Ev"), col = c(AdjEvCol, AdjEvCol, EvCol, 
                                                      RndEvCol), pch = c(19, 21, 20, 20), lty = c(1, 1, 1, 1))
  
  invisible(list(Retained = retained, AdjEv = AdjEv, Ev = Ev, 
                 RndEv = RndEv, Bias = Bias, SimEvs = SimEvs))
}