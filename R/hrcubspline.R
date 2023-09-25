.bootrcsHR_knot_more_than_3 <- function(data, idx , x , model , var1 , var2){
  form_str <- deparse(formula(model), width.cutoff = 500)
  has_vector_rcs <- "^.+ ~ .+\\s*\\*?\\s*rcs\\([^,]+, c\\(.*\\)\\).*$"
  matches_vector_rcs <- stringr::str_detect(form_str, has_vector_rcs)
  df <- data[idx,]
  mycall <- model$call
  mycall <- pryr::modify_call(mycall , list(data=quote(df)))
  if("cph" %in% class(model) | "lrm" %in% class(model)){
    myformula <- model$sformula
  } else {
    myformula <- model$formula
  }

  mymodel <- eval(mycall)
  coefMod <- coef(mymodel)

  if("cph" %in% class(model) | "lrm" %in% class(model)){
    separator <- " * "
    k <- mymodel$Design$parms[[var2]]
  } else {
    if (matches_vector_rcs) {
      matches <- stringr::str_match(form_str, "rcs\\((.*)\\)")
      rcs_content <- matches[, 2]
      split_str <- strsplit(rcs_content, ",")[[1]]
      var <- trimws(split_str[1])
      knots_str <- paste(split_str[-1], collapse = ",")
      knots <- eval(parse(text = paste("c(", knots_str, ")", sep="")))
      k <- knots
      separator <- ":"
      rcsTerm <- grep(":" , grep("rcs\\(" , attributes(mymodel$terms)$term.labels , value = TRUE) , invert = TRUE , value = TRUE)
      names(coefMod) <- gsub(rcsTerm , "" , names(coefMod) , fixed = TRUE)
    }   else {
      separator <- ":"
      # need to recreate the knot sequence for object of class coxph
      # remove the rcs part from the names of the variables in case of a class coxph model
      rcsTerm <- grep(":" , grep("rcs\\(" , attributes(model$terms)$term.labels , value = TRUE) , invert = TRUE , value = TRUE)
      names(coefMod) <- gsub(rcsTerm , "" , names(coefMod) , fixed = TRUE)
      indices <- length(grep(paste0("^", var2, "'*$"), names(coefMod)))+1
      k <- attributes(rms::rcs(df[[var2]], indices))$parms
    }
  }


  ### Get Iterative Fractions
  k_num = length(k)
  k_num_counter = length(k)

  fraction_list = list()
  j = 1
  while (k_num_counter >2) {
    temp_frac = (k[k_num] - k[j]) / (k[k_num]-k[k_num-1])
    temp_frac_2 = (k[k_num-1] - k[j]) / (k[k_num]-k[k_num-1])
    fraction_list = c(fraction_list,temp_frac)
    fraction_list = c(fraction_list,temp_frac_2)
    j = j + 1
    k_num_counter = k_num_counter - 1
  }

  ### New Variable Definition
  var_list <- c(var1, var2, paste(var1 , var2 , sep = separator),paste(var2 , var1 , sep = separator))

  for(i in 1:(length(k)-2)){
    var2_mod <- paste0(var2, paste(rep("'", i), collapse = ""))
    var_list <- c(var_list,
                  paste0(var2_mod),
                  paste(var1, var2_mod, sep = separator),
                  paste(var2_mod, var1, sep = separator))
  }


  myvars <- sort(intersect(var_list, names(coefMod)))

  mycoef <- coefMod[  myvars ]
  mycoefWhich <- sort(sapply( myvars , function(v) which( names(coefMod) %in% v )))

  num_ticks <- length(k)-2
  b <- mycoef[ var1 ]
  a <- NULL

  for (i in 0:num_ticks) {

    var_name <- paste0(var2, paste(rep("'", i), collapse = ""))

    if (var_name %in% names(mycoef)) {
      a <- c(a, mycoef[var_name])
    }
  }

  l <- mycoef[setdiff(setdiff(myvars, c(var1, names(a))), c(var1, names(a)))]

  ### Get Iterative Terms

  k_num = length(k)
  k_num_counter = length(k)
  term_list = list()
  j = 1
  frac_index = 1
  while (k_num_counter >2) {
    numer <- vapply(x , function(i) {
      max(i - k[j],0)^3  - (max(i - k[k_num-1],0)^3)*fraction_list[[frac_index]] + (max(i -             k[k_num],0)^3)*fraction_list[[frac_index+1]]
    } , numeric(1))
    denom <- (k[k_num] - k[1])^2
    numDem <- numer/denom

    term_list = c(term_list,numDem)

    k_num_counter = k_num_counter - 1
    frac_index = frac_index + 2
    j = j + 1
  }

  n <- length(term_list)
  elements_per_iteration <- length(x)

  subsets <- list()

  for (i in seq(1, n, by = elements_per_iteration)) {
    end <- min(i + elements_per_iteration - 1, n)
    subsets[[length(subsets) + 1]] <- term_list[i:end]
  }

  sp2 <- vapply(x , function(i) l[1]*i , numeric(1))

  for(i in 2:length(l)) {
    sp2 <- sp2 + l[[i]] * unlist(subsets[i-1])
  }


  HR <- unname(exp( b + sp2))
}


#' Restricted cubic spline interaction HR for more than 3 knots
#'
#' Generate HR values in a Cox model for a 1 unit increase in a variable at
#' specified points of another interacting variable splined with rcs(df >= 3)
#'
#' @param var2values numeric vector of var2 points to estimate
#' @param model model of class cph or coxph. If data is NULL, the function expects to find the data in model$x.
#' @param data data used in the model. If absent, we will attempt to recover the data from the model. Only used for bootstrap and coxph models
#' @param var1 variable that increases by 1 unit from 0.
#' @param var2 variable to spline. var2values belong to var2
#' @param ci calculate 95% CI?
#' @param conf confidence level. Default 0.95
#' @param ci.method confidence interval method. "delta" performs delta method. "bootstrap" performs bootstrapped CI (slower)
#' @param ci.boot.method one of the available bootstrap CI methods from \code{\link[boot]{boot.ci}}. Default percentile
#' @param R number of bootstrap samples if ci.method = "bootstrap". Default 100
#' @param parallel can take values "no", "multicore", "snow" if ci.method = "bootstrap". Default multicore
#' @param ... other parameters for boot
#' @examples
#' library(survival)
#' library(rms)
#' data(cancer)
#' myformula <- Surv(time, status) ~ ph.karno + ph.ecog + rcs(age,4)*sex
#' model <- cph(myformula , data = lung )
#' rcsHR( var2values = 40:80
#'        , model = model , data = lung , var1 ="sex", var2="age"
#'        , ci=TRUE , conf = 0.95 , ci.method = "delta")
#' @return if ci = FALSE, a dataframe with initial values and HR
#' , if ci = TRUE a dataframe with 5 columns, initial values, HR, lower CI, upper CI and SE
#' @importFrom rms cph rcs
#' @importFrom survival coxph
#' @importFrom pryr modify_call
#' @importFrom msm deltamethod
#' @importFrom boot boot boot.ci
#' @importFrom stats vcov coef as.formula qnorm sd formula
#' @importFrom stringr str_detect str_match
#' @importFrom utils tail
#' @export


rcsHR <- function(var2values , model , data=NULL , var1 , var2
                  , ci=TRUE , conf = 0.95 , ci.method = "delta"
                  , ci.boot.method = "perc" , R = 100 , parallel = "multicore" , ...) {

  ### Check correct class for model
  if( !any( c("cph","coxph") %in% class(model) ) ){
    stop("Cubic spline Cox model must be run with rms::cph or survival::coxph")
  }
  if(!is.numeric(var2values)){
    stop("var2values must be a numeric vector")
  }
  x <- var2values

  if(missing(data)){
    if(is.null(model$x)){
      if("bootstrap" %in% ci.method || class(model)[1] == "coxph"){
        stop("Missing data")
      }
    } else {
      data <- model$x
    }
  }

  if(!all(c(var1,var2) %in% colnames(data) )){
    stop("var1 or var2 not present in the data")
  }

  form_str <- deparse(formula(model), width.cutoff = 500)
  has_vector_rcs <- "^.+ ~ .+\\s*\\*?\\s*rcs\\([^,]+, c\\(.*\\)\\).*$"
  matches_vector_rcs <- stringr::str_detect(form_str, has_vector_rcs)

  coefMod <- coef(model)

  if("cph" %in% class(model) | "lrm" %in% class(model)){
    separator <- " * "
    k <- model$Design$parms[[var2]]
  } else {
    if (matches_vector_rcs) {
      matches <- stringr::str_match(form_str, "rcs\\((.*)\\)")
      rcs_content <- matches[, 2]
      split_str <- strsplit(rcs_content, ",")[[1]]
      var <- trimws(split_str[1])
      knots_str <- paste(split_str[-1], collapse = ",")
      knots <- eval(parse(text = paste("c(", knots_str, ")", sep="")))
      k <- knots
      separator <- ":"
      rcsTerm <- grep(":" , grep("rcs\\(" , attributes(model$terms)$term.labels , value = TRUE) , invert = TRUE , value = TRUE)
      names(coefMod) <- gsub(rcsTerm , "" , names(coefMod) , fixed = TRUE)
    }   else {
      separator <- ":"
      # need to recreate the knot sequence for object of class coxph
      # remove the rcs part from the names of the variables in case of a class coxph model
      rcsTerm <- grep(":" , grep("rcs\\(" , attributes(model$terms)$term.labels , value = TRUE) , invert = TRUE , value = TRUE)
      names(coefMod) <- gsub(rcsTerm , "" , names(coefMod) , fixed = TRUE)
      indices <- length(grep(paste0("^", var2, "'*$"), names(coefMod)))+1
      k <- attributes(rms::rcs(data[[var2]], indices))$parms
    }
  }

  ### Get Iterative Fractions
  k_num = length(k)
  k_num_counter = length(k)

  fraction_list = list()
  j = 1
  while (k_num_counter >2) {
    temp_frac = (k[k_num] - k[j]) / (k[k_num]-k[k_num-1])
    temp_frac_2 = (k[k_num-1] - k[j]) / (k[k_num]-k[k_num-1])
    fraction_list = c(fraction_list,temp_frac)
    fraction_list = c(fraction_list,temp_frac_2)
    j = j + 1
    k_num_counter = k_num_counter - 1
  }

  ### New Variable Definition
  var_list <- c(var1, var2, paste(var1 , var2 , sep = separator),paste(var2 , var1 , sep = separator))

  for(i in 1:(length(k)-2)){
    var2_mod <- paste0(var2, paste(rep("'", i), collapse = ""))
    var_list <- c(var_list,
                  paste0(var2_mod),
                  paste(var1, var2_mod, sep = separator),
                  paste(var2_mod, var1, sep = separator))
  }


  myvars <- sort(intersect(var_list, names(coefMod)))

  mycoef <- coefMod[  myvars ]
  mycoefWhich <- sort(sapply( myvars , function(v) which( names(coefMod) %in% v )))

  num_ticks <- length(k)-2
  b <- mycoef[ var1 ]
  a <- NULL

  for (i in 0:num_ticks) {

    var_name <- paste0(var2, paste(rep("'", i), collapse = ""))

    if (var_name %in% names(mycoef)) {
      a <- c(a, mycoef[var_name])
    }
  }

  l <- mycoef[setdiff(setdiff(myvars, c(var1, names(a))), c(var1, names(a)))]

  ### Get Iterative Terms

  k_num = length(k)
  k_num_counter = length(k)
  term_list = list()
  j = 1
  frac_index = 1
  while (k_num_counter >2) {
    numer <- vapply(x , function(i) {
      max(i - k[j],0)^3  - (max(i - k[k_num-1],0)^3)*fraction_list[[frac_index]] + (max(i -             k[k_num],0)^3)*fraction_list[[frac_index+1]]
    } , numeric(1))
    denom <- (k[k_num] - k[1])^2
    numDem <- numer/denom

    term_list = c(term_list,numDem)

    k_num_counter = k_num_counter - 1
    frac_index = frac_index + 2
    j = j + 1
  }

  n <- length(term_list)
  elements_per_iteration <- length(x)

  subsets <- list()

  for (i in seq(1, n, by = elements_per_iteration)) {
    end <- min(i + elements_per_iteration - 1, n)
    subsets[[length(subsets) + 1]] <- term_list[i:end]
  }

  sp2 <- vapply(x , function(i) l[1]*i , numeric(1))

  for(i in 2:length(l)) {
    sp2 <- sp2 + l[[i]] * unlist(subsets[i-1])
  }


  HR <- unname(exp( b + sp2))

  if(ci){
    alpha <- qnorm( 1 - (1-conf)/2)
    if(ci.method == "delta"){
      # This creates a vector like x1 , x2 , x3 , x7 , x8
      # that tells you the position of the regressor as it appears in the model
      xNum <- paste0("x" , mycoefWhich)
      last_terms = tail(xNum,k_num-1)
      vcovMod <- vcov(model)
      xNum_var1 <- paste0("x" , mycoefWhich[var1])
      formula_terms <- lapply(1:length(subsets) , function(k) {
        subsets[k]
      })
      names(formula_terms) <- paste0("formula_terms" , 1:length(subsets))
      # for (k in 1:length(subsets)) {
      #   assign(paste0("formula_terms", k), subsets[k], envir = .GlobalEnv)
      # }


      HRci <- t(vapply( seq_len(length(x)) , function(i) {
        x_i <- x[i]
        # numDems <- length(grep("formula_terms", ls(.GlobalEnv)))
        numDems <- length(formula_terms)

        formula_parts <- c(paste0("~(", xNum_var1, " + ", last_terms[1], "*(", x_i, ")"))

        for (j in 1:numDems) {
          if (j!=numDems){
            # numDem_i <- unlist(get(paste0("formula_terms", j)))[i]
            numDem_i <- unlist(formula_terms[paste0("formula_terms", j)])[i]
            formula_parts <- c(formula_parts, paste0("+ ", last_terms[j+1], "*(", numDem_i, ")"))}
          else{
            # numDem_i <- unlist(get(paste0("formula_terms", j)))[i]
            numDem_i <- unlist(formula_terms[paste0("formula_terms", j)])[i]
            formula_parts <- c(formula_parts, paste0("+ ", last_terms[j+1], "*(", numDem_i, "))"))
          }
        }
        formula <- paste(formula_parts, collapse = " ")

        SE <- NULL
        try(SE<-msm::deltamethod(as.formula(formula), coefMod, vcovMod), silent = TRUE)
        if(is.null(SE)){
          return(c(HR[i] , NA , NA , NA))
        }
        up<-exp(log(HR[i])+alpha*SE)
        lo<-exp(log(HR[i])-alpha*SE)
        c(HR[i] , lo , up , SE)
      } , numeric(4)))
      HRci <- cbind( Value = x , HRci)
      rownames(HRci) <- x
      colnames(HRci) <- c("Value" , "HR" , "CI_L" , "CI_U" , "SE")
      HRci <- as.data.frame(HRci)
      # class(HRci) <- c("HR" , class(HRci))
      return(HRci)
    } else if(ci.method == "bootstrap"){
      if(missing(parallel)){
        parallel <- "multicore"
      }
      if(missing(R)){
        R <- 100
      }
      myBoot <- boot::boot(data = data, statistic = .bootrcsHR_knot_more_than_3, x = x , model = model
                           , R = R , parallel = parallel, var1 = var1 , var2 = var2 )
      SE <- apply(myBoot$t , 2 , sd)
      HRci <- t(vapply( seq_len(length(x)) , function(idx) {
        bci <- boot::boot.ci(boot.out = myBoot,  index = idx , type = ci.boot.method , conf = conf)
        if(ci.boot.method %in% "norm"){
          c(bci$t0 , bci$normal[2] , bci$normal[3] , SE = SE[idx])
        } else {
          c(bci$t0 , bci[[4]][4] , bci[[4]][5] , SE = SE[idx])
        }
      } , numeric(4)))
      HRci <- cbind(x , HRci)
      colnames(HRci) <- c("Value" , "HR" , "CI_L" , "CI_U" , "SE")
      rownames(HRci) <- x
      HRci <- as.data.frame(HRci)
      # class(HRci) <- c("HR" , class(HRci))
      return(HRci)
    } else {
      stop("Only delta and bootstrap are valid CI methods")
    }
  } else {
    HR <- data.frame(Value = x , HR = HR)
    # class(HR) <- c("HR",class(HR))
    return(HR)
  }
}
