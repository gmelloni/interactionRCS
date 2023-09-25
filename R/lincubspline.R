.bootrcsLIN_more_than_3 <- function(data, idx , x , model , var1 , var2){
  form_str <- deparse(formula(model), width.cutoff = 500)
  has_vector_rcs <- "^.+ ~ .+\\s*\\*?\\s*rcs\\([^,]+, c\\(.*\\)\\).*$"
  matches_vector_rcs <- stringr::str_detect(form_str, has_vector_rcs)
  df <- data[idx, ]
  mycall <- model$call
  mycall <- pryr::modify_call(mycall , list(data=quote(df)))
  mymodel <- eval(mycall)
  coefMod <- coef(mymodel)


  if ("Glm" %in% class(mymodel)){
    k <- mymodel$Design$parms[[var2]]
    separator <- " * "
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
    }   else{
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


  LIN <- unname(b + sp2)

}



#' Restricted cubic spline interaction linear estimates for more than 3 knots
#'
#' Generate estimates in a linear model for a 1 unit increase in a variable at
#' specified points of another interacting variable splined with rcs(df >= 3)
#'
#' @param var2values numeric vector of var2 points to estimate
#' @param model model of class rms::Glm or stats::glm family gaussian. If data is NULL, the function expects to find the data in model$x.
#' @param data data used in the model. If absent, we will attempt to recover the data from the model. Only used for bootstrap
#' @param var1 variable that increases by 1 unit from 0
#' @param var2 variable to spline. var2values belong to var2
#' @param ci calculate 95% CI?
#' @param conf confidence level. Default 0.95
#' @param ci.method confidence interval method. "delta" performs delta method. "bootstrap" performs bootstrapped CI (slower)
#' @param ci.boot.method one of the available bootstrap CI methods from \code{\link[boot]{boot.ci}}. Default percentile
#' @param R number of bootstrap samples if ci.method = "bootstrap". Default 100
#' @param parallel can take values "no", "multicore", "snow" if ci.method = "bootstrap". Default multicore
#' @param ... other parameters for boot
#' @examples
#' library(rms)
#' library(mlbench)
#' data(PimaIndiansDiabetes)
#' # Recode diabetes as 0/1
#' PimaIndiansDiabetes$diabetes <- ifelse(PimaIndiansDiabetes$diabetes=="pos" , 1 , 0)
#' myformula <- glucose ~ mass + diabetes * rcs(age, 4)
#' model <- glm(myformula , data = PimaIndiansDiabetes , family="gaussian")
#' # Show the effect on glucose of being diabetic at age 20 to 80
#' rcsLIN( var2values = 20:80
#'        , model = model , data = PimaIndiansDiabetes , var1 ="diabetes", var2="age"
#'        , ci=TRUE , conf = 0.95 , ci.method = "delta")
#' @return if ci = FALSE, a dataframe with initial values and linear estimates
#' , if ci = TRUE a dataframe with 5 columns, initial values, linear estimates, lower CI, upper CI and SE
#' @importFrom rms Glm
#' @importFrom pryr modify_call
#' @importFrom msm deltamethod
#' @importFrom boot boot boot.ci
#' @importFrom stats vcov coef as.formula qnorm sd glm formula
#' @importFrom stringr str_detect str_match
#' @importFrom utils tail
#' @export

rcsLIN <- function(var2values , model , data=NULL , var1 , var2
                   , ci=TRUE , conf = 0.95 , ci.method = "delta"
                   , ci.boot.method = "perc" , R = 100 , parallel = "multicore" , ...) {
  # Check correct class for model
  if( !any( c("Glm","glm","lm") %in% class(model) ) ){
    stop("Cubic spline Logistic model must be run with rms::Glm or stats::glm or stats::lm")
  }
  # Check correct family
  if(!( "gaussian" %in% model$family$family | "quasi" %in% model$family$family | identical(class(model), "lm"))){
    stop("model of class glm but not family gaussian nor quasi")
  } else {
    if("glm" %in% class(model) && !"Glm" %in% class(model)){
      modelClass <- "glm"
    } else {
      modelClass <- "Glm"
    }
  }
  if(!is.numeric(var2values)){
    stop("var2values must be a numeric vector")
  }
  x <- var2values
  if(missing(data)){
    if(is.null(model$x)){
      stop("Missing data")
    } else {
      data <- model$x
    }
  }
  if(!all(c(var1,var2) %in% colnames(data) )){
    stop("var1 or var2 not present in the data")
  }
  # Check that var1 is a 0/1, if not check if the mean is 0
  # if(!all(data[[var1]] %in% c(0,1,NA))){
  #   if(!isTRUE(all.equal(mean(data[[var1]] , na.rm =TRUE),0))){
  #     warning("var1 is not centered on 0 nor a 0/1 variable, results are always reported for a 0 to 1 change in var1.")
  #   }
  # }

  coefMod <- coef(model)
  form_str <- deparse(formula(model), width.cutoff = 500)
  has_vector_rcs <- "^.+ ~ .+\\s*\\*?\\s*rcs\\([^,]+, c\\(.*\\)\\).*$"
  matches_vector_rcs <- stringr::str_detect(form_str, has_vector_rcs)

  if ("Glm" %in% class(model)){
    k <- model$Design$parms[[var2]]
    separator <- " * "
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
    }   else{
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


  LIN <- unname( b + sp2)

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


      LINci <- t(vapply( seq_len(length(x)) , function(i) {
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
          return(c(LIN[i] , NA , NA , NA))
        }
        up<-LIN[i]+alpha*SE
        lo<-LIN[i]-alpha*SE
        c(LIN[i] , lo , up , SE)
      } , numeric(4)))
      LINci <- cbind( Value = x , LINci)
      rownames(LINci) <- x
      colnames(LINci) <- c("Value" , "LIN" , "CI_L" , "CI_U" , "SE")
      LINci <- as.data.frame(LINci)
      # class(LINci) <- c("LIN" , class(LINci))
      return(LINci)
    } else if(ci.method == "bootstrap"){
      if(missing(parallel)){
        parallel <- "multicore"
      }
      if(missing(R)){
        R <- 100
      }
      myBoot <- boot::boot(data = data, statistic = .bootrcsLIN_more_than_3, x = x , model = model
                           , R = R , parallel = parallel, var1 = var1 , var2 = var2 )
      SE <- apply(myBoot$t , 2 , sd)
      LINci <- t(vapply( seq_len(length(x)) , function(idx) {
        bci <- boot::boot.ci(boot.out = myBoot,  index = idx , type = ci.boot.method , conf = conf)
        if(ci.boot.method %in% "norm"){
          c(bci$t0 , bci$normal[2] , bci$normal[3] , SE = SE[idx])
        } else {
          c(bci$t0 , bci[[4]][4] , bci[[4]][5] , SE = SE[idx])
        }
      } , numeric(4)))
      LINci <- cbind(x , LINci)
      colnames(LINci) <- c("Value" , "LIN" , "CI_L" , "CI_U" , "SE")
      rownames(LINci) <- x
      LINci <- as.data.frame(LINci)

      return(LINci)
    } else {
      stop("Only delta and bootstrap are valid CI methods")
    }
  } else {
    LIN <- data.frame(Value = x , LIN = LIN)
    return(LIN)
  }
}
