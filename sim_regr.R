## Startup: load packages
# suppressPackageStartupMessages({
  # # library(ggplot2)
# #  library(magrittr)
  # # library(dplyr)
  # # library(tidyr)
  # #library(dagitty)
  # #library(broom)
# #  library(DescTools)
# })

# Avoid using DescTools
OddsRatio <- function(x) {
  if (any(x == 0)) x <- x + 0.5
  lx <- log(x)
  or <- exp(lx[1, 1] + lx[2, 2] - lx[1, 2] - lx[2, 1])
  return(or)
}
RelRisk <- function(x) {
  x1 <- x[1, 1]
  x2 <- x[2, 1]
  n1 <- x[1, 1] + x[1, 2]
  n2 <- x[2, 1] + x[2, 2]
  rr <- (x[1, 1]/sum(x[1, ]))/(x[2, 1]/sum(x[2, ]))
  return(rr)
}
Rev <- function (x, margin, ...) {
  if (!is.array(x)) 
      stop("'x' is not an array")
  newdim <- rep("", length(dim(x)))
  newdim[margin] <- paste(dim(x), ":1", sep = "")[margin]
  z <- eval(parse(text = gettextf("x[%s, drop = FALSE]", paste(newdim, 
      sep = "", collapse = ","))))
  class(z) <- oldClass(x)
  return(z)
}       

## Simulate the data

sim_regr <- function(my_seed, 
                     generating_model, regression_model,
		     reps, N, 
		     p1_E, p1_C, 
		     sd_E, sd_C, sd_D,
		     ror_ED, ror_CD, ror_ECD) {

  set.seed(my_seed)

  bigN <- N * reps

  # Simulate E and C
  if (generating_model %in% c("odds", "risk")) {
    E <- rbinom(bigN, 1, p1_E)  
    C <- rbinom(bigN, 1, p1_C)
  } else if (generating_model %in% c("linear")) {
    E <- rnorm(bigN, p1_E, sd_E)
    C <- rnorm(bigN, p1_C, sd_C)
  }

  # Net risk ratio for D conditional on E and C
  ror_D <- ror_ED^E * ror_CD^C * ror_ECD^(E*C) 

  # Computing a probability is more difficult for RR: need an appropriate intercept
  p1_D <- switch(generating_model,
                 risk=exp(-(log(ror_ED) + log(ror_CD) + log(ror_ECD)) -log(2)) * ror_D,
                 odds=ror_D/(ror_D + 1))
  
  if (generating_model %in% c("odds", "risk")) {
    D <- rbinom(bigN, 1, p1_D)  # Sample from the distribution of D
  } else if (generating_model %in% c("linear")) {
    D <- rnorm(bigN, log(ror_ED) * E + log(ror_CD) * C + log(ror_ECD) * E * C, sd_D)
  }

  big_data_set <- data.frame(rep=rep(1:N, reps), E=E, C=C, D=D) 
  big_data_set <- big_data_set[order(big_data_set$rep),]

  ## Tabulate the data (only meaningful in binary case)

  if (generating_model %in% c("odds", "risk")) {
    orfr <- function(x) { 
              sprintf("Odds ratio: %5.3f, relative risk: %5.3f", 
                      OddsRatio(x), RelRisk(Rev(x)))}

    # sprintf("----\nStratum C==0\n") |> cat()
    # xtabs(~E+D, data=big_data_set |> filter(C==0)) %T>% print() |> orfr() |> print()
    # sprintf("----\nStratum C==1\n") |> cat()
    # xtabs(~E+D, data=big_data_set |> filter(C==1)) %T>% print() |> orfr() |> print()
    # sprintf("----\nFull data set\n") |> cat()
    # xtabs(~E+D, data=big_data_set) %T>% 
      # print() %T>% 
      # (function(x) print(orfr(x))) ->
      # c_table
    # Avoiding magrittr:
    sprintf("----\nStratum C==0\n") |> cat()
    c_table <- xtabs(~E+D, data=big_data_set[big_data_set$C==0, ])
    print(c_table)
    print(orfr(c_table))
    sprintf("----\nStratum C==1\n") |> cat()
    c_table <- xtabs(~E+D, data=big_data_set[big_data_set$C==1, ])
    print(c_table)
    print(orfr(c_table))
    sprintf("----\nFull data set\n") |> cat()
    c_table <- xtabs(~E+D, data=big_data_set)
    print(c_table)
    print(orfr(c_table))

  }

  switch(regression_model, 
           "risk"=RelRisk(Rev(c_table)),
           "odds"=OddsRatio(c_table),
           "linear"=ror_ED) ->
    ror_0

  ## Perform regressions
  reg_fam <- switch(regression_model, 
                    "risk"=binomial(link=log),
                    "odds"=binomial(link=logit),
                    "linear"=gaussian())

  regression_formulae <- list(Naive=formula(D ~ E), 
                              Adjusted_C=formula(D ~ E + C),
                              Adjusted_C_interact=formula(D ~ E*C))
  regrs <- regression_formulae
  for (i in 1:length(regrs)) {
    regrs[[i]] <- vector("list", length=reps)
  }

  for (i in 1:reps) {
    this_rep_data <- big_data_set[(i-1)*N + (1:N),]
    for (j in names(regrs)) {
      try(regrs[[j]][[i]] <- glm(regression_formulae[[j]], 
                                 family=reg_fam,
                                 data=this_rep_data))
      if (is.null(regrs[[j]][[i]])) {
      # Rerun the regression but with the logit link to establish 
      # initial values
        logitlink <- glm(regression_formulae[[j]], 
                           family=binomial(),
                           data=this_rep_data) 
        logitlinkcoef <- coef(logitlink)
        # Work on the principle that the coefficients representing 
        # ORs approximate those for RRs and then need the initial
        # intercept to give the log of a probability.  
        starting_intercept <- -max(logitlinkcoef[-1] %*% 
                                     t(model.matrix(logitlink, this_rep_data)[,-1])) - 0.1
        regrs[[j]][[i]] <- glm(regression_formulae[[j]], 
                               family=reg_fam,
                               data=this_rep_data,
                               start=c(starting_intercept, logitlinkcoef[-1]))
        sprintf("%s: %d", j, i) |> print()
      }
    }
  }

  # An alternative because using tidy() seems to be slow
  bindy <- function(z) { do.call("rbind", z) }
  estimates <- 
    lapply(regrs, function(y) { lapply(y,
                                       function(x) { data.frame(coef(summary(x)))["E",] }) |>
                                  bindy() } ) 

  for (i in names(regrs)) {
    # estimates[[i]] <-
      # estimates[[i]] |> 
      # mutate(Type=i) |> 
      # mutate(term="E", .before=everything()) 
    estimates[[i]]$Type <- i
    estimates[[i]]$term <- "E"
    estimates[[i]] <- estimates[[i]][, c(ncol(estimates[[i]]), 1:(ncol(estimates[[i]])-1))]
  }

  estimates <- do.call("rbind", estimates)
print(str(estimates))
  estimates$Type <- factor(estimates$Type, levels=names(regrs)) 
  estimates <- estimates[order(estimates$Type), ]
  estimates$rep <- 1:reps

  if (regression_model %in% c("odds", "risk")) { 
    # estimates <- 
      # estimates |> 
      # mutate(Estimate = exp(Estimate))
    estimates$Estimate <- exp(estimates$Estimate)
  }

  names(estimates) <- c("term", "estimate", "std.error", "statistic", "p.value", names(estimates)[-(1:5)])

  ## Plot individual estimates

  estimates |>
    # pivot_wider(id_cols=c(term, rep), names_from=Type, values_from=estimate:p.value) ->
    data.frame() |>
    reshape(direction="wide", 
            idvar=c("term", "rep"), 
            timevar="Type", 
            v.names=c("estimate", "std.error", "statistic", "p.value"), 
            sep="_") ->
    estimates_to_plot

  # ggplot(estimates_to_plot, aes(estimate_Naive, estimate_Adjusted_C)) +
    # geom_point() +
    # geom_abline(slope=1, intercept=0, colour="#FF000080", linewidth=4) ->
    # linear_estimates_coef
  plot(estimate_Adjusted_C ~ estimate_Naive, data=estimates_to_plot,
       xlab="Estimate (no adjustment)",
       ylab="Estimate (adjusting for covariate)")
  abline(0, 1, col="#FF000080", lwd=4)
  
  # linear_estimates_coef
  
  # ggplot(estimates_to_plot, aes(std.error_Naive, std.error_Adjusted_C)) +
    # geom_point() +
    # geom_abline(slope=1, intercept=0, colour="#FF000080", linewidth=4) ->
    # linear_estimates_se

  plot(std.error_Adjusted_C ~ std.error_Naive, estimates_to_plot,
       xlab="Standard error (no adjustment)",
       ylab="Standard error (adjusting for covariate)")
  abline(0, 1, col="#FF000080", lwd=4)

  # linear_estimates_se
}

# sim_regr(42, "odds", "odds", 400, 400, 0.5, 0.5, 1, 1, 1, 2, 2, 1)
