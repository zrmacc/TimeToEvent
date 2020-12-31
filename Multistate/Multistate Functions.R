library(cowplot)
library(ggplot2)

# -----------------------------------------------------------------------------
# Data Generation.
# -----------------------------------------------------------------------------

#' Simulate Single Subject.
#' 
#' @param id Subject index.
#' @param rate_c Censoring rate.
#' @param rate_01 Rate for healthy to ill transition.
#' @param rate_02 Rate for healthy to death transition.
#' @param rate_12 Rate for ill to death transition.
#' @param tau Truncation time.

SimSubj <- function(
  id,
  rate_c,
  rate_01,
  rate_02,
  rate_12,
  tau
) {
  
  # Censoring time.
  censor <- min(rexp(n = 1, rate = rate_c), tau)
  
  # First transition.
  total_rate_0 <- rate_01 + rate_02
  first_time <- rexp(n = 1, rate = total_rate_0)
  first_type <- sample(
    x = c(1, 2), 
    size = 1, 
    prob = c(rate_01, rate_02) / total_rate_0
  )
  
  # Record first transition.
  if(censor < first_time) {
    out <- data.frame(
      id = id,
      entry = 0,
      exit = censor,
      from = 0,
      to = 0,
      status = 0
    )
    return(out)
  } else {
    out <- data.frame(
      id = id,
      entry = 0,
      exit = first_time,
      from = 0,
      to = first_type,
      status = 1
    )
    if (first_type == 2) {return(out)}
  }
  
  # Simulate second transition.
  second_time <- first_time + rexp(n = 1, rate = rate_12)
  
  if(censor < second_time) {
    out2 <- data.frame(
      id = id,
      entry = first_time,
      exit = censor,
      from = 1,
      to = 1,
      status = 0
    )
    out <- rbind(out, out2)
    return(out)
  } else {
    out2 <- data.frame(
      id = id,
      entry = first_time,
      exit = second_time,
      from = 1,
      to = 3,
      status = 1
    )
    out <- rbind(out, out2)
    return(out)
  }
}

#' Simulate Data from Illness-Death Model
#' 
#' @param n Sample size.
#' @param rate_c Censoring rate.
#' @param rate_01 Rate for healthy to ill transition.
#' @param rate_02 Rate for healthy to death transition.
#' @param rate_12 Rate for ill to death transition.
#' @param tau Truncation time.
#' @return Data.frame with subject 'id', the 'entry' and 'exit' of the
#'   observation period, the initial and final states ('from' and 'to'), and the
#'   status, 1 for observed, 0 for censored. Death after illness is coded as
#'   3 rather than 2. Note that the column names correspomd to those expected by mvna.

SimData <- function(
  n,
  rate_c,
  rate_01,
  rate_02,
  rate_12,
  tau
) {
  
  aux <- function(i) {
    return(SimSubj(i, rate_c, rate_01, rate_02, rate_12, tau))
  }
  data <- lapply(seq_len(n), aux)
  data <- do.call(rbind, data)
  
  # Convert states to factors.
  data$from <- as.character(data$from)
  data$to <- as.character(data$to)
  data$to[data$status == 0] <- "cens"
  return(data)  
}


# -----------------------------------------------------------------------------
# Plotting.
# -----------------------------------------------------------------------------

#' Plot Hazard Curves
#' 
#' @param alpha Type I error level.
#' @param df Data.frame containing the outcome of `mvna`, together with
#' a factor 'trans' specifying the transition type.
#' @param title Title.
#' @param y_lab Y-axis label.

PlotHazardCurves <- function(
  alpha = 0.05, 
  df,
  title,
  y_lab
  ) {
  
  z <- qnorm(1 - alpha / 2)
  df$lower <- df$na - z * sqrt(df$var.aalen)
  df$lower <- pmax(df$lower, 0)
  df$upper <- df$na + z * sqrt(df$var.aalen)
  q <- ggplot(data = df) + 
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + 
    geom_ribbon(
      aes(x = time, ymin = lower, ymax = upper, fill = trans),
      alpha = 0.2
    ) + 
    geom_step(
      aes(x = time, y = na, color = trans)
    ) + 
    labs(
      x = "Time",
      y = y_lab,
      fill = "Transition",
      color = "Transition",
      title = title
    )
  return(q)
}


#' Plot Prob Curves
#' 
#' @param alpha Type I error.
#' @param color Color.
#' @param probs Probabilities.
#' @param ses Standard errors.
#' @param times Time.

PlotProbCurve <- function(
  alpha = 0.05,
  color = "#6385B8",
  probs,
  ses,
  times,
  title = NULL
) {
  
  df <- data.frame(
    probs = probs,
    ses = ses,
    times = times
  )
  z <- qnorm(1 - alpha / 2)
  df$lower <- df$probs - z * df$ses
  df$lower <- pmax(df$lower, 0)
  df$upper <- df$probs + z * df$ses
  df$upper <- pmin(df$upper, 1)
  
  q <- ggplot(data = df) + 
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + 
    geom_ribbon(
      aes(x = times, ymin = lower, ymax = upper),
      alpha = 0.2,
      fill = color
    ) + 
    geom_step(
      aes(x = times, y = probs),
      color = color
    ) + 
    labs(
      x = "Time",
      y = "Probability",
      title = title
    ) +
    ylim(-0.05, 1.05)
  return(q)
}

# -----------------------------------------------------------------------------
