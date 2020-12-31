# -----------------------------------------------------------------------------
# Data Generation
# -----------------------------------------------------------------------------

#' Generate Competing Risks Data
#' 
#' Assumes exponential cause-specific hazards and independent exponential
#' censoring.
#' 
#' @param n Sample size.
#' @param censor_rate Censoring rate.
#' @param event_rate Event rate.
#' @param death_rate Death rate.
#' @param tau Optional truncation time needed if censoring rate is set to 0.
#' @return Data.frame containing the observation 'time' and event 'status'.
#' @importFrom stats rexp rmultinom
#' @export 

GenData <- function(
  n,
  censor_rate,
  event_rate,
  death_rate,
  tau = NULL
) {
  
  # Censoring times.
  if (censor_rate > 0) {
    censor <- rexp(n = n, rate = censor_rate)
  } else {
    if (is.null(tau)) {stop("Tau is required if censoring is absent.")}
    censor <- rep(tau, n)
  }
  
  
  # Overall hazard.
  overall_hazard <- event_rate + death_rate
  
  # Arrivals.
  event_times <- rexp(n = n, rate = overall_hazard)
  
  # Event type.
  pi <- c(event_rate, death_rate) / overall_hazard
  status <- rmultinom(n = n, size = 1, prob = pi)
  status <- apply(status, 2, which.max)
  
  # Observed data.
  out <- data.frame(
    time = pmin(event_times, censor),
    status = (event_times <= censor) * status
  )
  
  # Add columns needed by etm.
  out$from <- 0
  out$to <- as.character(out$status)
  out$to[out$to == "0"] <- "cens"
  return(out)
}