library(tidyverse)
library(deSolve)
library(adaptivetau)
library(patchwork)

# SIRS model with:
# multiple hosts 
#   (inputted as n hosts and using vectors to easily add more hosts)
#   (host relative abundances drawn from geometric, log, or lognormal curve)
#   (lognormal curve doesn't drop as steeply for rarer species)
# Stochasticity via tau-leaping (rather than full Gillespie)
# Seasonal, pulsed births

# To be added: 
#   Variation in duration of infection: Two classes (one carrier, one variable through multiple compartments)
#   Variation in immunity: Two classes (one fast waner, one regular)
#   Variation in symmetry of transmission
#   Mortality rate differs among hosts and follows allometry (rarer species are larger and longer lived)

# Overall idea is to have knobs for:
#     n_hosts, popn_size, host abundance curves, 
#     birth synchrony (shape and offset), 

# Knobs:
#   n_hosts_in          : how many host species/populations to simulate
#   abundance_type      : "geometric", "log", or "lognormal"
#   total_N             : total population summed across all hosts, split
#                         according to the relative abundance curve
#   I0_per_host         : how many initial infections to seed per host
#                         (default: 1 in the most abundant host, 0 elsewhere)
#   beta_within         : a vector of length of n_hosts for within-host transmission
#                         rates (recycled/truncated if needed, or auto-generated)
#   cross_contact       : off-diagonal (between-host) transmission scaling
#   gamma, omega, mu    : vectors of length of n_hosts (recycled if shorter)
#   birth_peak_day_host1: day-of-year of peak births for host 1
#                         (scenario 1 default: 120)
#   birth_width         : von Mises pulse width in days, SAME for every host
#   birth_offset_days   : days between each host's peak births and the next
#                         Default: 0 => scenario 1: all hosts share identical birth seasonality
#                         If >0 => staggered peaks
#   birth_scale         : birth rate (a vector of length n_hosts) (recycled if shorter)

# First a function to generate host abundance curves
f.abundance_curver <- function(n_species, k = 0.5, alpha = 5, sigma = 1.5) {
  
  ranks <- 1:n_species
  N_internal <- 100  # arbitrary, doesn't matter since working in proportions
  
  # 1) Geometric series (Motomura)
  df.geometric <- tibble(rank = ranks) %>%
    mutate(raw_abundance = k * (1 - k)^(rank - 1),
           distribution.name.long = "Geometric series (Motomura)",
           distribution.name = "geometric")
  
  # 2) Fisher's log-series
  x <- N_internal / (N_internal + alpha)
  df.log <- tibble(rank = ranks) %>%
    mutate(raw_abundance = alpha * x^rank / rank,
           distribution.name.long = "Log-series (Fisher)",
           distribution.name = "log")
  
  # 3) Preston's lognormal
  df.lognormal <- tibble(rank = ranks) %>%
    mutate(percentile = (n_species - rank + 0.5) / n_species,
           raw_abundance = qlnorm(percentile, meanlog = 0, sdlog = sigma),
           distribution.name.long = "Lognormal (Preston)",
           distribution.name = "lognormal")
  
  # Combine, normalize to relative proportions within each distribution
  df.out <- bind_rows(df.geometric, df.log, df.lognormal) %>%
    group_by(distribution.name) %>%
    mutate(abundance.relative = raw_abundance / sum(raw_abundance)) %>%
    ungroup() %>%
    mutate(species = paste0("Host_Species_", rank)) %>%
    select(distribution.name, distribution.name.long, rank, species, abundance.relative)
  
  return(df.out)
}

# Second a function to generate birth pulses

# Inputs:
#   n_hosts       : number of hosts
#   peak_day_host1: day-of-year (0-365) of peak births for host 1
#   width         : concentration of the pulse, expressed as an
#                   approximate "days" width (converted internally to
#                   the von Mises kappa parameter). Same width used for
#                   every host
#   offset_days   : days between each host's peak and the next
#                   (host i peak = host1 peak + (i-1) * offset_days,
#                   wrapped into 0-365). offset_days = 0 (default)
#                   means every host shares identical birth seasonality.

generate_birth_pulses <- function(n_hosts, peak_day_host1 = 120,
                                  width = 15, offset_days = 0) {
  
  # Claude suggestion for wrapping around year
  kappa <- (365 / (2 * pi * width))^2
  
  peak_days <- (peak_day_host1 + (0:(n_hosts - 1)) * offset_days) %% 365
  
  tibble(
    host_index = 1:n_hosts,
    peak_day   = peak_days,
    kappa      = kappa
  )
}

# Circular (von Mises) seasonal factor, vectorized over hosts.
# day_of_year, peak_day, kappa can each be scalars or vectors of length n_hosts.
von_mises_seasonal_factor <- function(day_of_year, peak_day, kappa) {
  theta <- 2 * pi * day_of_year / 365
  mu    <- 2 * pi * peak_day / 365
  raw   <- exp(kappa * cos(theta - mu))
  # normalize so the max value (at theta == mu) is 1, matching the old
  # Gaussian's peak height of 1
  raw / exp(kappa)
}


# Third, compute rates

compute_rates <- function(state, parameters, time) {
  with(parameters, {
    
    n <- n_hosts
    S <- state[1:n]
    I <- state[(n+1):(2*n)]
    R <- state[(2*n+1):(3*n)]
    N <- S + I + R
    
    #Apply birth pulses
    day_of_year <- time %% 365
    seasonal_factor <- von_mises_seasonal_factor(day_of_year, birth_peak_day, birth_kappa)
    birth_rate <- mu * N * seasonal_factor * birth_scale
    
    foi <- as.vector(beta_matrix %*% (I / N))
    
    infection <- foi * S
    recovery  <- gamma * I
    waning    <- omega * R
    birth     <- birth_rate
    death_S   <- mu * S
    death_I   <- mu * I
    death_R   <- mu * R
    
    rates <- numeric(7 * n)
    for (i in 1:n) {
      idx <- ((i - 1) * 7 + 1):(i * 7)
      rates[idx] <- c(infection[i], recovery[i], waning[i], birth[i],
                      death_S[i], death_I[i], death_R[i])
    }
    
    rates[rates < 0] <- 0
    rates
  })
}

# Fourth, build_transitions

build_transitions <- function(n_hosts, state_names) {
  n <- n_hosts
  transitions <- list()
  
  for (i in 1:n) {
    zero_vec <- setNames(rep(0, length(state_names)), state_names)
    
    inf <- zero_vec; inf[paste0("S", i)] <- -1; inf[paste0("I", i)] <- +1
    rec <- zero_vec; rec[paste0("I", i)] <- -1; rec[paste0("R", i)] <- +1
    wan <- zero_vec; wan[paste0("R", i)] <- -1; wan[paste0("S", i)] <- +1
    bir <- zero_vec; bir[paste0("S", i)] <- +1
    dS  <- zero_vec; dS[paste0("S", i)]  <- -1
    dI  <- zero_vec; dI[paste0("I", i)]  <- -1
    dR  <- zero_vec; dR[paste0("R", i)]  <- -1
    
    transitions <- c(transitions, list(inf, rec, wan, bir, dS, dI, dR))
  }
  transitions
}

# Fifth, run tau leaps
run_tau_leap <- function(parameters, initial_state, t_max, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  transitions <- build_transitions(parameters$n_hosts, names(initial_state))
  
  result <- ssa.adaptivetau(
    init.values  = initial_state,
    transitions  = transitions,
    rateFunc     = compute_rates,
    params       = parameters,
    tf           = t_max
  )
  
  as.data.frame(result)
}


# Now build scenarios: turn n_hosts + abundance shape + birth pulse knobs into full parameters + initial_state

build_scenario <- function(n_hosts_in,
                           abundance_type = "geometric",
                           total_N = 20000,
                           I0_per_host = NULL,
                           cross_contact = 0.5,
                           beta_within = NULL, #default dummy of 0.35 used for now
                           gamma_rate = 1/20,
                           omega_rate = 1/365,
                           mu_days = c(5, 3, 8),
                           birth_peak_day_host1 = 120,
                           birth_width = 15,
                           birth_offset_days = 0,
                           birth_scale = 8,
                           abundance_k = 0.5,
                           abundance_alpha = 5,
                           abundance_sigma = 1.5) {
  
  # --- helper to recycle a param vector out to n_hosts_in ---
  recycle_to_n <- function(x, n) {
    if (is.null(x)) return(NULL)
    rep(x, length.out = n)
  }
  
  # --- 1) Get relative abundances from f.abundance_curver() ---
  df.abund <- f.abundance_curver(n_species = n_hosts_in,
                                 k = abundance_k,
                                 alpha = abundance_alpha,
                                 sigma = abundance_sigma) %>%
    filter(distribution.name == abundance_type) %>%
    arrange(rank)
  
  if (nrow(df.abund) != n_hosts_in) {
    stop("Abundance curve did not return one row per host - check n_hosts_in / abundance_type.")
  }
  
  host_pop_sizes <- round(total_N * df.abund$abundance.relative)
  
  # --- 2) Build initial_state from abundance-derived population sizes ---
  if (is.null(I0_per_host)) {
    I0_per_host <- rep(0, n_hosts_in)
    I0_per_host[1] <- 1   # seed infection in the most abundant host by default
  } else {
    I0_per_host <- recycle_to_n(I0_per_host, n_hosts_in)
  }
  
  S0 <- host_pop_sizes - I0_per_host
  R0 <- rep(0, n_hosts_in)
  
  initial_state <- c(
    setNames(S0, paste0("S", 1:n_hosts_in)),
    setNames(I0_per_host, paste0("I", 1:n_hosts_in)),
    setNames(R0, paste0("R", 1:n_hosts_in))
  )
  
  # --- 3) Build parameter list, recycling scalars/vectors to n_hosts_in ---
  if (is.null(beta_within)) {
    beta_within <- rep(0.35, n_hosts_in)   # sensible default if not supplied
  } else {
    beta_within <- recycle_to_n(beta_within, n_hosts_in)
  }
  
  gamma <- recycle_to_n(gamma_rate, n_hosts_in)
  omega <- recycle_to_n(omega_rate, n_hosts_in)
  mu    <- 1 / (recycle_to_n(mu_days, n_hosts_in) * 365)
  birth_scale <- recycle_to_n(birth_scale, n_hosts_in)
  
  beta_matrix <- matrix(cross_contact, nrow = n_hosts_in, ncol = n_hosts_in)
  diag(beta_matrix) <- 1
  beta_matrix <- sweep(beta_matrix, 2, beta_within, `*`)
  
  # --- 4) Birth pulses via dedicated generator ---
  # birth_offset_days = 0 (default) -> identical seasonality across hosts.
  # birth_offset_days > 0 -> host i peak shifts by (i-1)*offset_days days,
  # wrapping around the year. Width/shape stays identical across hosts.
  pulses <- generate_birth_pulses(
    n_hosts        = n_hosts_in,
    peak_day_host1 = birth_peak_day_host1,
    width          = birth_width,
    offset_days    = birth_offset_days
  )
  
  parameters <- list(
    n_hosts        = n_hosts_in,
    beta_matrix    = beta_matrix,
    gamma          = gamma,
    omega          = omega,
    mu             = mu,
    birth_peak_day = pulses$peak_day,
    birth_kappa    = pulses$kappa,
    birth_scale    = birth_scale
  )
  
  list(
    parameters      = parameters,
    initial_state   = initial_state,
    abundance_table = df.abund,
    host_pop_sizes  = host_pop_sizes,
    birth_pulses    = pulses
  )
}


# Example scenario: identical birth seasonality across all hosts (offset_days = 0)
scenario <- build_scenario(
  n_hosts_in           = 5,
  abundance_type       = "lognormal",   # "geometric", "log", or "lognormal"
  total_N              = 20000,
  cross_contact        = 0.5,
  birth_peak_day_host1 = 120,
  birth_width          = 15,
  birth_offset_days    = 0              # identical seasonality (scenario 1)
)

#Printing for convenience
scenario$abundance_table   # relative abundances used
scenario$initial_state     # derived S/I/R starting values
scenario$birth_pulses      # peak_day / kappa per host
scenario$parameters$beta_matrix

sim_tau <- run_tau_leap(
  parameters    = scenario$parameters,
  initial_state = scenario$initial_state,
  t_max         = 365 * 2,
  seed          = 1
)

#Making a vector of host names for convenience
host_names <- paste("Host", 1:scenario$parameters$n_hosts)

sim_tau_long <- sim_tau %>%
  mutate(time_years = time / 365) %>%
  pivot_longer(cols = -c(time, time_years), names_to = "var", values_to = "count") %>%
  mutate(
    state = substr(var, 1, 1),
    host_index = as.integer(sub("^[SIR]", "", var)),
    host = host_names[host_index]
  )

#Plotting
ggplot(sim_tau_long, aes(x = time_years, y = count, color = state)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ host, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c(S = "steelblue", I = "firebrick", R = "darkgreen")) +
  labs(
    title = paste0("Tau-Leap ", scenario$parameters$n_hosts, "-Host SIRS Model (",
                   unique(scenario$abundance_table$distribution.name), " abundance)"),
    x = "Time (years)", y = "Count", color = "Compartment"
  ) +
  theme_minimal(base_size = 14)


# Staggered birth seasonality scenario (offset_days > 0)
scenario_staggered <- build_scenario(
  n_hosts_in           = 5,
  abundance_type       = "lognormal",
  total_N              = 20000,
  cross_contact        = 0.5,
  birth_peak_day_host1 = 120,
  birth_width          = 15,
  birth_offset_days    = 60   # each subsequent host's peak shifts 60 days later
)

scenario_staggered$birth_pulses

sim_tau.staggered <- run_tau_leap(
  parameters    = scenario_staggered$parameters,
  initial_state = scenario_staggered$initial_state,
  t_max         = 365 * 2,
  seed          = 1
)

sim_tau_long.staggered <- sim_tau.staggered %>%
  mutate(time_years = time / 365) %>%
  pivot_longer(cols = -c(time, time_years), names_to = "var", values_to = "count") %>%
  mutate(
    state = substr(var, 1, 1),
    host_index = as.integer(sub("^[SIR]", "", var)),
    host = host_names[host_index]
  )

ggplot(sim_tau_long.staggered, aes(x = time_years, y = count, color = state)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ host, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c(S = "steelblue", I = "firebrick", R = "darkgreen")) +
  labs(
    title = paste0("Tau-Leap ", scenario_staggered$parameters$n_hosts, "-Host SIRS Model (",
                   unique(scenario_staggered$abundance_table$distribution.name), " abundance)"),
    x = "Time (years)", y = "Count", color = "Compartment"
  ) +
  theme_minimal(base_size = 14)






# Define the mortality function
mu_func <- function(w) {
  0.4 * w^(-0.26)
}

# Create a range of body masses (adjust range as needed)
mass_df <- data.frame(w = seq(0.01, 100, length.out = 500))
mass_df$mu <- mu_func(mass_df$w)

# Plot
ggplot(mass_df %>% filter(w > 0.5), aes(x = w, y = mu)) +
  geom_line(linewidth = 1, color = "steelblue") +
  labs(x = "Body mass, w (i)",
       y = expression(paste("Mortality rate, ", mu, "(i)")),
       title = "Host mortality rate as a function of body mass") +
  theme_minimal()

ggplot(mass_df, aes(x = w, y = mu)) +
  geom_line(linewidth = 1, color = "steelblue") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Body mass, w (i) [log scale]",
       y = expression(paste("Mortality rate, ", mu, "(i)", " [log scale]")),
       title = "Host mortality rate vs. body mass (log-log)") +
  theme_minimal()

