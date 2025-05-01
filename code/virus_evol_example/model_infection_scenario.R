args <- commandArgs(trailingOnly=T)

# where to save output
setwd(args[1])

# define parameters
R0 = 1.1 # basic reproduction number
inf_interval = 21 # average number of days between infection and transmission
N_inf_t0 = 3 # number of individuals infected on day 0

# START SIM ----
set.seed(1234567)
# create dataframe to record infection events
inf_events = as.data.frame(matrix(ncol=7, nrow=0)) # record infection events
colnames(inf_events) <- c("inf_step", "inf_id", "inf_source", "host_id", "infected_host_ids", "transmission_day", "overall_day")
# seed the infection sequence
H_idx <- N_inf_t0
H_inf0 <- paste0("H", 1:H_idx)
start_day <- rpois(1, inf_interval)
inf_events[1,] <- c(0, 0, NA, "seed", paste0(H_inf0, collapse=";"), start_day, start_day)
# continue infection sequence
for (x in 1:20) {
  N_infected <- rpois(length(H_inf0), R0) # determine how many hosts are infected by each infected host
  if (sum(N_infected) > 0) {  # at least one transmission event needs to happen in order to continue
    H_inf <- paste0("H", (H_idx+1):(H_idx+sum(N_infected))) # generate IDs for newly infected hosts
    for (y in seq_along(N_infected)) {
      if (N_infected[y] > 0) { # looping over all currently infected individuals, needs to transmit to at least one other individual
        for (j in 1:nrow(inf_events)){if (H_inf0[y] %in% unlist(strsplit(inf_events$infected_host_ids[j], split=";"))){break}} # find the row index corresponding to the source infection by searching through the dataframe of infection events
        inf_source <- inf_events[j, "inf_id"] # what is the inf_id of the infection that this came from
        # determine indices of the H_inf vector that correspond to each infected host 
        if (y == 1) { 
          i <- 1:N_infected[y] # if y = 1, it's just the start of the vector
        } else {
          i <- (sum(N_infected[1:(y-1)])+1):(sum(N_infected[1:(y-1)])+N_infected[y])
        }
        # determine start time of the current infection
        start_day <- as.numeric(inf_events[j, "overall_day"])
        # pick a transmission day for the current infection
        transmission_day <- rpois(1, inf_interval)
        inf_events[nrow(inf_events)+1,] <- c(x, nrow(inf_events), inf_source, H_inf0[y], paste0(H_inf[i], collapse=";"), transmission_day, start_day + transmission_day)
      }
    }
    # updated H_idx and H_inf0
    H_idx <- H_idx+sum(N_infected)
    H_inf0 <- H_inf
  }
}

write.csv(inf_events, "inf_sequence.csv", quote=F, row.names=F)
