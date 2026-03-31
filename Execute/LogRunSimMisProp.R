# Generate a log file name with the current date and time
log_file <- paste0("Log/ScriptOutput_", paste0("af_70_sim_mis_prop", format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".log"))

# Redirect console output and messages to the log file
sink(log_file)

# Track start time
start_time <- Sys.time()
cat("Start time:", format(start_time), "\n")

# Source the R script, logging all output
source("AnalysisCodes/sim_mis_prop.R", echo = T, max.deparse.length = 100000)

# Track stop time
end_time <- Sys.time()
cat("End time:", format(end_time), "\n")

# Stop logging and return output back to the console
sink()  # Close the standard output sink
