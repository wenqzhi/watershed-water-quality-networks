# Local data configuration for the WRR1 analysis scripts.
# Keep this repository code-only; store raw data and generated results outside it.
# Set WRR_DATA_ROOT to the folder that contains NTN_2023_MonthlyLoadTable.csv,
# NTN-station.txt, TE_network_results/, PNet_2011_2023_BlockTracing/, etc.

get_data_root <- function() {
  data_root <- Sys.getenv("WRR_DATA_ROOT", unset = "")
  if (!nzchar(data_root)) {
    stop(
      "Set the WRR_DATA_ROOT environment variable to your local WRR data directory.\n",
      "See README.md for setup instructions.",
      call. = FALSE
    )
  }
  if (!dir.exists(data_root)) {
    stop("WRR_DATA_ROOT does not exist: ", data_root, call. = FALSE)
  }
  normalizePath(data_root, winslash = "/", mustWork = TRUE)
}
