# Benchmark subject-level OpenMP likelihood evaluation in dynr
#
# This script builds the LinearSDE model from the dynr demo, expands the
# Oscillator data to multiple subjects, and compares serial and parallel fits.

library(dynr)

# Confirm that the installed dynr includes the OpenMP interface.
if (!"n_threads" %in% names(formals(dynr::dynr.cook))) {
  stop(
    "The installed dynr::dynr.cook() does not have an n_threads argument. ",
    "Reinstall the feature/openmp-likelihood branch."
  )
}

# ---------------------------------------------------------------------------
# Benchmark controls

thread_counts <- c(1L, 2L, 4L, 8L)

# The original Oscillator data contain one subject. Replication creates
# independent subject likelihood contributions for the OpenMP benchmark.
n_subjects <- 32L

# Repeat each fit to reduce timing noise. Increase after the first successful run.
n_repetitions <- 3L

# Set TRUE only when debugging. Parallel likelihood evaluation intentionally
# falls back to serial mode when verbose output is enabled.
verbose_fit <- FALSE

# ---------------------------------------------------------------------------
# Construct the LinearSDE model from demo("LinearSDE", package = "dynr")

measurement <- prep.measurement(
  values.load = matrix(c(1, 0), 1, 2),
  params.load = matrix(c("fixed", "fixed"), 1, 2),
  state.names = c("Position", "Velocity"),
  obs.names = "y1"
)

noise <- prep.noise(
  values.latent = diag(c(0, 1), 2),
  params.latent = diag(c("fixed", "dnoise"), 2),
  values.observed = diag(1.5, 1),
  params.observed = diag("mnoise", 1)
)

initial <- prep.initial(
  values.inistate = c(0, 1),
  params.inistate = c("inipos", "fixed"),
  values.inicov = diag(1, 2),
  params.inicov = diag("fixed", 2)
)

dynamics <- prep.matrixDynamics(
  values.dyn = matrix(c(0, -0.1, 1, -0.2), 2, 2),
  params.dyn = matrix(
    c("fixed", "spring", "fixed", "friction"),
    2,
    2
  ),
  isContinuousTime = TRUE
)

# ---------------------------------------------------------------------------
# Create a representative multi-subject data set

data("Oscillator", package = "dynr")

required_columns <- c("id", "times", "y1")
missing_columns <- setdiff(required_columns, names(Oscillator))
if (length(missing_columns) > 0L) {
  stop(
    "Oscillator is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

base_subject <- Oscillator
base_subject$id <- 1L

# Replicate the observed series across independent subjects. A tiny deterministic
# shift avoids perfectly identical observed trajectories while preserving the
# same model structure. The shift is deliberately small relative to the data.
benchmark_data <- do.call(
  rbind,
  lapply(seq_len(n_subjects), function(subject_id) {
    subject_data <- base_subject
    subject_data$id <- subject_id
    subject_data$y1 <- subject_data$y1 +
      (subject_id - 1L) * .Machine$double.eps
    subject_data
  })
)

row.names(benchmark_data) <- NULL

dynr_data <- dynr.data(
  benchmark_data,
  id = "id",
  time = "times",
  observed = "y1"
)

model <- dynr.model(
  dynamics = dynamics,
  measurement = measurement,
  noise = noise,
  initial = initial,
  data = dynr_data,
  outfile = "LinearSDE-openmp-benchmark.c"
)

model$ub <- c(
  friction = 101,
  spring = 100,
  inipos = 103,
  mnoise = 100,
  dnoise = 99
)

cat(
  "Benchmark data:",
  n_subjects,
  "subjects and",
  nrow(benchmark_data),
  "rows\n"
)
cat("Thread counts:", paste(thread_counts, collapse = ", "), "\n")
cat("Repetitions:", n_repetitions, "\n\n")

# ---------------------------------------------------------------------------
# Fit helper

fit_once <- function(n_threads) {
  fit <- NULL
  error_message <- NA_character_

  timing <- system.time({
    fit <- tryCatch(
      dynr.cook(
        dynrModel = model,
        n_threads = n_threads,
        verbose = verbose_fit,
        hessian_flag = FALSE
      ),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
  })

  if (is.null(fit)) {
    return(
      data.frame(
        repetition = NA_integer_,
        n_threads = n_threads,
        elapsed_seconds = unname(timing[["elapsed"]]),
        neg_log_likelihood = NA_real_,
        exitflag = NA_real_,
        error = error_message,
        stringsAsFactors = FALSE
      )
    )
  }

  data.frame(
    repetition = NA_integer_,
    n_threads = n_threads,
    elapsed_seconds = unname(timing[["elapsed"]]),
    neg_log_likelihood = fit@neg.log.likelihood,
    exitflag = fit@exitflag,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Warm-up
#
# The first call includes dynamic-library loading, generated-code compilation,
# and filesystem/cache effects. Exclude it from the reported timings.

cat("Running one-thread warm-up fit...\n")
warmup <- fit_once(1L)

if (!is.na(warmup$error)) {
  stop("Warm-up fit failed: ", warmup$error)
}

cat(
  "Warm-up complete. Negative log likelihood:",
  format(warmup$neg_log_likelihood, digits = 12),
  "\n\n"
)

# ---------------------------------------------------------------------------
# Timed benchmark
#
# Interleave thread counts by repetition rather than running all one-thread
# fits first. This reduces confounding from changing system load or temperature.

results <- vector(
  "list",
  length(thread_counts) * n_repetitions
)
result_index <- 1L

for (repetition in seq_len(n_repetitions)) {
  cat("Repetition", repetition, "of", n_repetitions, "\n")

  for (n_threads in thread_counts) {
    gc()

    cat("  n_threads =", n_threads, "... ")
    current <- fit_once(n_threads)
    current$repetition <- repetition

    if (is.na(current$error)) {
      cat(
        format(current$elapsed_seconds, digits = 5),
        "seconds; NLL =",
        format(current$neg_log_likelihood, digits = 12),
        "\n"
      )
    } else {
      cat("ERROR:", current$error, "\n")
    }

    results[[result_index]] <- current
    result_index <- result_index + 1L
  }
}

results <- do.call(rbind, results)

# ---------------------------------------------------------------------------
# Numerical-equivalence checks

successful <- is.na(results$error)
serial_nll <- results$neg_log_likelihood[
  successful & results$n_threads == 1L
]

if (length(serial_nll) > 0L) {
  reference_nll <- serial_nll[[1L]]
  results$nll_difference <- results$neg_log_likelihood - reference_nll
  results$nll_equivalent <- abs(results$nll_difference) <= 1e-7
} else {
  results$nll_difference <- NA_real_
  results$nll_equivalent <- NA
}

# ---------------------------------------------------------------------------
# Timing summary

successful_results <- results[successful, , drop = FALSE]

if (nrow(successful_results) > 0L) {
  timing_summary <- aggregate(
    elapsed_seconds ~ n_threads,
    data = successful_results,
    FUN = function(x) {
      c(
        median = median(x),
        minimum = min(x),
        maximum = max(x)
      )
    }
  )

  timing_summary <- data.frame(
    n_threads = timing_summary$n_threads,
    median_seconds = timing_summary$elapsed_seconds[, "median"],
    minimum_seconds = timing_summary$elapsed_seconds[, "minimum"],
    maximum_seconds = timing_summary$elapsed_seconds[, "maximum"],
    row.names = NULL
  )

  serial_median <- timing_summary$median_seconds[
    timing_summary$n_threads == 1L
  ]

  if (length(serial_median) == 1L) {
    timing_summary$speedup <- serial_median /
      timing_summary$median_seconds
    timing_summary$efficiency <- timing_summary$speedup /
      timing_summary$n_threads
  } else {
    timing_summary$speedup <- NA_real_
    timing_summary$efficiency <- NA_real_
  }
} else {
  timing_summary <- data.frame()
}

cat("\nIndividual fits\n")
print(results, row.names = FALSE)

cat("\nTiming summary\n")
print(timing_summary, row.names = FALSE)

# Fail clearly if any fit failed or if likelihoods disagree materially.
if (any(!successful)) {
  warning("At least one benchmark fit failed. Inspect results$error.")
}

if (any(results$nll_equivalent %in% FALSE, na.rm = TRUE)) {
  warning(
    "At least one parallel fit differs from the serial likelihood ",
    "by more than 1e-7."
  )
}

# ---------------------------------------------------------------------------
# Save reproducible outputs

write.csv(
  results,
  file = "openmp-benchmark-individual.csv",
  row.names = FALSE
)

write.csv(
  timing_summary,
  file = "openmp-benchmark-summary.csv",
  row.names = FALSE
)

saveRDS(
  list(
    controls = list(
      thread_counts = thread_counts,
      n_subjects = n_subjects,
      n_repetitions = n_repetitions
    ),
    warmup = warmup,
    individual = results,
    summary = timing_summary,
    session_info = sessionInfo()
  ),
  file = "openmp-benchmark-results.rds"
)

cat(
  "\nSaved:\n",
  "  openmp-benchmark-individual.csv\n",
  "  openmp-benchmark-summary.csv\n",
  "  openmp-benchmark-results.rds\n",
  sep = ""
)
