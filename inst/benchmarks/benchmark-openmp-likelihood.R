# Benchmark serial versus parallel likelihood optimization in dynr.
# Replace `model` and `infile` with a representative multi-subject dynr model.

thread_counts <- c(1L, 2L, 4L, 8L)
results <- lapply(thread_counts, function(n_threads) {
  gc()
  elapsed <- system.time({
    fit <- dynr.cook(
      dynrModel = model,
      infile = infile,
      n_threads = n_threads,
      verbose = FALSE,
      hessian_flag = FALSE
    )
  })[["elapsed"]]

  data.frame(
    n_threads = n_threads,
    elapsed_seconds = elapsed,
    neg_log_likelihood = fit@neg.log.likelihood,
    exitflag = fit@exitflag
  )
})

results <- do.call(rbind, results)
results$speedup <- results$elapsed_seconds[results$n_threads == 1L] / results$elapsed_seconds
results$efficiency <- results$speedup / results$n_threads
print(results)
