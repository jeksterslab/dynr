# OpenMP likelihood evaluation

This branch adds optional subject-level parallelism to the negative log-likelihood evaluations requested by NLopt SLSQP. The optimizer itself remains serial.

## Usage

```r
fit <- dynr.cook(
  model,
  infile = infile,
  n_threads = 4L,
  verbose = FALSE
)
```

The default is `n_threads = 1L`. When OpenMP is unavailable, the code compiles and executes serially. Parallel likelihood evaluation is disabled under `verbose = TRUE` because R console output is not thread-safe.

## Design

Each subject is evaluated independently using a private `Param` object and private filtering workspace. The subject contributions are summed with an OpenMP reduction. The final Extended Kim filter and smoother remain serial because they populate shared full-sample output arrays.

## Validation

Compare serial and parallel fits using numerical tolerances rather than bitwise identity because OpenMP reduction can change floating-point summation order.
