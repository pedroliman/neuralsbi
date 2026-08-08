# Shared across the benchmark scripts: resolve a --task flag to the matching
# nsbi_task object. Sourced by 01_generate_data.R, 03_run_neuralsbi.R and
# 04_compare.R so the three don't drift on which task names are supported.
benchmark_task <- function(task_name) {
  switch(task_name,
    gaussian_linear = task_gaussian_linear(),
    two_moons = task_two_moons(),
    slcp = task_slcp(),
    stop("unknown task: ", task_name)
  )
}
