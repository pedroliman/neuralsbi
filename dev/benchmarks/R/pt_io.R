# Reading the tensors that sbibm ships as PyTorch `.pt` / `.torch` files,
# without Python.
#
# Three of the tasks depend on frozen constants that sbibm stores as pickled
# torch tensors rather than as CSV: the Bernoulli GLM's input stimulus and
# design matrix, and the SLCP-distractor noise distribution. Reproducing those
# constants from their generating code is not possible in R, because they came
# out of numpy's RNG. So we read the files.
#
# The files use torch's *legacy* (pre-1.6, non-zip) format: five pickles
# followed by the storages in the order they appear in the fifth pickle, each
# storage written as an int64 element count followed by its raw little-endian
# payload. We do not parse the pickles; we locate storages by their known
# element counts, which is enough and keeps this to a few lines.

#' Read a raw storage from a legacy torch file, given its element count.
#'
#' @param path File path.
#' @param numel Number of elements in the storage.
#' @param what "double" (float64), "float" (float32) or "int64".
#' @param which Which matching storage to take when several share `numel`:
#'   "first" or "last".
read_torch_storage <- function(path, numel, what = c("float", "double", "int64"),
                               which = c("first", "last")) {
  what <- match.arg(what)
  which <- match.arg(which)
  size <- switch(what, float = 4L, double = 8L, int64 = 8L)
  n_bytes <- numel * size

  con <- file(path, "rb")
  on.exit(close(con))
  all_bytes <- readBin(con, "raw", n = file.info(path)$size)

  offsets <- find_int64(all_bytes, numel)
  if (length(offsets) == 0) {
    stop("No storage with ", numel, " elements found in ", path, call. = FALSE)
  }
  # A storage header is only plausible if the payload fits in the remaining file.
  offsets <- offsets[offsets + 8 + n_bytes - 1 <= length(all_bytes)]
  if (length(offsets) == 0) {
    stop("Storage header for ", numel, " elements does not fit in ", path,
         call. = FALSE)
  }
  start <- if (which == "first") offsets[1] else offsets[length(offsets)]
  payload <- all_bytes[(start + 8):(start + 8 + n_bytes - 1)]

  if (what == "int64") {
    # readBin caps integer size at 4 bytes; take the low word of each int64.
    idx <- seq(1, n_bytes, by = 8)
    lo <- vapply(idx, function(i) {
      as.numeric(readBin(payload[i:(i + 3)], "integer", n = 1, size = 4,
                         endian = "little"))
    }, numeric(1))
    return(lo)
  }
  readBin(payload, "numeric", n = numel, size = size, endian = "little")
}

#' Byte offsets (1-based) where a little-endian int64 equal to `value` starts.
find_int64 <- function(bytes, value) {
  target <- c(writeBin(as.integer(value), raw(), size = 4, endian = "little"),
              as.raw(c(0, 0, 0, 0)))
  cand <- which(bytes == target[1])
  cand <- cand[cand + 7 <= length(bytes)]
  keep <- vapply(cand, function(i) all(bytes[i:(i + 7)] == target), logical(1))
  cand[keep]
}
