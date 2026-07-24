% Hand-written (mirrors roxygen tags); regenerate with devtools::document()
\name{plot_coverage}
\alias{plot_coverage}
\title{Plot nominal vs. empirical credible-interval coverage}
\usage{
plot_coverage(sbc_result, levels = seq(0.05, 0.95, by = 0.05))
}
\arguments{
\item{sbc_result}{An \code{nsbi_sbc} object from \code{\link{sbc}()}.}

\item{levels}{Nominal credibility levels to evaluate.}
}
\value{
A \code{ggplot} object (also drawn as a side effect), invisibly.
}
\description{
Well-calibrated posteriors lie on the diagonal. Curves above the diagonal
mean the posterior is too wide (conservative); below means overconfident.
A shaded band shows the Monte-Carlo uncertainty from the finite number of
SBC trials.
}
