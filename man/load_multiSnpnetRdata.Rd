% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{load_multiSnpnetRdata}
\alias{load_multiSnpnetRdata}
\title{Given a results directory and lambda index (optional), load the corresponding R Data file}
\usage{
load_multiSnpnetRdata(results_dir, lambda_idx = NULL)
}
\arguments{
\item{results_dir}{The results directory}

\item{lambda_idx}{The lambda index. If not specified, we call find_prev_iter() and load the last lambda index available on the file system.}
}
\description{
Given a results directory and lambda index (optional), load the corresponding R Data file
}
