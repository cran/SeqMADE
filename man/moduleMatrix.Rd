\name{moduleMatrix}
\alias{moduleMatrix}
\title{Modulematrix Construction}
\description{
A function of constructing the modulematrix for the modules was used to indicate whether genes belong to a given module or not.
}
\usage{
moduleMatrix(exprs, networkModule)
}
\arguments{
  \item{exprs}{
exprs is a data frame or matrix for two groups or conditions, with rows as variables (genes) and columns as samples
}
  \item{networkModule}{
NetworkModule is the gene sets or modules in the biological network or metabolic pathway, with the 1th column as the module names and the 2th columnn as the gene symbols constituting the module}
}
\details{
 Modulematrix is a matrix, in which the indicator variables 1 or 0 represent whether a gene belong to a given module or not.
}
\author{
Mingli Lei, Jia Xu, Li-Ching Huang, Lily Wang, Jing Li
}
\examples{
data(exprs)
data(networkModule)
modulematrix <- moduleMatrix(exprs,networkModule)
}
\keyword{ modulematrix }
\keyword{ module }
