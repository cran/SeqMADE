\name{Factor}
\alias{Factor}
\title{
 Construction of Variable Factors
}
\description{
  A function of constructing the Group variables, Direction variables, and the Count variables.
  }
\usage{
Factor(exprs, case, control)
}
\arguments{
  \item{exprs}{
  exprs is a data frame or matrix for two groups or conditions, with rows as variables (genes) and columns as samples.
}
  \item{case}{
  case is the sample names in case groups.
}
  \item{control}{
  control is the sample names in control groups.
}
}
\details{
  Two indicator variables Group and Direction corresponding to the different groups and the direction of the gene expression changes in the context of an RNA-Seq experiment, respectively.
   And in this part, 1 represents that a gene belongs to case group or up-regulated and 0 represents a gene belongs to control group or down-regulated. Besides, Count variables are the 
   expression value in different samples for genes.
}
\value{
 \item{Count}{The gene expression count values.}
 \item{Group}{The indicator variables represent that whether a gene belongs to case group or not.}
 \item{Direction}{The indicator variables represent that a gene is up-regulated or down-regulated.}
}
\author{
Mingli Lei, Jia Xu, Li-Ching Huang, Lily Wang, Jing Li
}
\examples{
data(exprs)
case <- c("A1","A2","A3","A4","A5","A6","A7")
control <- c("B1","B2","B3","B4","B5","B6","B7")
factors <- Factor(exprs, case, control) 
}
\keyword{ Group }
\keyword{ Direction }
\keyword{ Count }
