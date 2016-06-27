\name{nbGLM}
\alias{nbGLM}
\title{Identify Differential Expression Modules Based on the Generalized Linear Model}
\description{
 The algorithm identify differential expression modules using Generalized Linear Model (GLM)
 for differential expression analysis in RNA-Seq data, and in the model two indicator variables Group and 
 Module are adopted to fit the GLM.
}
\usage{
nbGLM(factors, N, networkModule, modulematrix, distribution = c("poisson", "NB")[1])
}
\arguments{
  \item{factors}{Factors with three variables including Count, Group, Direction.}
  \item{N}{The total sample sizes.}
  \item{networkModule}{NetworkModule is the gene sets or modules in the biological network or metabolic pathway, with the 1th column as the module names and the 2th columnn as the gene symbol constituting the module.}
  \item{modulematrix}{Modulematrix is a matrix, in which the indicator variables 1 or 0 represent whether a gene belong to a given module or not.}
  \item{distribution}{a character string indicating the distribution of RNA-Seq count value, default is 'NB'.}
}
\details{
 The GLM method was determined by the distribution of RNA-Seq count value including Poisson and Negative Binomial distribution, and there are two indicator
 variables Group and Module, Module=1 when a gene belongs to the module and Module= 0 otherwise; Group=1 for case values and Group=0 for control values.
 Group * Module represents the interaction effects between Group and Module, and the significance of a module is decided by the interaction and adjusted
 p-values are calculated to correct for multiple testing.
}
\value{
 The nominal pvalue and FDR for the significance of each gene set or module.
}
\author{
Mingli Lei, Jia Xu, Li-Ching Huang, Lily Wang, Jing Li
}
\seealso{
glm()
}
\examples{
data(exprs)
data(networkModule)
case <- c("A1","A2","A3","A4","A5","A6","A7")
control <- c("B1","B2","B3","B4","B5","B6","B7")
factors <- Factor(exprs, case, control) 
modulematrix <- moduleMatrix(exprs,networkModule)
Result <- nbGLM(factors, 14, networkModule, modulematrix, distribution = "NB")
}
\keyword{ nbGLM }
