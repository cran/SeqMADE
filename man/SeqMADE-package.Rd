\name{SeqMADE-package}
\alias{SeqMADE-package}
\alias{SeqMADE}
\docType{package}
\title{Network Module-Based Model in the Differential Expression Analysis for RNA-Seq}
\description{A network module-based generalized linear model for differential expression analysis with the count-based sequence data from RNA-Seq.
}
\details{
\tabular{ll}{
Package: \tab SeqMADE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-06-27\cr
License: \tab GPL (>2)\cr
LazyLoad: \tab yes\cr
}
The main functions in this package are
\code{\link{Factor}}   A function of constructing the Group variables, Direction variables, and the Count variables,
\code{\link{moduleMatrix}} a function of constructing the modulematrix for all the modules,
\code{\link{nbGLM}} Identify differential expression modules based on the GLM method using Group and Module variables,
\code{\link{nbGLMdir}} Identify differential expression modules based on the Generalized Linear Model(GLM) using Group, Module and Direction variables, and
\code{\link{nbGLMdirperm}} Identify differential expression modules based on the GLM method by shuffling the phenotypic variables.
}
\author{
Mingli Lei, Jia Xu, Li-Ching Huang, Lily Wang, Jing Li
Maintainer: Mingli Lei<leimingli2013@sjtu.edu.cn>
}
\references{
Xu, J., Wang, L. and Li, J. (2014) Biological network module-based model for the analysis of differential expression in shotgun proteomics,
J Proteome Res, 13, 5743-5750.}
\keyword{ package }
\seealso{
glm(),lm()
}
\examples{
data(exprs)
data(networkModule)
case <- c("A1","A2","A3","A4","A5","A6","A7")
control <- c("B1","B2","B3","B4","B5","B6","B7")
factors <- Factor(exprs,case,control)
modulematrix <- moduleMatrix(exprs,networkModule)
Result1<- nbGLM(factors,14,networkModule,modulematrix,distribution="NB")
Result2<- nbGLMdir(factors,14,networkModule,modulematrix,distribution="NB")
Result3<- nbGLMdirperm(exprs,case,control,factors,networkModule,
                       modulematrix,10,distribution="NB")
}


