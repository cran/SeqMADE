\name{nbGLMdirperm}
\alias{nbGLMdirperm}
\title{Identify Differential Expression Modules Based on the GLM Method by Shuffling the Phenotypic Variables}
\description{
Identify differential expression modules based on the Generalized Linear Model(GLM), including Group, Module and Direction
variables, then generate the empirical null distribution for the statistic z-values and calculate a empirical estimate of p-value 
of each module in the permutation null distribution by shuffling the phenotypic variables.
}
\usage{
nbGLMdirperm(exprs, case, control, factors, 
             networkModule, modulematrix, N,
			 distribution = c("poisson", "NB")[1])
}
\arguments{
  \item{exprs}{exprs is a data frame or matrix for two groups or conditions, with rows as variables (genes) and columns as samples.}
  \item{case}{case is the sample names in case groups.}
  \item{control}{control is the sample names in control groups.}
  \item{factors}{Factors with three variables including Count, Group, Direction.}
  \item{networkModule}{NetworkModule is the gene sets or modules in the biological network or metabolic pathway, with the 1th column as the module names and the 2th columnn as the gene symbol constituting the module.}
  \item{modulematrix}{Modulematrix is a matrix, in which the indicator variables 1 or 0 represent whether a gene belong to a given module or not.}
  \item{N}{permutation times. If N>0, the permutation step will be implemented. The default value for N is 0.}
  \item{distribution}{a character string indicating the distribution of RNA-Seq count value, default is 'NB'.}

}

\details{
 The GLM method was determined by the distribution of RNA-Seq count value including poisson and Negative Binomial distribution, and there are 
 three indicator variables Group, Module and Direction, in which Module=1 when a gene belongs to the module and Module= 0 otherwise; Group=1 for case values and 
 Group=0 for control values;Direction=1 for up-regulated and Direction=-1 for down-regualted. We therefore construct the contrast vector to
 test the null hypothesis by fitting the GLM and then focus on the interaction term Group*Module*Direction. Then the samples between the two
 conditions will be disturbed and by shuffling the phenotypic variables, we can generate the empirical null distribution for each module. 
 Repeat the above process for N times. Pool all the z score together to form a null distribution of z-value. The corresponding statistical 
 significance (p-value) is estimated against null statistics.
}
\value{
The matrix for the sigificance of each module in differential expression analysis.
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
result <- nbGLMdirperm(exprs,case,control,factors,
                       networkModule, modulematrix,
					   5, distribution="NB")
}