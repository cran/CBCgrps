\name{CBCgrps2.4-package}
\alias{CBCgrps}
\docType{package}
\title{
Compare Baseline Characteristics Between Groups
}
\description{
The package aims to automate the process of comparing Baseline Characteristics between groups.
}
\details{

The DESCRIPTION file:
In clinical studies emplying electronic medical records, the variables to be investigated are usually large in number. It is sometimes cumbersome to compare these variables between two or more groups one by one. I design this package to automate the process of statistical description and bivariate statistical inference. 
}
\author{
Zhongheng Zhang
Department of emergency medicine, Sir Run-Run Shaw Hospital, Zhejiang University School of Medicine, Hangzhou, 310016, China.
\email{zh_zhang1984@zju.edu.cn}
Maintainer: Zhongheng Zhang
}
\references{
Zhang Z. Univariate description and bivariate statistical inference: the first step delving into data.Ann Transl Med. 2016 Mar;4(5):91.
}
\keyword{ bivariate analysis; statisitcal description }
\seealso{
No ther reference.
}
\examples{
data(df)	
a<-twogrps(df,"mort")
}
