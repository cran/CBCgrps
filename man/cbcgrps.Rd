\name{cbcgrps}
\alias{cbcgrps}
\title{
Compare Baseline Characteristics between two groups
}
\description{
The main function of the CBCgrps package
}
\usage{
cbcgrps(df, gvar, p.rd = 3, norm.rd = 2, sk.rd = 2, 
tabNA = "no", cat.rd = 2, maxfactorlevels = 30, 
minfactorlevels = 10, fisher.flg = "yes", workspace = 2e+05)
}
\arguments{
  \item{df}{
The data frame on which statistical descroption and inference are performed.
}
  \item{gvar}{
The group variable.
}
  \item{p.rd}{
Decimal space of p value to be displayed.
}
  \item{norm.rd}{
Decimal space of normally distributed numeric variables to be displayed.
}
  \item{sk.rd}{
Decimal space of skewed numeric variables to be displayed.
}
  \item{tabNA}{
Weather categorical variables with NA be displayed or omitted. "no" to be omitted, "ifany" to be displayed. The default value is "no".
}
  \item{cat.rd}{
Decimal space of categorical variables to be displayed.
}
  \item{maxfactorlevels}{
The maximum levels for factor variables, the default is 30. The argument is used to avoid treating date or time variables as factor variables. 
}
  \item{minfactorlevels}{
If a numeric variable has only several values, it is treated as categorical variable. The default value is 10.
}
  \item{fisher.flg}{
The default is to switch on the Fisher's exact test by "yes". "no" to switch off the test to save worspace.
}
  \item{workspace}{
If the fisher.test() fucntion requires more workspace, it can be defined here. 
}
}
\details{
It is common that some categorical variables contain numeric or integer values. For example, the gender variable may contain values 1 and 2, representing male and female respectively. Such a variable can be identified by counting the number of integer values. Thus, the minfactorlevels argument is used to define the minimum value for a variable to be considered as numeric variable.
}
\value{

 \item{table }{The table with string values. The mean and standard error are put in a single cell, and connected by plus and minus symbol.}
 \item{table.norm }{The table contains only normally distributed numeric variables. Mean and standard deviation are put in separate cells}
 \item{table.skew }{The table contains only skewed numeric variables. Median and interquartile range are put in separate cells}
 \item{table.cat }{The table contains only factor variables. Number and proportion are put in separate cells}
 \item{g1 }{level of group 1}
 \item{g2 }{level of group 2}

}
\references{
Zhang Z. Univariate description and bivariate statistical inference: the first step delving into data.Ann Transl Med. 2016 Mar;4(5):91. 
}
\author{
Zhongheng Zhang
Department of emergency medicine, Sir Run-Run Shaw Hospital, Zhejiang University School of Medicine, Hangzhou, 310016, China.
\email{zh_zhang1984@zju.edu.cn} 
}
\note{
No further notes 
}
\seealso{
No other reference
}
\examples{
data(df)	
a<-cbcgrps(df,"mort")
}
\keyword{Compare}
\keyword{baseline}
