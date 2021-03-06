\name{twogrps}
\alias{twogrps}
\title{
Compare Baseline Characteristics between two groups
}
\description{
The main function of the CBCgrps package. The function compares differences in categorical and continuous variables between two groups. The function automatically judges the distribution of the continuous variable and use appropriate description for them. Wilcoxon rank sum test is employed for non-normal data. 
}
\usage{
twogrps(df, gvar, varlist = NULL, p.rd = 3, 
    skewvar=NULL, norm.rd = 2, 
    sk.rd = 2, tabNA = "no", cat.rd = 0, pnormtest = 0.05,
    maxfactorlevels = 30, minfactorlevels = 10, sim = FALSE, 
    workspace = 2e+05,ShowStatistic = F,ExtractP = 0.05)
}
\arguments{
  \item{df}{
The data frame on which statistical descroption and inference are performed.
}
  \item{gvar}{
The group variable.
}
  \item{p.rd}{
Decimal space of p value to be displayed. If the p value is less than the minimum value of that decimal space, it will print less than that value. For instance, if p.rd = 3 and p = 0.00045, then it will print "<0.001" in the p column of the table. 
}
  \item{varlist}{
Specify a vector of variable names to be compared between groups (i.e. not all variables in the df will be compared and users can choose which variables to be compared). This argument also allows to specify the order in which the variables will appear in the table.
}
  \item{skewvar}{
Specify a vector of variable names which are considered to be not normally distributed. This function is useful for some large datasets where the normality test is too sensitive and users may want to specify skew variables by their own judgement. skewvar is NULL by default. If it is not NULL and skew data were specified by the users, the statistical test for normality is switched off. 
}
  \item{norm.rd}{
Decimal space of normally distributed numeric variables to be displayed.
}
  \item{sk.rd}{
Decimal space of skewed numeric variables to be displayed.
}
  \item{tabNA}{
Weather categorical variables with NA be displayed or not. "no" to be omitted, "ifany" to be displayed. The default value is "no".
}
  \item{cat.rd}{
Decimal space of categorical variables (the proportion) to be displayed.
}
  \item{pnormtest}{
Significance level for the normal test. It is 0.05 by convention (default). However, for some large datasets the test will be too sensitive that only a small deviation in magnitude from the normal distribution will give a p value less than 0.05. In this situation, users can specify smaller significance level. Note that the normality test will be switched off by specifying skewvar argument.
}
  \item{maxfactorlevels}{
The maximum levels for factor variables, the default is 30. The argument is used to avoid treating date or time variables as factor variables. 
}
  \item{minfactorlevels}{
If a numeric variable has only several values, it is treated as categorical variable. The default value is 10.
}
  \item{sim}{
A logical indicating whether to compute p-values by Monte Carlo simulation, in larger than 2 by 2 tables. The default is FALSE.
}
  \item{workspace}{
If the fisher.test() fucntion requires more workspace, it can be defined here. The default is workspace=2e+05.
}
  \item{ShowStatistic}{
logic value for whether showing statistics or not. The defaule is FALSE for not showing statistics. Statistics is used for statistical inference such as F value for Chi-square test and T value for student t test. No statistic will be shown because Fisher's exact test jumps past a test statistic and goes straight to a p-value. 
}
  \item{ExtractP}{
Some variables with p value less than a certain threshoud can be extracted for subsequent multivariate regression modeling. The parameter specifies the threshoud of p value for a variable to be extrtacted.
}
}
\details{
It is common that some categorical variables contain numeric or integer values. For example, the gender variable may contain values 1 and 2, representing male and female respectively. Such a variable can be identified by counting the number of integer values. Thus, the minfactorlevels argument is used to define the minimum value for a variable to be considered as numeric variable.
}
\value{
\item{Table }{A table with string values. The mean and standard error are put in a single cell, and connected by plus and minus symbol.}

\item{VarExtract }{A character vector containing variable names. These extracted variables have p value less than ExtractP in univariate analysis.}
}
\references{
Zhang Z. Univariate description and bivariate statistical inference: the first step delving into data.Ann Transl Med. 2016 Mar;4(5):91. 

Zhang Z, Gayle AA, Wang J, Zhang H, Cardinal-Fernandez P. Comparing baseline characteristics between groups: an introduction to the CBCgrps package. Ann Transl Med. 2017 Dec;5(24):484. doi: 10.21037/atm.2017.09.39.
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
a<-twogrps(df,"mort")
print(a,quote = TRUE)
# define skewed variables manually
print(twogrps(df,"mort",skewvar=c("crp","wbc")),
      quote = TRUE)
}
\keyword{Compare}
\keyword{baseline}
