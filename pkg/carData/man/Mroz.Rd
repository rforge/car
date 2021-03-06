\name{Mroz}
\alias{Mroz}
\docType{data}

\title{U.S. Women's Labor-Force Participation}
\usage{Mroz}
\description{
The \code{Mroz} data frame has 753 rows and 8 columns.
The observations, from the Panel Study of Income Dynamics (PSID),
are married women.

}
\format{
  This data frame contains the following columns:
  \describe{
    \item{lfp}{labor-force participation; a factor with levels:
        \code{no};
        \code{yes}.
    }
    \item{k5}{number of children 5 years old or younger.}
    \item{k618}{number of children 6 to 18 years old.}
    \item{age}{in years.}
    \item{wc}{wife's college attendance; a factor with levels:
        \code{no};
        \code{yes}.
    }
    \item{hc}{husband's college attendance; a factor with levels:
        \code{no};
        \code{yes}.
    }
    \item{lwg}{log expected wage rate; for women in the labor force, the actual
      wage rate; for women not in the labor force, an imputed value based on the
      regression of \code{lwg} on the other variables.}
    \item{inc}{family income exclusive of wife's income.}
  }
}
 
\source{
  Mroz, T. A. (1987)
  The sensitivity of an empirical model of married women's hours of work to
  economic and statistical assumptions. 
  \emph{Econometrica} \bold{55}, 765--799.
}

\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  

  Fox, J. (2000)
  \emph{Multiple and Generalized Nonparametric Regression.} Sage.
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
  
  Long. J. S. (1997)
  \emph{Regression Models for Categorical and Limited Dependent Variables.}
  Sage.
}
 
\keyword{datasets}
