\name{get_rep}
\alias{get_rep}
\title{Generates number of replications in a design}

\description{
The get_rep function finds the number of replications of the treaments in a given design.
}

\usage{
get_rep(design)
}

\arguments{
\item{design}{A design in matrix format}
}


\details{
This function reads a design (D) in matrix format and obtain the number of replications (r) of the treatments (v) of the given design where rows are treated as blocks.
}

\value{
It returns a scalar value denoting the number of replications of the treatments in the design.
}

\author{
1) Sayantani Karmakar
2) Md. Ashraful Haque
3) Cini Varghese
4) Seema Jaggi
5) Eldho Varghese
6) Mohd. Harun
}
