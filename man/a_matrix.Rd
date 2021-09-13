\name{a_matrix}
\alias{a_matrix}
\title{Generates the 'a_matrix' (A) of a design}

\description{
The a_matrix function creates an 'a_matrix' (A) with help of design D2 and 'n_matrix' (N) of design D1.
}

\usage{
a_matrix(design2, n_matrix)
}
\arguments{
\item{design2}{A design in matrix format}

\item{n_matrix}{An incidence matrix of a design}
}

\details{
This functions reads the second design (D2) and the incidence matrix or 'n_matrix' (N) of the first design (D1) as input in matrix format. The dimention of the 'a_matrix' (A) is same as the input incidence matrix (N).
}


\value{
It returns matrix (A) with the same dimension as of 'n_matrix' (N).
}

\author{
1) Sayantani Karmakar
2) Md. Ashraful Haque
3) Cini Varghese
4) Seema Jaggi
5) Eldho Varghese
6) Mohd. Harun
}
