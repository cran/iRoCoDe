\name{n_matrix}
\alias{n_matrix}
\title{Generates the Incidence Matrix (N) of a Design}

\description{
The n_matrix functions creates an incidence matrix/'n_matrix' (N) from a given design.
}

\usage{
n_matrix(design, treatments, blocks)
}
\arguments{
\item{design}{A design in matrix format}

\item{treatments}{Number of treatments in the given design}

\item{blocks}{Number of blocks in the given design}
}

\details{
This functions reads a design (D1) in matrix format, its no. of treatments (v1) and its no. of blocks (b1) as input. Here the rows of design D1 are treated as blocks. The generated incidence matrix (N) will have the dimension of (v1 x b1).
}

\value{
It returns a matrix (N) with the dimension of (v1 x b1)
}

\author{
1) Sayantani Karmakar
2) Md. Ashraful Haque
3) Cini Varghese
4) Seema Jaggi
5) Eldho Varghese
6) Mohd. Harun
}
