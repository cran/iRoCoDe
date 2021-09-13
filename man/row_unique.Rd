\name{row_unique}
\alias{row_unique}
\title{Checks uniqness of elements in the rows of a matrix}

\description{
The row_unique function checks whether each row contains unique treatments.
}

\usage{
row_unique(a_matrix, treatments, rows, cols)
}
\arguments{
\item{a_matrix}{An 'a_matrix' of the first design (D1)}
\item{treatments}{The treatments of the second design (D2)}
\item{rows}{Number of rows in the given 'a_matrix' }
\item{cols}{Number of columns in the given 'a_matrix'}
}

\details{
This function reads an 'a_matrix' (A) generated from two designs (D1 and D2) and checks whether each row contains all the treatments from D2 exactly once.
}


\value{
It returns boolean value (either True or False).
}

\author{
1) Sayantani Karmakar
2) Md. Ashraful Haque
3) Cini Varghese
4) Seema Jaggi
5) Eldho Varghese
6) Mohd. Harun
}
