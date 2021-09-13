\name{rearrange_mat}
\alias{rearrange_mat}
\title{Rearranges the given matrix}

\description{
The rearrange_mat function rearranges the column elements of a given 'a_matrix'.
}

\usage{
rearrange_mat(a_matrix, treatments)
}
\arguments{
\item{a_matrix}{A a_matrix of the first design (D1)}

\item{treatments}{The number of treatments of the second design (D2)}
}

\details{
This function reads an 'a_matrix' (A) generated from two designs (D1 and D2) and rearranges the column elements based on the treatments numbers in the second design (D2). It iteratively rearrange the items row-wise until the new A is generated with each row containing all the treatments of D2 exaclty once.
}


\value{
It returns new 'a_matrix' (A) with each row containing all the treatments of D2 exaclty once.
}

\author{
1) Sayantani Karmakar
2) Md. Ashraful Haque
3) Cini Varghese
4) Seema Jaggi
5) Eldho Varghese
6) Mohd. Harun
}
