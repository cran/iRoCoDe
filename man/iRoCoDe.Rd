\name{iRoCoDe}
\alias{iRoCoDe}
\title{Generate Incomplete Row-Column Designs}

\description{
The 'iRoCoDe' function uses an algorithmic approach for
constructing row-column designs with incomplete rows and columns by using two incomplete block designs. A wide range of incomplete block designs, viz., balanced incomplete block designs (BIBD)/ partially balanced incomplete block designs (PBIBD)/ t-designs can be used in this function to construct these designs.
}

\usage{
iRoCoDe(design1, design2)
}

\arguments{
\item{design1}{first design (D1) having parameters: no. of treaments = v1, No. of blocks = b1, block size = k1, and treatment-replications = r1}

\item{design2}{second design (D2) having parameters: no. of treaments = v2, No. of blocks = b2, block size = k2, and treatment-replications = r2}
}

\details{
The function takes two incomplete block designs D1 and D2 (having parameters v1, b1, k1, r1 and v2, b2, k2, r2) as input in matrix format. The input designs (D1 and D2) must satisfy the parametric relationships: v1 = r2, b1 = b2 and k1 = k2 for generating the new design using this function.

}

\value{
It returns a matrix with dimention v1 x v2.
}

\references{
McSorley, J.P., Phillips, N.C., Wallis, W.D. & Yucas, J.L. (2005). Double arrays, triple arrays and balanced grids. Designs, Codes and Cryptography, 35(1), 21 - 45
}
\author{
1) Sayantani Karmakar
2) Md. Ashraful Haque
3) Cini Varghese
4) Seema Jaggi
5) Eldho Varghese
6) Mohd. Harun
}

\examples{
library(iRoCoDe)
D1= matrix(c(1, 2, 5,
            2, 3, 6,
            3, 4, 1,
            4, 5, 6,
            5, 1, 3,
            6, 1, 4,
            2, 3, 4,
            3, 5, 6,
            6, 1, 2,
            2, 4, 5), nrow = 10, ncol = 3, byrow = TRUE)
D2 = matrix(c(1, 4, 5,
             2, 4, 5,
             1, 3, 4,
             1, 2, 5,
             2, 3, 5,
             1, 3, 5,
             1, 2, 4,
             1, 2, 3,
             2, 3, 4,
             3, 4, 5), nrow = 10, ncol = 3, byrow = TRUE)

D = iRoCoDe(D1,D2)
## Output:
#[1] "Incomplete Row-Column Design"
#     [,1] [,2] [,3] [,4] [,5]
#[1,]    6    9    3    1    5
#[2,]    7    2   10    9    1
#[3,]    8    7    5    3    2
#[4,]    3    4    6    7   10
#[5,]    1    5    8   10    4
#[6,]    4    8    9    2    6
}
