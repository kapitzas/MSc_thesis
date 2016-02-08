require(devtools)
install.packages("/Users/Simon/Downloads/demogR")

library(demogR)
options(digits = 3, scipen = 1)

data("goodman")
goodman
names(goodman)
ven <- with(goodman, life.table(x = age, nKx = ven.nKx, nDx = ven.nDx)) 
mad <- with(goodman, life.table(x = age, nKx = mad.nKx, nDx = mad.nDx)) 
usa <- with(goodman, life.table(x = age, nKx = usa.nKx, nDx = usa.nDx))

data("thar")
thar.lt <- with(thar, life.table(x = age, nDx = deaths, nKx = count, type = "cohort", iwidth = 1, width12 = c(1, 1)))

thar.lt

ven.mx <- with(goodman, ven.bx/ven.nKx)
mad.mx <- with(goodman, mad.bx/mad.nKx)
A <- leslie.matrix(lx = ven$nLx, mx = ven.mx)
B <- leslie.matrix(lx = mad$nLx, mx = mad.mx)

usa <- with(goodman, life.table(x = age, nKx = usa.nKx, nDx = usa.nDx))
usa.mx <- goodman$usa.bx/goodman$usa.nKx 
C <- leslie.matrix(lx = usa$nLx, mx = usa.mx)
C
no <- goodman$usa.nKx[3:11]
no <- c(sum(goodman$usa.nKx[1:2]), no)
no
tmax <- 20 
N <- project.leslie(A = C, no = no, tmax = tmax)
N
cols <- rgb(0, (10:1)/10, (1:10)/10)
plot(5 * (0:20), N[1, ]/100000, type = "l", xlab = "Years", ylab = "Population Size (x100,000)", ylim = c(16, 175), col = cols[1])
for (i in 2:10) lines(5 * (0:20), N[i, ]/100000, col = cols[i])