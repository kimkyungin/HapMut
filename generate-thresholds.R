hfile <- "./thresholds.h"

n <- 10:600
perr <- 10:60

lbnd <- matrix(, nrow=length(n), ncol=length(perr))
colnames(lbnd) <- sprintf("err%d", perr)
rownames(lbnd) <- n

th0 <- 10^(-perr/10)
th1 <- 0.5
arg <- commandArgs(TRUE)
decision <- as.numeric(arg[1])

for (i in 1:ncol(lbnd)) {
    e0 <- th0[i]
    e1 <- th1
    for (j in 1:nrow(lbnd)) {
        nj <- n[j]
        pp <- 1/(1 + (e0/e1)^(0:nj)*((1 - e0)/(1 - e1))^(nj-0:nj))
        lbnd[j,i] <- (0:nj)[which.max(pp > decision)]
    }
}

sbnd <- matrix(, nrow=length(n), ncol=length(perr))
colnames(sbnd) <- sprintf("err%d", perr)
rownames(sbnd) <- n

for (i in 1:ncol(sbnd)) {
    e0 <- th0[i]
    e1 <- 0.1  ## upper bound for sequencing error rate
    for (j in 1:nrow(sbnd)) {
        nj <- n[j]
        xj <- 0:nj
        num <- dbinom(xj, nj, e0)
        nxj <- length(xj)
        denum <- numeric(nxj)
        for (k in 1:nxj) {
            ## e1 <- runif(B, min=e0, max=.5)
            ## denum[k] <- mean(dbinom(xj[k], nj, e1))
            denum[k] <- dbinom(xj[k], nj, e1)
        }
        pp <- 1/(1 + num/denum)
        sbnd[j,i] <- xj[which.max(pp > decision)]
    }
}


acat <- function(x, ...) cat(x, file=hfile, append=TRUE, ...)

cat(file=hfile)

acat("const int germ_err[", length(perr), "][", length(n), "] = {\n", sep="")
for (i in 1:length(perr)) {
    acat("{")
    acat(lbnd[,i], sep=",")
    acat("},\n")
}
acat("};\n\n")

acat("const int bial_err[", length(perr), "][", length(n), "] = {\n", sep="")
for (i in 1:length(perr)) {
    acat("{")
    acat(n - lbnd[,i], sep=",")
    acat("},\n")
}
acat("};\n\n")

acat("const int seq_err[", length(perr), "][", length(n), "] = {\n", sep="")
for (i in 1:length(perr)) {
    acat("{")
    acat(sbnd[,i], sep=",")
    acat("},\n")
}
acat("};\n\n")

q <- 0:100
acat("const double phred2error[", length(q), "] = {", sep="")
acat(10^(-q/10), sep=",")
acat("};\n\n")

acat("const int read_deps[", length(n), "] = {", sep="")
acat(n, sep=",")
acat("};\n\n")



