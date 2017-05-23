library('mvtnorm')
library('mvnfast')
library('MCMCpack')


rinvchisq <- function (N, nu, tausq) { nu * tausq / rchisq(N, nu) }

chainest <- function (chain, nTrim = 200) {
    ## chain should be a list. Last dimension is the sample index
    trimmask <- 1:nTrim

    est <- Map( function (x) {
        if (!is.matrix(x)) {
            mean(x)
        } else {
            margin_max <- length(dim(as.matrix(x)))-1
            dimens <- 1:margin_max
            apply(x, dimens, mean)
        }
    }, chain )
    est
}

goutlie <- function (x, nsamp = 800,
                                 mu0 = 0, tausq0 = 2e-2,
                                 sigsq0 = 3^2, nu0 = 2,
                                 eps = 0.05, xi = 1e-4) {
    ## x is a vector of data points.
    n <- length(x)

    sigsq <- rep(-1e6, nsamp+1)
    sigsq[1] <- (max(x) - min(x))^2 / 4
    mu <- rep(-1e6, nsamp+1)
    mu[1] <- 0

    del <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    del[,1] <- 0
    A <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    A[,1] <- 0

    rmu_sigsq_y <- function (sigsq, y, tausq0, mu0) {
        n <- length(y)
        pressq <- tausq0 + n/sigsq
        rnorm( 1, (tausq0*mu0 + sum(y)/sigsq) / pressq, sqrt(pressq)^-1 )
    }

    rsigsq_mu_y <- function (mu, y, nu0, sigsq0) {
        n <- length(y)
        nu <- nu0+n
        rinvchisq( 1, nu, (nu0*sigsq0 + sum((y - mu)^2)) / nu )
    }

    for (r in 2:(nsamp+1)) {
        correctx <- x - c(del[,r-1] * A[,r-1])
        mu[r] <- rmu_sigsq_y(sigsq[r-1], correctx, tausq0, mu0)
        sigsq[r] <- rsigsq_mu_y(mu[r], correctx, nu0, sigsq0)

        outlieprob <- (tmp <- eps * dnorm(x, mu[r]+A[,r-1], sqrt(sigsq[r]))) /
            (tmp + (1-eps)*dnorm(x, mu[r], sqrt(sigsq[r])))
        del[,r] <- rbinom(n, 1, outlieprob)

        outlying <- del[,r] == 1
        for (i in which(outlying))
            A[i,r] <- rmu_sigsq_y(sigsq[r], x[i] - mu[r], xi, 0)
        A[!outlying,r] <- rnorm(sum(!outlying), 0, sqrt(1/xi))
    }

    list(del = del[,-1], sigsq = sigsq[-1], mu = mu[-1], A = A[,-1])
}

## With epsilon beta distributed.
goutlie.betaeps <- function (x, nsamp = 800,
                                mu0 = 0, tausq0 = 2e-2,
                                sigsq0 = 3^2, nu0 = 2,
                                epsmean = 0.05, xi = 1e-4) {
    ## x is a vector of data points.
    n <- length(x)

    sigsq <- rep(-1e6, nsamp+1)
    sigsq[1] <- (max(x) - min(x))^2 / 4
    mu <- rep(-1e6, nsamp+1)
    mu[1] <- 0

    del <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    del[,1] <- 0
    A <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    A[,1] <- 0

    eps <- rep(-1e6, nsamp+1)
    eps[1] <- epsmean

    ## Convert prior mean of epsilon to a beta distribution, with
    ## P(epsilon < 0.5) = 0.99
    alpha <- uniroot( function (alpha) {
        beta <- (1/epsmean - 1) * alpha
        qbeta(0.99, alpha, beta) - 0.5
    }, interval=c(1e-7, 1e7))$root
    beta <- (1/epsmean - 1) * alpha


    rmu_sigsq_y <- function (sigsq, y, tausq0, mu0) {
        n <- length(y)
        pressq <- tausq0 + n/sigsq
        rnorm( 1, (tausq0*mu0 + sum(y)/sigsq) / pressq, sqrt(pressq)^-1 )
    }

    rsigsq_mu_y <- function (mu, y, nu0, sigsq0) {
        n <- length(y)
        nu <- nu0+n
        rinvchisq( 1, nu, (nu0*sigsq0 + sum((y - mu)^2)) / nu )
    }

    for (r in 2:(nsamp+1)) {
        correctx <- x - c(del[,r-1] * A[,r-1])
        mu[r] <- rmu_sigsq_y(sigsq[r-1], correctx, tausq0, mu0)
        sigsq[r] <- rsigsq_mu_y(mu[r], correctx, nu0, sigsq0)

        outlieprob <- (tmp <- eps[r-1] * dnorm(x, mu[r]+A[,r-1], sqrt(sigsq[r]))) /
            (tmp + (1-eps[r-1])*dnorm(x, mu[r], sqrt(sigsq[r])))
        del[,r] <- rbinom(n, 1, outlieprob)

        outlying <- del[,r] == 1
        for (i in which(outlying))
            A[i,r] <- rmu_sigsq_y(sigsq[r], x[i] - mu[r], xi, 0)
        A[!outlying,r] <- rnorm(sum(!outlying), 0, sqrt(1/xi))

        n_out <- sum(outlying)
        eps[r] <- rbeta(1, alpha + n_out, beta + n - n_out)
    }

    list(del = del[,-1], sigsq = sigsq[-1], mu = mu[-1], A = A[,-1])
}


############## Normal mixture version

## With epsilon beta distributed.
goutlie.betaeps <- function (x, nsamp = 800,
                                mu0 = 0, tausq0 = 2e-2,
                                sigsq0 = 3^2, nu0 = 2,
                                epsmean = 0.05, xi = 1e-4) {
    ## x is a vector of data points.
    n <- length(x)

    sigsq <- rep(-1e6, nsamp+1)
    sigsq[1] <- (max(x) - min(x))^2 / 4
    mu <- rep(-1e6, nsamp+1)
    mu[1] <- 0

    del <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    del[,1] <- 0
    A <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    A[,1] <- 0

    eps <- rep(-1e6, nsamp+1)
    eps[1] <- epsmean

    ## Convert prior mean of epsilon to a beta distribution, with
    ## P(epsilon < 0.5) = 0.99
    alpha <- uniroot( function (alpha) {
        beta <- (1/epsmean - 1) * alpha
        qbeta(0.99, alpha, beta) - 0.5
    }, interval=c(1e-7, 1e7))$root
    beta <- (1/epsmean - 1) * alpha

    rmu_sigsq_y <- function (sigsq, y, tausq0, mu0) {
        n <- length(y)
        pressq <- tausq0 + n/sigsq
        rnorm( 1, (tausq0*mu0 + sum(y)/sigsq) / pressq, sqrt(pressq)^-1 )
    }

    rsigsq_mu_y <- function (mu, y, nu0, sigsq0) {
        n <- length(y)
        nu <- nu0+n
        rinvchisq( 1, nu, (nu0*sigsq0 + sum((y - mu)^2)) / nu )
    }

    for (r in 2:(nsamp+1)) {
        correctx <- x - c(del[,r-1] * A[,r-1])
        mu[r] <- rmu_sigsq_y(sigsq[r-1], correctx, tausq0, mu0)
        sigsq[r] <- rsigsq_mu_y(mu[r], correctx, nu0, sigsq0)

        outlieprob <- (tmp <- eps[r-1] * dnorm(x, mu[r]+A[,r-1], sqrt(sigsq[r]))) /
            (tmp + (1-eps[r-1])*dnorm(x, mu[r], sqrt(sigsq[r])))
        del[,r] <- rbinom(n, 1, outlieprob)

        outlying <- del[,r] == 1
        for (i in which(outlying))
            A[i,r] <- rmu_sigsq_y(sigsq[r], x[i] - mu[r], xi, 0)
        A[!outlying,r] <- rnorm(sum(!outlying), 0, sqrt(1/xi))

        n_out <- sum(outlying)
        eps[r] <- rbeta(1, alpha + n_out, beta + n - n_out)
    }

    list(del = del[,-1], sigsq = sigsq[-1], mu = mu[-1], A = A[,-1])
}

## Mixture of normals with outlier detection
gibmix.out <- function (y, k, nsamp = 1000,
                    mu0 = numeric(k), tausq0 = rep(1e-7,k),
                    nu0 = rep(4,k), sigsq0 = rep(var(y)/k, k),
                    alpha0 = rep(0.5,k),
                    epsmean = 0.05, xi = 1e-5) {
    n <- length(y)

    xisd <- sqrt(1/xi)                  # Standard dev. version of xi.
    rmu_sigsq_y <- function (sigsq, y, tausq0, mu0) {
        n <- length(y)
        if (n == 0)
            rnorm(1, mu0, sqrt(tausq0)^-1)
        else {
            pressq <- tausq0 + n/sigsq
            rnorm( 1, (tausq0*mu0 + sum(y)/sigsq) / pressq, sqrt(pressq)^-1 )
        }
    }

    rsigsq_mu_y <- function (mu, y, nu0, sigsq0) {
        n <- length(y)
        if (n == 0)
            rinvchisq(1, nu0, sigsq0)
        else {
            nu <- nu0 + n
            rinvchisq( 1, nu, (nu0*sigsq0 + sum((y - mu)^2)) / nu )
        }
    }
    rpi <- function (z) { (g <- rgamma(k, tabulate(z) + alpha0, 1)) / sum(g) }

    ## Allocate memory for the chain
    pi_chn <- cbind( rep(1/k,k), matrix(NA, k, nsamp) )
    z_chn <- cbind( sample(k, length(y), T), matrix(NA, length(y), nsamp) )
    mu_chn <- cbind( as.matrix(mu0), matrix(NA, k, nsamp) )
    sigsq_chn <- cbind( sigsq0, matrix(NA, k, nsamp) )

    del <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    del[,1] <- 0
    A <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    A[,1] <- 0

    eps <- rep(-1e6, nsamp+1)
    eps[1] <- epsmean

    ## Convert prior mean of epsilon to a beta distribution, constraint
    ## to P(epsilon < 0.5) = 0.99 and E(epsilon) = epsmean
    alpha <- uniroot( function (alpha) {
        beta <- (1/epsmean - 1) * alpha
        qbeta(0.99, alpha, beta) - 0.5
    }, interval=c(1e-7, 1e7))$root
    beta <- (1/epsmean - 1) * alpha

    for (r in 2:(nsamp+1)) {
        pi_chn[,r] <- pi_new <- rpi(z_chn[,r-1])

        clusprob <- sapply( 1:k, function (c) {
            pi_new[c] * dnorm(y, mu_chn[c,r-1], sqrt(sigsq_chn[c,r-1]))
        })
        if ( sum(is.na(clusprob)) > 0) {
            browser()
        }
        for (i in 1:n) {
            z_chn[i,r] <- sample(k, 1, T, clusprob[i,])
        }

        ## c'th row is the mask for the c'th cluster
        clustmasks <- t(matrix(z_chn[,r])[,rep(1, k)]) == 1:k

        correcty <- y - c(del[,r-1] * A[,r-1])
        mu_chn[,r] <- mu_new <- sapply(1:k, function (c) {
            rmu_sigsq_y(sigsq_chn[c,r-1], correcty[clustmasks[c,]],
                        mu0[c], tausq0[c])
        })
        sigsq_chn[,r] <- sapply(1:k, function (c) {
            rsigsq_mu_y(mu_new[c], correcty[clustmasks[c,]], nu0[c], sigsq0[c])
        })

        muobs <- mu_chn[z_chn[,r],r]
        sigsqobs <- sigsq_chn[z_chn[,r],r]
        sdobs <- sqrt(sigsqobs)
        outlieprob <- (tmp <- eps[r-1] * dnorm(y, muobs+A[,r-1], sdobs)) /
            (tmp + (1-eps[r-1]) * dnorm(y, muobs, sdobs) )
        del[,r] <- rbinom(n, 1, outlieprob)

        outlying <- del[,r] == 1
        n_out <- sum(outlying)

        for (i in which(outlying))
            A[i,r] <- rmu_sigsq_y( sigsqobs[i], y[i] - muobs[i], xi, 0)
        A[!outlying] <- rnorm(n - n_out, 0, xisd)

        eps[r] <- rbeta(1, alpha + n_out, beta + n - n_out)
    }
    list(pi = pi_chn, z = z_chn, mu = mu_chn, sigsq = sigsq_chn,
         eps = eps[-1], del = del[,-1], A = A[,-1])
}



## Usual Bayesian linear regression
blm <- function (x, y,
                 mu0 = numeric(ncol(x)+1),
                 omega0 = diag(rep(1e-7, ncol(x)+1)),
                 nu0 = 1e-7, sigsq0 = 10) {
    x1 <- cbind(1, as.matrix(x))        # Add constant terms to x
    y <- as.matrix(y)
    d <- t(x1) %*% x1
    omega_n <- d + omega0
    mu_n <- solve(omega_n) %*%
        ( d %*% solve(d) %*% t(x1) %*% y + omega0 %*% mu0 )
    nu_n <- nu0 + length(y)
    sigsq_n <- c( ( nu0 * sigsq0 + t(y) %*% y + t(mu0) %*% omega0 %*% mu0 -
                 t(mu_n) %*% omega_n %*% mu_n ) / nu_n )

    ## Return the parameters for the t distribution of beta
    t_sigma <- sigsq_n * solve(omega_n) # Sigma for multivariate t-distribution
    mu_n <- c(mu_n)
    colnames(t_sigma) <- rownames(t_sigma) <- names(mu_n) <- c('const', colnames(x))
    list(mu = mu_n, df = nu_n, sigma = t_sigma,
         omega_n = omega_n, sigsq_n= sigsq_n)
}


## Only one dimensional
blm.out <- function (x, y, nsamp = 1000,
                     mu0 = numeric(ncol(x)+1),
                     omega0 = diag(rep(1e-7, ncol(x)+1)),
                     nu0 = 1e-7, sigsq0 = 10,
                     epsmean = 0.05, xi = 1e-4) {
    x1 <- cbind(1, x)
    n <- nrow(x)

    del <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    del[,1] <- 0
    A <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    A[,1] <- 0

    eps <- rep(-1e6, nsamp+1)
    eps[1] <- epsmean

    ## Convert prior mean of epsilon to a beta distribution, with
    ## P(epsilon < 0.5) = 0.99
    alpha <- uniroot( function (alpha) {
        beta <- (1/epsmean - 1) * alpha
        qbeta(0.99, alpha, beta) - 0.5
    }, interval=c(1e-7, 1e7))$root
    beta <- (1/epsmean - 1) * alpha

    resid <- matrix(-1e6, nrow(x), nsamp+1)
    betcoef <- matrix(-1e6, ncol(x1), nsamp+1)

    rmu_sigsq_y <- function (sigsq, y, tausq0, mu0) {
        n <- length(y)
        pressq <- tausq0 + n/sigsq
        rnorm( 1, (tausq0*mu0 + sum(y)/sigsq) / pressq, sqrt(pressq)^-1 )
    }

    pred <- rep(1e-6, nrow(x))
    for (r in 2:(nsamp+1)) {
        inlying <- del[,r-1] == 0

        blmres <- blm(x[inlying,,drop=F], y[inlying,drop=F],
                      mu0 = mu0, omega0 = omega0, nu0 = nu0, sigsq0 = sigsq0)
        residvar <- rinvchisq(1, blmres$df, blmres$sigsq_n)
        betcoef[,r] <- c(rmvn(1, blmres$mu, residvar * solve(blmres$omega_n)))

        resid[,r] <- c(x1 %*% betcoef[,r]) - y
        ## resid[!inlying,r] <- rnorm(sum(!inlying), A[,r-1], sqrt(rinvchisq(1,nu0,sigsq0)))
        correctresid <- resid[,r] - del[,r-1] * A[,r-1]

        residsd <- sqrt(residvar)
        outlieprob <- pmin(1, exp((tmp <- log(eps[r-1]) + dnorm(resid[,r], A[,r-1], residsd, log=T)) -
            log(exp(tmp) + exp(log(1-eps[r-1]) + dnorm(resid[,r], 0, residsd, log=T)))))

        del[,r] <- rbinom(n, 1, outlieprob)

        outlying <- del[,r] == 1
        n_out <- sum(outlying)
        for (i in which(outlying))
            A[i,r] <- rmu_sigsq_y(residvar, resid[i,r], xi, 0)
        A[!outlying,r] <- rnorm(n - n_out, 0, sqrt(1/xi))

        eps[r] <- rbeta(1, alpha + n_out, beta + n - n_out)
    }

    list(bet = betcoef[,-1], del = del[,-1], A = A[,-1])
}


gmvmix.out <- function (x,
                        k = 2,
                        mu0 = matrix(0, nrow(x), k),
                        lamb0 = diag(1e-9, nrow(x)),
                        kappa0 = 1,
                        nu0 = nrow(x)+1,
                        alpha0 = rep(5, k),
                        epsmean = 0.05,
                        xi = diag(1e-4, nrow(x)),
                        nsamp = 1500,
                        detect.outlier = T) {
    ## x: an observation per column
    n <- ncol(x)
    p <- nrow(x)
    xivar <- solve(xi)

    del <- matrix(-1e6, nrow = n, ncol = nsamp + 1)
    del[,1] <- 0
    A <- array(-1e6, dim = c(p, n, nsamp + 1))
    A[,,1] <- 0

    eps <- rep(-1e6, nsamp+1)
    eps[1] <- epsmean

    ## Convert prior mean of epsilon to a beta distribution, with
    ## P(epsilon < 0.5) = 0.99
    alpha <- uniroot( function (alpha) {
        beta <- (1/epsmean - 1) * alpha
        qbeta(0.99, alpha, beta) - 0.5
    }, interval=c(1e-7, 1e7))$root
    beta <- (1/epsmean - 1) * alpha

    mu <- array(-1e6, dim = c(p, k, nsamp+1))
    for (c in 1:k)
        mu[,c,1] <- rmvn(1, mu0[,c], solve(lamb0)/kappa0)

    sig <- array(-1e6, dim = c(p, p, k, nsamp+1))
    for (c in 1:k)
        sig[,,c,1] <- solve(lamb0)
    z <- array(-1e6, dim = c(ncol(x), nsamp+1))
    z[,1] <- sample(k,n,T)
    w <- matrix(-1e6, k, nsamp+1)

    newsig <- sig[,,,1]
    newmu <- mu[,,1]
    newz <- z[,1]
    neww <- w[,1]
    newA <- A[,,1]
    newdel <- del[,1]

    nk <- rep(1e-6, k)     # Number of element in cluster
    allcorrectx <- array(1e-6, dim=dim(x))
    for (r in 2:(nsamp+1)) {
        for (c in 1:k) {
            clustmask <- newz == c
            nk[c] <- nc <- sum(clustmask)
            allcorrectx[,clustmask] <- correctx <-
                x[,clustmask] - t(newdel[clustmask] * t(newA[,clustmask]))

            if (nc == 0) {
                sig[,,c,r] <- newsig[,,c] <- solve(lamb0)
                mu[,c,r] <- newmu[,c] <- rmvn(1, mu0[,c], newsig[,,c]/kappa0)
                next()
            }
            if (nc == 1) {
                sig[,,c,r] <- newsig[,,c] <- solve(lamb0)
                mu[,c,r] <- newmu[,c] <- rmvn(1, c(correctx), newsig[,,c]/kappa0)
                next()
            }

            kappa_n <- kappa0 + nc

            correctxbar <- rowMeans(correctx)

            Q0 <- (nc - 1) * cov(t(correctx))
            sig[,,c,r] <- newsig[,,c] <-
                riwish( kappa_n,
                       lamb0 + Q0 +
                       kappa0*nc/(kappa0+nc) *
                       tcrossprod(correctxbar - mu0[,c]))

            mu[,c,r] <- newmu[,c] <-
                rmvn(1,
                     mu= (nc*correctxbar + kappa0*mu0[,c]) / (kappa0 + nc),
                     sigma= newsig[,,c]/kappa_n)
        }

        w[,r] <- neww <- rdirichlet(1, nk + alpha0)

        for (i in 1:n) {
            prob <- sapply(1:k, function (c) {
                neww[c] * dmvn(x[,i], newmu[,c], newsig[,,c])
            })
            z[i,r] <- newz[i] <- sample(k, 1, T, prob = prob)
        }

        sigx <- newsig[,,newz]
        mux <- newmu[,newz]
        outlieprob <- (tmp <- eps[r-1] * dmvn2(t(x), mux + newA, sigx)) /
            (tmp + (1-eps[r-1]) * dmvn2(t(x), mux, sigx))
        del[,r] <- newdel <- rbinom(n, 1, outlieprob)
        if (! detect.outlier) {
            ## Very dirty hack that works...
            del[,r] <- newdel <- rep(0, n)
        }

        outlying <- newdel == 1
        n_out <- sum(outlying)
        n_in <- n - n_out
        for (i in which(outlying))
            A[,i,r] <- rmvn(1, x[,i] - newmu[,newz[i]], solve(xi + solve(sigx[,,i])) )
        A[,!outlying,r] <- t(rmvn(n_in, rep(0,p), xivar))
        newA <- A[,,r]

        eps[r] <- rbeta(1, alpha + n_out, beta + n_in)
    }
    list(mu = mu[,,-1], sig = sig[,,,-1], z = z[,-1],
         w = w[,-1], del = del[,-1], A = A[,,-1], eps = eps[-1])
}

dmvn2 <- function (x, mu, sig, ...) {
    ## x is a vector per row to match the interface of dmvn
    n <- nrow(x)
    sapply(1:n, function (i) {
        c(dmvn(x[i,], mu[,i], sig[,,i]))
    })
}




