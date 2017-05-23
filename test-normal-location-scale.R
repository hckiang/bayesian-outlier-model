test_norm_outliers <- function () {
    library('lattice')
    library('latticeExtra')

    normvars <- c(rt(200, 19))
    normvars <- normvars[sample(length(normvars))]

    normres <- goutlie(normvars, 3000, xi = 0.01)
    normres_est <- chainest(normres)
    plot1 <- xyplot(normres_est$del ~ normvars, pch = 20,
           xlab = expression(X[i]), ylab = expression(E(delta[i])),
           main = 'Probability of being outlier')

    normres2 <- goutlie.betaeps(normvars, 3000, xi = 0.01)
    normres2_est <- chainest(normres)
    plot2 <- xyplot(normres2_est$del ~ normvars, pch = 20,
                    xlab = expression(X[i]), ylab = expression(E(delta[i])))


    plot.eps05 <- c('Known Epsilon (=0.05)' = plot1,
                    'Unknown Epsilon (E.V. at 0.05)' = plot2)

    normres <- goutlie(normvars, 3000, eps = 0.15, xi = 0.01)
    normres_est <- chainest(normres)
    plot3 <- xyplot(normres_est$del ~ normvars, pch = 20,
           xlab = expression(X[i]), ylab = expression(E(delta[i])),
           main = 'Probability of being outlier (epsilon=0.15)')

    normres2 <- goutlie.betaeps(normvars, 3000, epsmean = 0.15, xi = 0.01)
    normres2_est <- chainest(normres)
    plot4 <- xyplot(normres2_est$del ~ normvars, pch = 20,
                    xlab = expression(X[i]), ylab = expression(E(delta[i])))
    plot.eps15 <- c('Known Epsilon (=0.15)' = plot3,
                    'Unknown Epsilon (E.V. at 0.15)' = plot4)

    c(plot.eps05, plot.eps15, y.same = T, x.same = T)
    ## print( plot.eps05, split=c(1,1,1,2), more=T )
    ## print( plot.eps15, split=c(1,2,1,2))
}
test_norm_outliers()

rainfall_test <- function () {
    library('lattice')
    library('latticeExtra')
    rainfall <- read.csv('rainfall.dat', sep = '', header=F)[,1]

    gmixfit <- gibmix.out(rainfall, 2, nsamp = 800,
                          xi = 1000^(-2))

    gmixev <- list(mu = apply(gmixfit$mu[,-(1:100)], 1, mean),
                   sigsq = apply(gmixfit$sigsq[,-(1:50)], 1, mean),
                   pi = apply(gmixfit$pi[,-(1:100)], 1, mean),
                   eps = mean(gmixfit$eps[-(1:100)]),
                   del = apply(gmixfit$del[,-(1:100)], 1, mean),
                   A = apply(gmixfit$A[,-(1:100)], 1, mean))

    plot.fit <- histogram( rainfall, type='d', col=0, nint=52,
              main = 'Daily rain fall somewhere in Washinton',
              xlab = 'Rain fall (1/100 inch)',
              ylab = 'Density',
              panel = function (x, ...) {
                  allargs <- list(...)
                  xs <- seq(allargs$breaks[1], tail(allargs$breaks, 1), length.out = 2000)
                  panel.histogram(x, ...)
                  den1 <- dnorm(xs, gmixev$mu[1], sqrt(gmixev$sigsq[1]))
                  den2 <- dnorm(xs, gmixev$mu[2], sqrt(gmixev$sigsq[2]))
                  llines(xs, gmixev$pi %*% rbind(den1, den2), lty = 1)
                  panel.xyplot(rainfall, 0, pch = 20)
              }, key=list(space="top",
                          lines = list(col=c(4), lty=c(1)),
                          text=list(c('Mixture from Gibbs'))))

    plot.outly <- xyplot(gmixev$del ~ rainfall,
                         xlab = 'Rain fail (1/100 inch)',
                         ylab = expression(E(delta[i])),
                         main = 'Probability of observation being outlier')
    print(plot.fit, split = c(1,1,1,2), more = T)
    print(plot.outly, split = c(1,2,1,2))
}
rainfall_test()

many_mixture_test <- function () {
    library('lattice')
    library('latticeExtra')
    ##set.seed(1)
    rainfall <- c(rt(1000,2)+100, rt(400,2)*2.5 -500,
                  rt(900,2) +200, runif(15, -600, 300))

    gmixfit <- gibmix.out(rainfall, 3, nsamp = 1500,
                          xi = 1e-5, epsmean = 0.15)

    gmixev <- list(mu = apply(gmixfit$mu[,-(1:400)], 1, mean),
                   sigsq = apply(gmixfit$sigsq[,-(1:400)], 1, mean),
                   pi = apply(gmixfit$pi[,-(1:400)], 1, mean),
                   eps = mean(gmixfit$eps[-(1:400)]),
                   del = apply(gmixfit$del[,-(1:400)], 1, mean),
                   A = apply(gmixfit$A[,-(1:400)], 1, mean))

    plot.fit <- histogram( rainfall, type='d', col=0, nint=52,
              main = 'Daily rain fall somewhere in Washinton',
              xlab = 'Rain fall (1/100 inch)',
              ylab = 'Density',
              panel = function (x, ...) {
        allargs <- list(...)
        xs <- seq(allargs$breaks[1], tail(allargs$breaks, 1), length.out = 2000)
        den1 <- dnorm(xs, gmixev$mu[1], sqrt(gmixev$sigsq[1]))
        den2 <- dnorm(xs, gmixev$mu[2], sqrt(gmixev$sigsq[2]))
        den3 <- dnorm(xs, gmixev$mu[3], sqrt(gmixev$sigsq[3]))
        llines(xs, (1 - gmixev$eps) * gmixev$pi %*% rbind(den1, den2, den3),
               lty = 1)
        panel.xyplot(rainfall, 0, pch = 20)
    })

    plot.outly <- xyplot(gmixev$del ~ rainfall,
                         xlab = 'Rain fail (1/100 inch)',
                         ylab = expression(E(delta[i])),
                         main = 'Probability of observation being outlier')
    print(plot.fit, split = c(1,1,1,2), more = T)
    print(plot.outly, split = c(1,2,1,2))
}
many_mixture_test()



plotreg <- function () {
    library('lattice')
    library('latticeExtra')
    library('mvnfast')

    dispmpg <- t(mtcars[,c('mpg', 'disp')])
    dispmpg[1,] <- dispmpg[1,] / 1.60934

    x <- dispmpg
    ## Add artificial outlier
    x.err <- cbind(x, matrix(c(50, 150)), matrix(c(30, 103)))

    blmres <- blm(t(x['disp',,drop=F]), t(x[1,,drop=F]),
              sigsq0 = 5)
    nsample <- 5000
    sim <- rmvt(nsample, blmres$mu, sigma = blmres$sigma, blmres$df)
    colnames(sim) <- names(blmres$mu)

    blmres.err <- blm(t(x.err['disp',,drop=F]), t(x.err[1,,drop=F]),
              sigsq0 = 5)
    sim.err <- rmvt(nsample, blmres.err$mu, blmres.err$sigma, blmres.err$df)

    colnames(sim.err) <- names(blmres$mu)

    blmout <- blm.out(t(x.err['disp',,drop=F]), t(x.err[1,,drop=F]),
                      nsamp=5000, xi=1e-5)


    blmout$bet <- blmout$bet[,-(1:1000)]
    rownames(blmout$bet) <- names(blmres$mu)

    plot.orig <- xyplot(x.err['mpg',] ~ x.err['disp',],
                        col = c(rep(1, ncol(x)), 2,2),
                        main = 'Fuel efficiency of some cars vs. Engine displacement (mtcars data set in R)',
                        xlab = 'Engine displacement (Cubic inch)',
                        ylab = 'Kilometers per gallon',
                        panel = function (x, y, ...) {
                            allargs <- list(...)
                            minx <- min(x)
                            maxx <- max(x)
                            xs <- seq(min(x), max(x), length.out = 300)
                            predxs <- cbind(1, xs)

                            curves <- sapply( 1:nsample, function (r) t(sim[r,c('const', 'disp')]) %*% t(predxs) )
                            meancurve <- apply( curves, 1, mean )
                            intervals <- apply( curves, 1, function (f_est) {
                                quantile(f_est, probs = c(0.025, 0.975)) } )
                            llines(x = xs, y = meancurve, col = 1)
                            panel.polygon(x = c(xs, rev(xs)),
                                          y = c(intervals[1,], rev(intervals[2,])),
                                          col = 1, alpha = 0.1, border = 0)

                            curves.err <- sapply( 1:nsample, function (r) t(sim.err[r,c('const', 'disp')]) %*% t(predxs) )
                            meancurve.err <- apply( curves.err, 1, mean )
                            intervals.err <- apply( curves.err, 1, function (f_est) {
                                quantile(f_est, probs = c(0.025, 0.975)) } )
                            llines(x = xs, y = meancurve.err, col = 'red')
                            panel.polygon(x = c(xs, rev(xs)),
                                          y = c(intervals.err[1,], rev(intervals.err[2,])),
                                          col = 'red', alpha = 0.1, border = 0)

                            panel.points(x, y,...) })

    plot.corr <- xyplot(x.err['mpg',] ~ x.err['disp',],
                        col = c(rep(1, ncol(x)), 2,2),
                        panel = function (x, y, ...) {
                            allargs <- list(...)
                            minx <- min(x)
                            maxx <- max(x)
                            xs <- seq(min(x), max(x), length.out = 300)
                            predxs <- cbind(1, xs)

                            curves.err <- sapply( 1:nsample, function (r) t(sim.err[r,c('const', 'disp')]) %*% t(predxs) )
                            meancurve.err <- apply( curves.err, 1, mean )
                            intervals.err <- apply( curves.err, 1, function (f_est) {
                                quantile(f_est, probs = c(0.025, 0.975)) } )
                            llines(x = xs, y = meancurve.err, col = 'red')
                            panel.polygon(x = c(xs, rev(xs)),
                                          y = c(intervals.err[1,], rev(intervals.err[2,])),
                                          col = 'red', alpha = 0.1, border = 0)

                            curves.out <- sapply(1:ncol(blmout$bet),
                                                 function (r) blmout$bet[c('const', 'disp'),r] %*% t(predxs) )
                            meancurve.out <- apply( curves.out, 1, mean )
                            intervals.out <- apply( curves.out, 1, function (f_est) {
                                quantile(f_est, probs = c(0.025, 0.975)) } )
                            llines(x = xs, y = meancurve.out, col = '#006400')
                            panel.polygon(x = c(xs, rev(xs)),
                                          y = c(intervals.out[1,], rev(intervals.out[2,])),
                                          col = '#00BB00', alpha = 0.1, border = 0)
                            panel.points(x, y,...) })

    c('Manual outlier removal' = plot.orig, 'Our Method' = plot.corr)
}
plotreg()


##### Multivariate normal mixture
testgmvmix.out <- function () {
    library('lattice')
    library('latticeExtra')
    library('mixtools')

    testx <- t(rbind(rmvn( 200, mu= c(5, 7), sigma = diag(70,2)),
                     rmvn( 100, mu= rep(-20, 2), sigma = diag(50,2)),
                     rmvt( 2, mu=c(30,-70), df = 3, sigma = diag(1,2) ),
                     rmvt( 2, mu=c(0,80), df = 3, sigma = diag(1,2) )))

    xyplot(testx[2,] ~ testx[1,])

    res <- gmvmix.out(testx, k=2, nsamp = 3000)
    del <- apply(res$del[,-(1:1000)], c(1), mean)
    cols <- apply(res$z[,-(1:1000)], 1, function(x) which.max(tabulate(x)))
    cols[del > 0.6] <- 6

    mu <- apply(res$mu[,,-(1:1000)], c(1,2), mean)
    sig <- apply(res$sig[,,,-(1:1000)], c(1,2,3), mean)
    plot.out <- xyplot(testx[2,] ~ testx[1,], col = cols,
                       xlab = '', ylab = '',
                       main = 'Normal mixture with outlier detection',
                       key=list(space="top",
                                points = list(col=c(1,2,6), pch = 21),
                                text = list(c('Cluster 1', 'Cluster 2', 'P(Outlie) > 0.6'))),
                       panel = function (x, y, ...) {
                           elip <- mixtools::ellipse(mu[,1], sig[,,1],
                                                     alpha = 0.05, npoints = 400, draw = F )
                           panel.polygon(elip[,1], elip[,2], col = 1, alpha = 0.1)
                           elip <- mixtools::ellipse(mu[,2], sig[,,2],
                                                     alpha = 0.05, npoints = 400, draw = F )
                           panel.polygon(elip[,1], elip[,2], col = 2, alpha = 0.1)
                           panel.xyplot(x, y, ...)
                       })


    res.naive <- gmvmix.out(testx, k=2, nsamp = 3000, detect.outlier = F)
    del.naive <- apply(res.naive$del[,-(1:1000)], c(1), mean)
    cols.naive <- apply(res.naive$z[,-(1:1000)], 1, function(x) which.max(tabulate(x)))
    mu.naive <- apply(res.naive$mu[,,-(1:1000)], c(1,2), mean)
    sig.naive <- apply(res.naive$sig[,,,-(1:1000)], c(1,2,3), mean)

    plot.naive <- xyplot(testx[2,] ~ testx[1,], col = cols.naive,
                         xlab = '', ylab = '',
                         main = 'Naive normal mixture',
                         key=list(space="top",
                                  points = list(col=c(1,2), pch = 21),
                                  text = list(c('Cluster 1', 'Cluster 2'))),
                         panel = function (x, y, ...) {
                           elip <- mixtools::ellipse(mu.naive[,1], sig.naive[,,1],
                                                     alpha = 0.05, npoints = 400, draw = F )
                           panel.polygon(elip[,1], elip[,2], col = 1, alpha = 0.1)
                           elip <- mixtools::ellipse(mu.naive[,2], sig.naive[,,2],
                                                     alpha = 0.05, npoints = 400, draw = F )
                           panel.polygon(elip[,1], elip[,2], col = 2, alpha = 0.1)
                           panel.xyplot(x,y, ...)
                         })

    print(plot.naive, split = c(1,1,2,1), more = T)
    print(plot.out, split = c(2,1,2,1))
}

testgmvmix.out()

