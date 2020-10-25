##########################################
### Longitudinal Rate Regression Model ###
##########################################
### Using an independent covariance structure

require(compiler)

logmsg <- function(...) {
  message(Sys.time(), ': ', ...)
}

rdlm.ind <-
  function(DesignT, Outcome, Exposure, Base, Group, init.vals, max.it) {
    ## Variance and Mean Structure and Derivatives
    mean.str <- function(alphas, thetas, betas, DesignT, DesignE, DesignB) {
        mu <- DesignB %*% alphas + DesignT %*% betas %*% DesignE %*% thetas
        mu
    }

    dmean.dalphas <- function(DesignB, alpha.num) {
      mu.der <- DesignB[, alpha.num]
      mu.der
    }

    dmean.dthetas <- function(betas, DesignT, theta.val) {
      mu.der <- DesignT %*% betas * theta.val
      mu.der
    }

    dmean.dbetas <- function(thetas, DesignT, DesignE, beta.num) {
      mu.der <- DesignT[, beta.num] * DesignE %*% thetas
      mu.der
    }

    ## SCORE EQUATIONS ##
    score.eq <- function(dm.dp, mean.diff, variance) {
      sc <- (t(dm.dp) %*% mean.diff / exp(variance))
      sc
    }

    score.eq.sig <- function(variance, mean.diff, ni) {
      sc <- (-.5 * (
          ni / exp(variance) - t(mean.diff) %*% mean.diff / exp(2 * variance)
      ))
      sc
    }

    ## HESSIAN EQUATIONS  ##
    hess.eq <- function(dm.dp, variance) {
      hess <- (t(dm.dp) %*% dm.dp / exp(variance))
      hess
    }

    hess.eq.sig.sig <- function(variance, ni) {
      hess <- (.5 * ni / exp(2 * variance))
      hess
    }

    scorehess.mat.i <- function(alphas, thetas, betas, variance, ni, DesignT,
                                DesignE, DesignB, Outcome, p, r1, q1, p1,
                                vect.0) {
      logmsg("Starting setup")
      mean.diff <-
        (Outcome - mean.str(alphas, thetas, betas, DesignT, DesignE, DesignB))
      dm.da.mat <- matrix(sapply(c(1:r1),
                                 function(x) {
                                   dmean.dalphas(DesignB, x)
                                 }),
                          ncol = r1)
      dm.dt.mat <- matrix(sapply(c(1:q1),
                                 function(x) {
                                   dmean.dthetas(betas, DesignT,
                                                 DesignE[x + 1])
                                 }),
                          ncol = q1)
      dm.db.mat <- matrix(sapply(c(1:p1),
                                 function(x) {
                                   dmean.dbetas(thetas, DesignT,
                                                DesignE, x)
                                 }),
                          ncol = p1)
      dm.datb.mat <- cbind(dm.da.mat, dm.dt.mat, dm.db.mat)

      score.atb <- score.eq(dm.datb.mat, mean.diff, variance)
      vect.i <- c(score.atb, score.eq.sig(variance, mean.diff, ni))

      hess.atb <- hess.eq(dm.datb.mat, variance)
      mat <-
        rbind(cbind(hess.atb, vect.0), c(vect.0, hess.eq.sig.sig(variance, ni)))
      scorehess <- cbind(vect.i, mat)
      scorehess
    }
    require(compiler)
    scorehess.mat.i.new <- cmpfun(scorehess.mat.i)

    scorehess.mat <-
      function(alphas, thetas, betas, variance, M.m, ni.vect, DesignT.list,
               DesignE.list, DesignB.list, Outcome.list, p, r1, q1, p1, vect.0) {
        scorehesses <- lapply(1:M.m,
                              function(x) {
                                scorehess.mat.i.new(alphas, thetas, betas,
                                                    variance, ni.vect[x],
                                                    DesignT.list[[x]],
                                                    DesignE.list[[x]],
                                                    DesignB.list[[x]],
                                                    Outcome.list[[x]], p, r1,
                                                    q1, p1, vect.0)
                              })
        scorehess <- Reduce('+', scorehesses)
        return(list(score.vect = scorehess[, 1], hess.mat = scorehess[, -1]))
      }

    ## Likelihood Equation
    likelihood.1 <- function(alphas, thetas, betas, variance, ni, DesignT,
                             DesignE, DesignB, Outcome) {
      mean.diff <-
        (Outcome - mean.str(alphas, thetas, betas, DesignT, DesignE, DesignB))

      like <-
        (-.5 * (ni * variance + t(mean.diff) %*% mean.diff / exp(variance)))
      like
    }

    likelihood.n <- function(alphas, thetas, betas, variance, ni.vect,
                             DesignT.list, DesignE.list, DesignB.list,
                             Outcome.list) {
      likes <- lapply(1:M.m,
                      function(x) {
                        likelihood.1(alphas, thetas, betas, variance,
                                     ni.vect[x], DesignT.list[[x]],
                                     DesignE.list[[x]], DesignB.list[[x]],
                                     Outcome.list[[x]])
                      })
      like <- Reduce('+', likes)
      like
    }

    ## Setting up initial values
    p <- length(init.vals)
    p1 <- dim(as.matrix(DesignT))[2]
    q1 <- dim(as.matrix(Exposure))[2]
    r1 <- dim(as.matrix(Base))[2] + 1
    phi <- log(init.vals[p])
    parameters <- c(init.vals[-p], phi)
    ni.vect <- table(Group)
    id.vect <- sort(unique(Group))
    DesignT.list <-
      lapply(id.vect, function(x) {
        as.matrix(DesignT[Group == x, ])
      })
    DesignE.list <- lapply(id.vect,
                           function(x) {
                             cbind(1, as.matrix(Exposure)[Group == x, ])[1, ]
                           })
    DesignB.list <- lapply(id.vect,
                           function(x) {
                             cbind(1, as.matrix(Base)[Group == x, ])
                           })
    Outcome.list <-
      lapply(id.vect, function(x) {
        as.matrix(Outcome[Group == x])
      })
    M.m <- length(ni.vect)
    vect.0 <- rep(0, p - 1)

    ## Run NR iterations
    change <- rep(-999, p)
    count <- 0
    fudgeIND <- diag(rep(.0001, p))
    logmsg("Starting iterations")
    while (sum(abs(change)) > 0.001 & count < max.it) {
      iter <- scorehess.mat(parameters[1:r1],
                            c(1, parameters[(r1 + 1):(r1 + q1)]),
                            parameters[(r1 + q1 + 1):(r1 + q1 + p1)],
                            parameters[p], M.m, ni.vect, DesignT.list,
                            DesignE.list, DesignB.list, Outcome.list, p, r1, q1,
                            p1, vect.0)

      eig <- tryCatch(
        min(abs(eigen(iter$hess.mat)$value)),
        error = function(e) {
          return(9999)
        }
      )
      ifelse(
        eig < .00001,
        change <-
          (solve(iter$hess.mat + fudgeIND) %*% iter$score.vect),
        change <- (solve(iter$hess.mat) %*% iter$score.vect)
      )
      adjust <- ifelse(count < 4, .5, 1)
      parameters <- parameters + adjust * change
      count <- count + 1
      logmsg("Maximum step change after iteration", count)
      logmsg(max(abs(change)))
    }
    logmsg(count)
    logmsg("Finish Mean Iterations")
    return(parameters)
  }

#####################################################################
### Using random intercept and slope effect covariance structure

rdlm <- function(DesignT, Outcome, Exposure, Base, Group, init.vals) {
  ## Variance and Mean Structure and Derivatives
  V.str <- function(betas, variance, DesignT, one.ni, ni) {
    Var <- (
      exp(variance[1]) * one.ni %*% t(one.ni) + exp(variance[2]) * DesignT %*%
        betas %*% t(DesignT %*% betas) + (2 / pi) * atan(variance[3]) *
        exp((variance[1] + variance[2]) /
              2) * (
          one.ni %*% t(DesignT %*% betas) + DesignT %*% betas %*% t(one.ni)
              ) + exp(variance[4]) * diag(1, nrow = ni)
    )
    Var
  }

  dV.dbetas <- function(betas, variance, DesignT, one.ni, beta.num) {
    Var.der <-
      (
        exp(variance[2]) * (
          DesignT[, beta.num] %*% t(DesignT %*% betas) +
            DesignT %*% betas %*% t(DesignT[, beta.num])
        ) +
          (2 / pi) * atan(variance[3]) * exp((variance[1] + variance[2]) /
                                               2) *
          (one.ni %*% t(DesignT[, beta.num]) +
             DesignT[, beta.num] %*% t(one.ni))
      )
    Var.der
  }

  dV.dsig0 <- function(betas, variance, DesignT, one.ni) {
    Var.der <- (
      exp(variance[1]) * one.ni %*% t(one.ni) + (1 / pi) *
        atan(variance[3]) * exp((variance[1] + variance[2]) /
                                  2) *
        (
          one.ni %*% t(DesignT %*% betas) + DesignT %*% betas %*% t(one.ni)
        )
    )
    Var.der
  }

  dV.dsig1 <- function(betas, variance, DesignT, one.ni) {
    Var.der <- (
      exp(variance[2]) * DesignT %*% betas %*% t(DesignT %*% betas) +
        (1 / pi) * atan(variance[3]) * exp((variance[1] + variance[2]) /
                                             2) *
        (
          one.ni %*% t(DesignT %*% betas) + DesignT %*% betas %*% t(one.ni)
        )
    )
    Var.der
  }

  dV.drho <- function(betas, variance, DesignT, one.ni) {
    Var.der <-
      ((2 / (pi * (1 + variance[3] ^ 2))) * exp((variance[1] + variance[2]) /
                                                  2) *
         (
           one.ni %*% t(DesignT %*% betas) + DesignT %*% betas %*% t(one.ni)
         ))
    Var.der
  }

  dV.dsig <- function(variance, ni) {
    Var.der <- exp(variance[4]) * diag(1, nrow = ni)
    Var.der
  }

  mean.str <-
    function(alphas, thetas, betas, DesignT, DesignE, DesignB) {
      mu <- DesignB %*% alphas + DesignT %*% betas %*% DesignE %*% thetas
      mu
    }

  dmean.dalphas <- function(DesignB, alpha.num) {
    mu.der <- DesignB[, alpha.num]
    mu.der
  }

  dmean.dthetas <- function(betas, DesignT, theta.val) {
    mu.der <- DesignT %*% betas * theta.val
    mu.der
  }

  dmean.dbetas <- function(thetas, DesignT, DesignE, beta.num) {
    mu.der <- DesignT[, beta.num] * DesignE %*% thetas
    mu.der
  }

  ## SCORE EQUATIONS ##
  score.eq.mean <- function(dm.dp, mean.diff, inv.var) {
    sc <- (t(dm.dp) %*% inv.var %*% mean.diff)
    sc
  }

  score.eq.var <- function(dV.dv, mean.diff, inv.var) {
    sc <- (-.5 * (
      sum(diag(inv.var %*% dV.dv)) - t(mean.diff) %*% inv.var %*%
        dV.dv %*% inv.var %*% mean.diff
    ))
    sc
  }

  ## HESSIAN EQUATIONS  ##
  hess.eq.mean <- function(dm.dp, inv.var) {
    hess <- (t(dm.dp) %*% inv.var %*% dm.dp)
    hess
  }

  hess.eq.var <- function(dV.dv1, dV.dv2, inv.var) {
    hess <- (.5 * sum(diag(
      inv.var %*% dV.dv1 %*% inv.var %*% dV.dv2
    )))
    hess
  }

  # Score and Hessian Calculator function
  scorehess.mat.i <-
    function(alphas, thetas, betas, variance, ni, one.ni, DesignT, DesignE,
             DesignB, Outcome, p, r1, q1, p1, mat.0) {
      mean.diff <-
        (Outcome - mean.str(alphas, thetas, betas, DesignT, DesignE, DesignB))
      inv.var <- solve(V.str(betas, variance, DesignT, one.ni, ni))

      dm.da.mat <- matrix(sapply(c(1:r1),
                                 function(x) {
                                   dmean.dalphas(DesignB, x)
                                 }), ncol = r1)
      dm.dt.mat <- matrix(sapply(c(1:q1),
                                 function(x) {
                                   dmean.dthetas(betas, DesignT, DesignE[x + 1])
                                 }), ncol = q1)
      dm.db.mat <- matrix(sapply(c(1:p1),
                                 function(x) {
                                   dmean.dbetas(thetas, DesignT, DesignE, x)
                                 }),
                          ncol = p1)
      dm.datb.mat <- cbind(dm.da.mat, dm.dt.mat, dm.db.mat)

      dV.db.list <- lapply(c(1:p1),
                           function(x) {
                             dV.dbetas(betas, variance, DesignT, one.ni, x)
                           })
      dV.ds0 <- dV.dsig0(betas, variance, DesignT, one.ni)
      dV.ds1 <- dV.dsig1(betas, variance, DesignT, one.ni)
      dV.dr <- dV.drho(betas, variance, DesignT, one.ni)
      dV.ds <- dV.dsig(variance, ni)

      score.mean <- score.eq.mean(dm.datb.mat, mean.diff, inv.var)
      score.betas <- sapply(c(1:p1),
                            function(x) {
                              score.eq.var(dV.db.list[[x]], mean.diff, inv.var)
                            })
      score.mean[(r1 + q1 + 1):(r1 + q1 + p1)] <-
        score.mean[(r1 + q1 + 1):(r1 + q1 + p1)] + score.betas
      vect.i <- c(
        score.mean,
        score.eq.var(dV.ds0, mean.diff, inv.var),
        score.eq.var(dV.ds1, mean.diff, inv.var),
        score.eq.var(dV.dr, mean.diff, inv.var),
        score.eq.var(dV.ds, mean.diff, inv.var)
      )

      mat11 <- hess.eq.mean(dm.datb.mat, inv.var)
      hess.betas.betas <- sapply(c(1:p1),
                                 function(x) {
                                   sapply(c(1:p1), function(y) {
                                     hess.eq.var(dV.db.list[[x]], dV.db.list[[y]], inv.var)
                                   })
                                 })
      mat11[(r1 + q1 + 1):(r1 + q1 + p1), (r1 + q1 + 1):(r1 + q1 + p1)] <-
        mat11[(r1 + q1 + 1):(r1 + q1 + p1), (r1 + q1 + 1):(r1 + q1 + p1)] +
        hess.betas.betas

      hess.betas.sig0 <- sapply(c(1:p1),
                                function(x) {
                                  hess.eq.var(dV.db.list[[x]], dV.ds0, inv.var)
                                })
      hess.betas.sig1 <-
        sapply(c(1:p1), function(x) {
          hess.eq.var(dV.db.list[[x]], dV.ds1, inv.var)
        })
      hess.betas.rho <-
        sapply(c(1:p1), function(x) {
          hess.eq.var(dV.db.list[[x]], dV.dr, inv.var)
        })
      hess.betas.sig <-
        sapply(c(1:p1), function(x) {
          hess.eq.var(dV.db.list[[x]], dV.ds, inv.var)
        })
      mat12 <-
        rbind(hess.betas.sig0,
              hess.betas.sig1,
              hess.betas.rho,
              hess.betas.sig)

      hess.s0.s0 <- hess.eq.var(dV.ds0, dV.ds0, inv.var)
      hess.s0.s1 <- hess.eq.var(dV.ds0, dV.ds1, inv.var)
      hess.s0.r <- hess.eq.var(dV.ds0, dV.dr, inv.var)
      hess.s0.s <- hess.eq.var(dV.ds0, dV.ds, inv.var)
      hess.s1.s1 <- hess.eq.var(dV.ds1, dV.ds1, inv.var)
      hess.s1.r <- hess.eq.var(dV.ds1, dV.dr, inv.var)
      hess.s1.s <- hess.eq.var(dV.ds1, dV.ds, inv.var)
      hess.r.r <-  hess.eq.var(dV.dr, dV.dr, inv.var)
      hess.r.s <- hess.eq.var(dV.dr, dV.ds, inv.var)
      hess.s.s <- hess.eq.var(dV.ds, dV.ds, inv.var)
      mat22 <-
        matrix(c(hess.s0.s0, hess.s0.s1, hess.s0.r, hess.s0.s, hess.s0.s1,
                 hess.s1.s1, hess.s1.r, hess.s1.s, hess.s0.r, hess.s1.r, hess.r.r,
                 hess.r.s, hess.s0.s, hess.s1.s, hess.r.s, hess.s.s), 4, 4)

      mat <-
        rbind(cbind(mat11, rbind(t(mat.0), t(mat12))), cbind(mat.0, mat12, mat22))
      scorehess <- cbind(vect.i, mat)
      scorehess
    }
  scorehess.mat.i.new <- cmpfun(scorehess.mat.i)

  scorehess.mat <-
    function(alphas, thetas, betas, variance, M.m, ni.vect, one.list,
             DesignT.list, DesignE.list, DesignB.list, Outcome.list, p, r1, q1,
             p1, mat.0) {
      scorehesses <-
        lapply(1:M.m, function(x) {
          scorehess.mat.i.new(alphas, thetas, betas, variance, ni.vect[x],
                              one.list[[x]], DesignT.list[[x]],
                              DesignE.list[[x]], DesignB.list[[x]],
                              Outcome.list[[x]], p, r1, q1, p1, mat.0)
        })
      scorehess <- Reduce('+', scorehesses)
      return(list(score.vect = scorehess[, 1], hess.mat = scorehess[, -1]))
    }

  ## Likelihood Equation
  likelihood.1 <-
    function(alphas, thetas, betas, variance, ni, one.ni, DesignT, DesignE,
             DesignB, Outcome) {
      mean.diff <-
        (Outcome - mean.str(alphas, thetas, betas, DesignT, DesignE, DesignB))
      Var.str <- V.str(betas, variance, DesignT, one.ni, ni)

      like <-
        (-.5 * (log(det(Var.str)) + t(mean.diff) %*% solve(Var.str) %*% mean.diff))
      like
    }

  likelihood.n <-
    function(alphas, thetas, betas, variance, ni.vect, one.list, DesignT.list,
             DesignE.list, DesignB.list, Outcome.list) {
      likes <- lapply(1:M.m, function(x) {
        likelihood.1(alphas, thetas, betas, variance, ni.vect[x],
                     one.list[[x]], DesignT.list[[x]], DesignE.list[[x]],
                     DesignB.list[[x]], Outcome.list[[x]])
      })
      like <- Reduce('+', likes)

      like
    }

  ## Setting up initial values
  p1 <- dim(as.matrix(DesignT))[2]
  q1 <- dim(as.matrix(Exposure))[2]
  r1 <- dim(as.matrix(Base))[2] + 1
  p <- p1 + q1 + r1 + 4
  #values <- c(init.vals[1:(p-4)],sum(init.vals[c(p-3,p-2,p)]))
  #init.new <- rdlm.ind(DesignT,Outcome,Exposure,Base,Group,values,100)
  alpha <- log(init.vals[p - 3])
  delta <- log(init.vals[p - 2])
  psi <- tan(init.vals[p - 1] * pi / 2)
  phi <- log(init.vals[p])
  #parameters <- c(init.new[1:(p-4)],alpha,delta,psi,phi)
  parameters <- c(init.vals[1:(p - 4)], alpha, delta, psi, phi)
  ni.vect <- table(Group)
  id.vect <- sort(unique(Group))
  one.list <- lapply(ni.vect, function(x) { rep(1, x) })
  DesignT.list <-
    lapply(id.vect, function(x) { as.matrix(DesignT[Group == x, ]) })
  DesignE.list <-
    lapply(id.vect, function(x) {
      cbind(1, as.matrix(Exposure)[Group == x, ])[1, ]
    })
  DesignB.list <-
    lapply(id.vect, function(x) {
      cbind(1, as.matrix(Base)[Group == x, ])
    })
  Outcome.list <-
    lapply(id.vect, function(x) {
      as.matrix(Outcome[Group == x])
    })
  M.m <- length(ni.vect)
  mat.0 <- matrix(rep(0, (p - 4 - p1) * 4), nrow = 4)


  ## Run NR iterations
  change <- rep(-999, p)
  count <- 0
  fudgeME <- diag(rep(.0001, p))
  while (sum(abs(change)) > 0.000001 & count < 100) {
    iter <-
      scorehess.mat(parameters[1:r1], c(1, parameters[(r1 + 1):(r1 + q1)]),
                    parameters[(r1 + q1 + 1):(r1 + q1 + p1)],
                    parameters[(p - 3):p], M.m, ni.vect, one.list, DesignT.list,
                    DesignE.list, DesignB.list, Outcome.list, p, r1, q1, p1,
                    mat.0)

    eig <-
      tryCatch(
        min(abs(eigen(iter$hess.mat)$value)),
        error = function(e) {
          return(9999)
        }
      )
    ifelse(
      eig < .00001,
      change <- (solve(iter$hess.mat + fudgeME) %*% iter$score.vect),
      change <- (solve(iter$hess.mat) %*% iter$score.vect)
    )
    adjust <- ifelse(count < 4, .5, 1)
    parameters <- parameters + adjust * change
    count <- count + 1
    logmsg(paste("Total step change after iteration", count))
    logmsg(sum(abs(change)))
    logmsg("Parameter Values:")
    logmsg(parameters)
  }
  logmsg("Algorithm Convergence")
  hess.der <-
    scorehess.mat(parameters[1:r1], c(1, parameters[(r1 + 1):(r1 + q1)]),
                  parameters[(r1 + q1 + 1):(r1 + q1 + p1)],
                  parameters[(p - 3):p], M.m, ni.vect, one.list, DesignT.list,
                  DesignE.list, DesignB.list, Outcome.list, p, r1, q1, p1,
                  mat.0)$hess.mat
  logmsg("Final Step 1")
  sig0 <- exp(parameters[p - 3])
  sig1 <- exp(parameters[p - 2])
  rho <- (atan(parameters[p - 1]) * 2 / pi)
  sig <- exp(parameters[p])
  parameters[p - 3] <- sig0
  parameters[p - 2] <- sig1
  parameters[p - 1] <- rho
  parameters[p] <- sig
  logmsg("Final Step 2")

  hess <-
    (diag(c(rep(1, p - 4), parameters[p - 3] ^ (-1), parameters[p - 2] ^ (-1),
    (cos(parameters[p - 1] * pi / 2) ^ (-2) * pi / 2), parameters[p] ^ (-1))) %*% hess.der %*%
      diag(c(rep(1, p - 4), parameters[p - 3] ^ (-1), parameters[p - 2] ^ (-1),
      (cos(parameters[p - 1] * pi / 2) ^ (-2) * pi / 2), parameters[p] ^ (-1))))
  logmsg("Final Step 3")
  results <-
    (as.vector(parameters)[1:(p - 4)] + qnorm(0.975) *
       sqrt(diag(solve(hess)[1:(p - 4), 1:(p - 4)])) %o% c(0, -1, 1))
  logmsg("Final Step 4")

  return(list(
    coefficients = results,
    variance = parameters[(p - 3):p],
    Hessian = hess
  ))
}


################################################################################
########## Univariate Rate Regression with LDL decomposition
################################################################################

urrlm <-
  function(DesignT, Exposure, Base, Outcome, Group, init.mean, init.D, init.var) {
    ## Variance and Mean Structure and Derivatives
    V.str <- function(Dvar, var.e, Z.mat, ni) {
      Var <- (Z.mat %*% Dvar %*% t(Z.mat) + diag(exp(var.e), nrow = ni))
      Var <- (Var + t(Var)) / 2
      Var
    }

    dV.dbetas <- function(Dvar, Z.mat, Z.der) {
      Var.der <- (Z.der %*% Dvar %*% t(Z.mat) + Z.mat %*% Dvar %*% t(Z.der))
      Var.der
    }

    dV.dDmat <- function(Dvar.der, Z.mat) {
      Var.der <- (Z.mat %*% Dvar.der %*% t(Z.mat))
      Var.der
    }

    dV.dsig <- function(var.e, ni) {
      Var.der <- exp(var.e) * diag(1, nrow = ni)
      Var.der
    }

    mean.str <- function(alphas, thetas, betas, DesignT, DesignE, DesignB) {
      mu <- DesignB %*% alphas + DesignT %*% betas %*% DesignE %*% thetas
      mu
    }

    dmean.dalphas <- function(DesignB, alpha.num) {
      mu.der <- DesignB[, alpha.num]
      mu.der
    }

    dmean.dthetas <- function(betas, DesignT, theta.val) {
      mu.der <- DesignT %*% betas * theta.val
      mu.der
    }

    dmean.dbetas <- function(thetas, DesignT, DesignE, beta.num) {
      mu.der <- DesignT[, beta.num] * DesignE %*% thetas
      mu.der
    }

    ## SCORE EQUATIONS ##
    score.eq.mean <- function(dm.dp, mean.diff, inv.var) {
      sc <- (t(dm.dp) %*% inv.var %*% mean.diff)
      sc
    }

    score.eq.var <- function(dV.dv, mean.diff, inv.var) {
      sc <- (-.5 * (sum(diag(inv.var %*% dV.dv)) - t(mean.diff) %*%
                      inv.var %*% dV.dv %*% inv.var %*% mean.diff))
      sc
    }

    ## HESSIAN EQUATIONS  ##
    hess.eq.mean <- function(dm.dp, inv.var) {
      hess <- (t(dm.dp) %*% inv.var %*% dm.dp)
      hess
    }

    hess.eq.var <- function(dV.dv1, dV.dv2, inv.var) {
      hess <- (.5 * sum(diag(inv.var %*% dV.dv1 %*% inv.var %*% dV.dv2)))
      hess
    }


    # Score and Hessian Calculator function
    scorehess.mat.i <-
      function(alphas, thetas, betas, Dvar, Dvar.D11, Dvar.D22, Dvar.L, var.e,
               ni, one.ni, DesignT, DesignE, DesignB, Outcome, p, r1, q1, p1,
               mat.0) {
        Z.mat <- cbind(one.ni, DesignT %*% betas)
        mean.diff <-
          (Outcome - mean.str(alphas, thetas, betas, DesignT, DesignE, DesignB))
        inv.var <- solve(V.str(Dvar, var.e, Z.mat, ni))

        dm.da.mat <-
          matrix(sapply(c(1:r1), function(x) {
            dmean.dalphas(DesignB, x)
          }), ncol = r1)
        dm.dt.mat <-
          matrix(sapply(c(1:q1), function(x) {
            dmean.dthetas(betas, DesignT, DesignE[x + 1])
          }), ncol = q1)
        dm.db.mat <-
          matrix(sapply(c(1:p1), function(x) {
            dmean.dbetas(thetas, DesignT, DesignE, x)
          }), ncol = p1)
        dm.datb.mat <- cbind(dm.da.mat, dm.dt.mat, dm.db.mat)

        dV.db.list <-
          lapply(c(1:p1), function(x) {
            dV.dbetas(Dvar, Z.mat, cbind(one.ni - one.ni, DesignT[, x]))
          })
        dV.dD11 <- dV.dDmat(Dvar.D11, Z.mat)
        dV.dD22 <- dV.dDmat(Dvar.D22, Z.mat)
        dV.dL <- dV.dDmat(Dvar.L, Z.mat)
        dV.ds <- dV.dsig(var.e, ni)

        score.mean <- score.eq.mean(dm.datb.mat, mean.diff, inv.var)
        score.betas <-
          sapply(c(1:p1), function(x) {
            score.eq.var(dV.db.list[[x]], mean.diff, inv.var)
          })
        score.mean[(r1 + q1 + 1):(r1 + q1 + p1)] <-
          score.mean[(r1 + q1 + 1):(r1 + q1 + p1)] + score.betas
        vect.i <-
          c(
            score.mean,
            score.eq.var(dV.dD11, mean.diff, inv.var),
            score.eq.var(dV.dD22, mean.diff, inv.var),
            score.eq.var(dV.dL, mean.diff, inv.var),
            score.eq.var(dV.ds, mean.diff, inv.var)
          )

        mat11 <- hess.eq.mean(dm.datb.mat, inv.var)
        hess.betas.betas <-
          sapply(c(1:p1), function(x) {
            sapply(c(1:p1), function(y) {
              hess.eq.var(dV.db.list[[x]], dV.db.list[[y]], inv.var)
            })
          })
        mat11[(r1 + q1 + 1):(r1 + q1 + p1), (r1 + q1 + 1):(r1 + q1 + p1)] <-
          mat11[(r1 + q1 + 1):(r1 + q1 + p1), (r1 + q1 + 1):(r1 + q1 + p1)] + hess.betas.betas

        hess.betas.D11 <-
          sapply(c(1:p1), function(x) {
            hess.eq.var(dV.db.list[[x]], dV.dD11, inv.var)
          })
        hess.betas.D22 <-
          sapply(c(1:p1), function(x) {
            hess.eq.var(dV.db.list[[x]], dV.dD22, inv.var)
          })
        hess.betas.L <-
          sapply(c(1:p1), function(x) {
            hess.eq.var(dV.db.list[[x]], dV.dL, inv.var)
          })
        hess.betas.sig <-
          sapply(c(1:p1), function(x) {
            hess.eq.var(dV.db.list[[x]], dV.ds, inv.var)
          })
        mat12 <-
          rbind(hess.betas.D11,
                hess.betas.D22,
                hess.betas.L,
                hess.betas.sig)

        hess.d1.d1 <- hess.eq.var(dV.dD11, dV.dD11, inv.var)
        hess.d1.d2 <- hess.eq.var(dV.dD11, dV.dD22, inv.var)
        hess.d1.l <- hess.eq.var(dV.dD11, dV.dL, inv.var)
        hess.d1.s <- hess.eq.var(dV.dD11, dV.ds, inv.var)
        hess.d2.d2 <- hess.eq.var(dV.dD22, dV.dD22, inv.var)
        hess.d2.l <- hess.eq.var(dV.dD22, dV.dL, inv.var)
        hess.d2.s <- hess.eq.var(dV.dD22, dV.ds, inv.var)
        hess.l.l <-  hess.eq.var(dV.dL, dV.dL, inv.var)
        hess.l.s <- hess.eq.var(dV.dL, dV.ds, inv.var)
        hess.s.s <- hess.eq.var(dV.ds, dV.ds, inv.var)
        mat22 <-
          matrix(
            c(hess.d1.d1, hess.d1.d2, hess.d1.l, hess.d1.s, hess.d1.d2,
              hess.d2.d2, hess.d2.l, hess.d2.s, hess.d1.l, hess.d2.l, hess.l.l,
              hess.l.s, hess.d1.s, hess.d2.s, hess.l.s, hess.s.s),
            4,
            4
          )

        mat <-
          rbind(cbind(mat11, rbind(t(mat.0), t(mat12))), cbind(mat.0, mat12, mat22))
        scorehess <- cbind(vect.i, mat)
        scorehess
      }
    require(compiler)
    scorehess.mat.i.new <- cmpfun(scorehess.mat.i)

    scorehess.mat <-
      function(alphas, thetas, betas, var.D, var.L, var.e, M.m, ni.vect,
               one.list, DesignT.list, DesignE.list, DesignB.list, Outcome.list,
               p, r1, q1, p1, mat.0) {
        L.mat <- matrix(c(1, var.L, 0, 1), nrow = 2)
        D.mat <- diag(exp(var.D))
        Dvar <- L.mat %*% D.mat %*% t(L.mat)
        Dvar <- (Dvar + t(Dvar)) / 2
        Dvar.D11 <- L.mat %*% diag(c(exp(var.D[1]), 0)) %*% t(L.mat)
        Dvar.D11 <- (Dvar.D11 + t(Dvar.D11)) / 2
        Dvar.D22 <- L.mat %*% diag(c(0, exp(var.D[2]))) %*% t(L.mat)
        Dvar.D22 <- (Dvar.D22 + t(Dvar.D22)) / 2
        Dvar.L <- (matrix(c(0, 1, 0, 0), nrow = 2) %*% D.mat %*% t(L.mat) +
                     L.mat %*% D.mat %*% t(matrix(c(0, 1, 0, 0), nrow = 2)))
        Dvar.L <- (Dvar.L + t(Dvar.L)) / 2
        scorehesses <-
          lapply(1:M.m, function(x) {
            scorehess.mat.i.new(alphas, thetas, betas, Dvar, Dvar.D11, Dvar.D22,
                                Dvar.L, var.e, ni.vect[x], one.list[[x]],
                                DesignT.list[[x]], DesignE.list[[x]],
                                DesignB.list[[x]], Outcome.list[[x]], p, r1, q1,
                                p1, mat.0)
          })
        scorehess <- Reduce('+', scorehesses)
        return(list(score.vect = scorehess[, 1], hess.mat = scorehess[, -1]))
      }

    ## Likelihood Equation
    likelihood.1 <-
      function(alphas, thetas, betas, Dvar, var.e, ni, one.ni, DesignT, DesignE,
               DesignB, Outcome) {
        Z.mat <- cbind(one.ni, DesignT %*% betas)
        mean.diff <-
          (Outcome - mean.str(alphas, thetas, betas, DesignT, DesignE, DesignB))
        Var.str <- V.str(Dvar, var.e, Z.mat, ni)

        like <-
          (-.5 * (log(det(Var.str)) + t(mean.diff) %*% solve(Var.str) %*% mean.diff))
        like
      }

    likelihood.n <-
      function(alphas, thetas, betas, var.D, var.L, var.e, ni.vect, one.list,
               DesignT.list, DesignE.list, DesignB.list, Outcome.list) {
        L.mat <- matrix(c(1, var.L, 0, 1), nrow = 2)
        D.mat <- diag(exp(var.D))
        Dvar <- L.mat %*% D.mat %*% t(L.mat)
        Dvar <- (Dvar + t(Dvar)) / 2

        likes <-
          lapply(1:M.m, function(x) {
            likelihood.1(alphas, thetas, betas, Dvar, var.e, ni.vect[x],
                         one.list[[x]], DesignT.list[[x]], DesignE.list[[x]],
                         DesignB.list[[x]], Outcome.list[[x]])
          })
        like <- Reduce('+', likes)

        like
      }

    ## Setting up initial values
    p1 <- dim(as.matrix(DesignT))[2]
    q1 <- dim(as.matrix(Exposure))[2]
    r1 <- dim(as.matrix(Base))[2] + 1
    p <- p1 + q1 + r1 + 4

    D.mat.11 <- log(init.D[1, 1])
    L.mat.21 <- (exp(-D.mat.11) * init.D[2, 1])
    D.mat.22 <- log(init.D[2, 2] - L.mat.21 ^ 2 * exp(D.mat.11))

    parameters <- c(init.mean, D.mat.11, D.mat.22, L.mat.21, log(init.var))

    ni.vect <- table(Group)
    id.vect <- sort(unique(Group))
    one.list <- lapply(ni.vect, function(x) {
      rep(1, x)
    })
    DesignT.list <-
      lapply(id.vect, function(x) {
        as.matrix(DesignT[Group == x, ])
      })
    DesignE.list <-
      lapply(id.vect, function(x) {
        cbind(1, as.matrix(Exposure)[Group == x, ])[1, ]
      })
    DesignB.list <-
      lapply(id.vect, function(x) {
        cbind(1, as.matrix(Base)[Group == x, ])
      })
    Outcome.list <-
      lapply(id.vect, function(x) {
        as.matrix(Outcome[Group == x])
      })
    M.m <- length(ni.vect)
    mat.0 <- matrix(rep(0, (p - 4 - p1) * 4), nrow = 4)

    ## Run NR iterations
    change <- rep(-999, p)
    count <- 0
    fudgeME <- diag(rep(.0001, p))
    while (max(abs(change)) > 0.000001 & count < 100) {
      iter <-
        scorehess.mat(parameters[1:r1], c(1, parameters[(r1 + 1):(r1 + q1)]),
                      parameters[(r1 + q1 + 1):(r1 + q1 + p1)],
                      parameters[(p - 3):(p - 2)], parameters[(p - 1)],
                      parameters[p], M.m, ni.vect, one.list, DesignT.list,
                      DesignE.list, DesignB.list, Outcome.list, p, r1, q1, p1,
                      mat.0)

      eig <-
        tryCatch(
          min(abs(eigen(iter$hess.mat)$value)),
          error = function(e) {
            return(9999)
          }
        )
      ifelse(
        eig < .00001,
        change <- (solve(iter$hess.mat + fudgeME) %*% iter$score.vect),
        change <- (solve(iter$hess.mat) %*% iter$score.vect)
      )
      adjust <- ifelse(count < 4, .5, 1)
      parameters <- parameters + adjust * change
      count <- count + 1
    }
    L.mat <- matrix(c(1, parameters[(p - 1)], 0, 1), nrow = 2)
    D.mat <- diag(exp(parameters[(p - 3):(p - 2)]))
    return(list(
      coefficients = parameters[1:(p - 4)],
      Dvar = L.mat %*% D.mat %*% t(L.mat),
      variance = exp(parameters[p])
    ))
  }


run_analysis <- function(model_name, init.vals) {
  # Read person-level data, prepared in data frame wt_ht_data3
  load(paste0(model_name, '.RData'))

  # Calculating the run time for the LRR model
  fnc2 <- proc.time()

  LRRmodel.rslt <- rdlm(wt_ht_data3[,c("Basis1","Basis2","Basis3","Basis4")],
                        wt_ht_data3$wt,
                        wt_ht_data3[,c("abx_ind","sex.c", "has_asthma.c",
                                       "SterOral_count1.c", "SterOral_count2.c",
                                       "SterOral_count3.c", "SterOral_count4.c",
                                       "Infect_count1.c", "Infect_count2.c",
                                       "Infect_count3.c", "Infect_count4.c",
                                       "preterm.c", "ethnicity_black.c",
                                       "ethnicity_other.c", "ethnicity_white.c",
                                       "ethnicity_hispanic.c", "has_chronic.c")],
                        wt_ht_data3[,c("abx_ind","sex.n","has_asthma",
                                       "SterOral_count1","SterOral_count2",
                                       "SterOral_count3","SterOral_count4",
                                       "Infect_count1","Infect_count2",
                                       "Infect_count3","Infect_count4",
                                       "preterm","ethnicity_black",
                                       "ethnicity_other","ethnicity_white",
                                       "ethnicity_hispanic","has_chronic")],
                        wt_ht_data3$uniqueID, init.vals)

  # Calculating LRR run time
  logmsg('Model estimation required ', (proc.time() - fnc2) / 60 * 60, ' hours')

  ## Writing all estimated coefficients and confidence intervals to a .csv file
  write.csv(LRRmodel.rslt$coeff,file=paste0(model_name, '_coeff.csv'))

  ## Saving all model results
  save.image(paste0(model_name, '_model.RData'))

  invisible(LRRmodel.rslt)
}

## Example
##
## run_analysis('htwt_all',
##              init.vals = c(12.029980575,0.165520449,0.639503136, 0.104532019,
##                            0.026583171, 0.220233422, 0.174956190,-0.040834969,
##                            -0.163208194, -0.107124715, -0.197266061,
##                            -0.139317207, -0.578301315,0.560104259,0.619769606,
##                            0.467269052,0.752299977,-0.447815314, -0.030422789,
##                            -0.030716042, 0.091625298, 0.017754846,-0.012909890,
##                            0.048730112, -0.032194370, -0.042130199,
##                            -0.007286965, -0.016641516,-0.004910909,
##                            -0.018737812, 0.165001931,0.126840838,0.054675590,
##                            0.117270175,-0.016429072, 4.765527411, 6.720920495,
##                            10.895950039,8.251380747, 2.9291288, 0.1144845,
##                            0.4325124, 0.2688165))

