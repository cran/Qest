# REGOLA GENERALE: se c'? una Qest.family, e se X ? singolare,
# l'utente deve COMUNQUE passare un valore di "start" per ogni colonna di X.


# QUESTA ? la funzione che dovrebbe sostituire la "glm.family" in ingresso.
# (vedi come modificarla per darle un formato "simile" alle altre).
# Si potrebbe passare direttamente X e y, invece di "data" e "formula".
# Nel mio codice, "data" dovrebbe essere gi? il model.frame (altrimenti model.response non funziona),
# senza valori mancanti ma con possibili singolarit?.

# "Uniform" family. You can fit U(min = a0 + a1*x, max = b0 + b1*x).
# Use min = FALSE if you want to fit a U(0, b0 + b1*x).
# If "start" is supplied, it should be a vector (a,b) of size 2*q, if min = TRUE,
# and a vector b of size q, if min = FALSE.

# The argument "formula" may be without intercept.

# Qunif.family <- function(formula, data, start, min = TRUE){
#   y <- model.response(data)
#   X <- model.matrix(formula, data)
#   int <- (attr(X, "assign")[1] == 0) # Is there an intercept?
#   q <- ncol(X)
#
#   # X singularities
#
#   z <- qr(X)
#   Xok <- z$pivot[1:z$rank]  # Ho visto che tu avevi fatto in un modo diverso
#   X <- X[,Xok]
#   signX <- apply(X,2, function(x){all(x <= 0) | all(x >= 0)})
#   if(any(!signX) && !int){
#     stop("When some predictors include both negative and positive values,
#          and there is no intercept, the model is automatically ill-defined due to quantile crossing")}
#
#   # Starting points
#
# # E' cos? che farebbe glm? "start" della stessa dimensione di X PRIMA dell'eliminazione delle singolarit?;
# # Inoltre, possibilit? di start = NULL invece di missing.
#   if(!missing(start) && !is.null(start)){
#     if(any(is.na(start))){stop("NAs are not allowed in 'start'")}
#     if(length(start) != q*(1 + min)){stop("Wrong size of 'start'")}
#     if(min){start <- start[c(Xok, q + Xok)]}
#     else{start <- start[Xok]}
#
#     q <- ncol(X)
#     a <- (if(min) start[1:q] else NULL)
#     b <- start[(q*min + 1):(q*min + q)]
#     pred.a <- (if(min) X%*%a else 0)
#     pred.b <- X%*%b
#     if(any(pred.a > pred.b)){stop("At the provided starting points, quantile crossing occurs")}
#   }
#   else{
#
#     m <- lm.fit(X,y)
#     a <- (if(min) m$coef else NULL)
#     b <- m$coef
#     s <- sd(m$residuals)
#     start.ok <- FALSE
#
#     count <- 0
#     while(!start.ok){
#       if(int){
#        b[1] <- b[1] + 2*s
#        if(min){a[1] <- a[1] - 2*s}
#       }
#       else{b <- b*1.5}
#       pred.a <- (if(min) X%*%a else 0)
#       pred.b <- X%*%b
#       start.ok <- all(pred.a < pred.b)
#       if(count == 100){stop("Could not find starting points, due to quantile crossing.")}
#     }
#
#   }
#
#   #####
#
#   attr(X, "min") <- min
#   list(start = c(a, b), X = X) # VEDI TU se vuoi mettere la X come attributo
# }

# This function does not actually compute Q, which is never needed except in the loss.
# For the same reason, there is no "tau" argument
Qunif <- function(min = TRUE) {
  Q <- function(theta, X){
    q <- ncol(X)
    min <- attr(X, "min")

    a <- (if(min) theta[1:q] else NULL)
    b <- theta[(q*min + 1):(q*min + q)]
    pred.a <- (if(min) c(X%*%a) else 0)
    pred.b <- c(X%*%b)

    delta <- pred.b - pred.a
    list(pred.a = pred.a, delta = delta, crossing = any(delta < 0))
  }
  scale.unif <- function(X, y, z) {
    Stats <- list(X = X, y = y, z = z, mX = 0, sX = 1, sy = 1)
    list(Xc = X, yc = y, zc = z, Stats = Stats)
  }
  descale.unif <- function(theta, Stats) {
    theta
  }
  bfun.unif <- function(wtau = NULL, wtauoptions = list()){

    if(is.null(wtau)){
      WL <- function(tau){tau*(1 - .5*tau)}
      WR <- function(tau){.5*tau^2}
      WLbar <- 1/2 - 1/6
      WRbar <- 1/6
      return(list(WL = WL, WR = WR, WLbar = WLbar, WRbar = WRbar))
    }

    r <- 4999
    p <- (0:r)/r
    dp <- p[-1] - p[-r]

    w <- callwtau(wtau, p, wtauoptions)
    WL <- num.fun(dp,w*(1 - p))
    WR <- num.fun(dp,w*p)
    WLbar <- num.fun(dp,WL)[r]
    WRbar <- num.fun(dp,WR)[r]

    list(
      WL = splinefun(p,WL, method = "hyman"),
      WR = splinefun(p,WR, method = "hyman"),
      WLbar = WLbar, WRbar = WRbar
    )
  }
  fix.Fy.unif <- function(fit, theta, tau, y, X, Q){
    Q1 <- Q(theta, X)
    Q1 <- Q1$pred.a + tau$TAU*Q1$delta
    list(Fy = fit$Fy, delta = y - Q1)
  }
  initialize <- function(x, z, y, d, w, Q, start, tau, data, ok, Stats){
    min <- attr(x, "min")
    if(!missing(start)) {
      q <- length(ok)
      if(any(is.na(start))){stop("NAs are not allowed in 'start'")}
      if(length(start) != q*(1 + min)){stop("Wrong size of 'start'")}
      start <- if(min){start[rep(ok, 1 + min)]}else{start <- start[ok]}

      q <- ncol(x)
      a <- (if(min) start[1:q] else NULL)
      b <- start[(q*min + 1):(q*min + q)]
      pred.a <- (if(min) c(x %*% a) else 0)
      pred.b <- c(x %*% b)
      if(any(pred.a > pred.b)){stop("At the provided starting points, quantile crossing occurs")}
    }
    else{
      nms <- colnames(x)
      int <- "(Intercept)" %in% nms
      # m <- lm.fit(x, y)

      a <- NULL
      if(min){
        p1 <- 0.3; p2 <- 0.7
        m1 <- rq.fit.br2(x, y, tau = p1)$coef
        m2 <- rq.fit.br2(x, y, tau = p2)$coef
        # m1 <- rq.fit(x, y, tau = p1, method = "fn")$coef
        # m2 <- rq.fit(x, y, tau = p2, method = "fn")$coef
        a <- m1 + p1/(p1 - p2)*(m2 - m1)
        b <- m1 + (1 - p1)/(p1 - p2)*(m1 - m2)
      }
      else{b <- rq.fit.br2(x, y, tau = 0.5)$coef*2}
      # print(a); print(b)
      # else{b <- rq.fit(x, y, tau = 0.5, method = "fn")$coef*2}

      # a <- (if(min) m$coef else NULL)
      if(!is.null(a)) names(a) <- paste0("min: ", names(a))
      # b <- m$coef
      names(b) <- paste0("max: ", names(b))
      # s <- sd(m$residuals)
      # start.ok <- FALSE
      # count <- 0
      # while(!start.ok){
      #   # if(int){
      #   #   b[1] <- b[1] + 2*s
      #   #   if(min){a[1] <- a[1] - 2*s}
      #   # }
      #   # else{b <- b*1.5}
      #   b <- b*1.5
      #   pred.a <- (if(min) c(x %*% a) else 0)
      #   pred.b <- c(x %*% b)
      #   start.ok <- all(pred.a < pred.b)
      #   if(count == 100){stop("Could not find starting points, due to quantile crossing.")}
      # }
    }
    c(a, b)
  }
  structure(list(family = list(family = "uniform"), Q = Q, scale = scale.unif, descale = descale.unif,
                 bfun = bfun.unif, min = min, fix.Fy = fix.Fy.unif, initialize = initialize),
            class = "Qfamily")
}

##########################################################

# QUESTA VA MESSA COME ELEMENTO AGGIUNTIVO di "tau", da passare a newton-raphson.
# Credo che sia da chiamare dentro a build.tau.

# If you call Qunif.bfun() with no arguments, it is assumed that wtau = 1.
# Qunif.bfun <- function(wtau, wtauoptions){
#
#   if(missing(wtau)){
#     WL <- function(tau){tau*(1 - tau)}
#     WR <- function(tau){tau^2}
#     WLbar <- 1/2 - 1/3
#     WRbar <- 1/3
#     return(list(WL = WL, WR = WR, WLbar = WLbar, WRbar = WRbar))
#   }
#
#   r <- 4999
#   p <- (0:r)/r
#   dp <- p[-1] - p[-r]
#
#   w <- callwtau(wtau, p, wtauoptions)
#   WL <- num.fun(dp,w*(1 - p))
#   WR <- num.fun(dp,w*p)
#   WLbar <- num.fun(dp,WL)[r]
#   WRbar <- num.fun(dp,WR)[r]
#
#
#   list(
#     WL = splinefun(p,WL, method = "hyman"),
#     WR = splinefun(p,WR, method = "hyman"),
#     WLbar = WLbar, WRbar = WRbar
#   )
# }

QestUnif.ee.u <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  X <- attr(Q, "X")
  min <- attr(X, "min")
  bfun <- tau$bfun
  n <- length(y)

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, X) # this is not Q(tau | x), actually
    # if(Q1$crossing){return(list(g = Inf))}
    Fy <- (y - Q1$pred.a)/Q1$delta
    Fy <- pmin(1, pmax0(Fy))

    # All quantities below to be multiplied by X
    AyL <- (if(min) bfun$WL(Fy) else NULL) # A(F(y)), left
    AyR <- bfun$WR(Fy) # A(F(y)), right
    AA1L <- (if(min) bfun$WLbar else NULL) # AA(1), left
    AA1R <- bfun$WRbar # AA(1), right

    Ay <- cbind(AyL, AyR) # a matrix n*(min + 1)
    AA1 <- c(AA1L, AA1R) # a vector of (min + 1) elements
    AA1 <- t(matrix(AA1, min + 1, n))

    g.i <- (AA1 - Ay)
    g.i <- (if(min) cbind(g.i[,1]*X, g.i[,2]*X) else c(g.i)*X)
    g <- c(w %*% g.i)
    fy <- NULL
  }
  else{g <- EE$g; g.i <- EE$g.i; Q1 <- EE$Q1; Fy <- EE$Fy}

  # Hessian

  if(J){

    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy)) # w(F(y))
    fy <- 1/Q1$delta # f(y)
    Qtheta_FyL <- (if(min) 1 - Fy else NULL)
    Qtheta_FyR <- Fy
    Qtheta_Fy <- (if(min) cbind(Qtheta_FyL*X,  Qtheta_FyR*X) else Qtheta_FyR*X)

    # Assemble

    H1 <- -Qtheta_Fy*(w*fy)
    H2 <- -wy*Qtheta_Fy
    H <- crossprod(H1, H2)

    # Finish
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the hessian
  }

  list(g = g, g.i = g.i, J = J, Q1 = Q1, Fy = Fy, fy = fy)
}

# Quando esci da newton-raphson, devi fare:
# fit$Fy <- fix.Fy.unif(fit)
# dove:

# fix.Fy.unif <- function(fit, theta, tau, y, X){
#   Q1 <- Qunif(theta, X)
#   Q1 <- Q1$pred.a + tau$TAU*Q1$delta
#   list(Fy = fit$Fy, delta = y - Q1)
# }











