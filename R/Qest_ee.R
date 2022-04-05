# Qest estimating equations


# Estimating equation (= gradient of the loss) and Jacobian (= Hessian of the loss) for type = "u".
Qest.ee.u <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, tau$TAU, data) # Q(tau | x), a matrix n*ntau
    Qtheta <- derQtheta(theta, eps, Q, Q1, data, tau) # r matrices of size n*ntau
    Fy <- findp(y, tau, Q1) # F(y | x)

    ap <- lapply(Qtheta, function(x){x*tau$wtau_long}) # a_theta(p) = wtau(p)*Qtheta(p)
    Ap <- intA(ap, tau) # A_theta(p)
    Ay <- findAp(Fy, tau, Ap) # A_theta(F(y | x))
    AA1 <- intA(Ap, tau, ifun = FALSE) # AA_theta(1)

    g.i <- AA1 - Ay  # Note: no need of regularization on outR
    outL <- Fy$outL; g.i[outL,] <- AA1[outL,]
    g <- c(w %*% g.i)/eps
    fy <- NULL
  }
  else{g <- EE$g; g.i <- EE$g.i; Qtheta <- EE$Qtheta; Fy <- EE$Fy}

  # Hessian
  if(J){
    n <- length(y)
    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy$Fy))

    # f(y)
    FyL <- expit(logit(Fy$Fy) - 0.05)
    FyR <- expit(logit(Fy$Fy) + 0.05)
    fy <- 1/((Q(theta, FyR, data) - Q(theta, FyL, data))/(FyR - FyL))
    if(any(is.infinite(fy))) fy[is.infinite(fy)] <- 1e+200

    # A_theta_theta(p)
    Qthetatheta <- der2Qtheta(theta, eps, Q, Qtheta, data, tau)
    bp <- lapply(Qthetatheta, function(x){x*tau$wtau_long}) # a_theta_theta(p) = wtau(p)*Qthetatheta(p)
    Bp <- intA(bp, tau) # A_theta_theta(p)
    BB1 <- intA(Bp, tau, ifun = FALSE) # AA_theta_theta(1)
    By <- findAp(Fy, tau, Bp) # A_theta_theta(F(y | x))

    Fy$Fy[Fy$outL] <- tau$tau[2]; Fy$v[Fy$outL] <- 2 # See comment to findAp
    Qtheta_Fy <- findAp(Fy,tau,Qtheta)

    # Assemble

    npar <- length(theta)
    outL <- Fy$outL; outR <- Fy$outR

    h0 <- BB1 - By
    h0[outL,] <- BB1[outL,]
    h0 <- .colSums(w*h0, n, npar*(npar + 1)/2)
    H0 <- matrix(NA, npar, npar)
    H0[lower.tri(H0, diag = TRUE)] <- h0; H0 <- t(H0)
    H0[lower.tri(H0, diag = TRUE)] <- h0

    fy[c(outL, outR)] <- 0
    # if(sum(is.finite(fy)) > 0) fy[!is.finite(fy)] <- 0
    H1 <- -Qtheta_Fy*(w*fy)
    H2 <- -wy*Qtheta_Fy
    H1 <- crossprod(H1, H2)

    # Finish

    H <- H0 + H1
    H <- H/tcrossprod(eps) # each element is divided by both e1 and e2
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the hessian
  }

  list(g = g, g.i = g.i, J = J, Qtheta = Qtheta, Fy = Fy, fy = fy)
}

# Estimating equation and Jacobian for type = "c".
Qest.ee.c <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, tau$TAU, data) # Q(tau | x), a matrix n*ntau
    Qtheta <- derQtheta(theta, eps, Q, Q1, data, tau) # r matrices of size n*ntau
    Fy <- findp(y, tau, Q1) # F(y | x)
    Sy <- 1 - Fy$Fy # S(y | x)

    ap <- lapply(Qtheta, function(x){x*tau$wtau_long}) # a_theta(p) = wtau(p)*Qtheta(p)
    Ap <- intA(ap, tau) # A_theta(p)
    Ay <- findAp(Fy, tau, Ap) # A_theta(F(y | x))
    AAp <- intA(Ap, tau) # AA_theta(p)
    AAy <- findAp(Fy, tau, AAp) # AA_theta(F(y | x))
    AA1 <- sapply(AAp, function(x) x[,tau$ntau]) # AA_theta(1)

    g.i <- (AA1 - d*Ay) + (1 - d)/Sy*(AAy - AA1) # Note: no need of regularization on outR
    outL <- Fy$outL; g.i[outL,] <- d[outL]*AA1[outL,]
    g <- c(w %*% g.i)/eps
    fy <- NULL
  }
  else{
    g <- EE$g; g.i <- EE$g.i; Qtheta <- EE$Qtheta
    Fy <- EE$Fy; Sy <- EE$Sy; Ay <- EE$Ay; AAy <- EE$AAy; AA1 <- EE$AA1
  }

  # Hessian
  if(J){
    n <- length(y)
    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy$Fy))

    # f(y)
    FyL <- expit(logit(Fy$Fy) - 0.05)
    FyR <- expit(logit(Fy$Fy) + 0.05)
    fy <- 1/((Q(theta, FyR, data) - Q(theta, FyL, data))/(FyR - FyL))

    # A_theta_theta(p)
    Qthetatheta <- der2Qtheta(theta, eps, Q, Qtheta, data, tau)
    bp <- lapply(Qthetatheta, function(x){x*tau$wtau_long}) # a_theta_theta(p) = wtau(p)*Qthetatheta(p)
    Bp <- intA(bp, tau) # A_theta_theta(p)
    By <- findAp(Fy, tau, Bp) # A_theta_theta(F(y | x))
    BBp <- intA(Bp, tau) # AA_theta_theta(p)
    BBy <- findAp(Fy, tau, BBp) # AA_theta_theta(F(y | x))
    BB1 <- sapply(BBp, function(x) x[,tau$ntau]) # AA_theta_theta(1)

    Fy$Fy[Fy$outL] <- tau$tau[2]; Fy$v[Fy$outL] <- 2 # See comment to findAp
    Qtheta_Fy <- findAp(Fy,tau,Qtheta)

    # Assemble

    npar <- length(theta)
    outL <- Fy$outL; outR <- Fy$outR
    Uy <- (1 - d)/Sy

    h0 <- (BB1 - d*By) + Uy*(BBy - BB1)
    h0[outL,] <- d[outL]*BB1[outL,]
    h0 <- .colSums(w*h0, n, npar*(npar + 1)/2)

    H0 <- matrix(NA, npar, npar)
    H0[lower.tri(H0, diag = TRUE)] <- h0; H0 <- t(H0)
    H0[lower.tri(H0, diag = TRUE)] <- h0

    fy[c(outL, outR)] <- 0
    H1 <- -Qtheta_Fy*(w*fy)
    H2 <- -d*wy*Qtheta_Fy + Uy*(Ay + (AAy - AA1)/Sy)
    H1 <- crossprod(H1, H2)

    # Finish

    H <- H0 + H1
    H <- H/tcrossprod(eps) # each element is divided by both e1 and e2
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the hessian
  }

  list(g = g, g.i = g.i, J = J, Qtheta = Qtheta,
       Fy = Fy, Sy = Sy, fy = fy, Ay = Ay, AAy = AAy, AA1 = AA1)
}

# Estimating equation and Jacobian for type = "ct".
Qest.ee.ct <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, tau$TAU, data) # Q(tau | x), a matrix n*ntau
    Qtheta <- derQtheta(theta, eps, Q, Q1, data, tau) # r matrices of size n*ntau

    Fy <- findp(y, tau, Q1) # F(y | x)
    Sy <- 1 - Fy$Fy # S(y | x)
    Fz <- findp(z, tau, Q1) # F(z | x)
    Sz <- 1 - Fz$Fy # S(z | x)

    ap <- lapply(Qtheta, function(x){x*tau$wtau_long}) # a_theta(p) = wtau(p)*Qtheta(p)
    Ap <- intA(ap, tau) # A_theta(p)
    Ay <- findAp(Fy, tau, Ap) # A_theta(F(y | x))
    Az <- findAp(Fz, tau, Ap) # A_theta(F(z | x))
    AAp <- intA(Ap, tau) # AA_theta(p)
    AAy <- findAp(Fy, tau, AAp) # AA_theta(F(y | x))
    AAz <- findAp(Fz, tau, AAp) # AA_theta(F(z | x))
    AA1 <- sapply(AAp, function(x) x[,tau$ntau]) # AA_theta(1)


    gy.i <- -d*Ay + (1 - d)/Sy*(AAy - AA1)
    gz.i <- -1/Sz*(AAz - AA1)

    # Note: no need of regularization on outR
    outLy <- Fy$outL; outLz <- Fz$outL
    gy.i[outLy,] <- -(1 - d[outLy])*AA1[outLy,]
    gz.i[outLz,] <- AA1[outLz,]

    g.i <- gy.i + gz.i
    g <- c(w %*% g.i)/eps
    fy <- fz <- NULL
  }
  else{
    g <- EE$g; g.i <- EE$g.i; Qtheta <- EE$Qtheta
    Fy <- EE$Fy; Sy <- EE$Sy; Ay <- EE$Ay; AAy <- EE$AAy
    Fz <- EE$Fz; Sz <- EE$Sz; Az <- EE$Az; AAz <- EE$AAz
    AA1 <- EE$AA1
  }

  # Hessian
  if(J){
    n <- length(y)
    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy$Fy))

    # f(y)
    FyL <- expit(logit(Fy$Fy) - 0.05)
    FyR <- expit(logit(Fy$Fy) + 0.05)
    fy <- 1/((Q(theta, FyR, data) - Q(theta, FyL, data))/(FyR - FyL))

    # f(z)
    FzL <- expit(logit(Fz$Fy) - 0.05)
    FzR <- expit(logit(Fz$Fy) + 0.05)
    fz <- 1/((Q(theta, FzR, data) - Q(theta, FzL, data))/(FzR - FzL))

    # A_theta_theta(p)
    Qthetatheta <- der2Qtheta(theta, eps, Q, Qtheta, data, tau)
    bp <- lapply(Qthetatheta, function(x){x*tau$wtau_long}) # a_theta_theta(p) = wtau(p)*Qthetatheta(p)
    Bp <- intA(bp, tau) # A_theta_theta(p)
    By <- findAp(Fy, tau, Bp) # A_theta_theta(F(y | x))
    Bz <- findAp(Fz, tau, Bp) # A_theta_theta(F(z | x))
    BBp <- intA(Bp, tau) # AA_theta_theta(p)
    BBy <- findAp(Fy, tau, BBp) # AA_theta_theta(F(y | x))
    BBz <- findAp(Fz, tau, BBp) # AA_theta_theta(F(z | x))
    BB1 <- sapply(BBp, function(x) x[,tau$ntau]) # AA_theta_theta(1)

    Fy$Fy[Fy$outL] <- tau$tau[2]; Fy$v[Fy$outL] <- 2 # See comment to findAp
    Qtheta_Fy <- findAp(Fy,tau,Qtheta)
    Fz$Fy[Fz$outL] <- tau$tau[2]; Fz$v[Fz$outL] <- 2 # See comment to findAp
    Qtheta_Fz <- findAp(Fz,tau,Qtheta)

    # Assemble

    npar <- length(theta)
    outLy <- Fy$outL; outLz <- Fz$outL; outRy <- Fy$outR; outRz <- Fz$outR
    Uy <- (1 - d)/Sy
    Uz <- 1/Sz

    h0y <- -d*By + Uy*(BBy - BB1)
    h0y[outLy,] <- -(1 - d[outLy])*BB1[outLy,]
    h0z <- -Uz*(BBz - BB1)
    h0z[outLz,] <- BB1[outLz,]
    h0 <- h0z + h0y
    h0 <- .colSums(w*h0, n, npar*(npar + 1)/2)

    H0 <- matrix(NA, npar, npar)
    H0[lower.tri(H0, diag = TRUE)] <- h0; H0 <- t(H0)
    H0[lower.tri(H0, diag = TRUE)] <- h0

    fy[c(outLy, outRy)] <- 0
    fz[c(outLz, outRz)] <- 0
    H1y <- -Qtheta_Fy*(w*fy)
    H1z <- -Qtheta_Fz*(w*fz)

    H2y <- -d*wy*Qtheta_Fy + Uy*(Ay + (AAy - AA1)/Sy)
    H2z <- -Uz*(Az + (AAz - AA1)/Sz)

    H1y <- crossprod(H1y, H2y)
    H1z <- crossprod(H1z, H2z)

    # Finish

    H <- H0 + H1y + H1z
    H <- H/tcrossprod(eps) # each element is divided by both e1 and e2
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the hessian
  }

  list(g = g, g.i = g.i, J = J, Qtheta = Qtheta,
       Fy = Fy, Sy = Sy, fy = fy, Ay = Ay, AAy = AAy,
       Fz = Fz, Sz = Sz, fz = fz, Az = Az, AAz = AAz,
       AA1 = AA1)
}
