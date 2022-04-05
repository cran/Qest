
# Basis of a piecewise linear function.
plfcox <- function(y, knots, deriv = 0){ # knots incl. boundary

  if(length(knots) == 2){
    if(deriv == 0){return(cbind(y - knots[1]))}
    else{return(list(b = cbind(y - knots[1]), b1 = cbind(rep(1,length(y)))))}
  }
  intknots <- knots[2:(length(knots) - 1)]
  k <- length(intknots)
  ind1 <- 1
  ind2 <- NULL
  for (j in intknots){ind1 <- cbind(ind1,(y > j))}
  for (j in k:1) {ind2 <- cbind(ind2, ind1[, j] - ind1[, j + 1])}
  ind2 <- cbind(ind1[, k + 1], ind2)[, (k + 1):1, drop = FALSE]
  ind1 <- cbind(ind1, 0)

  b <- NULL
  for (j in 1:(k + 1)){b <- cbind(b, (y - knots[j])*ind2[, j] + (knots[j + 1] - knots[j])*ind1[, j + 1])}
  colnames(b) <- paste("b", 1:(k + 1), sep = "")
  attr(b, "knots") <- knots
  if(deriv == 0){return(b)}
  else{return(list(b = b, b1 = ind2))}
}

# Transform the scale of y and X in a Cox model.
scalevars.Qcoxph <- function(X,z,y,knots){

  mX <- colMeans(X)
  sX <- apply(X, 2, sd)
  X <- scale(X, center = mX, scale = sX)

  My <- max(y)
  y <- y/My*10
  z <- z/My*10
  knots <- knots/My*10
  Stats <- list(mX = mX, sX = sX, My = My)

  list(X = X, z = z, y = y, knots = knots, Stats = Stats)
}

descalecoef.Qcoxph <- function(theta, Stats){

  q <- length(Stats$mX)
  beta <- theta[1:q]
  log.gamma <- theta[(q + 1):length(theta)]

  beta <- beta/Stats$sX
  log.gamma <- log.gamma + log(10/Stats$My) - sum(beta*Stats$mX)

  c(beta, log.gamma)
}

#### check singularities
check.singularities <- function(X, scaleVars){
  QR <- qr(X); nok <- (abs(diag(qr.R(QR))) < 1e-10); ok <- !nok
  rank <- QR$rank
  QR <- structure(QR[c("qr", "qraux", "pivot", "rank")], class = "qr")
  X <- X[, ok, drop = FALSE]
  scaleVars$X <- scaleVars$X[,ok,drop=FALSE]
  scaleVars$Stats$mX <- scaleVars$Stats$mX[ok]
  scaleVars$Stats$sX <- scaleVars$Stats$sX[ok]
  return(list(X = X, scaleVars = scaleVars, nok = nok, ok = ok, rank = rank, QR = QR))
}

#### starting points
starting.points.Qcox <- function(X, Y, n, w, mf, knots){
  # ok <- (Y[,1] >= quantile(Y[,1], .1)) & (Y[,1] <= quantile(Y[,1], .9))
  ok <- (Y[,1] <= quantile(Y[,1], .95))
  # ok <- rep(TRUE, nrow(Y))
  m0 <- coxph.fit(X[ok,, drop = FALSE], Y[ok,], NULL, rep(0, sum(ok)), NULL,
    coxph.control(), w[ok], method = "efron", row.names(mf)[ok])
  Xs <- X[ok,, drop = FALSE]; X2s <- matrix(0, 1, ncol(X))
  risk <- drop(exp(Xs %*% m0$coefficients))
  newrisk <- drop(exp(X2s %*% m0$coefficients))
  a0 <- survfitcoxph.fit(Y[ok,], Xs, w[ok], X2s, risk, newrisk, rep(0, sum(ok)), TRUE, 3, 3, m0$var)
  H0 <- a0$cumhaz; t0 <- a0$time
# plot(t0, H0, pch=16)
#
# b0 <- plfcox(t0, knots); M0 <- lm.fit(b0, H0)
# lines(t0, M0$fitted, col=2)

  if(attr(knots, "user-defined")) {
    nknots <- length(knots) - 2
    Z <- matrix(t0, nrow=length(t0), ncol=nknots)
    PSI <- matrix(rep(knots[2:(nknots+1)], each=length(t0)), ncol=nknots)
    temp <- seg.lm.fit1(H0, t0, Z, PSI)$psi
    if(length(temp) == nknots) knots[2:(nknots+1)] <- temp
  }
  b0 <- plfcox(t0, knots); M0 <- lm.fit(b0, H0)
# lines(t0, M0$fitted, col=3)

  M0$coef <- adjust.coef(M0$coef)
  theta <- c(m0$coef, log(M0$coef))
  attr(theta, "knots") <- knots
  theta
}

adjust.coef <- function(theta){
  id <- which(theta <= 0)
  id1 <- which(theta > 0)
  k <- length(theta)
  if(length(id) != 0){
    for(i in rev(id)){
      theta[i] <- if(i == 1) theta[i+1]/k else median(theta[-unique(c(k:i, id))])/k
    }
  }
  theta
}

# Control parameters

Qcoxph.control <- function(tol = 1e-8, maxit, safeit, alpha0, display = FALSE){
  if(missing(maxit)){maxit <- NULL}
  if(missing(safeit)){safeit <- NULL}
  if(missing(alpha0)){alpha0 <- NULL}
  list(tol = tol, maxit = maxit, safeit = safeit, alpha0 = alpha0, display = display)
}

seg.lm.fit1<-function(y,XREG,Z,PSI,return.all.sol=FALSE){
  opz <- list(toll = 1e-3, h = .1, stop.if.error = NULL,
              dev0 = 100, visual = FALSE, it.max = 30, gap = FALSE,
              pow = c(1, 1), digits = NULL,
              conv.psi = FALSE, alpha = .02, fix.npsi = FALSE,
              min.step = .0001)

  useExp.k <- TRUE
  #-----------------
  XREG <- as.matrix(XREG)
  w <- rep(1, length(y))
  offs <- rep(0, length(y))
  fix.npsi <- TRUE #FALSE
  nomiOK <- paste("U",1:ncol(PSI),".x",sep="")
  id.psi.group <- rep(1, ncol(PSI))
  #-----------------
  est.k<-function(x1,y1,L0){
    ax<-log(x1)
    .x<-cbind(1,ax,ax^2)
    b<-drop(solve(crossprod(.x),crossprod(.x,y1)))
    const<-b[1]-L0
    DD<-sqrt(b[2]^2-4*const*b[3])
    kk<-exp((-b[2]+ DD) /(2*b[3]))
    return(round(kk))
  }
  #-----------------
  dpmax<-function(x,y,pow=1){
    #deriv pmax
    if(pow==1) -(x>y) #ifelse(x>y, -1, 0)
    else -pow*((x-y)*(x>y))^(pow-1)#-pow*pmax(x-y,0)^(pow-1)
  }
  #-----------
  mylm<-function(x,y,w,offs=rep(0,length(y))){
    x1<-x*sqrt(w)
    y<-y-offs
    y1<-y*sqrt(w)
    b<-drop(solve(crossprod(x1),crossprod(x1,y1)))
    fit<-drop(tcrossprod(x,t(b)))
    r<-y-fit
    o<-list(coefficients=b,fitted.values=fit,residuals=r, df.residual=length(y)-length(b))
    o
  }
  #-----------
  mylmADD<-function(invXtX, X, v, Xty, y){
    #v: new column to be added
    vtv<-sum(v^2)
    Xtv<-crossprod(X,v) #-colSums(X[v!=0,,drop=FALSE]) #oppure -.colSums(X[v!=0,,drop=FALSE],n,p)
    m<-invXtX %*% Xtv
    d<-drop(1/(vtv- t(Xtv) %*% m))
    r<- -d*m
    invF <- invXtX + d*tcrossprod(m)
    newINV<- rbind(cbind(invF, r), c(t(r), d))
    b<-crossprod(newINV, c(Xty, sum(v*y)))
    fit<- tcrossprod(cbind(X,v), t(b)) #cbind(X,v) %*% b
    r<-y-fit
    o<-list(coefficients=b,fitted.values=fit,residuals=r)
    o
  }
  #-----------
  in.psi<-function(LIM, PSI, ret.id=TRUE){
    #check if psi is inside the range
    a<-PSI[1,]<=LIM[1,]
    b<-PSI[1,]>=LIM[2,]
    is.ok<- !a & !b #TRUE se psi e' OK
    if(ret.id) return(is.ok)
    isOK<- all(is.ok) && all(!is.na(is.ok))
    isOK}
  #------------
  far.psi<-function(Z, PSI, id.psi.group, ret.id=TRUE, fc=.93) {
    #check if psi are far from the boundaries ..s
    #   returns TRUE, if fine.
    #id.far.ok<-sapply(unique(id.psi.group), function(.x) (table(rowSums(((Z>PSI)[,id.psi.group==.x,drop=FALSE])))>=2)[-1]) #[-1] esclude lo zero, x<psi[1]
    #id.far.ok<-sapply(unique(id.psi.group), function(.x) (tabulate(rowSums(((Z>PSI)[,id.psi.group==.x,drop=FALSE]))+1)>=2)[-1]) #[-1] esclude lo zero, x<psi[1]
    #16/01/20:
    # se un psi assume l'estremo superiore "Z>PSI" non se ne accorge, mentre Z>=PSI, sì.. Il contrario è vero con estremo inf e Z>PSI
    nSeg<-length(unique(id.psi.group))
    npsij<-tapply(id.psi.group,id.psi.group,length)
    nj<-sapply(unique(id.psi.group), function(.x) { tabulate(rowSums((Z>PSI)[,id.psi.group==.x,drop=FALSE])+1) }, simplify = FALSE)
    ff<-id.far.ok<-vector("list",length=nSeg)
    for(i in 1:nSeg){
      if(length(nj[[i]])!=npsij[i]+1) nj[[i]]<- tabulate(rowSums((Z>=PSI)[,id.psi.group==i,drop=FALSE])+1)
      id.ok<-(nj[[i]] >= 2)
      id.far.ok[[i]] <- id.ok[-length(id.ok)] & id.ok[-1] #esattamente uguale al numero di psi del gruppo i
      ff[[i]]<-ifelse(diff(nj[[i]])>0, 1/fc, fc)
    }
    id.far.ok<-unlist(id.far.ok)
    ff<-unlist(ff)
    if(!ret.id) {return(all(id.far.ok))
    } else {
      attr(id.far.ok,"factor") <- ff
      return(id.far.ok)
    }
    #if(ret.id) return(id.far.ok) else return(all(id.far.ok))
  } #end far.psi
  #-----------
  adj.psi<-function(psii, LIM) {pmin(pmax(LIM[1,],psii),LIM[2,])}
  #-----------
  n <- length(y)
  min.step <- opz$min.step
  rangeZ <- apply(Z, 2, range)
  alpha <- opz$alpha
  limZ <- apply(Z, 2, quantile, names=FALSE, probs=c(alpha, 1-alpha))
  psi <- PSI[1,]
  conv.psi <- opz$conv.psi
  h <- opz$h
  digits <- opz$digits
  pow <- opz$pow
  toll <- opz$toll
  gap <- opz$gap
  fix.npsi <- opz$fix.npsi
  dev.new <- opz$dev0
  visual <- opz$visual
  it.max <- old.it.max <- opz$it.max
  fc <- .95 #opz$fc
  names(psi) <- id.psi.group
  it <- 0
  epsilon <- 10
  k.values <- dev.values<- NULL
  psi.values <-list()
  psi.values[[length(psi.values) + 1]] <- NA

  #id.psi.ok<-rep(TRUE, length(psi))
  sel.col.XREG <- unique(sapply(colnames(XREG), function(x)match(x,colnames(XREG))))
  if(is.numeric(sel.col.XREG)) XREG <- XREG[,sel.col.XREG,drop=FALSE] #elimina le ripetizioni, ad es. le due intercette..
  #==================
  invXtX <- solve(crossprod(XREG))
  Xty <- drop(crossprod(XREG, y))
  #===================
  #browser()
  if(!in.psi(limZ, PSI, FALSE)) stop("starting psi out of the range", call.=FALSE)
  if(!far.psi(Z, PSI, id.psi.group, FALSE)) stop("psi values too close each other. Please change (decreases number of) starting values", call.=FALSE)
  n.psi1 <- ncol(Z)
  #==============================================
  U <- ((Z-PSI)*(Z>PSI)) #pmax((Z - PSI), 0)^pow[1]
  if(pow[1]!=1) U<-U^pow[1]
  obj0 <- mylm(cbind(XREG, U), y, w, offs) #lm.wfit(cbind(XREG, U), y, w, offs) #se 1 psi, si puo' usare la funz efficiente..
  L0 <- sum(obj0$residuals^2*w)
  n.intDev0<-nchar(strsplit(as.character(L0),"\\.")[[1]][1])
  dev.values[length(dev.values) + 1] <- opz$dev0 #del modello iniziale (senza psi)
  dev.values[length(dev.values) + 1] <- L0 #modello con psi iniziali
  psi.values[[length(psi.values) + 1]] <- psi #psi iniziali
  #==============================================
  if (visual) {
    cat(paste("iter = ", sprintf("%2.0f",0),
              "  dev = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), L0), #formatC(L1,width=8, digits=5,format="f"), #era format="fg"
              "  k = ", sprintf("%2.0f", NA),
              "  n.psi = ",formatC(length(unlist(psi)),digits=0,format="f"),
              "  ini.psi = ",paste(formatC(unlist(psi),digits=4,format="f"), collapse="  "), #sprintf('%.2f',x)
              sep=""), "\n")
  }
  #==============================================
  id.warn <- FALSE
  id.psi.changed <- rep(FALSE, it.max)
  id.fail.step <- FALSE
  ############################################# START WHILE
  while (abs(epsilon) > toll) {
    it <- it+1
    n.psi0 <- n.psi1
    n.psi1 <- ncol(Z)
    if(n.psi1!=n.psi0){
      U <- ((Z-PSI)*(Z>PSI)) #pmax((Z - PSI), 0)^pow[1]
      if(pow[1]!=1) U <- U^pow[1]
      obj0 <- mylm(cbind(XREG, U), y, w, offs)#lm.wfit(cbind(XREG, U), y, w, offs) #se 1 psi, si puo' usare la funz efficiente..
      L0 <- sum(obj0$residuals^2*w)
    }
    V <- dpmax(Z,PSI,pow=pow[2])# ifelse((Z > PSI), -1, 0)
    X <- cbind(XREG, U, V)
    rownames(X) <- NULL
    colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c(paste("U", 1:ncol(U), sep = ""), paste("V", 1:ncol(V), sep = ""))
    obj <- lm.wfit(x = X, y = y, w = w, offset = offs) #mylm(X, y, w, offs) #
    beta.c <- coef(obj)[paste("U", 1:ncol(U), sep = "")]
    gamma.c <- coef(obj)[paste("V", 1:ncol(V), sep = "")]

    if(any(is.na(c(beta.c, gamma.c)))){
      return(list(psi=psi.last))
    }
    psi.old <- psi
    psi <- psi.old + h*gamma.c/beta.c
    if(!is.null(digits)) psi <- round(psi, digits)
    PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
    #--modello con il nuovo psi
    U1 <- (Z-PSI)*(Z>PSI)
    if(pow[1]!=1) U1<-U1^pow[1]
    obj1 <- try(mylm(cbind(XREG, U1), y, w, offs), silent = TRUE) #lm.wfit(cbind(XREG, pmax(Z-PSI,0)), y, w, offs)
    L1<- if(class(obj1)[1]=="try-error") L0+10 else sum(obj1$residuals^2*w)

    use.k <- k <- 1
    L1.k <- NULL
    L1.k[length(L1.k)+1] <- L1

    id.fail.step <- FALSE
    psi.last <- psi.values[[it+1]] #old #<===== psi o psi.old???!?!??!?!????????!!!!!!!!!!
    while(L1 > L0){
      #if(it==5) browser()
      #ATTENZIONE: i gamma.c e beta.c vengono dal modello, ma poi dopo il modello (linee 152-167) viene fatto un controllo che puo' eliminare break e ridurre le colonne di Z.
      #Per cui puo' risultare ncol(PSI)>ncol(Z). Quindi o non si fanno i controlli ( potrebbe essere perche' tanto c'e' il try(..)) oppure semplicemente
      # si prendono le stime corrispondenti alle colonne "ok". psi <- psi.old[id.psi.ok] + (gamma.c[id.psi.ok]/beta.c[id.psi.ok])/(use.k*h)
      #browser()
      k<-k+1
      use.k <- if(useExp.k) 2^(k-1) else k
      psi <- psi.old + (gamma.c/beta.c)/(use.k*h)
      if(!is.null(digits)) psi <- round(psi, digits)
      PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
      #qui o si aggiusta psi per farlo rientrare nei limiti, o si elimina, oppure se obj1 sotto non funziona semplicemente continua..
      U1 <- (Z-PSI)*(Z>PSI)
      if(pow[1]!=1) U1 <- U1^pow[1]
      obj1 <- try(mylm(cbind(XREG, U1), y, w, offs), silent=TRUE) #lm.wfit(cbind(X,U1), y, w, offs)
      L1 <- if(class(obj1)[1]=="try-error") L0+10 else sum(obj1$residuals^2*w)
      L1.k[length(L1.k)+1] <- L1
      if(1/(use.k*h)<min.step){
        #        #warning("step halving too small")
        id.fail.step<-TRUE
        break
      }
    } #end while L0-L1
    if (visual) {
      cat(paste("iter = ", sprintf("%2.0f",it),
                "  dev = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), L1), #formatC(L1,width=8, digits=5,format="f"), #era format="fg"
                "  k = ", sprintf("%2.0f", k),
                "  n.psi = ",formatC(length(unlist(psi)),digits=0,format="f"),
                "  est.psi = ",paste(formatC(unlist(psi),digits=4,format="f"), collapse="  "), #sprintf('%.2f',x)
                sep=""), "\n")
    }
    #if step halving failed, come back to the "initial" (before the stephalving) values and stop with a warning message
    if(id.fail.step) {
      #browser()
      psi.old <- psi <- psi.last
      PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
      L1 <- L0
      U1 <- U
    }

    epsilon <- if(conv.psi) max(abs((psi - psi.old)/psi.old)) else (L0 - L1)/(abs(L0) + 0.1)
    L0 <- L1
    U <-U1

    k.values[length(k.values)+1] <- use.k
    psi.values[[length(psi.values) + 1]] <- psi
    dev.values[length(dev.values) + 1] <- L0

    #Mi sa che non servono i controlli.. soprattutto se non ha fatto step-halving
    #check if i psi ottenuti sono nel range o abbastanza lontani
    id.psi.far <- far.psi(Z, PSI, id.psi.group, TRUE, fc)
    id.psi.in <- in.psi(limZ, PSI, TRUE)
    id.psi.ok <- id.psi.in & id.psi.far

    if(!all(id.psi.ok)){
      return(list(psi=psi.last))
    }
    if (it >= it.max) {
      id.warn <- TRUE
      break
    }
  }
  ##=============================================================================
  id.psi.changed <- id.psi.changed[1:it]

  attr( psi.values, "dev") <- dev.values
  attr( psi.values, "k")<- k.values

  #ordina i breakpoints..
  psi <- unlist(tapply(psi, id.psi.group, sort))
  names(psi) <- id.psi.group
  names.coef <- names(obj$coefficients) #obj e' quello vecchio che include U1,.. V1,...

  PSI.old <- PSI
  PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))

  #U e V possono essere cambiati (rimozione/ordinamento psi.. ) per cui si deve ricalcolare il tutto, altrimenti sarebbe uguale a U1 e obj1
  #  browser()

  if(sd(PSI-PSI.old) > 0 || id.psi.changed[length(id.psi.changed)] || id.fail.step){
    U <- (Z-PSI)*(Z>PSI)
    colnames(U)<-paste("U", 1:ncol(U), sep = "")
    V <- -(Z>PSI)
    colnames(V)<-paste("V", 1:ncol(V), sep = "")
    obj <- lm.wfit(x = cbind(XREG, U), y = y, w = w, offset = offs)
    L1 <- sum(obj$residuals^2*w)
  }
  else {
    obj<-obj1
  }
  obj$coefficients<-c(obj$coefficients, rep(0,ncol(V)))
  names(obj$coefficients)<-names.coef
  obj$epsilon <- epsilon
  obj$it <- it

  obj<-list(obj=obj,it=it,psi=psi, psi.values=psi.values, U=U,V=V,rangeZ=rangeZ,
            epsilon=epsilon,nomiOK=nomiOK, SumSquares.no.gap=L1, id.psi.group=id.psi.group,
            id.warn=id.warn) #inserire id.psi.ok?
  return(obj)
}

