source("utils.R")



### Robust CP decomposition with regularization ###
RobRSVD_tnsr <- function(data, irobust=FALSE, rtype="huber", huberk=1.345, bisqk=3.88, inits, initu, initv, initw, lambdaRange, nlambda, degree=1,
                         niter=1000, tol=1e-5, istablize=FALSE, uspar=0, vspar=0, wspar=0, 
                         iugcv=FALSE, ivgcv=FALSE, iwgcv=FALSE,
                         iloguspar=TRUE, ilogvspar=TRUE, ilogwspar=TRUE, isspline=FALSE, verbose=FALSE) {
  
  # Robust 3-way functional prinicpal component analysis
  # data is 3-mode tensor (m x n x p)
  
  ugcvmat = vgcvmat = wgcvmat = as.numeric()
  
  # check if the parameters are appropriate #
  size_data = c(dim(data))
  m = size_data[1]
  n = size_data[2]
  p = size_data[3]
  if (istablize) {
    # The original data is centered by the median of the whole tensor
    myscale = 1.4785*median(abs(c(data)-median(data)))
    localdata = data / myscale
  } else {
    myscale = 1
    localdata = data
  }
  
  uold = initu
  vold = initv
  wold = initw
  sold = inits
  if (istablize) {
    sold = sold / myscale
  } 
  uold = sold * uold
  
  # use the initial estiamte of u, v, w vectors to reconstruct data tensor 
  Appold = uold %o% vold %o% wold 
  Rmat = localdata - Appold
  Rvec = c(Rmat)
  mysigma = median(abs(Rvec))/0.675
  iter = 1
  localdiff = 9999
  diffvec = as.numeric()
  nuspar = nlambda$u
  nvspar = nlambda$v
  nwspar = nlambda$w
  uspar_current = uspar
  vspar_current = vspar
  wspar_current = wspar
  if (isspline) { # spline smoother 
    Omega.u = ssmatls(m)$y
    Omega.v = ssmatls(n)$y
    Omega.w = ssmatls(p)$y
  } else { # difference smoother
    Omega.u = crossprod(makeDiffOper(degree=degree, dim=m))
    Omega.v = crossprod(makeDiffOper(degree=degree, dim=n))
    Omega.w = crossprod(makeDiffOper(degree=degree, dim=p))
  }
  
  ugcvscore = vgcvscore = wgcvscore = as.numeric()
  
  # main iteration
  while (localdiff > tol & iter < niter) {
    if (irobust) {
      if (rtype == "huber") {
        Wmat = huberWeightLS(Rmat/mysigma, huberk)
      } else {
        Wmat = bisqWeightLS(Rmat/mysigma, bisqk)
      }
    } else {
      Wmat = array(1, dim=c(m,n,p))
    }
    
    # update u given (v, w)
    if (!iugcv) {
      S.u.inv = diag(m) + uspar*Omega.u
      S.v.inv = diag(n) + vspar*Omega.v
      S.w.inv = diag(p) + wspar*Omega.w
      
      au = c(vold %o% wold) 
      Wmat.rtns = as.tensor(Wmat) # rTensor object
      Omega.u.vw = (c(t(vold) %*% S.v.inv %*% vold)*c(t(wold) %*% S.w.inv %*% wold)*S.u.inv - diag(sum(vold^2),m)*diag(sum(wold^2),m))*(2*mysigma^2)
      tAuWuAu = diag(colSums(diag(c(au^2)) %*%  t(rs_unfold(Wmat.rtns, m=1)@data))) 
      tAuWuYu = (rs_unfold(Wmat.rtns, m=1)@data * rs_unfold(as.tensor(localdata), m=1)@data) %*% au
      uterm1 = tAuWuAu + Omega.u.vw
      uterm2 = tAuWuYu
      unew = solve(uterm1, uterm2)

    } else {
      # choosing the penalty parameter by cross-validation according to GCV criteria
      if (nuspar<0)
        stop("number of smoothing parameter can not be negative")
      else {
        if (iloguspar) {
          usparvec = 10^seq(log10(lambdaRange$u[1]), log10(lambdaRange$u[2]), length.out=nuspar)
        } else {
          usparvec = seq(lambdaRange$u[1], lambdaRange$u[2], length.out=nuspar)
        } 
        ugcvvec = as.numeric()
        ugcvmat = as.numeric()
        
        au = c(vold %o% wold) 
        Wmat.rtns = as.tensor(Wmat)
        tAuWuAu = diag(colSums(diag(c(au^2)) %*%  t(rs_unfold(Wmat.rtns, m=1)@data))) 
        tAuWuYu = (rs_unfold(Wmat.rtns, m=1)@data * rs_unfold(as.tensor(localdata), m=1)@data) %*% au
        u_nsrobust = solve(tAuWuAu, tAuWuYu)
        
        S.v.inv = diag(n) + vspar_current*Omega.v
        S.w.inv = diag(p) + wspar_current*Omega.w
        
        for (iter_uspar in 1:nuspar) {
          S.u.inv = diag(m) + usparvec[iter_uspar]*Omega.u
          Omega.u.vw = (c(t(vold) %*% S.v.inv %*% vold)*c(t(wold) %*% S.w.inv %*% wold)*S.u.inv - diag(sum(vold^2),m)*diag(sum(wold^2),m))*(2*mysigma^2)
          usterm1 = tAuWuAu + Omega.u.vw
          usterm2 = tAuWuYu
          u_srobust = solve(usterm1, usterm2)
          smooth_u = solve(usterm1, tAuWuAu)
          gcv_ut = m*sum((u_nsrobust-u_srobust)^2)/(m - sum(diag(smooth_u)))^2
          ugcvvec[iter_uspar] = gcv_ut
          ugcvmat = cbind(ugcvmat, u_srobust/(sqrt(sum(u_srobust^2))))
        }
        uspar_current = usparvec[which.min(ugcvvec)]
        ugcvscore = cbind(usparvec, ugcvvec)
        
        S.u.inv = diag(m) + uspar_current*Omega.u
        Omega.u.vw = (c(t(vold) %*% S.v.inv %*% vold)*c(t(wold) %*% S.w.inv %*% wold)*S.u.inv - diag(sum(vold^2),m)*diag(sum(wold^2),m))*(2*mysigma^2)
        #tAuWuAu = diag(colSums(diag(c(au^2)) %*%  t(rs_unfold(Wmat.rtns, m=1)@data))) 
        #tAuWuYu = (rs_unfold(Wmat.rtns, m=1)@data * rs_unfold(as.tensor(localdata), m=1)@data) %*% au
        uterm1 = tAuWuAu + Omega.u.vw
        uterm2 = tAuWuYu
        unew = solve(uterm1, uterm2)
      }
    }
    
    # update v given (u, w)  
    if (!ivgcv) {
      S.u.inv = diag(m) + uspar*Omega.u
      S.v.inv = diag(n) + vspar*Omega.v
      S.w.inv = diag(p) + wspar*Omega.w
      
      av = c(unew %o% wold) 
      Wmat.rtns = as.tensor(Wmat) # rTensor object
      Omega.v.uw = (c(t(unew) %*% S.u.inv %*% unew)*c(t(wold) %*% S.w.inv %*% wold)*S.v.inv - diag(sum(unew^2),n)*diag(sum(wold^2),n))*(2*mysigma^2)
      tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
      tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
      vterm1 = tAvWvAv + Omega.v.uw
      vterm2 = tAvWvYv
      vnew = solve(vterm1, vterm2)
    } else {
      # choosing the penalty parameter by cross-validation according to GCV criteria
      if (nvspar<0)
        stop("number of smoothing parameter can not be negative")
      else {
        if (ilogvspar) {
          vsparvec = 10^seq(log10(lambdaRange$v[1]), log10(lambdaRange$v[2]), length.out=nvspar)
        } else {
          vsparvec = seq(lambdaRange$v[1], lambdaRange$v[2], length.out=nvspar)
        } 
        vgcvvec = as.numeric()
        vgcvmat = as.numeric()
        
        av = c(unew %o% wold) 
        Wmat.rtns = as.tensor(Wmat)
        tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
        tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
        v_nsrobust = solve(tAvWvAv, tAvWvYv)
        
        S.u.inv = diag(m) + uspar_current*Omega.u
        S.w.inv = diag(p) + wspar_current*Omega.w
        
        for (iter_vspar in 1:nvspar) {
          S.v.inv = diag(n) + vsparvec[iter_vspar]*Omega.v
          Omega.v.uw = (c(t(unew) %*% S.u.inv %*% unew)*c(t(wold) %*% S.w.inv %*% wold)*S.v.inv - diag(sum(unew^2),n)*diag(sum(wold^2),n))*(2*mysigma^2)
          vsterm1 = tAvWvAv + Omega.v.uw
          vsterm2 = tAvWvYv
          v_srobust = solve(vsterm1, vsterm2)
          smooth_v = solve(vsterm1, tAvWvAv)
          gcv_vt = n*sum((v_nsrobust-v_srobust)^2)/(n - sum(diag(smooth_v)))^2
          vgcvvec[iter_vspar] = gcv_vt
          vgcvmat = cbind(vgcvmat, v_srobust/(sqrt(sum(v_srobust^2))))
        }
        vspar_current = vsparvec[which.min(vgcvvec)]
        vgcvscore = cbind(vsparvec, vgcvvec)
        
        S.v.inv = diag(n) + vspar_current*Omega.v
        Omega.v.uw = (c(t(unew) %*% S.u.inv %*% unew)*c(t(wold) %*% S.w.inv %*% wold)*S.v.inv - diag(sum(unew^2),n)*diag(sum(wold^2),n))*(2*mysigma^2)
        #tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
        #tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
        vterm1 = tAvWvAv + Omega.v.uw
        vterm2 = tAvWvYv
        vnew = solve(vterm1, vterm2)
      }
    }
    
    # update w given (u, v)  
    if (!iwgcv) {
      S.u.inv = diag(m) + uspar*Omega.u
      S.v.inv = diag(n) + vspar*Omega.v
      S.w.inv = diag(p) + wspar*Omega.w
      
      aw = c(unew %o% vnew) 
      Wmat.rtns = as.tensor(Wmat) # rTensor object
      Omega.w.uv = (c(t(unew) %*% S.u.inv %*% unew)*c(t(vnew) %*% S.v.inv %*% vnew)*S.w.inv - diag(sum(unew^2),p)*diag(sum(vnew^2),p))*(2*mysigma^2)
      tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
      tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
      wterm1 = tAwWwAw + Omega.w.uv
      wterm2 = tAwWwYw
      wnew = solve(wterm1, wterm2)
    } else {
      # choosing the penalty parameter by cross-validation according to GCV criteria
      if (nwspar<0)
        stop("number of smoothing parameter can not be negative")
      else {
        if (ilogwspar) {
          wsparvec = 10^seq(log10(lambdaRange$w[1]), log10(lambdaRange$w[2]), length.out=nwspar)
        } else {
          wsparvec = seq(lambdaRange$w[1], lambdaRange$w[2], length.out=nwspar)
        } 
        wgcvvec = as.numeric()
        wgcvmat = as.numeric()
        
        aw = c(unew %o% vnew) 
        Wmat.rtns = as.tensor(Wmat)
        tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
        tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
        w_nsrobust = solve(tAwWwAw, tAwWwYw)
        
        S.u.inv = diag(m) + uspar_current*Omega.u
        S.v.inv = diag(n) + vspar_current*Omega.v
        
        for (iter_wspar in 1:nwspar) {
          S.w.inv = diag(p) + wsparvec[iter_wspar]*Omega.w
          Omega.w.uv = (c(t(unew) %*% S.u.inv %*% unew)*c(t(vnew) %*% S.v.inv %*% vnew)*S.w.inv - diag(sum(unew^2),p)*diag(sum(vnew^2),p))*(2*mysigma^2)
          wsterm1 = tAwWwAw + Omega.w.uv
          wsterm2 = tAwWwYw
          w_srobust = solve(wsterm1, wsterm2)
          smooth_w = solve(wsterm1, tAwWwAw)
          gcv_wt = p*sum((w_nsrobust-w_srobust)^2)/(p - sum(diag(smooth_w)))^2
          wgcvvec[iter_wspar] = gcv_wt
          wgcvmat = cbind(wgcvmat, w_srobust/(sqrt(sum(w_srobust^2))))
        }
        wspar_current = wsparvec[which.min(wgcvvec)]
        wgcvscore = cbind(wsparvec, wgcvvec)
        
        S.w.inv = diag(p) + wspar_current*Omega.w
        Omega.w.uv = (c(t(unew) %*% S.u.inv %*% unew)*c(t(vnew) %*% S.v.inv %*% vnew)*S.w.inv - diag(sum(unew^2),p)*diag(sum(vnew^2),p))*(2*mysigma^2)
        #tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
        #tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
        wterm1 = tAwWwAw + Omega.w.uv
        wterm2 = tAwWwYw
        wnew = solve(wterm1, wterm2)
      }
    }  
    
    Appnew = c(unew) %o% c(vnew) %o% c(wnew)
    Rmat = localdata - Appnew
    
    localdiff = max(abs(Appnew-Appold))
    Appold = Appnew
    uold = sqrt(sum(vnew^2))* sqrt(sum(wnew^2))*unew
    vold = vnew / sqrt(sum(vnew^2)) # v is normalized to length 1
    wold = wnew / sqrt(sum(wnew^2)) # w is normalized to length 1
    if (verbose) {
      print(paste("iter", iter, "fin"))
      print(paste("diff:",round(localdiff,4)))
    }
  
    iter = iter + 1
    diffvec = c(diffvec, localdiff)
    
  } # end of while
  
  v = vold
  w = wold
  s = myscale * sqrt(sum(uold^2))
  u = uold / sqrt(sum(uold^2))
  
  if (iugcv) {
    uspar = uspar_current
  }
  if (ivgcv) {
    vspar = vspar_current
  }
  if (iwgcv) {
    wspar = wspar_current
  }
  
  diagout = list(ugcvsore=ugcvscore, vgcvscore=vgcvscore, wgcvscore=wgcvscore, 
                 ugcvmat=ugcvmat, vgcvmat=vgcvmat, wgcvmat=wgcvmat, diffvec=diffvec)
  return(list(s=s, u=u, v=v, w=w, diagout=diagout))
}

wrap_RobRSVD_tnsr <- function(data, irobust=FALSE, rtype="huber", huberk=1.345, bisqk=3.88, R, inits, initu, initv, initw, lambdaRange, nlambda, degree=1,
                              niter=1000, tol=1e-5, istablize=FALSE, uspar=0, vspar=0, wspar=0, 
                              iugcv=FALSE, ivgcv=FALSE, iwgcv=FALSE,
                              iloguspar=TRUE, ilogvspar=TRUE, ilogwspar=TRUE, verbose=FALSE, isspline=FALSE) {
  
  s = rep(NA, R)
  U = matrix(NA, nrow=dim(data)[1], ncol=R)
  V = matrix(NA, nrow=dim(data)[2], ncol=R)
  W = matrix(NA, nrow=dim(data)[3], ncol=R)
  diag.list = vector("list", R)
  
  dataold = data
  
  for (r in seq_len(R)) {
    #cat("\nr=", r, "\n")
    res = RobRSVD_tnsr(data=dataold, irobust=irobust, rtype=rtype, istablize=istablize, inits=inits, initu=initu, initv=initv, initw=initw, 
                       lambdaRange = lambdaRange, nlambda=nlambda, degree=degree, iugcv=iugcv, ivgcv=ivgcv, iwgcv=iwgcv, isspline=isspline, 
                       verbose=verbose, niter=niter, tol=tol)
    
    s[r] = res$s
    U[,r] = res$u[,1]
    V[,r] = res$v[,1]
    W[,r] = res$w[,1]
    diag.list[[r]] = res$diagout
    dataold = dataold - res$s * res$u[,1] %o% res$v[,1] %o% res$w[,1]
    
    #### initialization using CP ####
    tmp.data = as.tensor(dataold)
    cpold = cp(tmp.data, num_components = 1, max_iter=100)
    #cp_s = c(cpold$lambdas)*normVector(cpold$U[[1]])*normVector(cpold$U[[2]])*normVector(cpold$U[[3]])
    cp_u1 = c(cpold$U[[1]])/normVector(cpold$U[[1]])
    cp_v1 = c(cpold$U[[2]])/normVector(cpold$U[[2]])
    cp_w1 = c(cpold$U[[3]])/normVector(cpold$U[[3]])
  
    inits = res$s
    initu = cp_u1
    initv = cp_v1
    initw = cp_w1
  }
  
  return(list(s=s, U=U, V=V, W=W, GCV_results=diag.list))
}

### ROBCP - additive penalty version ###
RobRSVD2_tnsr <- function(data, irobust=FALSE, rtype="huber", huberk=1.345, bisqk=3.88, inits, initu, initv, initw, lambdaRange, nlambda, degree=1,
                          niter=1000, tol=1e-5, istablize=FALSE, uspar=0, vspar=0, wspar=0, 
                          iugcv=FALSE, ivgcv=FALSE, iwgcv=FALSE,
                          iloguspar=TRUE, ilogvspar=TRUE, ilogwspar=TRUE, isspline=FALSE, verbose=FALSE) {
  
  # Robust 3-way functional prinicpal component analysis
  # data is 3-mode tensor (m x n x p)
  
  ugcvmat = vgcvmat = wgcvmat = as.numeric()
  
  # check if the parameters are appropriate #
  size_data = c(dim(data))
  m = size_data[1]
  n = size_data[2]
  p = size_data[3]
  if (istablize) {
    # The original data is centered by the median of the whole tensor
    myscale = 1.4785*median(abs(c(data)-median(data)))
    localdata = data / myscale
  } else {
    myscale = 1
    localdata = data
  }
  
  uold = initu
  vold = initv
  wold = initw
  sold = inits
  if (istablize) {
    sold = sold / myscale
  } 
  uold = sold * uold
  
  # use the initial estiamte of u, v, w vectors to reconstruct data tensor 
  Appold = uold %o% vold %o% wold 
  Rmat = localdata - Appold
  Rvec = c(Rmat)
  mysigma = median(abs(Rvec))/0.675
  iter = 1
  localdiff = 9999
  diffvec = as.numeric()
  nuspar = nlambda$u
  nvspar = nlambda$v
  nwspar = nlambda$w
  uspar_current = uspar
  vspar_current = vspar
  wspar_current = wspar
  if (isspline) { # spline smoother 
    Omega.u = ssmatls(m)$y
    Omega.v = ssmatls(n)$y
    Omega.w = ssmatls(p)$y
  } else { # difference smoother
    Omega.u = crossprod(makeDiffOper(degree=degree, dim=m))
    Omega.v = crossprod(makeDiffOper(degree=degree, dim=n))
    Omega.w = crossprod(makeDiffOper(degree=degree, dim=p))
  }
  
  ugcvscore = vgcvscore = wgcvscore = as.numeric()
  
  # main iteration
  while (localdiff > tol & iter < niter) {
    if (irobust) {
      if (rtype == "huber") {
        Wmat = huberWeightLS(Rmat/mysigma, huberk) # Huber's weight
      } else {
        Wmat = bisqWeightLS(Rmat/mysigma, bisqk) # bisquare weight
      }
    } else {
      Wmat = array(1, dim=c(m,n,p))
    }
    
    # update u given (v, w)
    if (!iugcv) {
      #S.u.inv = diag(m) + uspar*Omega.u
      #S.v.inv = diag(n) + vspar*Omega.v
      #S.w.inv = diag(p) + wspar*Omega.w
      
      au = c(vold %o% wold) 
      Wmat.rtns = as.tensor(Wmat) # rTensor object
      #Omega.u.vw = (c(t(vold) %*% S.v.inv %*% vold)*c(t(wold) %*% S.w.inv %*% wold)*S.u.inv - diag(sum(vold^2),m)*diag(sum(wold^2),m))*(2*mysigma^2)
      Omega.u.vw = uspar*Omega.u*(2*mysigma^2)
      tAuWuAu = diag(colSums(diag(c(au^2)) %*%  t(rs_unfold(Wmat.rtns, m=1)@data))) 
      tAuWuYu = (rs_unfold(Wmat.rtns, m=1)@data * rs_unfold(as.tensor(localdata), m=1)@data) %*% au
      uterm1 = tAuWuAu + Omega.u.vw
      uterm2 = tAuWuYu
      unew = solve(uterm1, uterm2)
    } else {
      # choosing the penalty parameter by cross-validation according to GCV criteria
      if (nuspar<0)
        stop("number of smoothing parameter can not be negative")
      else {
        if (iloguspar) {
          usparvec = 10^seq(log10(lambdaRange$u[1]), log10(lambdaRange$u[2]), length.out=nuspar)
        } else {
          usparvec = seq(lambdaRange$u[1], lambdaRange$u[2], length.out=nuspar)
        } 
        ugcvvec = as.numeric()
        ugcvmat = as.numeric()
        
        au = c(vold %o% wold) 
        Wmat.rtns = as.tensor(Wmat)
        tAuWuAu = diag(colSums(diag(c(au^2)) %*%  t(rs_unfold(Wmat.rtns, m=1)@data))) 
        tAuWuYu = (rs_unfold(Wmat.rtns, m=1)@data * rs_unfold(as.tensor(localdata), m=1)@data) %*% au
        u_nsrobust = solve(tAuWuAu, tAuWuYu)
        
        for (iter_uspar in 1:nuspar) {
          #S.u.inv = diag(m) + usparvec[iter_uspar]*Omega.u
          #S.v.inv = diag(n) + vspar_current*Omega.v
          #S.w.inv = diag(p) + wspar_current*Omega.w
          #Omega.u.vw = (c(t(vold) %*% S.v.inv %*% vold)*c(t(wold) %*% S.w.inv %*% wold)*S.u.inv - diag(sum(vold^2),m)*diag(sum(wold^2),m))*(2*mysigma^2)
          Omega.u.vw = usparvec[iter_uspar]*Omega.u*(2*mysigma^2)
          usterm1 = tAuWuAu + Omega.u.vw
          usterm2 = tAuWuYu
          u_srobust = solve(usterm1, usterm2)
          smooth_u = solve(usterm1, tAuWuAu)
          gcv_ut = m*sum((u_nsrobust-u_srobust)^2)/(m - sum(diag(smooth_u)))^2
          ugcvvec[iter_uspar] = gcv_ut
          ugcvmat = cbind(ugcvmat, u_srobust/(sqrt(sum(u_srobust^2))))
        }
        uspar_current = usparvec[which.min(ugcvvec)]
        ugcvscore = cbind(usparvec, ugcvvec)
        
        #S.u.inv = diag(m) + uspar_current*Omega.u
        #Omega.u.vw = (c(t(vold) %*% S.v.inv %*% vold)*c(t(wold) %*% S.w.inv %*% wold)*S.u.inv - diag(sum(vold^2),m)*diag(sum(wold^2),m))*(2*mysigma^2)
        Omega.u.vw = uspar_current*Omega.u*(2*mysigma^2)
        #tAuWuAu = diag(colSums(diag(c(au^2)) %*%  t(rs_unfold(Wmat.rtns, m=1)@data))) 
        #tAuWuYu = (rs_unfold(Wmat.rtns, m=1)@data * rs_unfold(as.tensor(localdata), m=1)@data) %*% au
        uterm1 = tAuWuAu + Omega.u.vw
        uterm2 = tAuWuYu
        unew = solve(uterm1, uterm2)
      }
    }
    
    # update v given (u, w)  
    if (!ivgcv) {
      #S.u.inv = diag(m) + uspar*Omega.u
      #S.v.inv = diag(n) + vspar*Omega.v
      #S.w.inv = diag(p) + wspar*Omega.w
      
      av = c(unew %o% wold) 
      Wmat.rtns = as.tensor(Wmat) # rTensor object
      #Omega.v.uw = (c(t(unew) %*% S.u.inv %*% unew)*c(t(wold) %*% S.w.inv %*% wold)*S.v.inv - diag(sum(unew^2),n)*diag(sum(wold^2),n))*(2*mysigma^2)
      Omega.v.uw = vspar*Omega.v*(2*mysigma^2)
      tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
      tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
      vterm1 = tAvWvAv + Omega.v.uw
      vterm2 = tAvWvYv
      vnew = solve(vterm1, vterm2)
    } else {
      # choosing the penalty parameter by cross-validation according to GCV criteria
      if (nvspar<0)
        stop("number of smoothing parameter can not be negative")
      else {
        if (ilogvspar) {
          vsparvec = 10^seq(log10(lambdaRange$v[1]), log10(lambdaRange$v[2]), length.out=nvspar)
        } else {
          vsparvec = seq(lambdaRange$v[1], lambdaRange$v[2], length.out=nvspar)
        } 
        vgcvvec = as.numeric()
        vgcvmat = as.numeric()
        
        av = c(unew %o% wold) 
        Wmat.rtns = as.tensor(Wmat)
        tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
        tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
        v_nsrobust = solve(tAvWvAv, tAvWvYv)
        
        for (iter_vspar in 1:nvspar) {
          #S.u.inv = diag(m) + uspar_current*Omega.u
          #S.v.inv = diag(n) + vsparvec[iter_vspar]*Omega.v
          #S.w.inv = diag(p) + wspar_current*Omega.w
          #Omega.v.uw = (c(t(unew) %*% S.u.inv %*% unew)*c(t(wold) %*% S.w.inv %*% wold)*S.v.inv - diag(sum(unew^2),n)*diag(sum(wold^2),n))*(2*mysigma^2)
          Omega.v.uw = vsparvec[iter_vspar]*Omega.v*(2*mysigma^2)
          vsterm1 = tAvWvAv + Omega.v.uw
          vsterm2 = tAvWvYv
          v_srobust = solve(vsterm1, vsterm2)
          smooth_v = solve(vsterm1, tAvWvAv)
          gcv_vt = n*sum((v_nsrobust-v_srobust)^2)/(n - sum(diag(smooth_v)))^2
          vgcvvec[iter_vspar] = gcv_vt
          vgcvmat = cbind(vgcvmat, v_srobust/(sqrt(sum(v_srobust^2))))
        }
        vspar_current = vsparvec[which.min(vgcvvec)]
        vgcvscore = cbind(vsparvec, vgcvvec)
        
        #S.v.inv = diag(n) + vspar_current*Omega.v
        #Omega.v.uw = (c(t(unew) %*% S.u.inv %*% unew)*c(t(wold) %*% S.w.inv %*% wold)*S.v.inv - diag(sum(unew^2),n)*diag(sum(wold^2),n))*(2*mysigma^2)
        Omega.v.uw = vspar_current*Omega.v*(2*mysigma^2)
        #tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
        #tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
        vterm1 = tAvWvAv + Omega.v.uw
        vterm2 = tAvWvYv
        vnew = solve(vterm1, vterm2)
      }
    }
    
    # update w given (u, v)  
    if (!iwgcv) {
      #S.u.inv = diag(m) + uspar*Omega.u
      #S.v.inv = diag(n) + vspar*Omega.v
      #S.w.inv = diag(p) + wspar*Omega.w
      
      aw = c(unew %o% vnew) 
      Wmat.rtns = as.tensor(Wmat) # rTensor object
      #Omega.w.uv = (c(t(unew) %*% S.u.inv %*% unew)*c(t(vnew) %*% S.v.inv %*% vnew)*S.w.inv - diag(sum(unew^2),p)*diag(sum(vnew^2),p))*(2*mysigma^2)
      Omega.w.uv = wspar*Omega.w*(2*mysigma^2)
      tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
      tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
      wterm1 = tAwWwAw + Omega.w.uv
      wterm2 = tAwWwYw
      wnew = solve(wterm1, wterm2)
    } else {
      # choosing the penalty parameter by cross-validation according to GCV criteria
      if (nwspar<0)
        stop("number of smoothing parameter can not be negative")
      else {
        if (ilogwspar) {
          wsparvec = 10^seq(log10(lambdaRange$w[1]), log10(lambdaRange$w[2]), length.out=nwspar)
        } else {
          wsparvec = seq(lambdaRange$w[1], lambdaRange$w[2], length.out=nwspar)
        } 
        wgcvvec = as.numeric()
        wgcvmat = as.numeric()
        
        aw = c(unew %o% vnew) 
        Wmat.rtns = as.tensor(Wmat)
        tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
        tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
        w_nsrobust = solve(tAwWwAw, tAwWwYw)
        
        for (iter_wspar in 1:nwspar) {
          #S.u.inv = diag(m) + uspar_current*Omega.u
          #S.v.inv = diag(n) + vspar_current*Omega.v
          #S.w.inv = diag(p) + wsparvec[iter_wspar]*Omega.w
          #Omega.w.uv = (c(t(unew) %*% S.u.inv %*% unew)*c(t(vnew) %*% S.v.inv %*% vnew)*S.w.inv - diag(sum(unew^2),p)*diag(sum(vnew^2),p))*(2*mysigma^2)
          Omega.w.uv = wsparvec[iter_wspar]*Omega.w*(2*mysigma^2)
          wsterm1 = tAwWwAw + Omega.w.uv
          wsterm2 = tAwWwYw
          w_srobust = solve(wsterm1, wsterm2)
          smooth_w = solve(wsterm1, tAwWwAw)
          gcv_wt = p*sum((w_nsrobust-w_srobust)^2)/(p - sum(diag(smooth_w)))^2
          wgcvvec[iter_wspar] = gcv_wt
          wgcvmat = cbind(wgcvmat, w_srobust/(sqrt(sum(w_srobust^2))))
        }
        wspar_current = wsparvec[which.min(wgcvvec)]
        wgcvscore = cbind(wsparvec, wgcvvec)
        
        #S.w.inv = diag(p) + wspar_current*Omega.w
        #Omega.w.uv = (c(t(unew) %*% S.u.inv %*% unew)*c(t(vnew) %*% S.v.inv %*% vnew)*S.w.inv - diag(sum(unew^2),p)*diag(sum(vnew^2),p))*(2*mysigma^2)
        Omega.w.uv = (wspar_current*Omega.w)*(2*mysigma^2)
        #tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
        #tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
        wterm1 = tAwWwAw + Omega.w.uv
        wterm2 = tAwWwYw
        wnew = solve(wterm1, wterm2)
      }
    }  
    
    Appnew = c(unew) %o% c(vnew) %o% c(wnew)
    Rmat = localdata - Appnew
    
    localdiff = max(abs(Appnew-Appold))
    Appold = Appnew
    uold = sqrt(sum(vnew^2))* sqrt(sum(wnew^2))*unew
    vold = vnew / sqrt(sum(vnew^2)) # v is normalized to length 1
    wold = wnew / sqrt(sum(wnew^2)) # w is normalized to length 1
    if (verbose) {
      print(paste("iter", iter, "fin"))
      print(paste("diff:",round(localdiff,4)))
    }
    
    iter = iter + 1
    diffvec = c(diffvec, localdiff)
    
  } # end of while
  
  v = vold
  w = wold
  s = myscale * sqrt(sum(uold^2))
  u = uold / sqrt(sum(uold^2))
  
  if (iugcv) {
    uspar = uspar_current
  }
  if (ivgcv) {
    vspar = vspar_current
  }
  if (iwgcv) {
    wspar = wspar_current
  }
  
  diagout = list(ugcvsore=ugcvscore, vgcvscore=vgcvscore, wgcvscore=wgcvscore, 
                 ugcvmat=ugcvmat, vgcvmat=vgcvmat, wgcvmat=wgcvmat, diffvec=diffvec)
  return(list(s=s, u=u, v=v, w=w, diagout=diagout))
}

wrap_RobRSVD2_tnsr <- function(data, irobust=FALSE, rtype="huber", huberk=1.345, bisqk=3.88, R, inits, initu, initv, initw, lambdaRange, nlambda, degree=1,
                              niter=1000, tol=1e-5, istablize=FALSE, uspar=0, vspar=0, wspar=0, 
                              iugcv=FALSE, ivgcv=FALSE, iwgcv=FALSE,
                              iloguspar=TRUE, ilogvspar=TRUE, ilogwspar=TRUE, verbose=FALSE, isspline=FALSE) {
  
  s = rep(NA, R)
  U = matrix(NA, nrow=dim(data)[1], ncol=R)
  V = matrix(NA, nrow=dim(data)[2], ncol=R)
  W = matrix(NA, nrow=dim(data)[3], ncol=R)
  diag.list = vector("list", R)
  
  dataold = data
  
  for (r in seq_len(R)) {
    #cat("\nr=", r, "\n")
    res = RobRSVD2_tnsr(data=dataold, irobust=irobust, rtype=rtype, inits=inits, initu=initu, initv=initv, initw=initw, lambdaRange = lambdaRange,
                       nlambda=nlambda, degree=degree, iugcv=iugcv, ivgcv=ivgcv, iwgcv=iwgcv, isspline=isspline, verbose=verbose, niter=niter, tol=tol)
    
    s[r] = res$s
    U[,r] = res$u[,1]
    V[,r] = res$v[,1]
    W[,r] = res$w[,1]
    diag.list[[r]] = res$diagout
    dataold = dataold - res$s * res$u[,1] %o% res$v[,1] %o% res$w[,1]
    
    #### initialization using CP ####
    tmp.data = as.tensor(dataold)
    cpold = cp(tmp.data, num_components = 1, max_iter=100)
    #cp_s = c(cpold$lambdas)*normVector(cpold$U[[1]])*normVector(cpold$U[[2]])*normVector(cpold$U[[3]])
    cp_u1 = c(cpold$U[[1]])/normVector(cpold$U[[1]])
    cp_v1 = c(cpold$U[[2]])/normVector(cpold$U[[2]])
    cp_w1 = c(cpold$U[[3]])/normVector(cpold$U[[3]])
  
    inits = res$s
    initu = cp_u1
    initv = cp_v1
    initw = cp_w1
  }
  
  return(list(s=s, U=U, V=V, W=W, GCV_results=diag.list))
}

### ROBCP-s - additive penalty version with smoothness constraint ###
RobRSVD2_sm_tnsr <- function(data, irobust=FALSE, rtype="huber", huberk=1.345, bisqk=3.88, inits, initu, initv, initw, lambdaRange, nlambda, degree=1,
                             niter=1000, tol=1e-5, istablize=FALSE, uspar=0.001, vspar=0.001, wspar=0.001, basistype="bspline",
                             iugcv=FALSE, ivgcv=FALSE, iwgcv=FALSE, nknots, norder=3, ismoothu=TRUE, ismoothv=TRUE, ismoothw=TRUE,
                             iloguspar=TRUE, ilogvspar=TRUE, ilogwspar=TRUE, isspline=FALSE, verbose=FALSE) {
  
  # Robust 3-way functional prinicpal component analysis
  # data is 3-mode tensor (m x n x p)
  
  ugcvmat = vgcvmat = wgcvmat = as.numeric()
  
  # check if the parameters are appropriate #
  size_data = c(dim(data))
  m = size_data[1]
  n = size_data[2]
  p = size_data[3]
  if (istablize) {
    # The original data is centered by the median of the whole tensor
    myscale = 1.4785*median(abs(c(data)-median(data)))
    localdata = data / myscale
  } else {
    myscale = 1
    localdata = data
  }
  
  # B-spline smoother matrix #
  ugrid = seq(0, 1, length.out=m)
  vgrid = seq(0, 1, length.out=n)
  wgrid = seq(0, 1, length.out=p)
  
  num.dfu = nknots$u + norder # cubic spline when norder = 3
  num.dfv = nknots$v + norder # cubic spline when norder = 3
  num.dfw = nknots$w + norder # cubic spline when norder = 3
  #Bu = splines::bs(x=ugrid, df=num.dfu) # m x (deg of u)
  #Bv = splines::bs(x=vgrid, df=num.dfv) # n x (deg of v)
  #Bw = splines::bs(x=wgrid, df=num.dfw) # p x (deg of w)
  
  fdabasis_u = fda::create.bspline.basis(rangeval = range(ugrid), nbasis = num.dfu)
  Bu = eval.basis(ugrid, fdabasis_u)
  
  fdabasis_v = create.bspline.basis(rangeval = range(vgrid), nbasis = num.dfv)
  Bv = eval.basis(vgrid, fdabasis_v)
  
  fdabasis_w = create.bspline.basis(rangeval = range(wgrid), nbasis = num.dfw)
  Bw = eval.basis(wgrid, fdabasis_w)
  
  if (basistype == "fourier") {
    fdabasis_u = create.fourier.basis(rangeval = range(ugrid), nbasis = num.dfu)
    Bu = eval.basis(ugrid, fdabasis_u)
    
    fdabasis_v = create.fourier.basis(rangeval = range(vgrid), nbasis = num.dfv)
    Bv = eval.basis(vgrid, fdabasis_v)
    
    fdabasis_w = create.fourier.basis(rangeval = range(wgrid), nbasis = num.dfw)
    Bw = eval.basis(wgrid, fdabasis_w)
  } 
  
  Du = (t(Bu) %*% Bu) * (1/(m-1))
  Dv = (t(Bv) %*% Bv) * (1/(n-1))
  Dw = (t(Bw) %*% Bw) * (1/(p-1))
  
  #phiu_old = solve(t(Bu) %*% Bu, t(Bu) %*% initu)
  #phiv_old = solve(t(Bv) %*% Bv, t(Bv) %*% initv)
  #phiw_old = solve(t(Bw) %*% Bw, t(Bw) %*% initw)
  
  uold = initu
  vold = initv
  wold = initw
  sold = inits
  if (istablize) {
    sold = sold / myscale
  } 
  uold = sold * uold
  
  # use the initial estiamte of u, v, w vectors to reconstruct data tensor 
  Appold = uold %o% vold %o% wold 
  Rmat = localdata - Appold
  Rvec = c(Rmat)
  mysigma = median(abs(Rvec))/0.675
  iter = 1
  localdiff = 9999
  diffvec = as.numeric()
  nuspar = nlambda$u
  nvspar = nlambda$v
  nwspar = nlambda$w
  uspar_current = uspar
  vspar_current = vspar
  wspar_current = wspar
  
  if (isspline) { # spline smoother 
    Omega.u = Du
    Omega.v = Dv
    Omega.w = Dw
  } else { # difference smoother
    Omega.u = crossprod(makeDiffOper(degree=degree, dim=num.dfu))
    Omega.v = crossprod(makeDiffOper(degree=degree, dim=num.dfv))
    Omega.w = crossprod(makeDiffOper(degree=degree, dim=num.dfw))
  } 
  
  if (!ismoothu) {
    Omega.u = diag(1, m)
  } 
  
  if (!ismoothv) {
    Omega.v = diag(1, n)
  } 
  
  if (!ismoothw) {
    Omega.w = diag(1, p)
  } 
  
  
  ugcvscore = vgcvscore = wgcvscore = as.numeric()
  
  # main iteration
  while (localdiff > tol & iter < niter) {
    if (irobust) {
      if (rtype == "huber") {
        Wmat = huberWeightLS(Rmat/mysigma, huberk)
      } else {
        Wmat = bisqWeightLS(Rmat/mysigma, bisqk)
      }
    } else {
      Wmat = array(1, dim=c(m,n,p))
    }
    
    # update u given (v, w)
    if (!iugcv) {
      #S.u.inv = diag(num.dfu) + uspar*Omega.u
      #S.v.inv = diag(num.dfv) + vspar*Omega.v
      #S.w.inv = diag(num.dfw) + wspar*Omega.w
      
      if (!ismoothu) {
        Bu = diag(1, m)
      }
      
      au = c(vold %o% wold) 
      Wmat.rtns = as.tensor(Wmat) # rTensor object
      Omega.u.vw = uspar*Omega.u*(2*mysigma^2)
      tAuWuAu = diag(colSums(diag(c(au^2)) %*%  t(rs_unfold(Wmat.rtns, m=1)@data))) 
      tAuWuYu = (rs_unfold(Wmat.rtns, m=1)@data * rs_unfold(as.tensor(localdata), m=1)@data) %*% au
      phiu_term1 = t(Bu) %*% tAuWuAu %*% Bu + Omega.u.vw
      phiu_term2 = t(Bu) %*% tAuWuYu
      phiu_new = solve(phiu_term1, phiu_term2)
      unew = Bu %*% phiu_new
    } else {
      # choosing the penalty parameter by cross-validation according to GCV criteria
      if (nuspar<0)
        stop("number of smoothing parameter can not be negative")
      else {
        if (iloguspar) {
          usparvec = 10^seq(log10(lambdaRange$u[1]), log10(lambdaRange$u[2]), length.out=nuspar)
        } else {
          usparvec = seq(lambdaRange$u[1], lambdaRange$u[2], length.out=nuspar)
        } 
        ugcvvec = as.numeric()
        ugcvmat = as.numeric()
        
        if (!ismoothu) {
          Bu = diag(1, m)
        }
        
        au = c(vold %o% wold) 
        Wmat.rtns = as.tensor(Wmat)
        tAuWuAu = diag(colSums(diag(c(au^2)) %*%  t(rs_unfold(Wmat.rtns, m=1)@data))) 
        tAuWuYu = (rs_unfold(Wmat.rtns, m=1)@data * rs_unfold(as.tensor(localdata), m=1)@data) %*% au
        phiu_nsrobust = solve(t(Bu) %*% tAuWuAu %*% Bu, t(Bu) %*% tAuWuYu)
        
        for (iter_uspar in 1:nuspar) {
          
          #S.u.inv = diag(num.dfu) + usparvec[iter_uspar]*Omega.u
          #S.v.inv = diag(num.dfv) + vspar_current*Omega.v
          #S.w.inv = diag(num.dfw) + wspar_current*Omega.w

          #Omega.u.vw = (c(t(vold) %*% S.v.inv %*% vold)*c(t(wold) %*% S.w.inv %*% wold)*S.u.inv - diag(sum(vold^2),m)*diag(sum(wold^2),m))*(2*mysigma^2)
          Omega.u.vw = usparvec[iter_uspar]*Omega.u*(2*mysigma^2)
          phius_term1 = t(Bu) %*% tAuWuAu %*% Bu + Omega.u.vw
          phius_term2 = t(Bu) %*% tAuWuYu
          phiu_srobust = solve(phius_term1, phius_term2)
          smooth_phiu = solve(phius_term1, t(Bu) %*% tAuWuAu %*% Bu)
          gcv_phiut = num.dfu*sum((phiu_nsrobust-phiu_srobust)^2)/(num.dfu - sum(diag(smooth_phiu)))^2
          ugcvvec[iter_uspar] = gcv_phiut
          ugcvmat = cbind(ugcvmat, phiu_srobust/(sqrt(sum(phiu_srobust^2))))
        }
        uspar_current = usparvec[which.min(ugcvvec)]
        ugcvscore = cbind(usparvec, ugcvvec)
        
        Omega.u.vw = uspar_current*Omega.u*(2*mysigma^2)
        phiu_term1 = t(Bu) %*% tAuWuAu %*% Bu + Omega.u.vw
        phiu_term2 = t(Bu) %*% tAuWuYu
        phiu_new = solve(phiu_term1, phiu_term2)
        unew = Bu %*% phiu_new
      }
    }
    
    # update v given (u, w)  
    if (!ivgcv) {
      #S.u.inv = diag(num.dfu) + uspar*Omega.u
      #S.v.inv = diag(num.dfv) + vspar*Omega.v
      #S.w.inv = diag(num.dfw) + wspar*Omega.w
      
      if (!ismoothv) {
        Bv = diag(1, n)
      } 
      #print(dim(Bv))
      
      
      av = c(unew %o% wold) 
      Wmat.rtns = as.tensor(Wmat) # rTensor object
      Omega.v.uw = vspar*Omega.v*(2*mysigma^2)
      tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
      tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
      
      phiv_term1 = t(Bv) %*% tAvWvAv %*% Bv + Omega.v.uw
      phiv_term2 = t(Bv) %*% tAvWvYv
      
      phiv_new = solve(phiv_term1, phiv_term2)
      vnew = Bv %*% phiv_new
    } else {
      # choosing the penalty parameter by cross-validation according to GCV criteria
      if (nvspar<0)
        stop("number of smoothing parameter can not be negative")
      else {
        if (ilogvspar) {
          vsparvec = 10^seq(log10(lambdaRange$v[1]), log10(lambdaRange$v[2]), length.out=nvspar)
        } else {
          vsparvec = seq(lambdaRange$v[1], lambdaRange$v[2], length.out=nvspar)
        } 
        vgcvvec = as.numeric()
        vgcvmat = as.numeric()
        
        if (!ismoothv) {
          Bv = diag(1, n)
        }
        
        av = c(unew %o% wold) 
        Wmat.rtns = as.tensor(Wmat)
        tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
        tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
        phiv_nsrobust = solve(t(Bv) %*% tAvWvAv %*% Bv, t(Bv) %*% tAvWvYv)
        
        for (iter_vspar in 1:nvspar) {
          #av = c(unew %o% wold) 
          #Wmat.rtns = as.tensor(Wmat)
          #tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
          #tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
          #phiv_nsrobust = solve(t(Bv) %*% tAvWvAv %*% Bv, t(Bv) %*% tAvWvYv)
          
          #S.u.inv = diag(m) + uspar_current*Omega.u
          #S.v.inv = diag(n) + vsparvec[iter_vspar]*Omega.v
          #S.w.inv = diag(p) + wspar_current*Omega.w
          #Omega.v.uw = (c(t(unew) %*% S.u.inv %*% unew)*c(t(wold) %*% S.w.inv %*% wold)*S.v.inv - diag(sum(unew^2),n)*diag(sum(wold^2),n))*(2*mysigma^2)
          Omega.v.uw = vsparvec[iter_vspar]*Omega.v*(2*mysigma^2)
          phivs_term1 = t(Bv) %*% tAvWvAv %*% Bv + Omega.v.uw
          phivs_term2 = t(Bv) %*% tAvWvYv
          phiv_srobust = solve(phivs_term1, phivs_term2)
          smooth_phiv = solve(phivs_term1, t(Bv) %*% tAvWvAv %*% Bv)
          gcv_phivt = num.dfv*sum((phiv_nsrobust-phiv_srobust)^2)/(num.dfv - sum(diag(smooth_phiv)))^2
          vgcvvec[iter_vspar] = gcv_phivt
          vgcvmat = cbind(vgcvmat, phiv_srobust/(sqrt(sum(phiv_srobust^2))))
        }
        vspar_current = vsparvec[which.min(vgcvvec)]
        vgcvscore = cbind(vsparvec, vgcvvec)
        
        #S.v.inv = diag(n) + vspar_current*Omega.v
        #Omega.v.uw = (c(t(unew) %*% S.u.inv %*% unew)*c(t(wold) %*% S.w.inv %*% wold)*S.v.inv - diag(sum(unew^2),n)*diag(sum(wold^2),n))*(2*mysigma^2)
        Omega.v.uw = vspar_current*Omega.v*(2*mysigma^2)
        #tAvWvAv = diag(colSums(diag(c(av^2)) %*%  t(rs_unfold(Wmat.rtns, m=2)@data))) 
        #tAvWvYv = (rs_unfold(Wmat.rtns, m=2)@data * rs_unfold(as.tensor(localdata), m=2)@data) %*% av
        phiv_term1 = t(Bv) %*% tAvWvAv %*% Bv + Omega.v.uw
        phiv_term2 = t(Bv) %*% tAvWvYv
        phiv_new = solve(phiv_term1, phiv_term2)
        vnew = Bv %*% phiv_new
      }
    }
    
    # update w given (u, v)  
    if (!iwgcv) {
      #S.u.inv = diag(num.dfu) + uspar*Omega.u
      #S.v.inv = diag(num.dfv) + vspar*Omega.v
      #S.w.inv = diag(num.dfw) + wspar*Omega.w
      
      if (!ismoothw) {
        Bw = diag(1, p)
      }
      
      aw = c(unew %o% vnew) 
      Wmat.rtns = as.tensor(Wmat) # rTensor object
      #Omega.w.uv = (c(t(unew) %*% S.u.inv %*% unew)*c(t(vnew) %*% S.v.inv %*% vnew)*S.w.inv - diag(sum(unew^2),p)*diag(sum(vnew^2),p))*(2*mysigma^2)
      Omega.w.uv = wspar*Omega.w*(2*mysigma^2)
      tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
      tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
      phiw_term1 = t(Bw) %*% tAwWwAw %*% Bw + Omega.w.uv
      phiw_term2 = t(Bw) %*% tAwWwYw
      phiw_new = solve(phiw_term1, phiw_term2)
      wnew = Bw %*% phiw_new
    } else {
      # choosing the penalty parameter by cross-validation according to GCV criteria
      if (nwspar<0)
        stop("number of smoothing parameter can not be negative")
      else {
        if (ilogwspar) {
          wsparvec = 10^seq(log10(lambdaRange$w[1]), log10(lambdaRange$w[2]), length.out=nwspar)
        } else {
          wsparvec = seq(lambdaRange$w[1], lambdaRange$w[2], length.out=nwspar)
        } 
        wgcvvec = as.numeric()
        wgcvmat = as.numeric()
        
        if (!ismoothw) {
          Bw = diag(1, p)
        }
        
        aw = c(unew %o% vnew) 
        Wmat.rtns = as.tensor(Wmat)
        tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
        tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
        phiw_nsrobust = solve(t(Bw) %*% tAwWwAw %*% Bw, t(Bw) %*% tAwWwYw)
        
        for (iter_wspar in 1:nwspar) {
          #aw = c(unew %o% vnew) 
          #Wmat.rtns = as.tensor(Wmat)
          #tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
          #tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
          #phiw_nsrobust = solve(t(Bw) %*% tAwWwAw %*% Bw, t(Bw) %*% tAwWwYw)
          
          #S.u.inv = diag(m) + uspar_current*Omega.u
          #S.v.inv = diag(n) + vspar_current*Omega.v
          #S.w.inv = diag(p) + wsparvec[iter_wspar]*Omega.w
          #Omega.w.uv = (c(t(unew) %*% S.u.inv %*% unew)*c(t(vnew) %*% S.v.inv %*% vnew)*S.w.inv - diag(sum(unew^2),p)*diag(sum(vnew^2),p))*(2*mysigma^2)
          Omega.w.uv = wsparvec[iter_wspar]*Omega.w*(2*mysigma^2)
          phiws_term1 = t(Bw) %*% tAwWwAw %*% Bw + Omega.w.uv
          phiws_term2 = t(Bw) %*% tAwWwYw
          phiw_srobust = solve(phiws_term1, phiws_term2)
          smooth_phiw = solve(phiws_term1, t(Bw) %*% tAwWwAw %*% Bw)
          gcv_phiwt = num.dfw*sum((phiw_nsrobust-phiw_srobust)^2)/(num.dfw - sum(diag(smooth_phiw)))^2
          wgcvvec[iter_wspar] = gcv_phiwt
          wgcvmat = cbind(wgcvmat, phiw_srobust/(sqrt(sum(phiw_srobust^2))))
        }
        wspar_current = wsparvec[which.min(wgcvvec)]
        wgcvscore = cbind(wsparvec, wgcvvec)
        
        #S.w.inv = diag(p) + wspar_current*Omega.w
        #Omega.w.uv = (c(t(unew) %*% S.u.inv %*% unew)*c(t(vnew) %*% S.v.inv %*% vnew)*S.w.inv - diag(sum(unew^2),p)*diag(sum(vnew^2),p))*(2*mysigma^2)
        Omega.w.uv = (wspar_current*Omega.w)*(2*mysigma^2)
        #tAwWwAw = diag(colSums(diag(c(aw^2)) %*%  t(rs_unfold(Wmat.rtns, m=3)@data))) 
        #tAwWwYw = (rs_unfold(Wmat.rtns, m=3)@data * rs_unfold(as.tensor(localdata), m=3)@data) %*% aw
        phiw_term1 = t(Bw) %*% tAwWwAw %*% Bw + Omega.w.uv
        phiw_term2 = t(Bw) %*% tAwWwYw
        phiw_new = solve(phiw_term1, phiw_term2)
        wnew = Bw %*% phiw_new
      }
    }  
    
    Appnew = c(unew) %o% c(vnew) %o% c(wnew)
    Rmat = localdata - Appnew
    
    localdiff = max(abs(Appnew-Appold))
    Appold = Appnew
    uold = sqrt(sum(vnew^2))* sqrt(sum(wnew^2))*unew
    vold = vnew / sqrt(sum(vnew^2)) # v is normalized to length 1
    wold = wnew / sqrt(sum(wnew^2)) # w is normalized to length 1
    if (verbose) {
      print(paste("iter", iter, "fin"))
      print(paste("diff:",round(localdiff,4)))
    }
    
    iter = iter + 1
    diffvec = c(diffvec, localdiff)
    
  } # end of while
  
  v = vold
  w = wold
  s = myscale * sqrt(sum(uold^2))
  u = uold / sqrt(sum(uold^2))
  
  if (iugcv) {
    uspar = uspar_current
  }
  if (ivgcv) {
    vspar = vspar_current
  }
  if (iwgcv) {
    wspar = wspar_current
  }
  
  diagout = list(ugcvsore=ugcvscore, vgcvscore=vgcvscore, wgcvscore=wgcvscore, 
                 ugcvmat=ugcvmat, vgcvmat=vgcvmat, wgcvmat=wgcvmat, diffvec=diffvec)
  return(list(s=s, u=u, v=v, w=w, diagout=diagout))
}

wrap_RobRSVD2_sm_tnsr <- function(data, R=2, irobust=FALSE, rtype="huber", huberk=1.345, bisqk=3.88, inits, initu, initv, initw, lambdaRange, nlambda, degree=1,
                                  niter=1000, tol=1e-5, istablize=FALSE, uspar=0.001, vspar=0.001, wspar=0.001, basistype="bspline",
                                  iugcv=FALSE, ivgcv=FALSE, iwgcv=FALSE, nknots, norder=3, ismoothu=TRUE, ismoothv=TRUE, ismoothw=TRUE,
                                  iloguspar=TRUE, ilogvspar=TRUE, ilogwspar=TRUE, isspline=FALSE, verbose=FALSE) {
  
  s = rep(NA, R)
  U = matrix(NA, nrow=dim(data)[1], ncol=R)
  V = matrix(NA, nrow=dim(data)[2], ncol=R)
  W = matrix(NA, nrow=dim(data)[3], ncol=R)
  diag.list = vector("list", R)
  
  dataold = data
  
  for (r in seq_len(R)) {
    #cat("\nr=", r, "\n")
    res = RobRSVD2_sm_tnsr(dataold, irobust=irobust, rtype=rtype, inits=inits, initu=initu, initv=initv, initw=initw, lambdaRange=lambdaRange, nlambda=nlambda, 
                           niter=niter, istablize=istablize, degree=degree, uspar=uspar, vspar=vspar, wspar=wspar, basistype=basistype,
                           iugcv=iugcv, ivgcv=ivgcv, iwgcv=iwgcv, nknots=nknots, norder=norder, ismoothu=ismoothu, ismoothv=ismoothv, ismoothw=ismoothw,
                           iloguspar=iloguspar, ilogvspar=ilogvspar, ilogwspar=ilogwspar, isspline=isspline, verbose=verbose, tol=tol)
    
    s[r] = res$s
    U[,r] = res$u[,1]
    V[,r] = res$v[,1]
    W[,r] = res$w[,1]
    #diag.list[[r]] = res$diagout
    dataold = dataold - res$s * res$u[,1] %o% res$v[,1] %o% res$w[,1]
    
    #### initialization using CP ####
    tmp.data = as.tensor(dataold)
    cpold = cp(tmp.data, num_components = 1, max_iter=100)
    #cp_s = c(cpold$lambdas)*normVector(cpold$U[[1]])*normVector(cpold$U[[2]])*normVector(cpold$U[[3]])
    cp_u1 = c(cpold$U[[1]])/normVector(cpold$U[[1]])
    cp_v1 = c(cpold$U[[2]])/normVector(cpold$U[[2]])
    cp_w1 = c(cpold$U[[3]])/normVector(cpold$U[[3]])

    inits = res$s
    initu = cp_u1
    initv = cp_v1
    initw = cp_w1
  }
  
  return(list(s=s, U=U, V=V, W=W))
}

###############################################
###############################################
###############################################

