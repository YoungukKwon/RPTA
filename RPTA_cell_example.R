source("pkg_load.R") 
source("utils.R")
source("robrcp_functions.R")


#### (Main article) rank R=1, normal noise add, random cell outliers ####
m = 100
n = p = 64

ugrid = seq(0,1, length.out=m)
u1 = 6*(ugrid-1/2)^2 - 2; u1 = u1/normVector(u1)

vgrid = seq(0, 1, length.out=n)
v1 = 1.5*sin(3*pi*2*(vgrid-1/3)); v1 = v1/normVector(v1)

wgrid = seq(0, 1, length.out=p)
set.seed(927)
w1 = gaussprocess(nt=p, K=function(s,t) {exp(-81*(s-t)^2)})$xt; w1 = w1/normVector(w1)
s1 = 1000

### generation 
X.sig = s1 * u1 %o% v1 %o% w1 
set.seed(1)
X.noise = array(rnorm(m*n*p, sd=sqrt(1)), dim=c(m,n,p)) # noise level(sd) 바꿔가면서 확인해보기 (1.0, 1.4, 1.8, ... ,3,0)
X = X.sig + X.noise
X.out = X

# replace random cells with outlying cells
out.rate = 0.04 # 4%, 8%, 12%, ... 
ncell.out = round(m*n*p*out.rate)
indx.out = sample(1:(m*n*p), ncell.out, replace = FALSE)
X.out[indx.out] = X.out[indx.out] + runif(ncell.out, min=max(X.sig), max=1.5*max(X.sig)) #add uniform noise

data = X.out

################################################################

#### M1: CP decomposition (L2 loss, no penalty) ####
tmp.data = as.tensor(data)
cpold = cp(tmp.data, num_components = 1, max_iter=100)
#cp_s = c(cpold$lambdas)*normVector(cpold$U[[1]])*normVector(cpold$U[[2]])*normVector(cpold$U[[3]])
cp_u1 = c(cpold$U[[1]])/normVector(cpold$U[[1]]); cp_u1 = sign_adjust(u1, cp_u1)
cp_v1 = c(cpold$U[[2]])/normVector(cpold$U[[2]]); cp_v1 = sign_adjust(v1, cp_v1)
cp_w1 = c(cpold$U[[3]])/normVector(cpold$U[[3]]); cp_w1 = sign_adjust(w1, cp_w1)

### initial (u,v,w) and s
#set.seed(1110)
init_s = 1000
init_u = cp_u1
init_v = cp_v1
init_w = cp_w1

## Penalty type 1: Allen (2013), Huang et al. (2009) 
## Penalty type 2: Hong and Lian (2013) <- additive penalty


#### M2: PCP (L2 loss, P1) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
nlambda = list(u=8, v=8, w=8)
rcp.out = RobRSVD_tnsr(data=data, irobust=FALSE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange, degree=1,
                            nlambda=nlambda, iugcv=TRUE, ivgcv=TRUE, iwgcv=TRUE, isspline=FALSE, verbose = TRUE, niter=20)
rcp_u1 = sign_adjust(u1, rcp.out$u[,1])
rcp_v1 = sign_adjust(v1, rcp.out$v[,1])
rcp_w1 = sign_adjust(w1, rcp.out$w[,1])


#### M3: RobCP (Huber loss, no penalty) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
nlambda = list(u=8, v=8, w=8)
robcp.out = RobRSVD_tnsr(data=data, irobust=TRUE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange, degree=1,
                         nlambda=nlambda, isspline=FALSE, verbose = TRUE, niter=20)
robcp_u1 = sign_adjust(u1, robcp.out$u[,1])
robcp_v1 = sign_adjust(v1, robcp.out$v[,1])
robcp_w1 = sign_adjust(w1, robcp.out$w[,1])


#### M4: RobPCP1 (Huber loss, P1) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
nlambda = list(u=8, v=8, w=8)
robrcp1.out = RobRSVD_tnsr(data=data, irobust=TRUE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange,
                          nlambda=nlambda, degree=3, iugcv=TRUE, ivgcv=TRUE, iwgcv=TRUE, isspline=FALSE, verbose = TRUE, niter=20)
robrcp1_u1 = sign_adjust(u1, robrcp1.out$u[,1])
robrcp1_v1 = sign_adjust(v1, robrcp1.out$v[,1])
robrcp1_w1 = sign_adjust(w1, robrcp1.out$w[,1])


#### M5: RobRCP2 (Huber loss, P2) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
nlambda = list(u=8, v=8, w=8)
robrcp2.out = RobRSVD2_tnsr(data=data, irobust=TRUE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange, degree=1,
                            nlambda=nlambda, iugcv=TRUE, ivgcv=TRUE, iwgcv=TRUE, isspline=FALSE, verbose = TRUE, niter=20)
robrcp2_u1 = sign_adjust(u1, robrcp2.out$u[,1])
robrcp2_v1 = sign_adjust(v1, robrcp2.out$v[,1])
robrcp2_w1 = sign_adjust(w1, robrcp2.out$w[,1])


#### M6: RobCCP (Huber loss, P2 with basis constraint) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
nlambda = list(u=8, v=8, w=8)
nknots = list(u=round(m*0.4), v=round(n*0.4), w=round(p*0.4))
robccp.out = RobRSVD2_sm_tnsr(data=data, irobust=TRUE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange, nknots=nknots, degree=3,
                                 nlambda=nlambda, iugcv=TRUE, ivgcv=TRUE, iwgcv=TRUE, uspar=0.001, vspar=0.001, wspar=0.001, isspline=FALSE, verbose = TRUE, niter=20)
robccp_u1 = sign_adjust(u1, robccp.out$u[,1])
robccp_v1 = sign_adjust(v1, robccp.out$v[,1])
robccp_w1 = sign_adjust(w1, robccp.out$w[,1])


### tensor reconstruction
cp.recon = cpold$est@data
rcp.recon =  rcp.out$s * rcp_u1 %o% rcp_v1 %o% rcp_w1
robcp.recon =  robcp.out$s * robcp_u1 %o% robcp_v1 %o% robcp_w1
robrcp1.recon =  robrcp1.out$s * robrcp1_u1 %o% robrcp1_v1 %o% robrcp1_w1
robrcp2.recon = robrcp2.out$s * robrcp2_u1 %o% robrcp2_v1 %o% robrcp2_w1
robccp.recon = robccp.out$s * robccp_u1 %o% robccp_v1 %o% robccp_w1
#robrcp_basis.recon = robrcp2_basis.out$s * robrcp2_basis_u1 %o% robrcp2_basis_v1 %o% robrcp2_basis_w1

### difference (RMSE, plot)

### singular value
data.frame(
  cp=MFPCA::ttv(cp.recon, v=list(cp_u1, cp_v1, cp_w1), dim=1:3), 
  rcp=rcp.out$s, 
  robcp=robcp.out$s,
  robrcp1=robrcp1.out$s,
  robrcp2=robrcp2.out$s, 
  robccp=robccp.out$s
) %>% signif(3)

# Comparison: u
data.frame(
  cp=sum((u1-cp_u1)^2), 
  rcp=sum((u1-rcp_u1)^2),
  robcp=sum((u1-robcp_u1)^2),
  robrcp1=sum((u1-robrcp1_u1)^2),
  robrcp2=sum((u1-robrcp2_u1)^2), 
  robccp=sum((u1-robccp_u1)^2)
) %>% signif(3)

par(mfrow=c(2,3), mar=c(2,2,2,2))
compare_plot(ugrid, u1, cp_u1, title="u1-CP")
compare_plot(ugrid, u1, rcp_u1, title="u1-PCP")
compare_plot(ugrid, u1, robcp_u1, title="u1-RobCP")
compare_plot(ugrid, u1, robrcp1_u1, title="u1-RobPCP1")
compare_plot(ugrid, u1, robrcp2_u1, title="u1-RobPCP2")
compare_plot(ugrid, u1, robccp_u1, title="u1-RobPCP2-sp")

# Comparison: v
data.frame(
  cp=sum((v1-cp_v1)^2), 
  rcp=sum((v1-rcp_v1)^2),
  robcp=sum((v1-robcp_v1)^2),
  robrcp1=sum((v1-robrcp1_v1)^2),
  robrcp2=sum((v1-robrcp2_v1)^2), 
  robccp=sum((v1-robccp_v1)^2)
) %>% signif(3)

par(mfrow=c(2,3), mar=c(2,2,2,2))
compare_plot(vgrid, v1, cp_v1, title="v1-CP")
compare_plot(vgrid, v1, rcp_v1, title="v1-PCP")
compare_plot(vgrid, v1, robcp_v1, title="v1-RobCP")
compare_plot(vgrid, v1, robrcp1_v1, title="v1-RobPCP1")
compare_plot(vgrid, v1, robrcp2_v1, title="v1-RobPCP2")
compare_plot(vgrid, v1, robccp_v1, title="v1-RobPCP2-sp")

# Comparison: w
data.frame(
  cp=sum((w1-cp_w1)^2), 
  rcp=sum((w1-rcp_w1)^2),
  robcp=sum((w1-robcp_w1)^2),
  robrcp1=sum((w1-robrcp1_w1)^2),
  robrcp2=sum((w1-robrcp2_w1)^2), 
  robccp=sum((w1-robccp_w1)^2)
) %>% signif(3)

par(mfrow=c(2,3), mar=c(2,2,2,2))
compare_plot(wgrid, w1, cp_w1, title="w1-CP")
compare_plot(wgrid, w1, rcp_w1, title="w1-PCP")
compare_plot(wgrid, w1, robcp_w1, title="w1-RobCP")
compare_plot(wgrid, w1, robrcp1_w1, title="w1-RobPCP1")
compare_plot(wgrid, w1, robrcp2_w1, title="w1-RobPCP2")
compare_plot(wgrid, w1, robccp_w1, title="w1-RobPCP2-sp")


# reconstructed tensor
data.frame(
  cp=mean(c((cp.recon - X.sig)^2)),
  rcp=mean(c((rcp.recon - X.sig)^2)),
  robcp=mean(c((robcp.recon - X.sig)^2)),
  robrcp1=mean(c((robrcp1.recon - X.sig)^2)),
  robrcp2=mean(c((robrcp2.recon - X.sig)^2)), 
  robccp=mean(c((robccp.recon - X.sig)^2))
) %>% signif(3)



############################################################
############################################################
############################################################


#### (Appendix) rank R=2, normal noise add, random cell outliers ####
m = 100
n = p = 64

ugrid = seq(0,1, length.out=m)
u1 = 6*(ugrid-0.25)^2 - 2 ; u1 = u1/normVector(u1)
u2 = -6*(ugrid-0.75)^2 - 2 ; u2 = u2/normVector(u2)

vgrid = seq(0, 1, length.out=n)
v1 = 1.5*sin(2*pi*2*(vgrid-1/3)); v1 = v1/normVector(v1)
v2 = 1.5*sin(2*pi*2*(vgrid-2/3)); v2 = v2/normVector(v2)

wgrid = seq(0, 1, length.out=p)
set.seed(927)
w1 = gaussprocess(nt=p, K=function(s,t) {exp(-81*(s-t)^2)})$xt; w1 = w1/normVector(w1)
w2 = gaussprocess(nt=p, K=function(s,t) {exp(-25*(s-t)^2)})$xt; w2 = w2/normVector(w2)

U = cbind(u1, u2)
gram_u = pracma::gramSchmidt(U)
u1 = gram_u$Q[,1]; u2 = gram_u$Q[,2]

V = cbind(v1, v2)
gram_v = pracma::gramSchmidt(V)
v1 = gram_v$Q[,1]; v2 = gram_v$Q[,2]

W = cbind(w1, w2)
gram_w = pracma::gramSchmidt(W)
w1 = gram_w$Q[,1]; w2 = gram_w$Q[,2]

s1 = 1000
s2 = 500

### generation 
X.sig = s1 * u1 %o% v1 %o% w1 + s2 * u2 %o% v2 %o% w2
set.seed(1)
X.noise = array(rnorm(m*n*p, sd=sqrt(1.5)), dim=c(m,n,p)) # noise level(sd) 바꿔가면서 확인해보기 (1.0, 1.4, 1.8, ... ,3,0)
X = X.sig + X.noise
X.out = X

# replace random cells with outlying cells
out.rate = 0.04 # 4%, 8%, 12%, ... 
ncell.out = round(m*n*p*out.rate)
indx.out = sample(1:(m*n*p), ncell.out, replace = FALSE)
X.out[indx.out] = X.out[indx.out] + runif(ncell.out, min=max(X.sig), max=1.5*max(X.sig)) #add uniform noise

data = X.out


#### M1: CP decomposition (L2 loss, no penalty) ####
tmp.data = as.tensor(data)
cpold = cp(tmp.data, num_components = 2, max_iter=100)
#cp_s = c(cpold$lambdas)*normVector(cpold$U[[1]])*normVector(cpold$U[[2]])*normVector(cpold$U[[3]])
first.ind = 1; second.ind = 2
if (cpold$lambdas[1] < cpold$lambdas[2]) {
  first.ind = 2; second.ind = 1
} 
cp_u1 = c(cpold$U[[1]][,first.ind])/normVector(cpold$U[[1]][,first.ind]); cp_u1 = sign_adjust(u1, cp_u1)
cp_v1 = c(cpold$U[[2]][,first.ind])/normVector(cpold$U[[2]][,first.ind]); cp_v1 = sign_adjust(v1, cp_v1)
cp_w1 = c(cpold$U[[3]][,first.ind])/normVector(cpold$U[[3]][,first.ind]); cp_w1 = sign_adjust(w1, cp_w1)
cp_u2 = c(cpold$U[[1]][,second.ind])/normVector(cpold$U[[1]][,second.ind]); cp_u2 = sign_adjust(u2, cp_u2)
cp_v2 = c(cpold$U[[2]][,second.ind])/normVector(cpold$U[[2]][,second.ind]); cp_v2 = sign_adjust(v2, cp_v2)
cp_w2 = c(cpold$U[[3]][,second.ind])/normVector(cpold$U[[3]][,second.ind]); cp_w2 = sign_adjust(w2, cp_w2)
### initial (u,v,w) and s
#set.seed(1110)
init_s = 800
init_u = cp_u1
init_v = cp_v1
init_w = cp_w1


#### M2: RCP (L2 loss, P2) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
rcp.out = wrap_RobRSVD_tnsr(data=data, R=2, irobust=FALSE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange,
                                 nlambda=nlambda, degree=3, iugcv=TRUE, ivgcv=TRUE, iwgcv=TRUE, isspline=FALSE, verbose = TRUE, niter=20, tol=5e-5)

rcp_u1 = sign_adjust(u1, rcp.out$U[,1])
rcp_u2 = sign_adjust(u2, rcp.out$U[,2])
rcp_v1 = sign_adjust(v1, rcp.out$V[,1])
rcp_v2 = sign_adjust(v2, rcp.out$V[,2])
rcp_w1 = sign_adjust(w1, rcp.out$W[,1])
rcp_w2 = sign_adjust(w2, rcp.out$W[,2])


#### M3: RobCP (Huber loss, no penalty) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
nlambda = list(u=8, v=8, w=8)
robcp.out = wrap_RobRSVD_tnsr(data=data, R=2, irobust=TRUE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange, degree=3,
                              nlambda=nlambda, isspline=FALSE, verbose = TRUE, niter=20, tol=5e-5)
robcp_u1 = sign_adjust(u1, robcp.out$U[,1])
robcp_u2 = sign_adjust(u2, robcp.out$U[,2])
robcp_v1 = sign_adjust(v1, robcp.out$V[,1])
robcp_v2 = sign_adjust(v2, robcp.out$V[,2])
robcp_w1 = sign_adjust(w1, robcp.out$W[,1])
robcp_w2 = sign_adjust(w2, robcp.out$W[,2])


#### M4: RobRCP1 (Huber loss, P1) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
nlambda = list(u=8, v=8, w=8)
robrcp1.out = wrap_RobRSVD_tnsr(data=data, R=2, irobust=TRUE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange,
                               nlambda=nlambda, degree=3, iugcv=TRUE, ivgcv=TRUE, iwgcv=TRUE, isspline=FALSE, verbose = TRUE, niter=20, tol=5e-5)
robrcp1_u1 = sign_adjust(u1, robrcp1.out$U[,1])
robrcp1_u2 = sign_adjust(u2, robrcp1.out$U[,2])
robrcp1_v1 = sign_adjust(v1, robrcp1.out$V[,1])
robrcp1_v2 = sign_adjust(v2, robrcp1.out$V[,2])
robrcp1_w1 = sign_adjust(w1, robrcp1.out$W[,1])
robrcp1_w2 = sign_adjust(w2, robrcp1.out$W[,2])


#### M5: RobRCP2 (Huber loss, P2) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
nlambda = list(u=8, v=8, w=8)
robrcp2.out = wrap_RobRSVD2_tnsr(data=data, R=2, irobust=TRUE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange,
                                 nlambda=nlambda, degree=3, iugcv=TRUE, ivgcv=TRUE, iwgcv=TRUE, isspline=FALSE, verbose = TRUE, niter=20, tol=5e-5)
robrcp2_u1 = sign_adjust(u1, robrcp2.out$U[,1])
robrcp2_u2 = sign_adjust(u2, robrcp2.out$U[,2])
robrcp2_v1 = sign_adjust(v1, robrcp2.out$V[,1])
robrcp2_v2 = sign_adjust(v2, robrcp2.out$V[,2])
robrcp2_w1 = sign_adjust(w1, robrcp2.out$W[,1])
robrcp2_w2 = sign_adjust(w2, robrcp2.out$W[,2])


#### M6: RobCCP (Huber loss, P2 with basis constraint) ####
lambdaRange = list(u=c(1e-6, 1e2), v=c(1e-6, 1e2), w=c(1e-6, 1e2))
nlambda = list(u=8, v=8, w=8)
nknots = list(u=round(m*0.4), v=round(n*0.4), w=round(p*0.4))
robccp.out = wrap_RobRSVD2_sm_tnsr(data=data, R=2, irobust=TRUE, inits=init_s, initu=init_u, initv=init_v, initw=init_w, lambdaRange = lambdaRange, nknots=nknots,
                                       nlambda=nlambda, degree=3, iugcv=TRUE, ivgcv=TRUE, iwgcv=TRUE, uspar=0.001, vspar=0.001, wspar=0.001, isspline=FALSE, verbose = TRUE, niter=20, tol=5e-5)
robccp_u1 = sign_adjust(u1, robccp.out$U[,1])
robccp_u2 = sign_adjust(u2, robccp.out$U[,2])
robccp_v1 = sign_adjust(v1, robccp.out$V[,1])
robccp_v2 = sign_adjust(v2, robccp.out$V[,2])
robccp_w1 = sign_adjust(w1, robccp.out$W[,1])
robccp_w2 = sign_adjust(w2, robccp.out$W[,2])


### tensor reconstruction
cp.recon = cpold$est@data
rcp.recon =  rcp.out$s[1] * rcp_u1 %o% rcp_v1 %o% rcp_w1 + rcp.out$s[2] * rcp_u2 %o% rcp_v2 %o% rcp_w2
robcp.recon =  robcp.out$s[1] * robcp_u1 %o% robcp_v1 %o% robcp_w1 + robcp.out$s[2] * robcp_u2 %o% robcp_v2 %o% robcp_w2
robrcp1.recon = robrcp1.out$s[1] * robrcp1_u1 %o% robrcp1_v1 %o% robrcp1_w1 +
  robrcp1.out$s[2] * robrcp1_u2 %o% robrcp1_v2 %o% robrcp1_w2
robrcp2.recon = robrcp2.out$s[1] * robrcp2_u1 %o% robrcp2_v1 %o% robrcp2_w1 +
  robrcp2.out$s[2] * robrcp2_u2 %o% robrcp2_v2 %o% robrcp2_w2
robccp.recon = robccp.out$s[1] * robccp_u1 %o% robccp_v1 %o% robccp_w1 +
  robccp.out$s[2] * robccp_u2 %o% robccp_v2 %o% robccp_w2

### difference 

### singular value
data.frame(
  rcp=rcp.out$s, 
  robcp=robcp.out$s, 
  robrcp1=robrcp1.out$s,
  robrcp2=robrcp2.out$s, 
  robccp=robccp.out$s
) %>% signif(3)

# u
data.frame(
  cp=c("u1"=sum((u1-cp_u1)^2), "u2"=sum((u2-cp_u2)^2)),
  rcp=c("u1"=sum((u1-rcp_u1)^2), "u2"=sum((u2-rcp_u2)^2)),
  robcp=c("u1"=sum((u1-robcp_u1)^2), "u2"=sum((u2-robcp_u2)^2)),
  robrcp1=c("u1"=sum((u1-robrcp1_u1)^2), "u2"=sum((u2-robrcp1_u2)^2)),
  robrcp2=c("u1"=sum((u1-robrcp2_u1)^2), "u2"=sum((u2-robrcp2_u2)^2)),
  robccp=c("u1"=sum((u1-robccp_u1)^2), "u2"=sum((u2-robccp_u2)^2))
) %>% signif(3)

data.frame(
  cp=c("u1"=sqrt(1-sum(u1*cp_u1)^2), "u2"=sqrt(1-sum(u2*cp_u2)^2)),
  rcp=c("u1"=sqrt(1-sum(u1*rcp_u1)^2), "u2"=sqrt(1-sum(u2*rcp_u2)^2)),
  robcp=c("u1"=sqrt(1-sum(u1*robcp_u1)^2), "u2"=sqrt(1-sum(u2*robcp_u2)^2)),
  robrcp1=c("u1"=sqrt(1-sum(u1*robrcp1_u1)^2), "u2"=sqrt(1-sum(u2*robrcp1_u2)^2)),
  robrcp2=c("u1"=sqrt(1-sum(u1*robrcp2_u1)^2), "u2"=sqrt(1-sum(u2*robrcp2_u2)^2)),
  robccp=c("u1"=sqrt(1-sum(u1*robccp_u1)^2), "u2"=sqrt(1-sum(u2*robccp_u2)^2))
) %>% signif(3)


par(mfrow=c(2,3), mar=c(2,2,2,2))
compare_plot(ugrid, u1, cp_u1, title="u1-CP")
compare_plot(ugrid, u1, rcp_u1, title="u1-PCP")
compare_plot(ugrid, u1, robcp_u1, title="u1-RobCP")
compare_plot(ugrid, u1, robrcp1_u1, title="u1-RobPCP1")
compare_plot(ugrid, u1, robrcp2_u1, title="u1-RobPCP2")
compare_plot(ugrid, u1, robccp_u1, title="u1-RobPCP2-sp")

par(mfrow=c(2,3), mar=c(2,2,2,2))
compare_plot(ugrid, u2, cp_u2, title="u2-CP")
compare_plot(ugrid, u2, rcp_u2, title="u2-PCP")
compare_plot(ugrid, u2, robcp_u2, title="u2-RobCP")
compare_plot(ugrid, u2, robrcp1_u2, title="u2-RobPCP1")
compare_plot(ugrid, u2, robrcp2_u2, title="u2-RobPCP2")
compare_plot(ugrid, u2, robccp_u2, title="u2-RobPCP2-sp")


# v
data.frame(
  cp=c("v1"=sum((v1-cp_v1)^2), "v2"=sum((v2-cp_v2)^2)),
  rcp=c("v1"=sum((v1-rcp_v1)^2), "v2"=sum((v2-rcp_v2)^2)),
  robcp=c("v1"=sum((v1-robcp_v1)^2), "v2"=sum((v2-robcp_v2)^2)),
  robrcp1=c("v1"=sum((v1-robrcp1_v1)^2), "v2"=sum((v2-robrcp1_v2)^2)),
  robrcp2=c("v1"=sum((v1-robrcp2_v1)^2), "v2"=sum((v2-robrcp2_v2)^2)),
  robccp=c("v1"=sum((v1-robccp_v1)^2), "v2"=sum((v2-robccp_v2)^2))
) %>% signif(3)

data.frame(
  cp=c("v1"=sqrt(1-sum(v1*cp_v1)^2), "v2"=sqrt(1-sum(v2*cp_v2)^2)),
  rcp=c("v1"=sqrt(1-sum(v1*rcp_v1)^2), "v2"=sqrt(1-sum(v2*rcp_v2)^2)),
  robcp=c("v1"=sqrt(1-sum(v1*robcp_v1)^2), "v2"=sqrt(1-sum(v2*robcp_v2)^2)),
  robrcp1=c("v1"=sqrt(1-sum(v1*robrcp1_v1)^2), "v2"=sqrt(1-sum(v2*robrcp1_v2)^2)),
  robrcp2=c("v1"=sqrt(1-sum(v1*robrcp2_v1)^2), "v2"=sqrt(1-sum(v2*robrcp2_v2)^2)),
  robccp=c("v1"=sqrt(1-sum(v1*robccp_v1)^2), "v2"=sqrt(1-sum(v2*robccp_v2)^2))
) %>% signif(3)

par(mfrow=c(2,3), mar=c(2,2,2,2))
compare_plot(vgrid, v1, cp_v1, title="v1-CP")
compare_plot(vgrid, v1, rcp_v1, title="v1-PCP")
compare_plot(vgrid, v1, robcp_v1, title="v1-RobCP")
compare_plot(vgrid, v1, robrcp1_v1, title="v1-RobPCP1")
compare_plot(vgrid, v1, robrcp2_v1, title="v1-RobPCP2")
compare_plot(vgrid, v1, robccp_v1, title="v1-RobPCP2-sp")

par(mfrow=c(2,3), mar=c(2,2,2,2))
compare_plot(vgrid, v2, cp_v2, title="v2-CP")
compare_plot(vgrid, v2, rcp_v2, title="v2-PCP")
compare_plot(vgrid, v2, robcp_v2, title="v2-RobCP")
compare_plot(vgrid, v2, robrcp1_v2, title="v2-RobPCP1")
compare_plot(vgrid, v2, robrcp2_v2, title="v2-RobPCP2")
compare_plot(vgrid, v2, robccp_v2, title="v2-RobPCP2-sp")


# w
data.frame(
  cp=c("w1"=sum((w1-cp_w1)^2), "w2"=sum((w2-cp_w2)^2)),
  rcp=c("w1"=sum((w1-rcp_w1)^2), "w2"=sum((w2-rcp_w2)^2)),
  robcp=c("w1"=sum((w1-robcp_w1)^2), "w2"=sum((w2-robcp_w2)^2)),
  robrcp1=c("w1"=sum((w1-robrcp1_w1)^2), "w2"=sum((w2-robrcp1_w2)^2)),
  robrcp2=c("w1"=sum((w1-robrcp2_w1)^2), "w2"=sum((w2-robrcp2_w2)^2)),
  robccp=c("w1"=sum((w1-robccp_w1)^2), "w2"=sum((w2-robccp_w2)^2))
) %>% signif(3)

data.frame(
  cp=c("w1"=sqrt(1-sum(w1*cp_w1)^2), "w2"=sqrt(1-sum(w2*cp_w2)^2)),
  rcp=c("w1"=sqrt(1-sum(w1*rcp_w1)^2), "w2"=sqrt(1-sum(w2*rcp_w2)^2)),
  robcp=c("w1"=sqrt(1-sum(w1*robcp_w1)^2), "w2"=sqrt(1-sum(w2*robcp_w2)^2)),
  robrcp1=c("w1"=sqrt(1-sum(w1*robrcp1_w1)^2), "w2"=sqrt(1-sum(w2*robrcp1_w2)^2)),
  robrcp2=c("w1"=sqrt(1-sum(w1*robrcp2_w1)^2), "w2"=sqrt(1-sum(w2*robrcp2_w2)^2)),
  robccp=c("w1"=sqrt(1-sum(w1*robccp_w1)^2), "w2"=sqrt(1-sum(w2*robccp_w2)^2))
) %>% signif(3)


par(mfrow=c(2,3), mar=c(2,2,2,2))
compare_plot(wgrid, w1, cp_w1, title="w1-CP")
compare_plot(wgrid, w1, rcp_w1, title="w1-PCP")
compare_plot(wgrid, w1, robcp_w1, title="w1-RobCP")
compare_plot(wgrid, w1, robrcp1_w1, title="w1-RobPCP1")
compare_plot(wgrid, w1, robrcp2_w1, title="w1-RobPCP2")
compare_plot(wgrid, w1, robccp_w1, title="w1-RobPCP2-sp")

par(mfrow=c(2,3), mar=c(2,2,2,2))
compare_plot(wgrid, w2, cp_w2, title="w2-CP")
compare_plot(wgrid, w2, rcp_w2, title="w2-PCP")
compare_plot(wgrid, w2, robcp_w2, title="w2-RobCP")
compare_plot(wgrid, w2, robrcp1_w2, title="w2-RobPCP1")
compare_plot(wgrid, w2, robrcp2_w2, title="w2-RobPCP2")
compare_plot(wgrid, w2, robccp_w2, title="w2-RobPCP2-sp")


# reconstructed tensor
data.frame(
  cp=mean(c((cp.recon - X.sig)^2)),
  rcp=mean(c((rcp.recon - X.sig)^2)),
  robcp=mean(c((robcp.recon - X.sig)^2)),
  robrcp1=mean(c((robrcp1.recon - X.sig)^2)),
  robrcp2=mean(c((robrcp2.recon - X.sig)^2)), 
  robccp=mean(c((robccp.recon - X.sig)^2))
) %>% signif(3)

