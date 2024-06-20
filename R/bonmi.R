#' BONMI: Block-wise Overlapping Noisy Matrix Integration
#'
#' @param W a list of PPMI matrices, with rownames and colnames being the features
#' @param r rank
#' @param weights The weight vector for the PPMI matrices. Default NULL, the weights will be estimated from data.
#' @param rsvd.use Bool. Default FALSE. If TRUE, we will use the 'rsvd' function to calculate the svd, which is much faster than the 'svd' function  when r is small.
#' @return An embedding matrix.
#' @examples
#' set.seed(1)
#' N = 3000
#' r = 10
#' m = 5
#' p = 0.1
#' X0 = matrix(rnorm(N*r),nrow=N)
#' W0 = X0 %*% t(X0)
#' rownames(W0) = colnames(W0) = paste0('code',1:N)
#'
#' W = list()
#' for(s in 1:m){
#'   ids = which(runif(N)<p)
#'   Ns = length(ids)
#'   Es = matrix(rnorm(Ns*Ns,sd=s*0.01),nrow=Ns)
#'   Es = Es + t(Es)
#'   Ws = W0[ids,ids] + Es
#'   W[[s]] = Ws
#' }
#'
#' Xhat <- bonmi(W,r,weights=NULL,rsvd.use=FALSE)
#' #Xhat <- bonmi(W,r,weights=NULL,rsvd.use=TRUE)
#' codes = rownames(Xhat)
#' What = Xhat%*%t(Xhat)
#'
#' id = match(codes,rownames(W0))
#' Wstar = W0[id,id]
#'
#' #bonmi's result
#' norm(What-Wstar,'F')/norm(Wstar,'F')
#'
#' @export
bonmi <- function(W, r, weights=NULL, rsvd.use=FALSE){

  m = length(W)
  codes = NULL
  for(s in 1:m){
    codes = union(codes, rownames(W[[s]]))
  }
  codes = sort(codes)

  if(is.null(weights)){
    weights = sapply(1:m, function(s){
      if(rsvd.use){
        set.seed(1)
        fit = rsvd::rsvd(W[[s]], r+1)
      }else{
        fit = svd(W[[s]], r+1, r+1)
      }
      w = fit$d[r+1] / sqrt(nrow(W[[s]]))
      return(w)
    })
  }

  Wc = matrix(0, nrow=length(codes), ncol=length(codes))
  C = matrix(0, nrow=length(codes), ncol=length(codes))
  Weis = matrix(0, nrow=length(codes), ncol=length(codes))

  rownames(Wc) = colnames(Wc) = codes
  rownames(C) = colnames(C) = codes
  rownames(Weis) = colnames(Weis) = codes

  for(s in 1:m){
    id = match(rownames(W[[s]]), codes)
    Wc[id, id] = Wc[id, id] + weights[s] * W[[s]]
    C[id, id] = C[id, id]+1
    Weis[id,id] = Weis[id,id]+weights[s]
  }
  Wc[C>0] = Wc[C>0]/Weis[C>0]
  Wo = Wc
  Wo[C==0] = NA

  W.new = list()
  for(s in 1:m){
    Is = match(rownames(W[[s]]),rownames(Wc))
    W.new[[s]] = Wc[Is, Is]
  }

  fits = lapply(1:m, function(s){
    if(rsvd.use){
      set.seed(1)
      return(rsvd::rsvd(W.new[[s]],r))
    }else{
      return(svd(W.new[[s]],r,r))
    }
  })

  Xs = lapply(1:m, function(s){
    U = embedding(fits[[s]], r)
    rownames(U) = rownames(W[[s]])
    U
  })

  Wm = matrix(0, nrow=nrow(Wc), ncol=ncol(Wc))
  M = matrix(0, nrow=nrow(Wc), ncol=ncol(Wc))
  rownames(Wm) = colnames(Wm) = rownames(Wc)
  rownames(M) = colnames(M) = rownames(Wc)

  for(s in 1:(m-1)){
    for(k in (s+1):m){
      name12 = intersect(rownames(W[[s]]), rownames(W[[k]]))
      ids = match(name12,rownames(W[[s]]))
      idk = match(name12,rownames(W[[k]]))
      Osk = Procrustes(Xs[[s]][ids, ], Xs[[k]][idk, ])
      Wsk = Xs[[s]][-ids,] %*% t(Osk) %*% t(Xs[[k]][-idk,])
      id1 = match(rownames(W[[s]])[-ids], rownames(Wm))
      id2 = match(rownames(W[[k]])[-idk], rownames(Wm))
      Wm[id1,id2] = Wm[id1,id2] + Wsk
      Wm[id2,id1] =  Wm[id2,id1] + t(Wsk)
      M[id1,id2] = M[id1,id2] + 1
      M[id2,id1] = M[id2,id1] + 1
    }
  }
  Wm[M>0] = Wm[M>0]/M[M>0]
  Wm[C>0] = Wc[C>0]


  if(rsvd.use){
    set.seed(1)
    fit.W = rsvd::rsvd(Wm, r)
  }else{
    fit.W = svd(Wm, r, r)
  }

  X = embedding(fit.W, r)
  rownames(X) = rownames(Wm)
  return(X)
}


#' Learn the orthogonal transformation matrix O
#' @param X1 matrix
#' @param X2 matrix
#' @noRd
Procrustes <- function(X1, X2){
  #Return Omega = arg min||X1 - X2 Omgea||_F
  H = t(X2) %*% X1
  mod = svd(H)
  return(mod$u %*% t(mod$v))
}

#' Learn the embedding
#' @param fit fit
#' @param r rank
#' @noRd
embedding <- function(fit, r){
  embed <- fit$u[,1:r] %*% diag(sqrt(fit$d[1:r]))
  embed
}



