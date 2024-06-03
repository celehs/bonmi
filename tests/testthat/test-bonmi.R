test_that("bonmi works", {
  set.seed(1)
  N = 3000
  r = 10
  m = 5
  p = 0.1
  X0 = matrix(rnorm(N*r),nrow=N)
  W0 = X0 %*% t(X0)
  rownames(W0) = colnames(W0) = paste0('code', 1:N)

  W = list()
  for(s in 1:m){
    ids = which(runif(N) < p)
    Ns = length(ids)
    Es = matrix(rnorm(Ns * Ns, sd = s * 0.01), nrow = Ns)
    Es = Es + t(Es)
    Ws = W0[ids,ids] + Es
    W[[s]] = Ws
  }

  Xhat <- bonmi(W, r, weights=NULL, rsvd.use=FALSE)
  #Xhat <- BONMI(W, r, weights=NULL, rsvd.use=TRUE)
  codes = rownames(Xhat)
  What = Xhat %*% t(Xhat)

  id = match(codes,rownames(W0))
  Wstar = W0[id, id]

  #BONMI's result
  result = norm(What - Wstar, 'F')/norm(Wstar, 'F')
  print(result)
  expect_equal(result, 0.0037467037)
})
