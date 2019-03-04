fuzzyROC = function(
  XX, # continuous marker assumes low levels of the marker predict class 0 = "disease"
  YY, # binary class
  grayMax = .2 # maximum tolerated indeterminacy zone (gamma in paper)
)
{
  XXu = sort(unique(XX))
  Nx = length(XXu)
  XXmid = ( XXu[1:(Nx-1)] + XXu[2:Nx] ) / 2
  Nx = length(XXmid)
  F0star = F0dagger = F1star = F1dagger = cLBest = cHBest = grayBest = rep(NA,Nx)
  for ( jj in 1:Nx ) {
    XXl = rev ( XXu[ XXu < XXmid[jj]] )
    XXr = XXu[ XXu > XXmid[jj]]
    if (length(XXl) > length(XXr)) XXr = c(XXr,rep(max(XXr),length(XXl) - length(XXr)))
    if (length(XXl) < length(XXr)) XXl = c(XXl,rep(min(XXl),length(XXr) - length(XXl)))
    UU = rep(NA,length(XXl))
    for (ii in 1:length(XXl)) {
      if ( mean(XX > XXl[ii] & XX < XXr[ii]) > grayMax ) {
        break
      }
      else{
        XXX = XX[ ! (XX > XXl[ii] & XX < XXr[ii]) ]
        YYY = YY[ ! (XX > XXl[ii] & XX < XXr[ii]) ]
        n0 = sum(YYY==0)
        n1 = sum(YYY==1)
        UU[ii] = wilcox.test(XXX[YYY==1],XXX[YYY==0])$statistic / (n1*n0)
      }
    }
#    if(jj==25) browser()
    UUBest = max( UU, na.rm = TRUE )
    iiBest =  min( (1:length(XXl))[ UU == UUBest ], na.rm = TRUE )
    cLBest[jj] = XXl[iiBest]
    cHBest[jj] = XXr[iiBest]
    grayBest[jj] = mean(XX > cLBest[jj] & XX < cHBest[jj])
    if (grayBest[jj]==0){
      XXstar = XXdagger = XX
    }
    else{
      XXstar = XXdagger = XX
      XXstar [ XX > cLBest[jj] & XX < cHBest[jj] & YY==1 ] = cHBest[jj]
      XXstar [ XX > cLBest[jj] & XX < cHBest[jj] & YY==0 ] = cLBest[jj]
      XXdagger [ XX > cLBest[jj] & XX < cHBest[jj] & YY==1 ] = cLBest[jj]
      XXdagger [ XX > cLBest[jj] & XX < cHBest[jj] & YY==0 ] = cHBest[jj]
    }
    F0star[jj] = sum( XXstar[YY==0] <= XXmid[jj] ) / sum( YY==0 )
    F1star[jj] = sum( XXstar[YY==1] <= XXmid[jj]  ) / sum( YY==1 )
    F0dagger[jj] = sum( XXdagger[YY==0] <= XXmid[jj] ) / sum( YY==0 )
    F1dagger[jj] = sum( XXdagger[YY==1] <= XXmid[jj] ) / sum( YY==1 )
    #  if(grayBest[jj]>0) browser()
  }
  return(list(XX=XX,
              YY=YY,
              grayMax=grayMax,
              cMid=XXmid,
              cLstar=cLBest, # vector of lower limits
              cHstar=cHBest, # vector of optimal upper limits
              grayBest=grayBest, # proportion of points in the optimal gray area
              F0star=F0star,
              F0dagger=F0dagger,
              F1star=F1star,
              F1dagger=F1dagger))
}
