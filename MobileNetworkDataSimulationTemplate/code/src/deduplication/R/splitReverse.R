#
splitReverse <-function(N, M){
  k=0
  forward <-(0:(M-1))%%M
  backward<-((M-1):0)%%M
  index<-c()
  while(length(index)<N){
    if(k%%2==0)
      index<-c(index,forward)
    else
      index<-c(index,backward)
    k=k+1
  }
  index<-index[1:N]
  ichunks<-list()
  for(i in 1:M) {
    ichunks[[i]]<-which(index==(i-1))
  }
  return(ichunks) 
}
