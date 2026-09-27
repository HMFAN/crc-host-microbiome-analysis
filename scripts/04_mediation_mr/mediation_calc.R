calc<-function(X,d,M,y,int,t0){
 ma<-lm.fit(cbind(X,d),M);if(ma$rank!=ncol(X)+1)stop('mediator rank')
 V<-if(int)cbind(X,d,M,dxM=d*M) else cbind(X,d,M)
 f<-suppressWarnings(glm.fit(V,y,family=binomial()));if(!f$converged || f$rank!=ncol(V) || any(!is.finite(f$coefficients)))stop('outcome fit')
 k<-ncol(X);mu<-drop(X%*%ma$coefficients[seq_len(k)]);b<-f$coefficients;eta<-drop(X%*%b[seq_len(k)])
 pred<-function(t,mm)plogis(eta+b[k+1]*t+b[k+2]*mm+if(int)b[k+3]*t*mm else 0)
 m0<-mu+ma$coefficients[k+1]*t0;m1<-mu+ma$coefficients[k+1]*(t0+1)
 p00<-pred(t0,m0);p01<-pred(t0,m1);p10<-pred(t0+1,m0);p11<-pred(t0+1,m1)
 c(ACME=mean(c(p01-p00,p11-p10)),ADE=mean(c(p10-p00,p11-p01)),total=mean(p11-p00))
}
