library(stats) 
library(ggplot2)
library(invgamma)
library(MASS)
library(data.table)
library(readxl)
library(writexl)


Data_gfg <- read_excel("C:/Users/cheny/OneDrive/桌面/碩二/新版model/data/standard_data 刑除人口 + 性別比 - 複製.xlsx")
Data_gfg
Data_C <- read_excel("C:/Users/cheny/OneDrive/桌面/碩二/新版model/C完整版_高雄.xlsx")
Data_C


set.seed(2023)
m=50000
beta  = array(2, dim=c(5,m))
alpha = array(1, dim=c(4,m))
alpha[1,1]=-1
alpha[3,1]=-1
Fi    = array(0, dim=c(8,m))
Fi0   = array(0, dim=c(1,m))
sig   = array(1, dim=c(8,m))
sig0  = array(1, dim=c(1,m))
sigR  = array(1, dim=c(1,m))
Eta   = array(1, dim=c(38,m))
z     = array(0, dim=c(38,m))
u     = array(0, dim=c(38,m))
delta = array(0, dim=c(38,m,8))
delta0= array(0, dim=c(38,m))
Q     = array(0, dim=c(38,38,9,m))


f=10
ee=10
#beta[,1] = rnorm(4,0,ee)     從零??????
#alpha[,1]= rnorm(4,0,f)      從零??????
#Fi[,1]   = runif(8,0,1)      ???0.5??????   
#Fi0[,1]  = runif(1,0,1)      ???0.5??????
a=4
b=13.4
aR=4
bR=13.4
#sig[,1]  = rinvgamma(8, shape = a, rate = b, scale = 1/rate)
#sig0[,1] = rinvgamma(1, shape = a, rate = b, scale = 1/rate)
#sigR[,1] = rinvgamma(1, shape = aR, rate = bR, scale = 1/rate)
d=1
#Eta[,1]  = rnorm(19,0,d)     從零??????
z[,1]    = rexp(38, 1 )
u[,1]    = rnorm(38,0,1)
Tow = 0.5





X     = array(0, dim=c(5,38,8))
X[1,,]=1
for (t in 0:7){
  X[2:4,,t+1]=t(Data_gfg[(1+(t*38)):(38+(t*38)),3:5])
  X[5,,]     =t(Data_gfg[(1+(t*38)):(38+(t*38)),9])
}
X[2,,]=(X[2,,]-min(X[2,,]))/(max(X[2,,])-min(X[2,,]))
X[3,,]=(X[3,,]-min(X[3,,]))/(max(X[3,,])-min(X[3,,]))
X[5,,]=(X[5,,]-min(X[5,,]))/(max(X[5,,])-min(X[5,,]))



W     = array(0, dim=c(4,38))
W[1,]=1
W[2:4,]=t(Data_gfg[1:38,6:8])
W[4,]=(W[4,]-min(W[4,]))/(max(W[4,])-min(W[4,]))



Y     = array(0, dim=c(38,8))
for (t in 0:7){
  print(t)
  Y[,t+1]=t(Data_gfg[(1+(38*t)):(38+(38*t)),2])
}


#E1=sum(Y)
#E2=0
#for (r in 1:8){
#  E2=0
#  for(c in 1:38){
#    if (Y[c,r]==0){
#      E2=E2+1
#    }
#  }
#  print(E2)
#}
E2=6

Y = t(Y)


C     = array(0, dim=c(38,38))
for (t in 0:37){
  C[1:38,t+1]=t(Data_C[1:38,t+1])
}

C_plus= matrix(0, nrow = 38, ncol = 38)
I19   = diag(38)
for (i in 1:38){
  C_plus[i,i]=sum(C[,i])
}


Fi_max=1/max(eigen(C)$values)
Fi0[,1]=0.5*Fi_max
#Fi0[,1]=runif(1,0,Fi_max)
Fi[,1]=0.5*Fi_max
#Fi[,1]=runif(8,0,Fi_max)


Q1 = -Fi[1,1] * C + I19
Q[,,1,1]=Q1

Q2 = -Fi[2,1] * C + I19
Q[,,2,1]=Q2

Q3 = -Fi[3,1] *C + I19
Q[,,3,1]=Q3

Q4 = -Fi[4,1] * C+ I19
Q[,,4,1]=Q4

Q5 = -Fi[5,1] * C + I19
Q[,,5,1]=Q5

Q6 = -Fi[6,1] * C + I19
Q[,,6,1]=Q6

Q7 = -Fi[7,1] * C + I19
Q[,,7,1]=Q7

Q8 = -Fi[8,1] * C + I19
Q[,,8,1]=Q8

Q0 = -Fi0[,1] * C + I19
Q[,,9,1]=Q0



zero19=Eta[,1]

delta0[,1] = Eta[,1]   + mvrnorm(n = 1, mu = zero19, Sigma = I19)%*% chol(sig0[,1]*solve(Q0))

delta[,1,1]= delta0[,1]+ mvrnorm(n = 1, mu = zero19, Sigma = I19)%*% chol(sig[1,1]*solve(Q1))

delta[,1,2]= delta0[,1]+ mvrnorm(n = 1, mu = zero19, Sigma = I19)%*% chol(sig[2,1]*solve(Q2))

delta[,1,3]= delta0[,1]+ mvrnorm(n = 1, mu = zero19, Sigma = I19)%*% chol(sig[3,1]*solve(Q3))

delta[,1,4]= delta0[,1]+ mvrnorm(n = 1, mu = zero19, Sigma = I19)%*% chol(sig[4,1]*solve(Q4))

delta[,1,5]= delta0[,1]+ mvrnorm(n = 1, mu = zero19, Sigma = I19)%*% chol(sig[5,1]*solve(Q5))

delta[,1,6]= delta0[,1]+ mvrnorm(n = 1, mu = zero19, Sigma = I19)%*% chol(sig[6,1]*solve(Q6))

delta[,1,7]= delta0[,1]+ mvrnorm(n = 1, mu = zero19, Sigma = I19)%*% chol(sig[7,1]*solve(Q7))

delta[,1,8]= delta0[,1]+ mvrnorm(n = 1, mu = zero19, Sigma = I19)%*% chol(sig[8,1]*solve(Q8))


check_beta    = array(0, dim=c(8))



Fi_full= function(Q,delta0,Eta,C,C_plus,sig0,Fi0_next){
  ((det(solve(Q)))^(-1/2)*exp(-((delta0-Eta)%*% (Fi0_next*(C_plus-C-I19)) %*%(delta0-Eta))/(2*sig0)))
}

Y_full_ui=function(Y,a,zi,ui,uinext,si,xi,bi,di,W,E1,E2){
  if (Y==0){
    if ((-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*uinext))+exp(xi%*%bi+di)))>700 || (-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))>700){
      if((-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*uinext))))>700 || (-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))>700){
        0
      }
      else{
        log((E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui)))))/
              (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*uinext))))))
      }
    }
    else{
      if((-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*uinext))))>700 || (-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))>700){
        log((E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*uinext))+exp(xi%*%bi+di))))/
              (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))))
      }
      else{
        log(((E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*uinext))+exp(xi%*%bi+di))))/
               (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))))/
              ((E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*uinext)))))/
                 (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui)))))))
      }
    }
    
    
  }
  else{
    
    log((exp((W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))    *E2+1)/(exp((W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*uinext))))*E2+1))
  }
}



Y_full_zi=function(Y,a,zi,zinext,ui,si,xi,bi,di,W,E1,E2){
  if (Y==0){
    
    if ((-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zinext+(2/(Tow*(1-Tow)))*(zinext^(1/2)*ui))+exp(xi%*%bi+di)))>700 || (-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))>700){
      if ((-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zinext+(2/(Tow*(1-Tow)))*(zinext^(1/2)*ui))))>700 || (-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))>700){
        0
      }
      else{
        log((E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui)))))/
              (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zinext+(2/(Tow*(1-Tow)))*(zinext^(1/2)*ui))))))
      }
    }
    else{
      if ((-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zinext+(2/(Tow*(1-Tow)))*(zinext^(1/2)*ui))))>700 || (-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))>700){
        log((E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zinext+(2/(Tow*(1-Tow)))*(zinext^(1/2)*ui))+exp(xi%*%bi+di))))/
              (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))))
      }
      else{
        log(((E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zinext+(2/(Tow*(1-Tow)))*(zinext^(1/2)*ui))+exp(xi%*%bi+di))))/
               (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))))/
              ((E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zinext+(2/(Tow*(1-Tow)))*(zinext^(1/2)*ui)))))/
                 (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui)))))))
      }
    }
    
    
  }
  else{
    
    log((exp((W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi    +(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))    *E2+1)/(exp((W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zinext+(2/(Tow*(1-Tow)))*(zinext^(1/2)*ui))))*E2+1))
    
  }
}


Y_full_sigR=function(Y,a,zi,ui,si,sinext,xi,bi,di,W,E1,E2){
  if (Y==0){
    if ((-(W%*%(a)+sinext*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))>700 || (-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))>700){
      if ((-(W%*%(a)+sinext*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))>700 || (-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))>700){
        0
      }
      else{
        log((E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui)))))/
              (E2+exp(-(W%*%(a)+sinext*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))))
        
      }
    }
    else{
      if((-(W%*%(a)+sinext*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))>700 || (-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))>700){
        log((E2+exp(-(W%*%(a)+sinext*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di))))/
              (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))))
      }
      else{
        log(((E2+exp(-(W%*%(a)+sinext*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di))))/
               (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))))/
              ((E2+exp(-(W%*%(a)+sinext*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui)))))/
                 (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui)))))))
      }
    }
    
    
  }
  else{
    
    log((exp((W%*%(a)+si    *(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))*E2+1)/(exp((W%*%(a)+sinext*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))*E2+1))
  }
}

Y_full_alpha=function(Y,a,anext,zi,ui,si,xi,bi,di,W,E1,E2){
  if (Y==0){
    
    
    log(((E2+exp(-(W%*%(anext)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di))))/
           (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))))/
          ((E2+exp(-(W%*%(anext)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui)))))/
             (E2+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui)))))))
    
  }
  else{
    
    log((exp((W%*%(a)    +si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))*E2+1)/(exp((W%*%(anext)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))))*E2+1))
    
  }
}


Y_full_beta=function(Y,a,zi,ui,si,xi,bi,binext,di,W,E1,E2){
  if (Y==0){
    
    log((1+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%binext+di)))/E2)/
          (1+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))/E2))
    
  }
  else{
    
    (-(exp(xi%*%binext+di)))+Y*(xi%*%binext+di)-(-(exp(xi%*%bi+di)))-Y*(xi%*%bi+di)
    
  }
}

Y_full_deltati=function(Y,a,zi,ui,si,xi,bi,di,dinext,W,E1,E2){
  if (Y==0){
    
    log((1+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+dinext)))/E2)/
          (1+exp(-(W%*%(a)+si*(((1-2*Tow)/(Tow*(1-Tow)))*zi+(2/(Tow*(1-Tow)))*(zi^(1/2)*ui))+exp(xi%*%bi+di)))/E2))
    
    
    
  }
  else{
    -(exp(xi%*%bi+dinext)-exp(xi%*%bi+di))+Y*(xi%*%bi+dinext)-Y*(xi%*%bi+di)
  }
}




D    =function(i,b,a,f,f0,cp){
  D=0
  for (t in 1:8){
    D=D+(B*(f[t]*cp[i,i]+(1-f[t])))/b[t]+(B*(f0*cp[i,i]+(1-f0)))/a
  }
  D
}

E =function(i,C,delta,delta0){
  E=0
  for (j in 1:38){
    E=E+C[i,j]*(delta[j]-delta0[j])
  }
  E
}


Ft =function(i,C,delta,delta0,Fi){
  Ft=0
  for (t in 1:8){
    Ft=Ft+(Fi[t]*C_plus[i,i]+1-Fi[t])*(delta[,t]-(Fi[t]*E(i,C,delta[,t],delta0)/(Fi[t]*C_plus[i,i]+1-Fi[t])))/sig[t,r]
  }
  Ft
}


rate_p    = array(0, dim=c(3))
r=2
rate_0 = (b+0.5*t(delta0[,r-1]-Eta[,r-1])%*% Q[,,4,r-1] %*%(delta0[,r-1]-Eta[,r-1]))
for (t in 1:3){
  rate_p[t]=(b+0.5*(delta[,r-1,t]-delta0[,r-1])%*% Q[,,t,r-1] %*%(delta[,r-1,t]-delta0[,r-1]))
}





start_time <- Sys.time()

for (r in 2:m){
  #print('----------------------------------------------------------------------------------------------------------------------------------------------------------------')
  print(r)
  
  #ui
  #print('ui')
  for (i in 1:38){
    uni = runif(1,0,1)
    ui_next=rnorm(1,u[i,r-1],sqrt(1))
    
    buf=0
    #print(i)
    for (t in 1:8){
      #print(Y[t,i])
      buf = buf+(Y_full_ui(Y[t,i],alpha[,r-1],z[i,r-1],u[i,r-1],ui_next,sigR[,r-1],X[,i,t],beta[,r-1],delta[i,r-1,t],W[,i],E1,E2))
      
    }
    buf = buf-(ui_next^2-u[i,r-1]^2)/2
    alpha_ui=min(0,buf)
    if (log(uni)<alpha_ui){ 
      u[i,r]=ui_next
    }
    else{
      u[i,r]=u[i,r-1]
    }
  }
  
  #zi
  #print('zi')
  for (i in 1:38){
    uni = runif(1,0,1)
    zi_next=rexp(1, 1/z[i,r-1] )
    while (zi_next<1e-200){
      zi_next=rexp(1, 1/z[i,r-1] )
    }
    
    buf=0
    
    for (t in 1:8){
      buf = buf+(Y_full_zi(Y[t,i],alpha[,r-1],z[i,r-1],zi_next,u[i,r],sigR[,r-1],X[,i,t],beta[,r-1],delta[i,r-1,t],W[,i],E1,E2))
    }
    buf = buf +z[i,r-1]-zi_next+log(z[i,r-1]/zi_next)+zi_next/z[i,r-1]-z[i,r-1]/zi_next
    alpha_zi=min(0,buf)
    if (log(uni)<alpha_zi){ 
      z[i,r]=zi_next
    }
    else{
      z[i,r]=z[i,r-1]
    }
  }
  
  #Eta i
  #print('Eta')
  for (i in 1:38){
    
    Eta[i,r]=rnorm(1,(delta0[i,r-1]-Fi0[,r-1]*E(i,C,delta0[,r-1],Eta[,r-1]))/(sig0[,r-1]^2*(1/d+1/sig0[,r-1]^2)),(1/(1/d+1/sig0[,r-1]^2))**0.5)
    
  }
  
  
  #delta0 i 
  #print('delta0')
  A=0
  for (t in 1:8){
    A=A+(1/sig[t,r-1])
  }
  B=sig0[,r]
  for (t in 1:8){
    B=B*sig[t,r-1]
  }
  
  for (i in 1:38){
    AAA=0
    for (t in 1:8){
      AAA=AAA+(delta[i,r-1,t]-Fi[t,r-1]*E(i,C,delta[,r-1,t],delta0[,r-1]))/sig[t,r-1]
    }
    AAA=AAA+(Eta[i,r]+Fi0[,r-1]*E(i,C,delta0[,r-1],Eta[,r]))/sig0[,r-1]
    
    delta0[i,r]=rnorm(1,AAA/A,(1/A)**0.5)
  } 
  
  
  #delta t i 
  #print('delta')
  for (i in 1:38){
    for (t in 1:8){
      
      
      uni = runif(1,0,1)
      delta_next=rnorm(1,delta[i,r-1,t],sqrt(1))
      
      buf=Y_full_deltati(Y[t,i],alpha[,r-1],z[i,r],u[i,r],sigR[,r-1],X[,i,t],beta[,r-1],delta[i,r-1,t],delta_next,W[,i],E1,E2)-
        ((delta_next-delta[i,r-1,t])*(delta_next+delta[i,r-1,t]-2*(delta0[i,r]+Fi[t,r-1]*E(i,C,delta[,r-1,t],delta0[,r]))))/2*sig[t,r-1]^2
      
      alpha_delta=min(0,buf)
      if (log(uni)<alpha_delta){ 
        delta[i,r,t]=delta_next
      }
      else{
        delta[i,r,t]=delta[i,r-1,t]
      }
    }
  }
  
  #sigma0
  #print('sigma0')
  if ((b+0.5*t(delta0[,r]-Eta[,r])%*% Q[,,9,r-1] %*%(delta0[,r]-Eta[,r]))>0){
    rate_0=(b+0.5*t(delta0[,r]-Eta[,r])%*% Q[,,9,r-1] %*%(delta0[,r]-Eta[,r]))
  }
  sig0[,r] = rinvgamma(1, shape = a+38/2, rate = 1/rate_0,  scale= 1/rate)
  
  
  #sigma
  #print('sigma')
  for (t in 1:8){
    
    if ((b+0.5*t(delta[,r,t]-delta0[,r])%*% Q[,,t,r-1] %*%(delta[,r,t]-delta0[,r])) >0){
      rate_p[t]=(b+0.5*t(delta[,r,t]-delta0[,r])%*% Q[,,t,r-1] %*%(delta[,r,t]-delta0[,r]))
    }
    sig[t,r]  = rinvgamma(1, shape = a+38/2,rate  = 1/rate_p[t],  scale= 1/rate)
  }
  
  
  #sigR
  #print('sigR')
  uni = runif(1,0,1)
  var_sigR=0.2
  sigR_next=rinvgamma(1, shape = 2+sigR[,r-1]/var_sigR,rate  = (1+sigR[,r-1]/var_sigR)*sigR[,r-1], scale = 1/rate)
  while (sigR_next < 1e-20){
    sigR_next=rinvgamma(1, shape = 2+sigR[,r-1]/var_sigR,rate  = (1+sigR[,r-1]/var_sigR)*sigR[,r-1], scale = 1/rate)
  }
  
  buf =0
  for (i in 1:38){
    for (t in 1:8){
      
      buf = buf+Y_full_sigR(Y[t,i],alpha[,r-1],z[i,r],u[i,r],sigR[,r-1],sigR_next,X[,i,t],beta[,r-1],delta[i,r,t],W[,i],E1,E2)
      +(aR+1)*log(sigR[,r-1]/sigR_next)-(bR/sigR_next)+(bR/sigR[,r-1])
      +log((dinvgamma( sigR_next, shape = 2+sigR[,r-1]/var_sigR,rate  = (1+sigR[,r-1]/var_sigR)*sigR[,r-1], scale = 1/rate))/(dinvgamma(sigR[,r-1], shape = 2+sigR_next/var_sigR ,rate  = (1+sigR_next/var_sigR)*sigR_next , scale = 1/rate)))
    }
  }
  
  alpha_sigR=min(0,buf)
  if (log(uni)<alpha_sigR){ 
    sigR[,r]=sigR_next
  }
  else{
    sigR[,r]=sigR[,r-1]
  }
  
  #alpha l
  #print('alpha')
  for (l in 1:4){
    #print(l)
    uni = runif(1,0,1)
    alpha_next=runif(1,alpha[l,r-1]-0.75,alpha[l,r-1]+0.75) #rnorm(1,alpha[l,r-1],0.1)
    
    af_old=alpha[,r-1]
    af_new=alpha[,r-1]
    if (l != 1){
      af_old[1:(l-1)]=alpha[1:(l-1),r]
      af_new[1:(l-1)]=alpha[1:(l-1),r]
    }
    af_new[l]=alpha_next
    buf =0
    
    for (i in 1:38){
      for (t in 1:8){
        
        buf=buf+Y_full_alpha(Y[t,i],af_old,af_new,z[i,r],u[i,r],sigR[,r],X[,i,t],beta[,r-1],delta[i,r,t],W[,i],E1,E2)
        
      }
    }
    
    alpha_alpha=min(0,buf)
    if (log(uni)<alpha_alpha){ 
      alpha[l,r]=alpha_next
    }
    else{
      alpha[l,r]=alpha[l,r-1]
    }
    
  }
  
  
  #beta k
  #print('beta')
  for (k in 1:5){
    #print('---------------------------------------------------------------------------------')
    #print(k)
    uni = runif(1,0,1)
    #beta_next=rnorm(1,beta[k,r-1],sqrt(1))
    if ( k==4 ){
      beta_next=rnorm(1,beta[k,r-1],sqrt(0.1))
    }
    else if(k==1){
      beta_next=rnorm(1,beta[k,r-1],sqrt(0.01))
    }
    else{
      beta_next=rnorm(1,beta[k,r-1],sqrt(1))
    }
    
    bt_old=beta[,r-1]
    bt_new=beta[,r-1]
    if (k != 1){
      bt_old[1:(k-1)]=beta[1:(k-1),r]
      bt_new[1:(k-1)]=beta[1:(k-1),r]
    }
    bt_new[k]=beta_next
    buf=0
    for (i in 1:38){
      for (t in 1:8){
        
        buf=buf+Y_full_beta(Y[t,i],alpha[,r],z[i,r],u[i,r],sigR[,r],X[,i,t],bt_old,bt_new,delta[i,r,t],W[,i],E1,E2)
        
      }
    }
    buf = buf-(beta_next^2-beta[k,r-1]^2)/(2*ee)
    alpha_beta=min(0,buf)
    if (log(uni)<alpha_beta){ 
      beta[k,r]=beta_next
    }
    else{
      beta[k,r]=beta[k,r-1]
    }
  }
  
  
  #Fi0
  #print('fi0')
  uni=runif(1,0,1)
  Fi0_next=runif(1,Fi0[,r-1]-0.005,Fi0[,r-1]+0.005)
  while (Fi0_next>Fi_max || Fi0_next<0){
    Fi0_next=runif(1,Fi0[,r-1]-0.005,Fi0[,r-1]+0.005)
  }
  
  Q_next= -Fi0_next*C + I19
  
  buf= (1/2)*(log(det(solve(Q[,,9,r-1])))-log(det(solve(Q_next)))+(Fi0[,r-1]-Fi0_next)*(t(delta0[,r]-Eta[,r])%*% C %*%(delta0[,r]-Eta[,r]))/sig0[,r])
  
  alpha_Fi0=min(0,buf)
  if (log(uni)<alpha_Fi0){
    Fi0[,r] =Fi0_next
    Q[,,9,r]=Q_next
  }
  else{
    Fi0[,r] =Fi0[,r-1]
    Q[,,9,r]=Q[,,9,r-1]
  }
  
  
  #Fi
  #print('fi')
  for (t in 1:8){
    uni=runif(1,0,1)
    Fi_next=runif(1,Fi[t,r-1]-0.0025,Fi[t,r-1]+0.0025)
    while (Fi_next>Fi_max || Fi_next<0){
      Fi_next=runif(1,Fi[t,r-1]-0.0025,Fi[t,r-1]+0.0025)
    }
    
    Q_next= -Fi_next*C + I19
    
    buf = (1/2)*(log(det(solve(Q[,,t,r-1])))-log(det(solve(Q_next)))+(Fi[t,r-1]-Fi_next)*(t(delta[,r,t]-delta0[,r])%*% C %*%(delta[,r,t]-delta0[,r]))/sig[t,r])
    
    alpha_Fi=min(0,buf)
    if (log(uni)<alpha_Fi){ 
      Fi[t,r]=Fi_next
      Q[,,t,r]=Q_next
    }
    else{
      Fi[t,r]=Fi[t,r-1]
      Q[,,t,r]=Q[,,t,r-1]
    }
  } 
  
  
}





change_rate = array(0, dim=c(76,11))

print('---------------------------------------------------------------------------------')
print('alpha????????????')
for (l in 1:4){
  rate =0
  
  for (i in 1:(m-1)){
    if (alpha[l,i]!=alpha[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l,1]=rate/(m-1)
  print(rate/(m-1))
}

print('---------------------------------------------------------------------------------')
print('beta????????????')
for (l in 1:5){
  rate =0
  for (i in 1:(m-1)){
    if (beta[l,i]!=beta[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l+4,1]=rate/(m-1)
  print(rate/(m-1))
}

print('---------------------------------------------------------------------------------')
print('Eta????????????')
for (l in 1:38){
  rate =0
  for (i in 1:(m-1)){
    if (Eta[l,i]!=Eta[l,i+1]){
      rate=rate+1
    }
  }
  
  print(rate/(m-1))
}

print('---------------------------------------------------------------------------------')
print('delta0????????????')
for (l in 1:38){
  rate =0
  for (i in 1:(m-1)){
    if (delta0[l,i]!=delta0[l,i+1]){
      rate=rate+1
    }
  }
  
  print(rate/(m-1))
}

print('---------------------------------------------------------------------------------')

print('u????????????')
for (l in 1:38){
  rate =0
  for (i in 1:(m-1)){
    if (u[l,i]!=u[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l,10]=rate/(m-1)
  print(rate/(m-1))
}

print('---------------------------------------------------------------------------------')

print('z????????????')
for (l in 1:38){
  rate =0
  for (i in 1:(m-1)){
    if (z[l,i]!=z[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l,11]=rate/(m-1)
  print(rate/(m-1))
}

print('---------------------------------------------------------------------------------')

rate =0
for (i in 1:(m-1)){
  if (sigR[1,i]!=sigR[1,i+1]){
    rate=rate+1
  }
}
change_rate[10,1]=rate/(m-1)
print('sigR????????????')
print(rate/(m-1))

print('---------------------------------------------------------------------------------')

rate =0
for (i in 1:(m-1)){
  if (Fi0[1,i]!=Fi0[1,i+1]){
    rate=rate+1
  }
}
change_rate[11,1]=rate/(m-1)
print('Fi0????????????')
print(rate/(m-1))

print('---------------------------------------------------------------------------------')

print('Fi????????????')
for (l in 1:8){
  
  rate =0
  for (i in 1:(m-1)){
    if (Fi[l,i]!=Fi[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l+11,1]=rate/(m-1)
  print(rate/(m-1))
}

print('---------------------------------------------------------------------------------')

for (t in 1:8){
  print('---------------------------------------------------------------------------------')
  print(t)
  print('delta????????????')
  
  for (l in 1:38){
    rate =0
    for (i in 1:(m-1)){
      if (delta[l,i,t]!=delta[l,i+1,t]){
        rate=rate+1
      }
    }
    change_rate[l,1+t]=rate/(m-1)
    print(rate/(m-1))
  }
}

print('---------------------------------------------------------------------------------')




















r=20000
print('---------------------------------------------------------------------------------')
print('alpha????????????')
for (l in 1:4){
  rate =0
  
  for (i in 30000:(m-1)){
    if (alpha[l,i]!=alpha[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l+38,1]=rate/r
  print(rate/r)
}

print('---------------------------------------------------------------------------------')
print('beta????????????')
for (l in 1:5){
  rate =0
  for (i in 30000:(m-1)){
    if (beta[l,i]!=beta[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l+4+38,1]=rate/r
  print(rate/r)
}

print('---------------------------------------------------------------------------------')
print('Eta????????????')
for (l in 1:38){
  rate =0
  for (i in 30000:(m-1)){
    if (Eta[l,i]!=Eta[l,i+1]){
      rate=rate+1
    }
  }
  
  print(rate/r)
}




print('---------------------------------------------------------------------------------')
print('delta0????????????')
for (l in 1:38){
  rate =0
  for (i in 30000:(m-1)){
    if (delta0[l,i]!=delta0[l,i+1]){
      rate=rate+1
    }
  }
  
  print(rate/r)
}

print('---------------------------------------------------------------------------------')

print('u????????????')
for (l in 1:38){
  rate =0
  for (i in 30000:(m-1)){
    if (u[l,i]!=u[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l+38,10]=rate/r
  print(rate/r)
}

print('---------------------------------------------------------------------------------')

print('z????????????')
for (l in 1:38){
  rate =0
  for (i in 30000:(m-1)){
    if (z[l,i]!=z[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l+38,11]=rate/r
  print(rate/r)
}

print('---------------------------------------------------------------------------------')

rate =0
for (i in 30000:(m-1)){
  if (sigR[1,i]!=sigR[1,i+1]){
    rate=rate+1
  }
}
change_rate[10+38,1]=rate/r
print('sigR????????????')
print(rate/r)

print('---------------------------------------------------------------------------------')

rate =0
for (i in 30000:(m-1)){
  if (Fi0[1,i]!=Fi0[1,i+1]){
    rate=rate+1
  }
}
change_rate[11+38,1]=rate/r
print('Fi0????????????')
print(rate/r)

print('---------------------------------------------------------------------------------')

print('Fi????????????')
for (l in 1:8){
  
  rate =0
  for (i in 30000:(m-1)){
    if (Fi[l,i]!=Fi[l,i+1]){
      rate=rate+1
    }
  }
  change_rate[l+11+38,1]=rate/r
  print(rate/r)
}

print('---------------------------------------------------------------------------------')

for (t in 1:8){
  print('---------------------------------------------------------------------------------')
  print(t)
  print('delta????????????')
  
  for (l in 1:38){
    rate =0
    for (i in 30000:(m-1)){
      if (delta[l,i,t]!=delta[l,i+1,t]){
        rate=rate+1
      }
    }
    change_rate[l+38,1+t]=rate/r
    print(rate/r)
  }
}
df <- data.frame(change_rate)
write_xlsx(df, path  = "change rate.xlsx")
end_time <- Sys.time()
print(end_time - start_time)


































  
  
  
  
  #trace plot
  
  for (i in 1:5){
    f=paste("beta",as.character(i-1),'.png', sep = "")
    png( 
      filename = f, # ????????????
      width = 480,height = 480,units = "px",bg = "white",res = 72)
    # 2. 绘图
    plot(beta[i,])
    # 3. ????????????
    dev.off()
  }
  
  
  for (i in 1:4){
    f=paste("alpha",as.character(i-1),'.png', sep = "")
    png( 
      filename = f, # ????????????
      width = 480,height = 480,units = "px",bg = "white",res = 72)
    # 2. 绘图
    plot(alpha[i,])
    # 3. ????????????
    dev.off()
  }
  
  
  png( 
    filename = "sigmaR.png", # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(sigR[,])
  # 3. ????????????
  dev.off()
  png( 
    filename = "sigma0.png", # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(sig0[1,50:50000])
  # 3. ????????????
  dev.off()
  png( 
    filename = "Fi0.png", # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(Fi0[1,])
  # 3. ????????????
  dev.off()
  
  
  for (i in 1:8){
    f=paste("sigma",as.character(i),'.png', sep = "")
    png( 
      filename = f, # ????????????
      width = 480,height = 480,units = "px",bg = "white",res = 72)
    # 2. 绘图
    plot(sig[i,50:50000])
    # 3. ????????????
    dev.off()
  }
  
for (i in 1:8){
  f=paste("Fi",as.character(i),'.png', sep = "")
  png( 
    filename = f, # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(Fi[i,])
   # 3. ????????????
   dev.off()
 }
  
  
  







#trace plot

for (i in 1:5){
  f=paste("30000beta",as.character(i-1),'.png', sep = "")
  png( 
    filename = f, # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(beta[i,30000:50000])
  # 3. ????????????
  dev.off()
}
for (i in 1:4){
  f=paste("30000alpha",as.character(i-1),'.png', sep = "")
  png( 
    filename = f, # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(alpha[i,30000:50000])
  # 3. ????????????
  dev.off()
}
png( 
  filename = "30000sigmaR.png", # ????????????
  width = 480,height = 480,units = "px",bg = "white",res = 72)
# 2. 绘图
plot(sigR[,30000:50000])
# 3. ????????????
dev.off()
png( 
  filename = "30000sigma0.png", # ????????????
  width = 480,height = 480,units = "px",bg = "white",res = 72)
# 2. 绘图
plot(sig0[1,30000:50000])
# 3. ????????????
dev.off()
png( 
  filename = "30000Fi0.png", # ????????????
  width = 480,height = 480,units = "px",bg = "white",res = 72)
# 2. 绘图
plot(Fi0[1,30000:50000])
# 3. ????????????
dev.off()
for (i in 1:8){
  f=paste("30000sigma",as.character(i),'.png', sep = "")
  png( 
    filename = f, # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(sig[i,30000:50000])
  # 3. ????????????
  dev.off()
}
for (i in 1:8){
  f=paste("30000Fi",as.character(i),'.png', sep = "")
  png( 
    filename = f, # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(Fi[i,30000:50000])
  # 3. ????????????
  dev.off()
}




















#acf plot

for (i in 1:5){
  f=paste("acf_bata",as.character(i-1),'.png', sep = "")
  png( 
    filename = f, # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(acf(beta[i,30000:50000],lag.max = 250))
  # 3. ????????????
  dev.off()
}

f=paste("acf_bata",'.png', sep = "")
png( 
  filename = f, # ????????????
  width = 705,height = 832,units = "px",bg = "white",res = 90)
par(mai=c(0.2, 0.5, 0.45,0.2), mfcol=c(5, 1))
(acf(beta[1,30000:50000],lag.max = 250))
(acf(beta[2,30000:50000],lag.max = 250))
(acf(beta[3,30000:50000],lag.max = 250))
(acf(beta[4,30000:50000],lag.max = 250))
(acf(beta[5,30000:50000],lag.max = 250))
dev.off()

for (i in 1:4){
  f=paste("acf_alpha",as.character(i-1),'.png', sep = "")
  png( 
    filename = f, # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(acf(alpha[i,30000:50000],lag.max = 500))
  # 3. ????????????
  dev.off()
}

f=paste("acf_Fi0.png", sep = "")
png( 
  filename = f, # ????????????
  width = 480,height = 480,units = "px",bg = "white",res = 72)
# 2. 绘图
plot(acf(Fi0[,30000:50000],lag.max = 500))
# 3. ????????????
dev.off()
f=paste("acf_sigR.png", sep = "")
png( 
  filename = f, # ????????????
  width = 480,height = 480,units = "px",bg = "white",res = 72)
# 2. 绘图
plot(acf(sigR[,30000:50000],lag.max = 250))
# 3. ????????????
dev.off()
f=paste("acf_sig0.png", sep = "")
png( 
  filename = f, # ????????????
  width = 480,height = 480,units = "px",bg = "white",res = 72)
# 2. 绘图
plot(acf(sig0[,30000:50000],lag.max = 250))
# 3. ????????????
dev.off()

for (i in 1:8){
  f=paste("acf_sigma",as.character(i),'.png', sep = "")
  png( 
    filename = f, # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(acf(sig[i,30000:50000],lag.max = 250))
  # 3. ????????????
  dev.off()
}

for (i in 1:8){
  f=paste("acf_Fi",as.character(i),'.png', sep = "")
  png( 
    filename = f, # ????????????
    width = 480,height = 480,units = "px",bg = "white",res = 72)
  # 2. 绘图
  plot(acf(Fi[i,30000:50000],lag.max = 500))
  # 3. ????????????
  dev.off()
}










#f=paste("acf_Fi.png", sep = "")
#png( 
#  filename = f, # ????????????
#  width = 480,height = 480,units = "px",bg = "white",res = 72)
# 2. 绘图
#plot(acf(Fi[5,30000:50000],lag.max = 1000))
# 3. ????????????
#dev.off()



















parameter_result = array(0, dim=c(28,28))

step=50
print(step)
step_b=400
save_data_steplength200= array(0, dim=c(28,step_b))


for (t in 1:5){
  for (i in 1:step_b){
    save_data_steplength200[t,i]=beta[t,30000+step*i]
  }
}
for (t in 1:4){
  for (i in 1:step_b){
    save_data_steplength200[5+t,i]=alpha[t,30000+step*i]
  }
}
for (t in 1:1){
  for (i in 1:step_b){
    save_data_steplength200[9+t,i]=sigR[t,30000+step*i]
    save_data_steplength200[10+t,i]=sig0[t,30000+step*i]
    save_data_steplength200[19+t,i]=Fi0[t,30000+step*i]
  }
}
for (t in 1:8){
  for (i in 1:step_b){
    save_data_steplength200[11+t,i]=sig[t,30000+step*i]
    save_data_steplength200[20+t,i]=Fi[t,30000+step*i]
  }
}


for (i in 1:28){
  print(mean(save_data_steplength200[i,]))
  parameter_result[i,2]=(mean(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(var(save_data_steplength200[i,]))
  parameter_result[i,3]=(var(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(median(save_data_steplength200[i,]))
  parameter_result[i,4]=(median(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(quantile(save_data_steplength200[i,], probs = c(.025)))
  parameter_result[i,5]=(quantile(save_data_steplength200[i,], probs = c(.025)))
}

for (i in 1:28){
  print(quantile(save_data_steplength200[i,], probs = c(.975)))
  parameter_result[i,6]=(quantile(save_data_steplength200[i,], probs = c(.975)))
}
parameter_result[,1]=''
parameter_result[1,1] ='beta'
parameter_result[6,1] ='alpha'
parameter_result[10,1] ='sigR'
parameter_result[11,1]='sig0'
parameter_result[12,1]='sigt'
parameter_result[20,1]='Fi0'
parameter_result[21,1]='Fi'
parameter_result[,7]  =''
for (i in 1:9){
  if (as.numeric(parameter_result[i,5])*as.numeric(parameter_result[i,6])>0){
    parameter_result[i,7]='*'
  }
}





step=100
print(step)
step_b=200
save_data_steplength200= array(0, dim=c(28,step_b))


for (t in 1:5){
  for (i in 1:step_b){
    save_data_steplength200[t,i]=beta[t,30000+step*i]
  }
}
for (t in 1:4){
  for (i in 1:step_b){
    save_data_steplength200[5+t,i]=alpha[t,30000+step*i]
  }
}
for (t in 1:1){
  for (i in 1:step_b){
    save_data_steplength200[9+t,i]=sigR[t,30000+step*i]
    save_data_steplength200[10+t,i]=sig0[t,30000+step*i]
    save_data_steplength200[19+t,i]=Fi0[t,30000+step*i]
  }
}
for (t in 1:8){
  for (i in 1:step_b){
    save_data_steplength200[11+t,i]=sig[t,30000+step*i]
    save_data_steplength200[20+t,i]=Fi[t,30000+step*i]
  }
}


for (i in 1:28){
  print(mean(save_data_steplength200[i,]))
  parameter_result[i,9]=(mean(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(var(save_data_steplength200[i,]))
  parameter_result[i,10]=(var(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(median(save_data_steplength200[i,]))
  parameter_result[i,11]=(median(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(quantile(save_data_steplength200[i,], probs = c(.025)))
  parameter_result[i,12]=(quantile(save_data_steplength200[i,], probs = c(.025)))
}

for (i in 1:28){
  print(quantile(save_data_steplength200[i,], probs = c(.975)))
  parameter_result[i,13]=(quantile(save_data_steplength200[i,], probs = c(.975)))
}
parameter_result[,8]  =''
parameter_result[1,8] ='beta'
parameter_result[6,8] ='alpha'
parameter_result[10,8] ='sigR'
parameter_result[11,8]='sig0'
parameter_result[12,8]='sigt'
parameter_result[20,8]='Fi0'
parameter_result[21,8]='Fi'
parameter_result[,14] =''
for (i in 1:9){
  if (as.numeric(parameter_result[i,13])*as.numeric(parameter_result[i,12])>0){
    parameter_result[i,14]='*'
  }
}





step=150
print(step)
step_b=133
save_data_steplength200= array(0, dim=c(28,step_b))


for (t in 1:5){
  for (i in 1:step_b){
    save_data_steplength200[t,i]=beta[t,30000+step*i]
  }
}
for (t in 1:4){
  for (i in 1:step_b){
    save_data_steplength200[5+t,i]=alpha[t,30000+step*i]
  }
}
for (t in 1:1){
  for (i in 1:step_b){
    save_data_steplength200[9+t,i]=sigR[t,30000+step*i]
    save_data_steplength200[10+t,i]=sig0[t,30000+step*i]
    save_data_steplength200[19+t,i]=Fi0[t,30000+step*i]
  }
}
for (t in 1:8){
  for (i in 1:step_b){
    save_data_steplength200[11+t,i]=sig[t,30000+step*i]
    save_data_steplength200[20+t,i]=Fi[t,30000+step*i]
  }
}


for (i in 1:28){
  print(mean(save_data_steplength200[i,]))
  parameter_result[i,16]=(mean(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(var(save_data_steplength200[i,]))
  parameter_result[i,17]=(var(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(median(save_data_steplength200[i,]))
  parameter_result[i,18]=(median(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(quantile(save_data_steplength200[i,], probs = c(.025)))
  parameter_result[i,19]=(quantile(save_data_steplength200[i,], probs = c(.025)))
}

for (i in 1:28){
  print(quantile(save_data_steplength200[i,], probs = c(.975)))
  parameter_result[i,20]=(quantile(save_data_steplength200[i,], probs = c(.975)))
}
parameter_result[,15]  =''
parameter_result[1,15] ='beta'
parameter_result[6,15] ='alpha'
parameter_result[10,15] ='sigR'
parameter_result[11,15]='sig0'
parameter_result[12,15]='sigt'
parameter_result[20,15]='Fi0'
parameter_result[21,15]='Fi'
parameter_result[,21] =''
for (i in 1:9){
  if (as.numeric(parameter_result[i,20])*as.numeric(parameter_result[i,19])>0){
    parameter_result[i,21]='*'
  }
}





step=200
print(step)
step_b=100
save_data_steplength200= array(0, dim=c(28,step_b))


for (t in 1:5){
  for (i in 1:step_b){
    save_data_steplength200[t,i]=beta[t,30000+step*i]
  }
}
for (t in 1:4){
  for (i in 1:step_b){
    save_data_steplength200[5+t,i]=alpha[t,30000+step*i]
  }
}
for (t in 1:1){
  for (i in 1:step_b){
    save_data_steplength200[9+t,i]=sigR[t,30000+step*i]
    save_data_steplength200[10+t,i]=sig0[t,30000+step*i]
    save_data_steplength200[19+t,i]=Fi0[t,30000+step*i]
  }
}
for (t in 1:8){
  for (i in 1:step_b){
    save_data_steplength200[11+t,i]=sig[t,30000+step*i]
    save_data_steplength200[20+t,i]=Fi[t,30000+step*i]
  }
}


for (i in 1:28){
  print(mean(save_data_steplength200[i,]))
  parameter_result[i,23]=(mean(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(var(save_data_steplength200[i,]))
  parameter_result[i,24]=(var(save_data_steplength200[i,]))
}

for (i in 1:28){
  print(median(save_data_steplength200[i,]))
  parameter_result[i,25]=(median(save_data_steplength200[i,]))
}
med_beta = parameter_result[1:5,25]
med_alpha= parameter_result[6:9,25]
med_sigR = parameter_result[10,25]
for (i in 1:28){
  print(quantile(save_data_steplength200[i,], probs = c(.025)))
  parameter_result[i,26]=(quantile(save_data_steplength200[i,], probs = c(.025)))
}

for (i in 1:28){
  print(quantile(save_data_steplength200[i,], probs = c(.975)))
  parameter_result[i,27]=(quantile(save_data_steplength200[i,], probs = c(.975)))
}
parameter_result[,22]  =''
parameter_result[1,22] ='beta'
parameter_result[6,22] ='alpha'
parameter_result[10,22] ='sigR'
parameter_result[11,22]='sig0'
parameter_result[12,22]='sigt'
parameter_result[20,22]='Fi0'
parameter_result[21,22]='Fi'
parameter_result[,28] =''
for (i in 1:9){
  if (as.numeric(parameter_result[i,20])*as.numeric(parameter_result[i,19])>0){
    parameter_result[i,28]='*'
  }
}



df <- data.frame(parameter_result)
write_xlsx(df, path  = "參數估計結果.xlsx")





parameter_delta = array(0, dim=c(38,8,20000/200))
parameter_z = array(0, dim=c(38,20000/200))
parameter_u = array(0, dim=c(38,20000/200))

med_delta= array(0, dim=c(38,8))
med_z    = array(0, dim=c(38))
med_u    = array(0, dim=c(38))

for (i in 1:38){
  for (t in 1:8){
    for (r in 1:100){
      parameter_delta[i,t,r]=delta[i,30000+r*200,t]
    }
    med_delta[i,t]=median(parameter_delta[i,t,])
  }
  for (r in 1:100){
    parameter_z[i,r]=z[i,30000+r*200]
    parameter_u[i,r]=u[i,30000+r*200]
  }
  med_z[i]=median(parameter_z[i,])
  med_u[i]=median(parameter_u[i,])
}

Y_hat = array(0, dim=c(8,38))
for (i in 1:38){
  for (t in 1:8){
    Y_hat[t,i]=(1-(1/(exp(-(W[,i]%*%med_alpha+med_sigR*(((1-2*Tow)/(Tow*(1-Tow)))*med_z[i]+(2/(Tow*(1-Tow)))*(med_z[i]^(1/2)*med_u[i]))))/E2+1)))*
      exp(X[,i,t]%*%med_beta+med_delta[i,t])
  }
}

sum((Y_hat-Y)*(Y_hat-Y))

for (i in 1:38){
  print(Y_hat[,i])
  print('----------------')
}






