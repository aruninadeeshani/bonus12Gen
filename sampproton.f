      Subroutine SAMPPROTON(gamma,yd,pT)

      IMPLICIT NONE
      real*8 e,ep,theta,epmin,epmax,tmin,tmax,sig,rc,mp,mp2
      real*8 q2,w2,qv,pp,cosp,radcon/57.2958/
      real*8 xd,yd,y,pt,pt2,gamma,fyint(4),fyoff(4),prob,ymax,ptmax
      real*8 fymax,norm,fyt,test,ran3
      INTEGER wfn,ibar,imom,rel,i,j,seed/84389723/
      logical firstr
      external ran3
      
CCCC  By default use WBArelativistic prescription with CDBONN potential  CCCC      
      
      wfn = 3      !!! 3 = CDBONN, 4 = AV18 
      ibar = 0     !!! Baryon # conservation.  Leave at 0  !!!
      rel = 1      !!! 1 = Relativistic WBA
      
     
      mp = 0.938272
      mp2 = mp*mp

CCCC  First normalize distribution for given gamma   CCCC
CCCC  For the moment assume the probability peaks at CCCC
CCCC  yD = 0.5, PT = 4 MeV.
      

c      write(6,*) "Here - sampproton"

      prob = 0.0
      yd = 0.0
      fymax = 0.0
      ymax = 1.0
      pTmax = 1000
       
      yd = 0.5
      pT = 4.0
      
      call dn(wfn,ibar,rel,yd,gamma,pT,FYINT,FYOFF)
      norm = FYINT(3)
      fyt = 0.0
      test = ran3(seed)
      
      dowhile(fyt.LT.test)
       yd = 1.0*ran3(seed)
       pT = pTmax*ran3(seed)
       call dn(wfn,ibar,rel,yd,gamma,pT,FYINT,FYOFF)
       fyt = FYINT(3)/norm
      enddo
      
c      write(6,*) gamma,yd,pT,fyt,test
            
      
 1001   format(1i8,5f11.4,2f11.4)  

      end


     
