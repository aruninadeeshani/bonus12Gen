      Subroutine SAMPPROTON(gamma,yd,pT)

      IMPLICIT NONE
      real*8 e,ep,theta,epmin,epmax,tmin,tmax,sig,rc,mp,mp2
      real*8 q2,w2,qv,pp,cosp,radcon/57.2958/
      real*8 xd,yd,y,pt,pt2,gamma,fyint(4),fyoff(4),prob,ymax,ptmax
      real*8 ptmin,fymax,norm,fyt,test,ran3
      INTEGER wfn,ibar,imom,rel,i,j,seed/84389723/
      logical firstr
      external ran3
      integer*4 numCount
      
CCCC  By default use WBArelativistic prescription with CDBONN potential  CCCC      
      
      wfn = 3      !!! 3 = CDBONN, 4 = AV18 
      ibar = 0     !!! Baryon # conservation.  Leave at 0  !!!
      rel = 1      !!! 1 = Relativistic WBA
      
     
      mp = 0.938272
      mp2 = mp*mp

CCCC  First normalize distribution for given gamma   CCCC
CCCC  For the moment assume the probability peaks at CCCC
CCCC  yD = 0.5, PT = 2 MeV.
      

c      write(6,*) "Here - sampproton"

      prob = 0.0
      yd = 0.0
      fymax = 0.0
      ymax = 1.0
      pTmax = 1000.0
      pTmin = 1.0
      numCount = 0
      pT = 26.0
      yd = 0.5
    
      call dn(wfn,ibar,rel,yd,gamma,pT,FYINT,FYOFF)
      norm = FYINT(3)*pT
      fyt = 0.0
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dowhile(fyt.LT.test)
        test = ran3(seed)
        yd = 5.0E-4+(5.0-1.0E-4)*ran3(seed)
        pT = pTmin+(pTmax-pTmin)*ran3(seed)
        call dn(wfn,ibar,rel,yd,gamma,pT,FYINT,FYOFF)
        fyt = (FYINT(3)*pT)/norm
c       write(*,*) numCount, (FYINT(3)*pT), norm, test, fyt
        numCount = numCount + 1
      enddo
      write(27, *) norm, test, gamma, yD, pT, fyt, FYINT(3)

      write(*,*) "********* num Counts = ", numCount
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
c      if(fyt.GT.1.0) write(30,*) gamma,yd,pT,fyt,test
c      write(6,*) gamma,yd,pT,fyt,test
      
 1001   format(1i8,5f11.4,2f11.4)  

      end


     
