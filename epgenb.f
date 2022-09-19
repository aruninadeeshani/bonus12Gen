      SUBROUTINE epgenb(first,e,epmin,epmax,tmin,tmax,type,seed,ep,th,
     &                                                 phi,sigint,sig)

C+______________________________________________________________________________
C-______________________________________________________________________________

CCCC   Version 1.1.  Written by M.E. Christy on 05/08/21      CCCC
CCCC   Samples inclusive e-p cross section using radiated     CCCC
CCCC   inelastic model from M.E. Christy.                     CCCC 
CCCC                                                          CCCC
CCCC   Requires input epmin,epmax,tmin,tmax,first             CCCC
CCCC   and returns sampled kinematics ep, theta               CCCC
CCCC   and the cross section, sig.                            CCCC
CCCC   tmin (tmax) is the min (max) angle in degrees          CCCC
CCCC                                                          CCCC
CCCC   Returns sampled ep, theta, and sig                     CCCC
CCCC                                                          CCCC
CCCC   Assumes cross section model in units of ub/Sr/Gev.     CCCC
CCCC   type = 1 (radiated inelastics), 2 (radiated elastics), CCCC
CCCC   3 (radiated QE)                                        CCCC
      
      implicit none

      real*8 e,epmin,epmax,tmin,tmax,ep,th,phi,sig
      real*8 pi,domega,ddcos,cost,sigtest,sigint
      real*8 ran3,rc,radcon,t1,mp,mp2,q2,w2
      integer i,j,nmax,seed,type
      logical samp,first
      common/mod/sigmax,eentries,exttab,firstr
      real*8 sigmax,exttab(100000,11)
      integer eentries
      logical firstr
       

c       if(first) then
c         call rc_modb(ep,th,type,sig)
c         write(6,*) "First: ", first
c       endif
c       first = .false.
          
        nmax = 200000
        pi = 3.14159
c        seed = -45647991

        mp = 0.938272
        mp2 = mp*mp
        radcon = 180./pi

        t1 = ran3(seed)

        ddcos =  cos(tmin/radcon)-cos(tmax/radcon)
        domega = 2.*pi*ddcos !!! Assume full phi coverage

        
CCCC

        if(first) then     !!!  Calculate integrated cross section
                           !!!  over generation limits using mean 
                           !!!  value theorem and find maximum.
         firstr=.true.     !!!  initialize radiated sigma routine
         sigint = 0.0
         sigmax = 0.0

      
c!$OMP PARALLEL PRIVATE (e)
c!$OMP DO 

         
         do i=1,nmax
c          write(6,*) i,sigint  
          ep = epmin + ran3(seed)*(epmax-epmin) !!! random sampling of E' 
          cost = cos(tmin/radcon)-ran3(seed)*ddcos  !!!  random sampling of cos(theta)
          th = radcon*acos(cost)
c          if(i.GT.1 ) first = .false.
          call rc_modb(ep,th/radcon,type,sig)

c          write(6,*) i,ep,theta,firstr,sig,sigmax
          
          firstr = .false.   
          sigint = sigint+sig
c          write(6,*) cos(tmin/radcon),cos(tmax/radcon),acos(cost)            

CCCC        find largest cross section in generation limits

          if(sig.GT.sigmax) sigmax = sig


          
       enddo

c!$OMP END DO 
c!$OMP END PARALLEL


       
         sigint = sigint/nmax*domega*(epmax-epmin)  !!! integrated cross section

        endif              !!!  Only on first call  


CCCC   Now sample the cross section distribution for particular event CCCC
        
        samp = .true.

c        do i=1,nmax
        dowhile(samp)     !!! keep sampling until condition met

          cost = cos(tmin/radcon)-ran3(seed)*ddcos  !!!  random sampling of cos(theta)
          ep = epmin + ran3(seed)*(epmax-epmin) !!! random sampling of E' 
          th = acos(cost)*radcon
          q2 = 4.*e*ep*(sin(th/radcon/2.))**2.0
          w2 = mp2+2.*mp*(e-ep)-q2
c          if(q2.GE.0..and.w2.GE.0) samp = .false.
          call rc_modb(ep,th/radcon,type,sig)
          sigtest = ran3(seed)*sigmax
           
          if(sigtest.LE.sig) samp = .false.
          
        enddo

        phi = 360.0*ran3(seed)

c        write(6,*) "epgen :",phi
        
    
        write(28,1001) ep,th,w2,q2,sig
        
        
 1001   format(5f10.3)

        return
        
        end

********************************************************************************

      function ran3(idum)

*     On first call idum should be less than zero

      INTEGER idum
      INTEGER MBIG,MSEED,MZ
      REAL*8 ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END

*******************************************************************************
