      subroutine rc_modb(xin,theta,type,sig)
      implicit none
      character*80 infile
      integer*4 i,j,k,type
      real*8 xin,rce,theta,thspect,thcentdeg,sig,sigb,sigrel,sigrin
      real*8 rc(15),rc_cent
      real*8 radtab_temp(11),thetadeg,thcent,thetalow,thetahigh
      real*8 thetatab,test,p/0.2/
      real*8 mp,mp2,radcon,thetarad
      real*8 xtab,xtab_next
      logical endof,extrap1,extrap2
 
      common/mod/sigmax,eentries,exttab,firstr
      real*8 sigmax,exttab(100000,11)
      integer eentries
      logical firstr

      
      infile = 'bonusrc.dat'

      if(type.EQ.1) then         !!! radiated inelastics  !!!
        k = 11
      elseif(type.EQ.2) then     !!! radiated elastics    !!!
        k = 9
      elseif(type.EQ.3) then     !!! radiated QE          !!!
        k = 10
      endif          
      

      if(firstr) open(unit=34,file=infile,status='old')    

      extrap1 = .false.           !!!  initialize to false       !!!
      extrap2 = .false.           !!!  if true then extrapolate  !!!

      radcon = 180./3.141593
      thetadeg = theta*radcon
      thcentdeg = thspect*radcon

   
      do i=1,15
       rc(i) = 0.
      enddo
   

CCCCCC              read in radcor table              CCCCCC

c      write(6,*)"here",firstr

c      firstr = .true.

      if (firstr) then 
       i = 1
       eentries = 0 
       endof = .false.
       dowhile(.not.endof)
        read(34,*,END=1001) radtab_temp
        do j=1,11
         exttab(i,j) = radtab_temp(j)
        enddo 
        eentries = eentries + 1
        i = i + 1 
       enddo

      endif

 1001 endof = .true.

c      write(6,*) "Nentries in radcor table is:  ",eentries

      close(34) 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    
CCCCCC       Calculate radiative correction and model by doing    CCCCCC
CCCCCC       linear interpolation in theta and xin                CCCCCC



        thetalow = 0.2*int(thetadeg/p)   !!! find angle below !!! 
        thetahigh = thetalow+0.2              !!! find angle above !!!

        
CCCCCC     do search for rcs to interpolate in theta and xin.     CCCCCC 
CCCCCC     thetahigh is the integer theta above the               CCCCCC
CCCCCC     central theta.                                         CCCCCC
 
        do j=1,eentries
         thetatab = exttab(j,3)
         xtab = exttab(j,2)
         xtab_next = exttab(j+1,2)  
         if(abs(thetatab-thetalow).LT.1.E-5) then 
          if(xin.GE.xtab) then
           if(exttab(j-1,3).NE.thetatab) then  !!!  extrapolate  !!!
            extrap1 = .true.          
            rc(1) = exttab(j,k)
            rc(2) = exttab(j+1,k)   
            rc(3) = (rc(2)-rc(1))/(xtab_next-xtab)*(xin-xtab)+rc(1)

c            write(6,*) thetalow,xin,xtab,xtab_next,rc(1),rc(2),rc(3)

           endif
          endif
          if(xin.LE.xtab.and.xin.GE.xtab_next) then !!!  interpolate  !!!
            rc(1) = exttab(j,k)
            rc(2) = exttab(j+1,k)
            rc(3) = ((xin-xtab_next)*rc(1)+(xtab-xin)*rc(2))/
     &         (xtab-xtab_next) 
 
c             write(6,*) xin,xtab,xtab_next,rc(1),rc(2),rc(3)

          endif
       
c           write(6,*) thetalow,xin,rc(1),rc(2),rc(3)    

         endif

         if(abs(thetatab-thetahigh).LT.1.E-5) then

          if(xin.GE.xtab) then
           if(exttab(j-1,3).NE.thetatab) then  !!!  extrapolate  !!!
            extrap2 = .true.
            rc(4) = exttab(j,k)
            rc(5) = exttab(j+1,k)
            rc(6) = (rc(5)-rc(4))/(xtab_next-xtab)*(xin-xtab)+rc(4)

c            write(6,*) thetahigh,xin,xtab,xtab_next,rc(4),rc(5),rc(6)
           endif
          endif

          if(xin.LE.xtab.and.xin.GE.xtab_next) then !!!  interpolate  !!!
           rc(4) = exttab(j,k)
           rc(5) = exttab(j+1,k)
           rc(6) = ((xin-xtab_next)*rc(4)+(xtab-xin)*rc(5))/
     &          (xtab-xtab_next)
          endif

         endif


CCCCCC                          End search                            CCCCCC

   
CCCCCC             Now do interpolation in theta                      CCCCCC

        sig = ((thetadeg-thetahigh)*rc(3) + 
     &      (thetalow-thetadeg)*rc(6))/(thetalow-thetahigh)  
   


c        sig = rce
      enddo 

c      write(6,*) "IN rc_mod:  ",xin,theta,sig,rc(3)
 
 8000 format(a80) 

      return

      end





















