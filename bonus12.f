      Program bonus12

CCCC  type = 1 (radiated inclusive inelastics), 2 (radiated elastics), 3 (radiated QE),
CCCC  type = 1, tagged = 1 (radiated neutron + recoil proton from e-d, BONuS process)
      
      implicit  none

      real*8 e,ep,theta,epmin,epmax,tmin,tmax,sig,rc,mp,mp2,eel
      real*8 thetap,q2,w2,nu,nuel,sin2,qv,pp,cosp,radcon/57.2958/
      real*8 lx,ly,lz,px,py,pz,phi,vx,vy,vz,ran3,lt/-1.0/
      real*8 me/0.000511/,epro,sigint,t1,testp,phip,gamma,pT,yd
      real*8 pgz,costhpq,qx,qy,qz,phiq,pphi,cospphi,sinpphi,pxy,pgzc
      real*8 cosqz,sinqz,cosqx,sinqx
      integer i,type,tagged,seed/280081277/,one,two,zero,pid,pidp,np
      real*8 mN/938.91897/,epsD/-2.224775/,MD,pi/3.14159/
      integer nevents
      logical first/.true./,firstr/.true./,check/.true./
      character*40 outfile
      external ran3

      MD = 2*mN + epsD          !!! deuteron mass in MeV/c     

      read(5,*) seed, outfile

      open(unit=26,file=outfile,status='new')
c      open(unit=29,file='check.dat' ,status='new')
      
c      seed = -180081277
      nevents = 200000
      
      type = 1                  !!! radiated inelastics
c       type = 2                  !!! radiated elastics
c      tagged = 1
       tagged = 0
       

CCCC      First define generation phase space   CCCC

c        tmin = 20.        !!!  min angle in degrees
c        tmax = 25.        !!!  max  "    "    "
c        epmin = 0.8       !!!  min momentum in Gev
c        epmax = 1.2       !!!  max  "       "  "

      e = 2.144
      epmin = 0.4
      epmax = 2.143
      tmin = 5.
      tmax = 30.0
      
     

      mp = 0.938272
      mp2 = mp*mp

      zero = 0
      one = 1
      two = 2
      np = 2
      pid = 11
      pidp = 2212
          
c      write(6,*) "test1"

     
      
      call epgenb(first,e,epmin,epmax,tmin,tmax,type,
     &         seed,ep,theta,phi,sigint,sig)         !!! initialize  !!
      first = .false.

      
      do i=1,nevents
c       write(6,*) i,"test2" 
       call epgenb(first,e,epmin,epmax,tmin,tmax,type,
     &         seed,ep,theta,phi,t1,sig)
       
       sin2 = sin(theta/radcon/2.)
       sin2 = sin2*sin2
       q2 = 4*e*ep*sin2
       nu = e-ep
       w2 = mp2 + 2*mp*nu-q2
       
       gamma = sqrt(1.0+nu*nu/q2)

c       write(6,*) "test2",nu,q2,gamma
       
       lx = ep*sin(theta/radcon)*cos(phi/radcon)
       ly = ep*sin(theta/radcon)*sin(phi/radcon)
       lz = ep*cos(theta/radcon)                  

CCCC  Next get components of photon q vector  CCCC       
       
       qx = -1.0D0*lx
       qy = -1.0D0*ly
       qv = sqrt(q2+nu*nu) 
       qz = e - lz                          !!!  looks correct to here.
       phiq = phi-180.D0
       cosqz = qz/qv
       sinqz = sqrt(1.0-cosqz*cosqz)
       cosqx = qx/sqrt(qx*qx+qy*qy)
       sinqx = qy/sqrt(qx*qx+qy*qy)
c       phiq = 180.0/pi*acos(qx/sqrt(qx*qx+qy*qy))

c       write(6,*) "check:  ", phiq+180.0, phi
       
       if(type.EQ.2) then
         eel = mp*ep/(mp-2.0*ep*sin2)
         nuel = eel-ep
         pp = sqrt(nuel*nuel+q2) 
         px = qx
         py = qy
         pz = eel-lz
         thetap = radcon*acos(pz/sqrt(pz*pz+px*px+py*py))
       elseif(type.EQ.1) then  
         thetap = 180*ran3(seed) !!! For now assume uniform distribution  !!!
                                  !!!  need to put in a reasonable proton momentum distribution sampling
         check = .true.
         do while(check)
           pp = 0.01+0.4*ran3(seed)   !!! generate only between 1-400 MeV/c  !!!
           testp = 0.01/pp            !!! pretend propability distribution   !!!         
           if(pp.LE.testp) check = .false. 
         enddo
         phip = 360.0*ran3(seed)
         px = pp*sin(thetap/radcon)*cos(phip/radcon)
         py = pp*sin(thetap/radcon)*sin(phip/radcon)
         pz = pp*cos(thetap/radcon)
       endif

       if(type.EQ.1.AND.tagged.EQ.1) then
c         write(6,*) "test3"
          call sampproton(gamma,yd,pT)  !!!  pT returned in MeV/c.  Component perpendicular to photon direction  !!!
         IF (gamma.GT.1.D0) THEN
	   pgz = ( DSQRT( (1.D0-yD)**2*MD**2              !!! Component of proton momentum along photon direction  !!!
     &           + (gamma**2-1.D0)*(pT**2 + mN**2) )
     &		 - MD*(1.D0-yD)*gamma )
     &           / (gamma**2-1.D0)
         ELSE
           pgz = (pT**2 + mN**2 - (1.D0-yD)**2*MD**2) / (2*(1.D0-yD)*MD)
         ENDIF
        
         pp = sqrt(pT*pT+pgz*pgz)         !!! total proton momentum
         costhpq = pgz/pp                !!! cosine of angle between proton and qvec  !!!

         pphi = 360.0*ran3(seed) !!! angle between pT and x-y plane, rotating about qv 
         cospphi = cos(pphi/radcon)
         sinpphi = sin(pphi/radcon)
         
c         write(6,*) "test  ", pphi,cospphi

         pz = pgz*cosqz+pT*sinpphi*sinqz !!! proton momentum along z-axis (Looks ok)
         pxy = sqrt(pp*pp-pz*pz)  !!!  magnitude of proton momentum in x-y plane    

         px = pT*(cospphi*sinqx-sinpphi*cosqz*cosqx) + pgz*qx/qv
         
         py = -1.0*pT*(cospphi*cosqx+sinpphi*cosqz*sinqx)+pgz*qy/qv     
         
c         px = pxy*cospphi                  !!! proton momentum along x-axis
c         py = sqrt(pp*pp-px*px-pz*pz)      !!! proton momentum along y-axiy
         
         pgzc = px*qx/qv+py*qy/qv+pz*qz/qv  !!! check of proton momentum along qv
         
c         write(6,1003) e,ep,theta,pgz,pgzc,pp,pxy,pz,sqrt(px*px+py*py)
       endif 
      
       vx = 0.0           !!!  Vertex position   !!!
       vy = 0.0
       vz = 20.0*(1.0-2.0*ran3(seed))
       epro = sqrt(mp2+pp*pp)   !!!  Proton total energy
 
       write(26,1002) np, one,one, zero,zero,pid,e,pidp,zero,one
       write(26,1001) one,lt,one,pid,zero,zero,
     &      lx,ly,lz,ep,me,vx,vy,vz

       write(26,1001) two,-1.0*lt,one,pidp,zero,zero,
     &      px,py,pz,epro,mp,vx,vy,vz

c       write(29,*) e,ep,theta,thetap,pp
       
c       write(6,*) sig
       
c       if(abs(theta-5.5).LT.0.2) write(6,*) i,theta,eel,ep,thetap,phi
       
c       if(w2.GE.0.8.AND.ep.LT.0.6) write(6,*) i,e,ep,theta,pp


       
c       pp = 

c     write(6,1001) i,e,ep,theta,w2,q2,pp,cosp,sig
       write(6,*) i,nevents
      enddo
      write(6,*) "total cross section: ", sigint

      
 1001 format(1i5,1f7.2,4i6,8f10.4)
 1002 format(6i5,1f7.4,4i6,3i5)    
 1003 format(9f10.4)
      
      end


     
