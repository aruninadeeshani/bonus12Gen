      Program bonus12

CCCC  type = 1 (radiated inclusive inelastics), 2 (radiated elastics), 3 (radiated QE),
CCCC  type = 1, tagged = 1 (radiated neutron + recoil proton from e-d, BONuS process)
      use cl_option
      implicit  none

      real*8 e,ep,theta,epmin,epmax,tmin,tmax,sig,rc,mp,mp2,eel
      real*8 thetap,q2,w2,nu,nuel,q2el,sin2,qv,pp,cosp,radcon/57.2958/
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
      integer*4 numOptions
      real*8 tempDump, p2, bjorkenX, gamma2, theta_pq
      real*8 y0, y_test

      MD = 2*mN + epsD          !!! deuteron mass in MeV/c     
      numOptions = iargc()
      call read_commandPara(numOptions)
      seed = cl_seed 
      outfile = cl_outfile
      nevents = cl_numEvents
      write(*,*) "Input parameter: ", seed, outfile
c      read(5,*) seed, outfile

      open(unit=26,file=outfile,status='new')
      open(unit=29,file='junk.dat', action='write', status='replace')
      write(29, *) "e   e'   q2   nu    w2   xb  yD  pT   p  theta_pq"
      open(unit=27,file='chkPT3.dat', action='write', status='replace')
      write(27, *) "In sampproton.f:  norm  test  gamma   yD   pT   fyt
     &   fyint"
      
c      seed = -180081277
c      nevents = 5
      
      type = cl_type                  !!! radiated inelastics
c      type = 1                  !!! radiated inelastics
c       type = 2                  !!! radiated elastics
      tagged = 1
c       tagged = 0
       

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
       w2 = mp2 + 2*mp*nu-q2 !// only for scattered proton? may not correct for the deuteron.
       
       gamma = sqrt(1.0+q2/(nu*nu))
       bjorkenX = q2/(2.0*MD*nu/1000.0)
c       gamma2 = sqrt(1 + 4*bjorkenX*bjorkenX*MD*MD/1000000.0/q2)
       gamma2 = sqrt(1.0+q2/nu/nu)
       write(*,*) gamma, gamma2, gamma-gamma2
       if (gamma2 .lt. 1.0) then
          write(*,*) "gamma is less than 1"
       endif

c       write(6,*) "test2",nu,q2,gamma
       
       lx = ep*sin(theta/radcon)*cos(phi/radcon)
       ly = ep*sin(theta/radcon)*sin(phi/radcon)
       lz = ep*cos(theta/radcon)                  

CCCC  Next get components of photon q vector  CCCC       
       
       qx = -1.0D0*lx
       qy = -1.0D0*ly
       qv = sqrt(q2+nu*nu) !// qv = qx^2+qy^2+qz^2
       qz = e - lz                          !!!  looks correct to here.
       phiq = phi-180.D0
       cosqz = qz/qv
       sinqz = sqrt(1.0-cosqz*cosqz)
       cosqx = qx/sqrt(qx*qx+qy*qy)
       sinqx = qy/sqrt(qx*qx+qy*qy)
c       phiq = 180.0/pi*acos(qx/sqrt(qx*qx+qy*qy))

c       write(6,*) "check:  ", phiq+180.0, phi
       pp = 0.0
       p2 = 0.0
       if(type.EQ.2) then
         eel = mp*ep/(mp-2.0*ep*sin2) 
         nuel = eel-ep
         q2el = 4.0*eel*ep*sin2
         pp = sqrt(nuel*nuel+q2el)!sqrt(nuel*nuel+q2) 
         px = qx
         py = qy
         pz = eel-lz
         p2 = sqrt(px*px+py*py+pz*pz)
         thetap = radcon*acos(pz/sqrt(pz*pz+px*px+py*py))
         tempDump = abs(pp - p2)
         write (6,*) "type ", type, ": pp = ", pp,
     &      " p2 = ", p2, " diff p = ", tempDump
c         write (29,*) pp, p2, tempDump
       elseif(type.EQ.1) then  
         thetap = 180*ran3(seed) !!! For now assume uniform distribution  !!!
                                  !!!  need to put in a reasonable proton momentum distribution sampling
         check = .true.
         do while(check)
           pp = 0.01+(0.4-0.01)*ran3(seed)   !!! generate only between 1-400 MeV/c  !!!
           testp = 0.01/pp            !!! pretend propability distribution   !!!         
           if(pp.LE.testp) check = .false. 
         enddo
         phip = 360.0*ran3(seed)
         px = pp*sin(thetap/radcon)*cos(phip/radcon)
         py = pp*sin(thetap/radcon)*sin(phip/radcon)
         pz = pp*cos(thetap/radcon)
       endif

       if(type.EQ.1.AND.tagged.EQ.1) then
          call sampproton(gamma,yd,pT)  !!!  pT returned in MeV/c.  Component perpendicular to photon direction  !!!
c          write(*, *) "yd = ", yD, " pT = ", pT
         IF (gamma.GT.1.D0) THEN
	   pgz = ( DSQRT( (1.D0-yD)**2*MD**2              !!! Component of proton momentum along photon direction  !!!
     &           + (gamma**2-1.D0)*(pT**2 + mN**2) )
     &		 - MD*(1.D0-yD)*gamma )
     &           / (gamma**2-1.D0)
         ELSE
           pgz = (pT**2 + mN**2 - (1.D0-yD)**2*MD**2) / (2*(1.D0-yD)*MD)
         ENDIF
        
         pp = sqrt(pT*pT+pgz*pgz)         !!! total proton momentum in MeV/c
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

         px = px/1000.0           !!!  put in GeV/c
         py = py/1000.0
         pz = pz/1000.0
         pp = pp/1000.0   
         
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
c       write(6,*) "writeInFile",-1.0*lt,one,pidp,zero,zero,
c     &      px,py,pz,epro,mp,vx,vy,vz

c       write(29,*) e,ep,theta,thetap,pp

c       *****check*****         
       p2 = sqrt(px*px + py*py + pz*pz)
       tempDump = (px*qx + py*qy + pz*qz)/(p2*qv)
       theta_pq = (180.0/pi)*acos(tempDump)
       write(*,*) theta_pq, px, py, pz

          write(29, *) e, ep, q2, nu, w2, bjorkenX, yD, pT, p2, theta_pq
       
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


     
