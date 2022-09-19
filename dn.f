C *********************************************************************
      SUBROUTINE DN (wfn,ibar,rel,yD,gamma,pT,FYINT,FYOFF)
C
C     Integrand of WBA and AQV nucleon momentum distributions
C     (for nonrelativistic and relativistic kinematics).
C     
C     Includes also the p^2 weigted (e.g., for the KP off-shell model)
C     
C     Expressions valid only in the deuteron rest frame.
*
*
*     INPUT:
*
*     wfn      = (i) what wave function, see below for explanation
*     ibar     = (i) Baryon number normalization: 0=M/p0  1=const
*                    Note: inly for AQV 
*     rel      = (i) relativ. corrections 0=WBA  1=WBAREL  2=AQV
*     yD       = (dp) nucleon fractional momentum
*                     NOTE: yD=p.q/PD.q with PD the unrescaled Deuteron momentum
*     gamma    = (dp) gamma factor (NOTE: not gamma^2)
*     pT [MeV] = (dp) nucleon transverse momentum (unrescaled)
*
*     OUTPUT:
* 
*     FYINT    = (dp) Integrand for the smearing function
*     FYOFF    = (dp) (p^2-mN^2)/mN^2 * FYINT
*
C *********************************************************************
      IMPLICIT NONE
      INTEGER wfn,ibar,imom,i
      REAL*8 yD,gamma,pT,FYINT(4),FYOFF(4)
      REAL*8 yDmax,pz,pv,pv_cut,Es,p0,p2,eps
     &     ,Jac,flux,rhoN,rhoNT,corr,bnorm
      REAL*8 U,W,VS,VT,rho,CC,prefact
      REAL*8 pi,hc,mN,MD,epsD
      integer rel
      COMMON /con/ pi,hc
      COMMON /mas/ mN,MD,epsD
      integer*8 ipv
      common/flag/ ipv
      
C...Constants
      pi = 4*DATAN(1.D0)
      hc = 197.327D0 		! conversion factor (MeV.fm)
      mN = 938.91897D0          ! nucleon mass
      epsD = -2.224575D0	! deuteron binding energy
      MD = 2*mN + epsD          ! deuteron mass

      pv_cut = 1000.0        ! pT<pv_cut)

      
      do i = 1, 4
         FYINT(i) = 0.D0
         FYOFF(i) = 0.D0
      end do

C...Kinematics...........................................................

      IF (rel.eq.0) THEN	! non-relativistic (sqrt in pz > 0)

         yDmax = (1.D0+gamma**2/2.0 + epsD/mN - pT**2/
     &       (2.0*mN**2))*(mN/MD)
         pz = mN * ( gamma - DSQRT(2*MD/mN*(yDmax-yD)) )

         
      ELSE                      ! relativistic 

         yDmax = 1.0d10            ! proxy for yDmax = \infty
         IF (gamma.GT.1.D0) THEN
	    pz = ( DSQRT( (1.D0-yD)**2*MD**2
     &           + (gamma**2-1.D0)*(pT**2 + mN**2) )
     &		 - MD*(1.D0-yD)*gamma )
     &           / (gamma**2-1.D0)
         ELSE
	    pz = (pT**2 + mN**2 - (1.D0-yD)**2*MD**2) / (2*(1.D0-yD)*MD)
         ENDIF

c         write(6,*) pT,mN,yD,MD,pz    !!! ok to here
         
      ENDIF


      IF (yDmax.LE.0.D0 .OR. yD.GE.(yDmax-0.000001)) RETURN ! to avoid the singular point

c        write(6,*) "pv = ",pv
      
      
C...Total nucleon 3-momentum
	pv = DSQRT(pT**2 + pz**2)
C...Cut-off integration at pv=pv_cut
C...Added command line flag to vary pv_cut values--SL, 05.2020


c        write(6,*) pv,pv_cut
        
	IF (pv.GT.pv_cut) RETURN

C...On-shell spectator nucleon energy
        if (rel.eq.0) then
           Es = mN + pv**2/(2*mN)
        else 
           Es = DSQRT(mN**2 + pv**2)
        end if

C...Separation energy
        eps = MD - Es - mN  ! valid for both rel and non-rel options
        
C...Off-shell interacting nucleon energy
	p0 = MD - Es
!	IF (p0.LE.0.D0) RETURN		!!! for MST only !!!

C...Nucleon virtuality
	p2 = p0**2 - pv**2


c        write(6,*) eps,p0,p2
        
!C...light-cone variable --- NO LONGER NEEDED [AA: 25 June 2020]
!C...solution of y = y0/2 (1 + gamma - (gamma-1)/y0^2 (p2+pT2)/M^2
!	IF (p2+pT**2.LT.-y**2*MD**2/(gamma**2-1.D0)) RETURN ! checks sign of sqrt argument below
!	y0 = y				!!! corrected 4/21/11
!     &	   * ( 1.D0
!     &	     + DSQRT( 1.D0 + (gamma**2-1.D0)*(p2+pT**2)/MD**2/y**2 ) )
!     &     / (1.D0+gamma) 


C...Deuteron wavefunction................................................

!...  Calculates (u^2+v^2+w^2) = |psi|^2 / (4\pi)
!...(--> the 1/(4\pi) term is included in the Jacobian below)  
!...NR:  1=Paris, 3=CDBonn, 4=AV18
!...Rel: 12=WJC-1, 13=WJC-2
!...[high-precision ones are: CDBonn, AV18, WJC-1, WJC-2]
        
      U = 0.D0                  ! S-wave
      W = 0.D0                  ! D-wave
      VS = 0.D0                 ! Singlet P-wave
      VT = 0.D0                 ! Triplet P-wave

      ! momenta in MeV

      !relativistic wave functions
c      if (wfn.EQ.12) then
c         imom = 0                       ! wfns normalized to ~105%
c         CALL WJC (1,pv/hc,U,W,VS,VT)   ! -- 5% in V' term
c      else IF (wfn.EQ.13) then 
c         imom = 0                       ! wfns normalized to ~102%
c         CALL WJC (2,pv/hc,U,W,VS,VT)   ! -- 2% in V' term 

      !Nonrelativistic wave functions (VT=VS=0)
c      else IF (wfn.EQ.1) then 
c         imom = 0
c         CALL PARIS (1,pv/hc,U,W)
      IF (wfn.EQ.3) then 
         imom = 0
         CALL CDBONN (pv/hc,U,W)
      else IF (wfn.EQ.4) then 
         imom = 1
         CALL AV18 (pv/hc,rho)
c         write(6,*) "TEST:  ", rho
      end if

      !Output wavefunctions in fm^3/2 => MeV^-3/2
      if (imom.eq.0) then   
         U = U / hc**1.5D0
         W = W / hc**1.5D0
         VS = VS / hc**1.5D0
         VT = VT / hc**1.5D0
         CC = (U**2 + W**2 + VS**2 + VT**2)
      else                   !Wave fns in terms of momentum distributions
         CC = rho / hc**3.D0
      end if


C...Nucleon distributions in yD=p.q/pA.q 
C...(building blocks for smearing functions)

C    ...flux factors and finite Q^2 corrections
      IF (rel.eq.2) THEN        ! AQV
C       ...Baryon number normalization factor
         if (ibar.eq.0) then    
            bnorm = mN/p0 ! -> for the future: onsider mN/(p0+gamma*pz) = 1/y, instead
         else                   ! taken care by normalizing for gamma=1
            bnorm = 1.0d0         ! in the calling routine
         end if
         flux = yD*MD/mN * bnorm 
         corr = (gamma**2-1.D0)/(yD*MD/mN)**2
     &	     * ( (1.D0 + eps/mN)**2 + (pv**2-3*pz**2)/2.0d0/mN**2 )
         ! i.e., corr = (gamma^2-1)/y^2 * (2p^2+3pT^2)/(2mN^2)
      elseif (rel.eq.1) then   ! WBAREL
         if (ibar.ne.0) then
            write(*,*) 'ERROR(SMEARFNS): ibar out of bounds for rel=1:'
     &           , ibar
            stop
         endif
         ! bnorm=mN/p0 already included in flux factor
         flux = (1.D0 + gamma*pz/mN)    ! neglects O(eps_D/mN) compared to AQV
         corr = (gamma**2-1.D0)/(yD*MD/mN)**2
     &          * ( (1.D0 + eps/mN)**2 + (pv**2-3*pz**2)/2d0/mN**2 )
         ! corr is same as AQV 
      ELSE                      ! WBA
         if (ibar.ne.0) then
            write(*,*) 'ERROR(SMEARFNS): ibar out of bounds for rel=0:'
     &           , ibar
            stop
         end if
         flux = (1.D0 + gamma*pz/mN)  ! same as in WBAREL
         corr = (gamma**2-1.D0)/(yD*MD/mN)**2
     &        * (1.D0 + 2.0*eps/mN + (pv**2-3.0*pz**2)/(2.0*mN**2))
         ! corr is same as WBAREL, but expanded to O(eps/mN)=O(pvec^2/mN^2)
      ENDIF

C...  "Jacobian" .....................................
!...  d3p -> dy*dpT2 & normalization of |u^2+v^2| read in from file
      if (rel.eq.0) then   ! WBA
C         prefact = gamma * mN*MD / 4.D0/ (gamma**2*mN + MD*(1-yD) - Es ) !: AA
          prefact = MD/4.0/(gamma-pz/mN )!  equivalent to formula above

      else  ! relativistic (WBAREL and AQV)
         prefact =  gamma*Es / (44.0*((1.0-yD)+(gamma**2-1.0)*Es/MD))
      endif

C...  FINAL CALCULATION


C...On-shell
      FYINT(1) = prefact        ! diagonal for 2xF1, xFL
     &         * flux * CC
      FYINT(2) = prefact        ! offdiagonal for 2xF1, xFL
     &	       * flux * CC
     &         * (gamma**2-1d0)/(yD*MD/mN)**2 * (pT/mN)**2
      FYINT(3) = prefact        ! F2 smearing function (diagonal for F2)
     &         * flux * CC
     &         * (1.D0 + corr)/gamma**2  
      FYINT(4) = prefact	! xF3 smearing function (diagonal for xF3)
     &         * (1.D0 + pz/mN/gamma) * CC
               ! ^^^^^^^^^^^^^^^^^^^^ CHECK implementation of bnorm here!! [AA 25 June 2020]


c      write(6,*) "In DN: ",gamma,yd,pT,fyint(1),prefact,bnorm,eps,es,Md
      
C...Off-shell weighted part
      do i = 1, 4
         FYOFF(i) = FYINT(i) * (p2-mN**2)/mN**2
      end do
         
      RETURN
      END
