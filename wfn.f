C ***********************************************************************
	SUBROUTINE AV18 (p,rho)
C
C  Deuteron momentum distribution rho = u^2 + w^2, where u and w are
C  S- and D-state wave functions in momentum space, normalized s.t.
C  int dp p^2 rho(p) = 1.
C  Ref: Fantoni and Pandharipande, Nucl. Phys. A427 (1984) 473
C
C  Input file has rho normalized s.t. 4*pi int dp p^2 rho(p) = 1,
C  so that for output, rho needs to be multiplied by 4*pi.
C
C  p in 1/fm, rho in 1/fm^3
C
C.. Uses IMSL interpolation routine DQDVAL.
C.. For compilation on jlabs1:
C.. > use imsl
C.. > f77 -o objectfile file.f -R/site/vni/lib/lib.solaris
C..       -L/site/vni/lib/lib.solaris -limsl -lsocket -lnsl
C.. IMSL decommissioned 8/30/11
C
C ***********************************************************************
        IMPLICIT NONE
        INTEGER	ip,np
        PARAMETER (np=200)
        REAL*8  p,rho
        REAL*8  Parr(np),RHOarr(np),rho_int,dum
        LOGICAL readin /.FALSE./
	REAL*8	pi
	SAVE

C...Value of pi
        pi = 4*DATAN(1.D0)

	rho = 0.D0

        IF (readin) GO TO 123
C...Read data from file
c        OPEN (10,FILE='/u/home/wmelnitc/Work/EMC/D/Wfn/av18.dat',
        OPEN (10,FILE='av18.dat',
     &		FORM='FORMATTED')
	DO ip=1,9
          READ (10,*)
	ENDDO
        DO ip=1,np
	  READ (10,*) Parr(ip), RHOarr(ip)
        ENDDO
        CLOSE (10)
        readin = .TRUE.
c        print *, '... AV18 data read ...'

 123	IF (p.LE.Parr(1)) rho = 4*pi * RHOarr(1)
	IF (p.GT.Parr(1) .AND. p.LE.Parr(np)) THEN
	  CALL Pinterp (Parr,RHOarr,np,p,rho_int,dum,2)
	  rho = 4*pi * rho_int
c     &    rho = 4*pi * DQDVAL(p,np,Parr,RHOarr,.FALSE.)
	ENDIF

c	print *, 'Parr(1),Parr(np)=',Parr(1),Parr(np)
c	print *, 'p,rho=',P, RHO

        RETURN
        END


C **********************************************************************
	SUBROUTINE CDBONN (q,u0,u2)
C
C  Deuteron wave function from CD-Bonn NN potential model.
C  q in 1/fm, u0,u2 in fm^3/2.
C
C  Normalization \int dq q^2 (u0^2+u2^2) = 1.
C
C  Sent by Charlotte Elster, April 8, 2009.
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=95)
        REAL*8  q,u0,u2,
     &		qgrid(nq),uqgrid(nq),wqgrid(nq),weight(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hc
	SAVE

        pi = 4*DATAN(1.D0)
        hc = 197.327D0		! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
c     &		FILE='/u/home/wmelnitc/Work/EMC/D/Wfn/cdbn.qwave',
     &	       FILE='cdbn.qwave',
     &	       STATUS='OLD')
c        READ (10,100)
        READ (10,*)
        READ (10,*)
C...Momentum space [qgrid in MeV, uqgrid in MeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id), weight(id), uqgrid(id), wqgrid(id)
c          write(*,*) qgrid(id), weight(id), uqgrid(id), wqgrid(id)
          qgrid(id) = qgrid(id) / hc            ! MeV => 1/fm
          uqgrid(id) = uqgrid(id) * hc**1.5D0   ! MeV^-3/2 => fm^3/2
          wqgrid(id) = wqgrid(id) * hc**1.5D0
        ENDDO
 100    FORMAT (2(1X/))
 101    FORMAT (3X,D13.6,3X,D20.13,3X,D20.13)
c	PRINT *, '... CD-Bonn model read...'
	init = .TRUE.

C...Evaluate wave function
c 999	u0 = DQDVAL (q,nq,qgrid,uqgrid,.FALSE.)
c	u2 = DQDVAL (q,nq,qgrid,wqgrid,.FALSE.)
 999	CALL Pinterp (qgrid,uqgrid,nq,q,u0,dum,2)
	CALL Pinterp (qgrid,wqgrid,nq,q,u2,dum,2)

        RETURN
        END


C **********************************************************************
	SUBROUTINE WJC1 (q,u,w,vt,vs)
C
C  Deuteron wavefunction from WJC-1 (Gross et al.) NN potential model.
C  Gross & Stadler, Phys. Rev. C78, 014005 (2008).
C
C  Note: q in 1/fm, wfns in fm^3/2.
C
C  Wave function normalization
C    \int dq q^2 (u^2+w^2+vt^2+vs^2) + V' term = 1
C    so that wave functions themselves normalized to ~105%
C    => renormalize to 1 for structure functions
C
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=60)
        REAL*8  q,u,w,vt,vs,DQDVAL,
     &		qgrid(nq),ugrid(nq),wgrid(nq),vtgrid(nq),vsgrid(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hcM,hcG
	SAVE

        pi = 4*DATAN(1.D0)
        hcM = 197.327D0		! GeV.fm conversion factor
        hcG = 0.197327D0	! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
c     &	      FILE='/u/home/wmelnitc/Work/EMC/D/Wfn/wjc-1.dat',
     &	      FILE='wjc-1.dat',
     &	      STATUS='OLD')

C...Momentum space [qgrid in MeV, ugrid in GeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id),
     &		      ugrid(id), wgrid(id), vtgrid(id), vsgrid(id)
          qgrid(id) = qgrid(id) / hcM		! MeV => 1/fm
          ugrid(id)  = ugrid(id)  * hcG**1.5D0	! GeV^-3/2 => fm^3/2
          wgrid(id)  = wgrid(id)  * hcG**1.5D0
          vtgrid(id) = vtgrid(id) * hcG**1.5D0
          vsgrid(id) = vsgrid(id) * hcG**1.5D0
        ENDDO
c	PRINT *, '... WJC-1 model read...'
	init = .TRUE.

C...Evaluate wavefunction
c 999	u  = DQDVAL (q,nq,qgrid,ugrid,.FALSE.)
c	w  = DQDVAL (q,nq,qgrid,wgrid,.FALSE.)
c	vt = DQDVAL (q,nq,qgrid,vtgrid,.FALSE.)
c	vs = DQDVAL (q,nq,qgrid,vsgrid,.FALSE.)
 999	CALL Pinterp (qgrid,ugrid,nq,q,u,dum,1)
	CALL Pinterp (qgrid,wgrid,nq,q,w,dum,1)
	CALL Pinterp (qgrid,vtgrid,nq,q,vt,dum,1)
	CALL Pinterp (qgrid,vsgrid,nq,q,vs,dum,1)

        RETURN
        END


C **********************************************************************
	SUBROUTINE WJC2 (q,u,w,vt,vs)
C
C  Deuteron wavefunction from WJC-2 (Gross et al.) NN potential model.
C  Gross & Stadler, Phys. Rev. C78, 014005 (2008).
C
C  Note: q in 1/fm, wfns in fm^3/2.
C
C  Wave function normalization
C    \int dq q^2 (u^2+w^2+vt^2+vs^2) + V' term = 1
C    so that wave functions themselves normalized to ~102%
C    => renormalize to 1 for structure functions
C
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=60)
        REAL*8  q,u,w,vt,vs,
     &		qgrid(nq),ugrid(nq),wgrid(nq),vtgrid(nq),vsgrid(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hcM,hcG
	SAVE

        pi = 4*DATAN(1.D0)
        hcM = 197.327D0		! GeV.fm conversion factor
        hcG = 0.197327D0	! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
     &	      FILE='wjc-2.dat',
     &	      STATUS='OLD')

C...Momentum space [qgrid in MeV, ugrid in GeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id),
     &		      ugrid(id), wgrid(id), vtgrid(id), vsgrid(id)
          qgrid(id) = qgrid(id) / hcM		! MeV => 1/fm
          ugrid(id)  = ugrid(id)  * hcG**1.5D0	! GeV^-3/2 => fm^3/2
          wgrid(id)  = wgrid(id)  * hcG**1.5D0
          vtgrid(id) = vtgrid(id) * hcG**1.5D0
          vsgrid(id) = vsgrid(id) * hcG**1.5D0
        ENDDO
c	PRINT *, '... WJC-2 model read...'
	init = .TRUE.

C...Evaluate wavefunction
c 999	u  = DQDVAL (q,nq,qgrid,ugrid,.FALSE.)
c	w  = DQDVAL (q,nq,qgrid,wgrid,.FALSE.)
c	vt = DQDVAL (q,nq,qgrid,vtgrid,.FALSE.)
c	vs = DQDVAL (q,nq,qgrid,vsgrid,.FALSE.)
 999	CALL Pinterp (qgrid,ugrid,nq,q,u,dum,1)
	CALL Pinterp (qgrid,wgrid,nq,q,w,dum,1)
	CALL Pinterp (qgrid,vtgrid,nq,q,vt,dum,1)
	CALL Pinterp (qgrid,vsgrid,nq,q,vs,dum,1)

        RETURN
        END


