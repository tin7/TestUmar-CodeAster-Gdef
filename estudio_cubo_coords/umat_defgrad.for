C===========================================================================
C  UMAT Neo-Hooke compresible (autocontenida, sin fichero externo)
C
C  Psi = (mu/2)*(Ibar1 - 3) + (kappa/2)*(J - 1)^2
C
C  PROPS(1) = E, PROPS(2) = nu
C     mu = E/(2*(1+nu)),  kappa = E/(3*(1-2*nu))
C
C  F se calcula internamente como F = I + sym(STRAN + DSTRAN)
C  (aproximacion valida para DEFORMATION='PETIT')
C
C  Cauchy:  sigma = (mu/J)*dev(bbar) + kappa*(J-1)*I
C           bbar  = J^{-2/3}*b,  b = F*F^T
C
C  STATEV(1:3) = coords GP actuales
C  STATEV(4:6) = coords GP anteriores
C===========================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     3 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     4 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     5 KSPT,KSTEP,KINC)

      IMPLICIT NONE

      CHARACTER*80 CMNAME
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,
     1        KSPT,KSTEP,KINC

      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),SSE,SPD,SCD,RPL,DDSDDT(NTENS),
     2 DRPLDE(NTENS),DRPLDT,STRAN(NTENS),DSTRAN(NTENS),
     3 TIME(2),DTIME,TEMP,DTEMP,PREDEF(*),DPRED(*),
     4 PROPS(NPROPS),COORDS(*),DROT(3,3),PNEWDT,CELENT,
     5 DFGRD0(3,3),DFGRD1(3,3)

C---- Local variables
      INTEGER I,J,K
      DOUBLE PRECISION EE,ANU,AMU,AKAPPA,ALAM
      DOUBLE PRECISION F(3,3),DETJ,DETJ23
      DOUBLE PRECISION BB(3,3),BBAR(3,3),TRBBAR
      DOUBLE PRECISION ETOT(6)

      DOUBLE PRECISION HALF,ZERO,ONE,TWO,THREE,THIRD,TWOTHD
      PARAMETER (ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0)
      PARAMETER (TWO=2.0D0, THREE=3.0D0)
      PARAMETER (THIRD=ONE/THREE, TWOTHD=TWO/THREE)

C==== Inicializar salidas =================================================
      DO I=1,NTENS
         DDSDDT(I) = ZERO
         DRPLDE(I) = ZERO
      END DO
      DRPLDT = ZERO
      SPD    = ZERO
      SCD    = ZERO
      RPL    = ZERO
      SSE    = ZERO
      PNEWDT = ONE

C==== Propiedades =========================================================
      EE     = PROPS(1)
      ANU    = PROPS(2)
      AMU    = EE / (TWO*(ONE + ANU))
      AKAPPA = EE / (THREE*(ONE - TWO*ANU))
      ALAM   = AKAPPA - TWOTHD*AMU

C==== Deformacion total (STRAN + DSTRAN) ==================================
C     Voigt: e11, e22, e33, g12, g13, g23  (g = 2*e para corte)
      DO I=1,NTENS
         ETOT(I) = STRAN(I) + DSTRAN(I)
      END DO

C==== F = I + sym(epsilon_total) ==========================================
C     Aproximacion para DEFORMATION='PETIT' (ignora rotacion)
      F(1,1) = ONE + ETOT(1)
      F(2,2) = ONE + ETOT(2)
      F(3,3) = ONE + ETOT(3)
      F(1,2) = HALF*ETOT(4)
      F(2,1) = HALF*ETOT(4)
      F(1,3) = ZERO
      F(2,3) = ZERO
      F(3,1) = ZERO
      F(3,2) = ZERO
      IF (NTENS .GE. 5) THEN
         F(1,3) = HALF*ETOT(5)
         F(3,1) = HALF*ETOT(5)
      END IF
      IF (NTENS .GE. 6) THEN
         F(2,3) = HALF*ETOT(6)
         F(3,2) = HALF*ETOT(6)
      END IF

C==== det(F) ==============================================================
      DETJ = F(1,1)*(F(2,2)*F(3,3) - F(2,3)*F(3,2))
     1     - F(1,2)*(F(2,1)*F(3,3) - F(2,3)*F(3,1))
     2     + F(1,3)*(F(2,1)*F(3,2) - F(2,2)*F(3,1))

C---- Proteccion contra inversion de elemento -----------------------------
      IF (DETJ .LT. 1.0D-2) THEN
         PNEWDT = 0.25D0
         DO I=1,NTENS
            DO J=1,NTENS
               DDSDDE(I,J) = ZERO
            END DO
         END DO
         DO I=1,NDI
            DO J=1,NDI
               DDSDDE(I,J) = ALAM
            END DO
            DDSDDE(I,I) = ALAM + TWO*AMU
         END DO
         DO I=NDI+1,NTENS
            DDSDDE(I,I) = AMU
         END DO
         GOTO 500
      END IF

C==== Left Cauchy-Green: b = F*F^T ========================================
      DO I=1,3
         DO J=1,3
            BB(I,J) = ZERO
            DO K=1,3
               BB(I,J) = BB(I,J) + F(I,K)*F(J,K)
            END DO
         END DO
      END DO

C==== Isochoric bbar = J^{-2/3} b ========================================
      DETJ23 = DETJ**(-TWOTHD)
      DO I=1,3
         DO J=1,3
            BBAR(I,J) = DETJ23 * BB(I,J)
         END DO
      END DO
      TRBBAR = BBAR(1,1) + BBAR(2,2) + BBAR(3,3)

C==== Cauchy: sigma = (mu/J)*dev(bbar) + kappa*(J-1)*I ====================
C     Voigt: S11 S22 S33 S12 S13 S23
      STRESS(1) = (AMU/DETJ)*(BBAR(1,1) - THIRD*TRBBAR)
     1          + AKAPPA*(DETJ - ONE)
      STRESS(2) = (AMU/DETJ)*(BBAR(2,2) - THIRD*TRBBAR)
     1          + AKAPPA*(DETJ - ONE)
      STRESS(3) = (AMU/DETJ)*(BBAR(3,3) - THIRD*TRBBAR)
     1          + AKAPPA*(DETJ - ONE)
      STRESS(4) = (AMU/DETJ)*BBAR(1,2)
      IF (NTENS .GE. 5) STRESS(5) = (AMU/DETJ)*BBAR(1,3)
      IF (NTENS .GE. 6) STRESS(6) = (AMU/DETJ)*BBAR(2,3)

C==== Tangente Neo-Hooke analitico ========================================
C     c_ijkl = kappa*(2J-1)*d_ij*d_kl
C            + (mu/J)*(d_ik*bbar_jl + d_il*bbar_jk)/2          [simetria]
C            - (2*mu)/(3*J)*(bbar_ij*d_kl + d_ij*bbar_kl)
C            + (2*mu)/(9*J)*tr(bbar)*d_ij*d_kl
C
C     Nota: estos son los componentes de Voigt con la convencion
C     de que los indices de corte van con factor 1 (no 2).
      DO I=1,NTENS
         DO J=1,NTENS
            DDSDDE(I,J) = ZERO
         END DO
      END DO

C     Bloque diag-diag (i,j = 1..3): c_iijj
C     = kappa*(2J-1) + (2*mu)/(9*J)*tr(bbar)
C       - (2*mu)/(3*J)*(bbar_ii + bbar_jj)   para i != j
C       + (mu/J)*bbar_ii                      para i == j (extra)

C     Simplificacion: usar formula cerrada por componente
C     c_1111 = kappa*(2J-1) + (2mu)/(9J)*trBB - (4mu)/(3J)*BB11 + (mu/J)*BB11
C            = kappa*(2J-1) + (2mu)/(9J)*trBB - (mu/J)*BB11/3
C     Hmm esto se complica. Uso la forma directa.

C     Vol + Dev contribution:
      DDSDDE(1,1) = AKAPPA*(TWO*DETJ - ONE)
     1  + (AMU/DETJ)*(TWOTHD*BBAR(1,1) + TWOTHD*TRBBAR*THIRD)
      DDSDDE(2,2) = AKAPPA*(TWO*DETJ - ONE)
     1  + (AMU/DETJ)*(TWOTHD*BBAR(2,2) + TWOTHD*TRBBAR*THIRD)
      DDSDDE(3,3) = AKAPPA*(TWO*DETJ - ONE)
     1  + (AMU/DETJ)*(TWOTHD*BBAR(3,3) + TWOTHD*TRBBAR*THIRD)

C     Off-diagonal normal: c_iijj (i!=j)
      DDSDDE(1,2) = AKAPPA*(TWO*DETJ - ONE)
     1  - (AMU/DETJ)*THIRD*(BBAR(1,1) + BBAR(2,2))
     2  + (AMU/DETJ)*TWOTHD*TRBBAR*THIRD*THIRD
      DDSDDE(1,3) = AKAPPA*(TWO*DETJ - ONE)
     1  - (AMU/DETJ)*THIRD*(BBAR(1,1) + BBAR(3,3))
     2  + (AMU/DETJ)*TWOTHD*TRBBAR*THIRD*THIRD
      DDSDDE(2,3) = AKAPPA*(TWO*DETJ - ONE)
     1  - (AMU/DETJ)*THIRD*(BBAR(2,2) + BBAR(3,3))
     2  + (AMU/DETJ)*TWOTHD*TRBBAR*THIRD*THIRD

      DDSDDE(2,1) = DDSDDE(1,2)
      DDSDDE(3,1) = DDSDDE(1,3)
      DDSDDE(3,2) = DDSDDE(2,3)

C     Corte: c_1212 = (mu/J) * (bbar_11+bbar_22)/2
C     Pero en Voigt: DDSDDE(4,4) corresponde a d(S12)/d(g12)
C     Como g12 = 2*e12, y sigma depende de b que depende de F^2,
C     el coeficiente tiene un factor 1/2 para el shear:
      DDSDDE(4,4) = (AMU/DETJ) * HALF*(BBAR(1,1)+BBAR(2,2))
      IF (NTENS .GE. 5)
     1  DDSDDE(5,5) = (AMU/DETJ) * HALF*(BBAR(1,1)+BBAR(3,3))
      IF (NTENS .GE. 6)
     1  DDSDDE(6,6) = (AMU/DETJ) * HALF*(BBAR(2,2)+BBAR(3,3))

 500  CONTINUE

C==== Log (solo EL=20, GP=1 y GP=5) ======================================
      IF (NOEL.EQ.20 .AND. (NPT.EQ.1 .OR. NPT.EQ.5)) THEN
         OPEN(UNIT=99,FILE='/work/umat_defgrad.log',
     1        POSITION='APPEND',STATUS='UNKNOWN')
         WRITE(99,'(A,I6,A,I3,A,F8.4,A,I4)')
     1     ' EL=',NOEL,' GP=',NPT,' t=',TIME(2),' it=',KINC
         WRITE(99,'(A,6F12.6)') '   STRAN =',
     1     (STRAN(I),I=1,NTENS)
         WRITE(99,'(A,6F12.6)') '   DSTRAN=',
     1     (DSTRAN(I),I=1,NTENS)
         WRITE(99,'(A,3F12.6)') '   F diag=',F(1,1),F(2,2),F(3,3)
         WRITE(99,'(A,F12.6)')  '   det(F)=',DETJ
         WRITE(99,'(A,3F12.6)') '   DFGRD1 r1=',
     1     DFGRD1(1,1),DFGRD1(1,2),DFGRD1(1,3)
         WRITE(99,'(A,3F12.6)') '   DFGRD1 r2=',
     1     DFGRD1(2,1),DFGRD1(2,2),DFGRD1(2,3)
         WRITE(99,'(A,3F12.6)') '   DFGRD1 r3=',
     1     DFGRD1(3,1),DFGRD1(3,2),DFGRD1(3,3)
         WRITE(99,'(A,6F12.4)') '   sigma=',
     1     (STRESS(I),I=1,NTENS)
         CLOSE(99)
      END IF

      RETURN
      END
