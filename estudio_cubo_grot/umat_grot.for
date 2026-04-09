C===========================================================================
C  UMAT Hencky hiperelastico para DEFORMATION='GDEF_LOG'
C
C  Con GDEF_LOG, Code_Aster pasa deformacion logaritmica en STRAN/DSTRAN.
C  La ley constitutiva es formalmente identica a elasticidad lineal:
C
C    sigma = lambda*tr(eps)*I + 2*mu*eps
C
C  pero eps = deformacion logaritmica (Hencky), lo que produce
C  un modelo hiperelastico objetivo para grandes deformaciones.
C
C  PROPS(1) = E, PROPS(2) = nu
C
C  TEST PRINCIPAL: verificar si DFGRD0/DFGRD1 son poblados por CA.
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
      INTEGER I,J,A
      DOUBLE PRECISION EE,ANU,AMU,AKAPPA,ALAM
      DOUBLE PRECISION ETOT(6),TRACE
      DOUBLE PRECISION F1TF1(3,3)

      DOUBLE PRECISION ZERO,ONE,TWO,THREE,TWOTHD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      PARAMETER (TWO=2.0D0, THREE=3.0D0)
      PARAMETER (TWOTHD=TWO/THREE)

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

C==== Deformacion total logaritmica (Voigt) ===============================
C     Voigt: eps11, eps22, eps33, 2*eps12, 2*eps13, 2*eps23
      DO I=1,NTENS
         ETOT(I) = STRAN(I) + DSTRAN(I)
      END DO

C==== Cauchy stress: sigma = lambda*tr(eps)*I + 2*mu*eps ==================
      TRACE = ETOT(1) + ETOT(2) + ETOT(3)

C     Componentes normales
      DO I=1,NDI
         STRESS(I) = ALAM*TRACE + TWO*AMU*ETOT(I)
      END DO
C     Componentes de corte: sigma_ij = 2*mu*eps_ij = mu*gamma_ij
      DO I=NDI+1,NTENS
         STRESS(I) = AMU*ETOT(I)
      END DO

C==== Tangente: d(sigma)/d(eps) = elasticidad isotropa ====================
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

C==== Log (solo EL=20, GP=1 y GP=5) ======================================
      IF (NOEL.EQ.20 .AND. (NPT.EQ.1 .OR. NPT.EQ.5)) THEN
         OPEN(UNIT=99,FILE='/work/umat_grot.log',
     1        POSITION='APPEND',STATUS='UNKNOWN')
         WRITE(99,'(A,I6,A,I3,A,F8.4,A,I4)')
     1     ' EL=',NOEL,' GP=',NPT,' t=',TIME(2),' it=',KINC
         WRITE(99,'(A,6F12.6)') '   STRAN =',
     1     (STRAN(I),I=1,NTENS)
         WRITE(99,'(A,6F12.6)') '   DSTRAN=',
     1     (DSTRAN(I),I=1,NTENS)
         WRITE(99,'(A,F12.6)')  '   tr(e) =',TRACE
         WRITE(99,'(A,3F12.6)') '   DFGRD0 r1=',
     1     DFGRD0(1,1),DFGRD0(1,2),DFGRD0(1,3)
         WRITE(99,'(A,3F12.6)') '   DFGRD0 r2=',
     1     DFGRD0(2,1),DFGRD0(2,2),DFGRD0(2,3)
         WRITE(99,'(A,3F12.6)') '   DFGRD0 r3=',
     1     DFGRD0(3,1),DFGRD0(3,2),DFGRD0(3,3)
         WRITE(99,'(A,3F12.6)') '   DFGRD1 r1=',
     1     DFGRD1(1,1),DFGRD1(1,2),DFGRD1(1,3)
         WRITE(99,'(A,3F12.6)') '   DFGRD1 r2=',
     1     DFGRD1(2,1),DFGRD1(2,2),DFGRD1(2,3)
         WRITE(99,'(A,3F12.6)') '   DFGRD1 r3=',
     1     DFGRD1(3,1),DFGRD1(3,2),DFGRD1(3,3)
C        Verificacion: F1^T*F1 diagonal (si DFGRD1 poblado)
         DO I=1,3
            DO J=1,3
               F1TF1(I,J) = ZERO
               DO A=1,3
                  F1TF1(I,J) = F1TF1(I,J)+DFGRD1(A,I)*DFGRD1(A,J)
               END DO
            END DO
         END DO
         WRITE(99,'(A,3F12.6)') '   F1TF1 d=',
     1     F1TF1(1,1),F1TF1(2,2),F1TF1(3,3)
         WRITE(99,'(A,6F12.4)') '   sigma=',
     1     (STRESS(I),I=1,NTENS)
         CLOSE(99)
      END IF

      RETURN
      END
