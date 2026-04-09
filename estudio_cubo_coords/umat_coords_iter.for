C===========================================================================
C UMAT elastica con coordenadas inyectadas via STATEV.
C
C STATEV layout (NSTATV = 6):
C   V1:V3   -> Coords Gauss ACTUALES  (gestionadas por .comm, NO tocar)
C   V4:V6   -> Coords Gauss ANTERIORES (gestionadas por .comm, NO tocar)
C
C Las coordenadas se inyectan desde el .comm mediante la pipeline
C DISC -> FORMULE -> EVAL -> ASSE -> ETAT_INIT (patron de PT_test.comm).
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
     1        KSPT,KSTEP,KINC,I,J

      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),SSE,SPD,SCD,RPL,DDSDDT(NTENS),
     2 DRPLDE(NTENS),DRPLDT,STRAN(NTENS),DSTRAN(NTENS),
     3 TIME(2),DTIME,TEMP,DTEMP,PREDEF(*),DPRED(*),
     4 PROPS(NPROPS),COORDS(*),DROT(3,3),PNEWDT,CELENT,
     5 DFGRD0(3,3),DFGRD1(3,3)

      DOUBLE PRECISION E,NU,MU,LAM,ONE,TWO,ZERO
      DOUBLE PRECISION DS(6)

      PARAMETER (ONE=1.0D0, TWO=2.0D0, ZERO=0.0D0)

C---- Inicializacion -------------------------------------------------------
      DO I=1,NTENS
         DDSDDT(I) = ZERO
         DRPLDE(I) = ZERO
         DS(I)     = ZERO
      END DO

      DRPLDT = ZERO
      SPD    = ZERO
      SCD    = ZERO
      RPL    = ZERO
      SSE    = ZERO
      PNEWDT = ONE

C---- Elasticidad isotropa lineal -----------------------------------------
      E  = PROPS(1)
      NU = PROPS(2)

      MU  = E / (TWO*(ONE+NU))
      LAM = E*NU / ((ONE+NU)*(ONE-TWO*NU))

      DO I=1,NTENS
         DO J=1,NTENS
            DDSDDE(I,J) = ZERO
         END DO
      END DO

      DDSDDE(1,1) = LAM + TWO*MU
      DDSDDE(1,2) = LAM
      DDSDDE(1,3) = LAM

      DDSDDE(2,1) = LAM
      DDSDDE(2,2) = LAM + TWO*MU
      DDSDDE(2,3) = LAM

      DDSDDE(3,1) = LAM
      DDSDDE(3,2) = LAM
      DDSDDE(3,3) = LAM + TWO*MU

      DDSDDE(4,4) = MU
      DDSDDE(5,5) = MU
      DDSDDE(6,6) = MU

      DO I=1,NTENS
         DO J=1,NTENS
            DS(I) = DS(I) + DDSDDE(I,J)*DSTRAN(J)
         END DO
      END DO

      DO I=1,NTENS
         STRESS(I) = STRESS(I) + DS(I)
      END DO

C---- Imprimir coordenadas para dos GPs: uno cerca del piso, otro arriba --
C     STATEV(1:3) = coords actuales, STATEV(4:6) = coords anteriores

      IF (NPT.EQ.1 .OR. NPT.EQ.5) THEN
         OPEN(UNIT=99,FILE='/work/umat_coords.log',
     1        POSITION='APPEND',STATUS='UNKNOWN')
         WRITE(99,'(A,I6,A,I3,A,F8.4,A,3F12.6,A,3F12.6)')
     1     ' Elemento=',NOEL,' GP=',NPT,' t=',TIME(2),
     2     ' Coord_actual=',STATEV(1),STATEV(2),STATEV(3),
     3     ' Coord_prev=',STATEV(4),STATEV(5),STATEV(6)
         CLOSE(99)
      END IF

      RETURN
      END
