C=========================================================================================
C UMAT elástica - Martin Sebastián Armoa 2025. Lenguaje fortran
C-----------------------------------------------------------------------------------------
C Material elástico lineal isotrópico 3D (pequeñas deformaciones).
C PROPS(1) = E   (Módulo de Young)
C PROPS(2) = NU  (Coeficiente de Poisson)
C
C Flujo:
C   1) Construye DDSDDE = C (6x6, Voigt con cortantes ingenieriles).
C   2) DS = DDSDDE * DSTRAN
C   3) STRESS := STRESS + DS
C
C Compilar: gfortran -c elastica.for \
c  -O3 -ffixed-form -march=native -funroll-loops \
c  -ftree-vectorize -fno-protect-parens \
c  -fopt-info-vec-optimized -fopt-info-vec-missed

c  -O3 activa vectorización; -ftree-vectorize la fuerza explícitamente.
c  -march=native habilita las instrucciones SIMD del CPU local.
c  -fopt-info-vec-* te imprime qué bucles vectorizó o no.
c  -fno-protect-parens permite re-asociación (útil para FMA), ojo si necesitás aritmética estricta.
c  ! PARA ENLAZAR libumat.so PARA CODE ASTER
c  ! gfortran -shared -o libumat.so  elastica.o
C=========================================================================================
C Tensor de constantes elásticas C en notación de Voigt (6x6)
C Convención Abaqus: Voigt con cortantes ingenieriles
C   ε = [ ε11 ε22 ε33 γ12 γ23 γ13 ]^T
C   σ = [ σ11 σ22 σ33 σ12 σ23 σ13 ]^T
C
C        | λ+2μ   λ      λ      0      0      0 |
C        |  λ   λ+2μ     λ      0      0      0 |
C  C  =  |  λ     λ    λ+2μ     0      0      0 |
C        |  0     0      0      μ      0      0 |
C        |  0     0      0      0      μ      0 |
C        |  0     0      0      0      0      μ |
C
C con:
C   μ  = E / [2(1+ν)]        (módulo de corte)
C   λ  = Eν / [(1+ν)(1-2ν)]  (constante de Lamé)
C--------------------------------------------------------------------------- 

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     3 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     4 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     5 KSPT,KSTEP,KINC)

      IMPLICIT NONE

C---- Tipos y argumentos ----------------------------------------------------
      CHARACTER*80 CMNAME
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,
     1        KSPT,KSTEP,KINC,I,J
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),SSE,SPD,SCD,RPL,DDSDDT(NTENS),
     2 DRPLDE(NTENS),DRPLDT,STRAN(NTENS),DSTRAN(NTENS),
     3 TIME(2),DTIME,TEMP,DTEMP,PREDEF(*),DPRED(*),
     4 PROPS(NPROPS),COORDS(*),DROT(3,3),PNEWDT,CELENT,
     5 DFGRD0(3,3),DFGRD1(3,3)

C---- Constantes y auxiliares ----------------------------------------------
      DOUBLE PRECISION E,NU,MU,LAM,ONE,TWO,ZERO
      DOUBLE PRECISION DS(6)
      PARAMETER (ONE=1.0D0, TWO=2.0D0, ZERO=0.0D0)

C---- Inicializaciones ------------------------------------------------------
      DO I=1,6
         DDSDDT(I) = ZERO
         DRPLDE(I) = ZERO
         DS(I)     = ZERO
      END DO
      DRPLDT = ZERO
      SPD    = ZERO
      SCD    = ZERO
      RPL    = ZERO
      PNEWDT = ONE

C---- Leer propiedades y constantes de Lamé --------------------------------
      E  = PROPS(1)
      NU = PROPS(2)
      MU  = E / (TWO*(ONE+NU))
      LAM = E*NU / ((ONE+NU)*(ONE-TWO*NU))

C---- Limpiar y construir DDSDDE (isótropa 3D en Voigt) --------------------
      DO I=1,6
         DO J=1,6
            DDSDDE(I,J) = ZERO
         END DO
      END DO

C
      DDSDDE(1,1) = LAM + TWO*MU
      DDSDDE(1,2) = LAM
      DDSDDE(1,3) = LAM

      DDSDDE(2,1) = LAM
      DDSDDE(2,2) = LAM + TWO*MU
      DDSDDE(2,3) = LAM

      DDSDDE(3,1) = LAM
      DDSDDE(3,2) = LAM
      DDSDDE(3,3) = LAM + TWO*MU

C
      DDSDDE(4,4) = MU
      DDSDDE(5,5) = MU
      DDSDDE(6,6) = MU

C---- DS = DDSDDE * DSTRAN --------------------------------------------------
      DO I=1,6
         DO J=1,6
            DS(I) = DS(I) + DDSDDE(I,J)*DSTRAN(J)
         END DO
      END DO

C---- Actualización de tensiones -------------------------------------------
      DO I=1,6
         STRESS(I) = STRESS(I) + DS(I)
      END DO

C---- (Opcional) Energía elástica: SSE = SSE + 0.5*DS:DE -------------------
C     Descomentá si querés acumular energía:
C     DO I=1,3
C        SSE = SSE + 0.5D0*(STRESS(I))*DSTRAN(I)
C     END DO
C     DO I=4,6
C        SSE = SSE + 0.5D0*(STRESS(I))*DSTRAN(I)
C     END DO

      RETURN
      END
