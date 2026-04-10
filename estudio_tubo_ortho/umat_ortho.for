C===========================================================================
C UMAT ortotropia cilindrica para el tubo con presion interna.
C
C Material: 9 propiedades en ejes cilindricos (r, theta, z):
C   PROPS(1) = E_r       modulo de Young radial
C   PROPS(2) = E_th      modulo de Young circunferencial
C   PROPS(3) = E_z       modulo de Young axial
C   PROPS(4) = nu_rth    coef. de Poisson r-theta
C   PROPS(5) = nu_rz     coef. de Poisson r-z
C   PROPS(6) = nu_thz    coef. de Poisson theta-z
C   PROPS(7) = G_rth     modulo de corte r-theta
C   PROPS(8) = G_rz      modulo de corte r-z
C   PROPS(9) = G_thz     modulo de corte theta-z
C
C Estrategia de rotacion (igual que UMAT4001.f90):
C   Se construye R_CYL cuyas FILAS son los vectores base locales expresados
C   en el sistema global cartesiano:
C         Fila 1: e_r   = ( x/r,  y/r, 0 )
C         Fila 2: e_th  = (-y/r,  x/r, 0 )
C         Fila 3: e_z   = ( 0,    0,   1 )
C   Rotacion de tensor 2do orden:
C         Global -> Local:  T_loc = R . T_glob . R^T
C         Local  -> Global: T_glob = R^T . T_loc . R
C   Las deformaciones de ingenieria (gamma) se convierten a tensoriales
C   antes de rotar y se reconvierten al extraer.
C   La rigidez se rota columna a columna (procedimiento de UMAT4001).
C
C   Las coordenadas actuales se leen de STATEV(1:3) —  Code_Aster 17.4
C   no puebla COORDS() en la interfaz UMAT (vease FIX5 de UMAT4001.f90).
C
C STATEV layout (NSTATV = 6):
C   V1:V3  -> Coords Gauss ACTUALES  (inyectadas por .comm)
C   V4:V6  -> Coords Gauss ANTERIORES (inyectadas por .comm)
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

C---- Variables locales ---------------------------------------------------
      DOUBLE PRECISION ER,ETH,EZ,NURTH,NURZ,NUTHZ
      DOUBLE PRECISION GRTH,GRZ,GTHZ
      DOUBLE PRECISION S11,S12,S13,S22,S23,S33,DET3
      DOUBLE PRECISION CCYL(6,6)
      DOUBLE PRECISION R(3,3)
      DOUBLE PRECISION XP,YP,RAD
      DOUBLE PRECISION STR_LOC(6),DST_LOC(6),DS_LOC(6)
      DOUBLE PRECISION CCART(6,6)

      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)

C---- Inicializacion -------------------------------------------------------
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

C---- Propiedades del material ---------------------------------------------
      ER    = PROPS(1)
      ETH   = PROPS(2)
      EZ    = PROPS(3)
      NURTH = PROPS(4)
      NURZ  = PROPS(5)
      NUTHZ = PROPS(6)
      GRTH  = PROPS(7)
      GRZ   = PROPS(8)
      GTHZ  = PROPS(9)

C---- Rigidez en el marco cilindrico: C_cyl = inv(S_cyl) ------------------
C     Voigt cilindrico: (rr, thth, zz, r-th, r-z, th-z)
C     Compliance (simetria de Maxwell: nu_ij/E_i = nu_ji/E_j):
      S11 = ONE/ER
      S22 = ONE/ETH
      S33 = ONE/EZ
      S12 = -NURTH/ER
      S13 = -NURZ/ER
      S23 = -NUTHZ/ETH

      DET3 = S11*(S22*S33 - S23*S23)
     1      -S12*(S12*S33 - S23*S13)
     2      +S13*(S12*S23 - S22*S13)

      DO I=1,6
         DO J=1,6
            CCYL(I,J) = ZERO
         END DO
      END DO

      CCYL(1,1) = (S22*S33 - S23*S23)/DET3
      CCYL(2,2) = (S11*S33 - S13*S13)/DET3
      CCYL(3,3) = (S11*S22 - S12*S12)/DET3
      CCYL(1,2) = (S13*S23 - S12*S33)/DET3
      CCYL(2,1) = CCYL(1,2)
      CCYL(1,3) = (S12*S23 - S13*S22)/DET3
      CCYL(3,1) = CCYL(1,3)
      CCYL(2,3) = (S12*S13 - S11*S23)/DET3
      CCYL(3,2) = CCYL(2,3)
      CCYL(4,4) = GRTH
      CCYL(5,5) = GRZ
      CCYL(6,6) = GTHZ

C---- Matriz de rotacion R_CYL (FILAS = ejes locales en global) -----------
C     Code_Aster 17.4 no puebla COORDS(); se usan STATEV(1:3)
C     (mismo patron que UMAT4001.f90, seccion 2)
      XP  = STATEV(1)
      YP  = STATEV(2)
      RAD = DSQRT(XP*XP + YP*YP)

      DO I=1,3
         DO J=1,3
            R(I,J) = ZERO
         END DO
      END DO

      IF (RAD .GT. 1.0D-9) THEN
C        Fila 1: direccion radial   e_r  = ( x/r,  y/r,  0 )
         R(1,1) =  XP/RAD
         R(1,2) =  YP/RAD
         R(1,3) =  ZERO
C        Fila 2: direccion hoop     e_th = (-y/r,  x/r,  0 )
         R(2,1) = -YP/RAD
         R(2,2) =  XP/RAD
         R(2,3) =  ZERO
C        Fila 3: direccion axial    e_z  = ( 0,     0,    1 )
         R(3,1) =  ZERO
         R(3,2) =  ZERO
         R(3,3) =  ONE
      ELSE
C        En el eje (r=0): identidad como fallback
         R(1,1) = ONE
         R(2,2) = ONE
         R(3,3) = ONE
      END IF

C---- Girar deformacion incremental: Global -> Local (G2L) ----------------
C     T_loc = R . T_tens . R^T   (T_tens con shear tensorial = gamma/2)
C     En ingenieria: V_local(4:6) = 2 * T_loc(shear)
      CALL ROT_STRAIN_G2L6(DSTRAN, DST_LOC, R)

C---- Girar tension actual: Global -> Local --------------------------------
      CALL ROT_STRESS_G2L6(STRESS, STR_LOC, R)

C---- Incremento de tension en marco local: DS_LOC = CCYL * DST_LOC ------
      DO I=1,6
         DS_LOC(I) = ZERO
         DO J=1,6
            DS_LOC(I) = DS_LOC(I) + CCYL(I,J)*DST_LOC(J)
         END DO
         STR_LOC(I) = STR_LOC(I) + DS_LOC(I)
      END DO

C---- Girar tension actualizada: Local -> Global --------------------------
C     T_glob = R^T . T_loc . R
      CALL ROT_STRESS_L2G6(STR_LOC, STRESS, R)

C---- Girar rigidez: Local -> Global (columna a columna) ------------------
C     Mismo procedimiento que ROT_STIFFNESS_ENG_L2G de UMAT4001.f90
      CALL ROT_STIFF_L2G6(CCYL, CCART, R)

      DO I=1,NTENS
         DO J=1,NTENS
            DDSDDE(I,J) = CCART(I,J)
         END DO
      END DO

C---- Registro de coordenadas: primer GP de cada elemento -----------------
      IF (NPT.EQ.1) THEN
         OPEN(UNIT=99,FILE='/work/umat_ortho.log',
     1        POSITION='APPEND',STATUS='UNKNOWN')
         WRITE(99,'(A,I6,A,I3,A,F8.4,A,3F12.6,A,3F12.6)')
     1     ' Elemento=',NOEL,' GP=',NPT,' t=',TIME(2),
     2     ' Coord_actual=',STATEV(1),STATEV(2),STATEV(3),
     3     ' Coord_prev=',STATEV(4),STATEV(5),STATEV(6)
         CLOSE(99)
      END IF

      RETURN
      END

C===========================================================================
C ROT_STRESS_G2L6: rota tensor de TENSION  Global -> Local
C    T_loc = R . T_glob . R^T
C    Voigt: (11, 22, 33, 12, 13, 23)   shear COMPLETO (no factor 1/2)
C===========================================================================
      SUBROUTINE ROT_STRESS_G2L6(VG, VL, R)
      IMPLICIT NONE
      DOUBLE PRECISION VG(6),VL(6),R(3,3)
      DOUBLE PRECISION T(3,3),RT(3,3),TR(3,3)
      INTEGER I,J,K

      T(1,1)=VG(1); T(2,2)=VG(2); T(3,3)=VG(3)
      T(1,2)=VG(4); T(2,1)=VG(4)
      T(1,3)=VG(5); T(3,1)=VG(5)
      T(2,3)=VG(6); T(3,2)=VG(6)

C     RT = R . T
      DO I=1,3
         DO J=1,3
            RT(I,J)=0.0D0
            DO K=1,3
               RT(I,J)=RT(I,J)+R(I,K)*T(K,J)
            END DO
         END DO
      END DO
C     TR = RT . R^T  = R . T . R^T
      DO I=1,3
         DO J=1,3
            TR(I,J)=0.0D0
            DO K=1,3
               TR(I,J)=TR(I,J)+RT(I,K)*R(J,K)
            END DO
         END DO
      END DO

      VL(1)=TR(1,1); VL(2)=TR(2,2); VL(3)=TR(3,3)
      VL(4)=TR(1,2); VL(5)=TR(1,3); VL(6)=TR(2,3)
      RETURN
      END

C===========================================================================
C ROT_STRESS_L2G6: rota tensor de TENSION Local -> Global
C    T_glob = R^T . T_loc . R
C===========================================================================
      SUBROUTINE ROT_STRESS_L2G6(VL, VG, R)
      IMPLICIT NONE
      DOUBLE PRECISION VL(6),VG(6),R(3,3)
      DOUBLE PRECISION T(3,3),RT(3,3),TR(3,3)
      INTEGER I,J,K

      T(1,1)=VL(1); T(2,2)=VL(2); T(3,3)=VL(3)
      T(1,2)=VL(4); T(2,1)=VL(4)
      T(1,3)=VL(5); T(3,1)=VL(5)
      T(2,3)=VL(6); T(3,2)=VL(6)

C     RT = R^T . T   (RT(i,j) = sum_k R(k,i)*T(k,j))
      DO I=1,3
         DO J=1,3
            RT(I,J)=0.0D0
            DO K=1,3
               RT(I,J)=RT(I,J)+R(K,I)*T(K,J)
            END DO
         END DO
      END DO
C     TR = RT . R = R^T . T . R
      DO I=1,3
         DO J=1,3
            TR(I,J)=0.0D0
            DO K=1,3
               TR(I,J)=TR(I,J)+RT(I,K)*R(K,J)
            END DO
         END DO
      END DO

      VG(1)=TR(1,1); VG(2)=TR(2,2); VG(3)=TR(3,3)
      VG(4)=TR(1,2); VG(5)=TR(1,3); VG(6)=TR(2,3)
      RETURN
      END

C===========================================================================
C ROT_STRAIN_G2L6: rota tensor de DEFORMACION  Global -> Local
C    Convierte gamma (Voigt ingenieria) <-> epsilon tensorial antes/despues
C    T_loc = R . T_tens . R^T   luego extrae gamma = 2*epsilon_shear
C===========================================================================
      SUBROUTINE ROT_STRAIN_G2L6(VG, VL, R)
      IMPLICIT NONE
      DOUBLE PRECISION VG(6),VL(6),R(3,3)
      DOUBLE PRECISION T(3,3),RT(3,3),TR(3,3)
      INTEGER I,J,K

C     Convertir a notacion tensorial (gamma -> epsilon = gamma/2)
      T(1,1)=VG(1);       T(2,2)=VG(2);       T(3,3)=VG(3)
      T(1,2)=0.5D0*VG(4); T(2,1)=T(1,2)
      T(1,3)=0.5D0*VG(5); T(3,1)=T(1,3)
      T(2,3)=0.5D0*VG(6); T(3,2)=T(2,3)

      DO I=1,3
         DO J=1,3
            RT(I,J)=0.0D0
            DO K=1,3
               RT(I,J)=RT(I,J)+R(I,K)*T(K,J)
            END DO
         END DO
      END DO
      DO I=1,3
         DO J=1,3
            TR(I,J)=0.0D0
            DO K=1,3
               TR(I,J)=TR(I,J)+RT(I,K)*R(J,K)
            END DO
         END DO
      END DO

C     Reconvertir a notacion ingenieria (epsilon -> gamma = 2*epsilon)
      VL(1)=TR(1,1); VL(2)=TR(2,2); VL(3)=TR(3,3)
      VL(4)=2.0D0*TR(1,2)
      VL(5)=2.0D0*TR(1,3)
      VL(6)=2.0D0*TR(2,3)
      RETURN
      END

C===========================================================================
C ROT_STIFF_L2G6: rota rigidez 6x6 Local -> Global (columna a columna)
C    Para cada columna j: aplica una deformacion unitaria global E_j,
C    la lleva al marco local, multiplica por DL, rota la tension resultante
C    de vuelta al global. Identico a ROT_STIFFNESS_ENG_L2G de UMAT4001.f90.
C===========================================================================
      SUBROUTINE ROT_STIFF_L2G6(DL, DG, R)
      IMPLICIT NONE
      DOUBLE PRECISION DL(6,6),DG(6,6),R(3,3)
      DOUBLE PRECISION EPS_G(6),EPS_L(6),SIG_L(6),SIG_G(6)
      INTEGER I,J,K

      DO J=1,6
         DO I=1,6
            EPS_G(I)=0.0D0
         END DO
         EPS_G(J)=1.0D0

         CALL ROT_STRAIN_G2L6(EPS_G, EPS_L, R)

         DO I=1,6
            SIG_L(I)=0.0D0
            DO K=1,6
               SIG_L(I)=SIG_L(I)+DL(I,K)*EPS_L(K)
            END DO
         END DO

         CALL ROT_STRESS_L2G6(SIG_L, SIG_G, R)

         DO I=1,6
            DG(I,J)=SIG_G(I)
         END DO
      END DO
      RETURN
      END
