!=======================================================================
!  UMAT4001.f90 — CRNL-4003 (Rev. 1) para Zr-2.5Nb / Tubos de Presión CANDU
!  Método: Integración Implícita (Backward Euler) + Newton-Raphson
!  Base: ch_implicit.f90 (M. Armoa, IFIR, Feb. 2026)
!
!  CORRECCIONES APLICADAS RESPECTO A ch_implicit.f90:
!  ---------------------------------------------------
!  [FIX1] ROT_STRESS_G2L, ROT_STRESS_L2G y ROT_STRAIN_G2L tenían las
!         fórmulas de rotación de tensor de 2do orden INTERCAMBIADAS.
!         G→L correcto: R·T·Rᵀ    L→G correcto: Rᵀ·T·R
!         (donde las FILAS de R son los ejes locales en global).
!
!  [FIX2] Tasas de creep por cizallamiento RC(4:6) = 0 en todo el cálculo.
!         Se implementan correctamente usando la derivada de la función de
!         Hill con respecto a los esfuerzos de cizallamiento:
!         ∂s_eff/∂τ = 2·L·τ/s_eff  →  ε̇_τ = 2·L·τ·K*(exponencial)
!         Los parámetros de cizallamiento se leen de PROPS(46:57) o se
!         usan los valores estándar de CRNL-4003 Rev.1 si NPROPS < 57.
!
!  [FIX3] Acumulación de fluencia con unidades incorrectas.
!         Error: STATEV(1) += flux[n/m²/s] * dt_h[h]  (unidades n·h/m²·s)
!         Correcto: STATEV(1) += flux[n/m²/s] * dt_h[h] * 3600[s/h]  [n/m²]
!
!  [FIX4] MAT_INV_6X6 sin pivoteo parcial puede ser numéricamente inestable.
!         Se agrega intercambio de filas (pivoteo parcial) en Gauss-Jordan.
!
!  NOTA: la ecuación de crecimiento (growth) no tiene componentes de
!  cizallamiento (stress-independent), RC(4:6) solo afectan creep.
!  Los parámetros Hill de cizallamiento TH1, TH2, CREEP(L,M,N) usan
!  valores por defecto de CRNL-4003 Rev.1 si no se proveen en PROPS.
!  Este UMAT es exclusively para PT (Pressure Tube). CT no implementado.
!
!  PROPS(1:9)   : E_a, E_t, E_r, nu_at, nu_ar, nu_tr, G_at, G_ar, G_tr
!  PROPS(10:37) : K1,K2,K3,Q1,Q3,F1,G1,F2,G2,Kc,Q4,K5,K41,K42,Fb,Ff,Gb,Gf,
!                 Kg,Q6,K61,K62,Cg,Bg,Ga_b,Ga_f,Gt_b,Gt_f
!  PROPS(38)    : flux_SI [n/m²/s].  Si < 0: |val| = pico, perfil coseno²
!  PROPS(39)    : factor tiempo (ej: 1.0 si DTIME ya en horas)
!  PROPS(40)    : x_scale_m (factor escala coordenada axial a metros)
!  PROPS(41)    : Ltube [m]
!  PROPS(42)    : AX (eje axial global: 1=X, 2=Y, 3=Z)
!  PROPS(43)    : x0 (origen axial)
!  PROPS(44)    : T_offset (para conversión °C → K, típicamente 273.15)
!  PROPS(45)    : flag creep (1=activo, 0=solo elástico)
!  PROPS(46:57) : [OPCIONAL] parámetros shear (TH1L,TH1M,TH1N,TH2L,TH2M,TH2N,
!                 CREEPLB,CREEPLF,CREEPMB,CREEPMF,CREEPNB,CREEPNF)
!                 Por defecto: valores estándar CRNL-4003 Rev.1
!
!  Autor original: M. Armoa - IFIR
!  Correcciones:   J.A. (IFIR) — Marzo 2026
!=======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
     & RPL,DDSDDT,DRPLDE,DRPLDT, &
     & STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPREDEF, &
     & CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS, &
     & DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT, &
     & LAYER,KSPT,KSTEP,KINC)

      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE

      ! --- Variables de Interfaz Standard UMAT ---
      INTEGER, INTENT(IN) :: NTENS, NSTATV, NPROPS, NDI, NSHR
      INTEGER, INTENT(IN) :: NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      CHARACTER(LEN=8), INTENT(IN) :: CMNAME

      DOUBLE PRECISION, INTENT(INOUT) :: STRESS(NTENS), STATEV(NSTATV)
      DOUBLE PRECISION, INTENT(OUT)   :: DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION, INTENT(INOUT) :: SSE, SPD, SCD, RPL
      DOUBLE PRECISION, INTENT(OUT)   :: DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT
      DOUBLE PRECISION, INTENT(IN)    :: STRAN(NTENS), DSTRAN(NTENS)
      DOUBLE PRECISION, INTENT(IN)    :: TIME(2), DTIME, TEMP, DTEMP
      DOUBLE PRECISION, INTENT(IN)    :: PREDEF(1), DPREDEF(1)
      DOUBLE PRECISION, INTENT(IN)    :: PROPS(NPROPS), COORDS(3)
      DOUBLE PRECISION, INTENT(IN)    :: DROT(3,3), CELENT
      DOUBLE PRECISION, INTENT(INOUT) :: PNEWDT
      DOUBLE PRECISION, INTENT(IN)    :: DFGRD0(3,3), DFGRD1(3,3)

      ! --- Variables Locales del Algoritmo ---
      INTEGER :: I, J, K, ITER, AX
      INTEGER, PARAMETER :: MAX_ITER = 25
      DOUBLE PRECISION, PARAMETER :: TOL_RES  = 1.0D-5  ! Tolerancia residual [MPa]
      DOUBLE PRECISION, PARAMETER :: S3600    = 3600.0D0 ! [FIX3] segundos/hora

      ! Geometría y Estado
      DOUBLE PRECISION :: R_CYL(3,3)              ! Rotación Global→Local
      DOUBLE PRECISION :: STR_LOC(NTENS)           ! Stress local
      DOUBLE PRECISION :: DST_LOC(NTENS)           ! dStrain local
      DOUBLE PRECISION :: D_EL_LOC(NTENS,NTENS)   ! Rigidez elástica local
      DOUBLE PRECISION :: D_CONS_LOC(NTENS,NTENS) ! Tangente consistente local
      DOUBLE PRECISION :: my_X, my_Y, my_Z, axial_pos

      ! Newton-Raphson
      DOUBLE PRECISION :: SIG_TRIAL(NTENS)
      DOUBLE PRECISION :: SIG_CURR(NTENS)
      DOUBLE PRECISION :: RESID(NTENS)
      DOUBLE PRECISION :: JAC(NTENS,NTENS)
      DOUBLE PRECISION :: INV_JAC(NTENS,NTENS)
      DOUBLE PRECISION :: DELTA_SIG(NTENS)

      ! Perturbación Numérica para Jacobiano
      DOUBLE PRECISION :: SIG_PERT(NTENS), RATE_PERT(NTENS)
      DOUBLE PRECISION :: RATE_CREEP(NTENS), RATE_GROWTH(NTENS)
      DOUBLE PRECISION :: DST_CREEP(NTENS), DST_GROWTH(NTENS)
      DOUBLE PRECISION :: H_PERT, VAL_PERT
      DOUBLE PRECISION, PARAMETER :: EPS_PERT = 1.0D-8

      ! Variables Físicas
      DOUBLE PRECISION :: flux_SI, dt_h, temp_K, fluence
      DOUBLE PRECISION :: x_m, xbar, Ltube, x0, x_scale_m
      DOUBLE PRECISION :: rad

      ! Propiedades Elásticas
      DOUBLE PRECISION :: E_a, E_t, E_r, nu_at, nu_ar, nu_tr, G_at, G_ar, G_tr

      ! ----------------------------------------------------------------
      ! 1. INICIALIZACIÓN Y CONVERSIÓN DE UNIDADES
      ! ----------------------------------------------------------------
      DDSDDT = 0.0D0; DRPLDE = 0.0D0; DRPLDT = 0.0D0
      PNEWDT = 1.0D0

      dt_h      = DTIME * PROPS(39)    ! Paso de tiempo en horas
      flux_SI   = PROPS(38)            ! Flujo [n/m²/s]
      temp_K    = TEMP + PROPS(44)     ! Temperatura en Kelvin
      IF (.NOT. IEEE_IS_FINITE(temp_K)) temp_K = PROPS(44)
      x_scale_m = PROPS(40)
      Ltube     = PROPS(41)
      AX        = NINT(PROPS(42))
      x0        = PROPS(43)
      fluence   = STATEV(1)            ! Fluencia acumulada [n/m²]

      IF (ABS(Ltube)     < 1.0D-12) Ltube     = 1.0D0
      IF (ABS(x_scale_m) < 1.0D-12) x_scale_m = 1.0D0

      ! Llamada inicial (DTIME≈0): devolver matriz elástica ortótropa
      IF (ABS(DTIME) < 1.0D-12) THEN
         E_a = PROPS(1); E_t = PROPS(2); E_r = PROPS(3)
         nu_at = PROPS(4); nu_ar = PROPS(5); nu_tr = PROPS(6)
         G_at = PROPS(7); G_ar = PROPS(8); G_tr = PROPS(9)
         CALL GET_ELASTIC_STIFFNESS_ORTO(E_a,E_t,E_r,nu_at,nu_ar,nu_tr, &
                                         G_at,G_ar,G_tr, DDSDDE, NTENS)
         STRESS = 0.0D0
         RETURN
      ENDIF

      ! ----------------------------------------------------------------
      ! 2. MATRIZ DE ROTACIÓN (GLOBAL → LOCAL CILÍNDRICO)
      !    R_CYL: FILAS = vectores base locales en global
      !    Transformación 2do orden G→L: T_loc = R · T_glob · Rᵀ  [FIX1]
      !    [FIX5] Code_Aster 17.4 NO puebla COORDS(3) en su interfaz UMAT.
      !    Coordenadas se pre-cargan vía ETAT_INIT en STATEV(NSTATV-2:NSTATV).
      !    Priorizar STATEV backup cuando contiene datos > 0.
      ! ----------------------------------------------------------------
      IF (NSTATV >= 3) THEN
         IF (ABS(STATEV(NSTATV-2)) + ABS(STATEV(NSTATV-1)) + &
             ABS(STATEV(NSTATV)) > 1.0D-12) THEN
            my_X = STATEV(NSTATV-2)
            my_Y = STATEV(NSTATV-1)
            my_Z = STATEV(NSTATV)
         ELSE
            my_X = COORDS(1)
            my_Y = COORDS(2)
            my_Z = COORDS(3)
         ENDIF
      ELSE
         my_X = COORDS(1)
         my_Y = COORDS(2)
         my_Z = COORDS(3)
      ENDIF

      AX = NINT(PROPS(42))
      IF (AX .EQ. 1) THEN
         axial_pos = my_X
      ELSE IF (AX .EQ. 2) THEN
         axial_pos = my_Y
      ELSE
         axial_pos = my_Z
      ENDIF

      x_m  = (axial_pos - x0) * x_scale_m
      xbar = x_m / Ltube
      IF (xbar < 0.0D0) xbar = 0.0D0
      IF (xbar > 1.0D0) xbar = 1.0D0

      ! --- Perfil axial de flujo cos² cuando PROPS(38) < 0 ---
      !     Pico = |PROPS(38)|, ceros cerca de los extremos del tubo.
      !     Factor 0.92 → flujo ~1.5 % del pico en xbar=0 y xbar=1.
      IF (flux_SI < 0.0D0) THEN
         flux_SI = ABS(flux_SI) * COS(3.14159265358979D0 &
                   * (xbar - 0.5D0) * 0.92D0)**2
      ENDIF

      R_CYL = 0.0D0
      IF (AX .EQ. 3) THEN
         rad = SQRT(my_X**2 + my_Y**2)
         IF (rad > 1.0D-9) THEN
            R_CYL(1,3) = 1.0D0                         ! Axial  = ê_z (global)
            R_CYL(2,1) = -my_Y/rad; R_CYL(2,2) = my_X/rad  ! Hoop
            R_CYL(3,1) =  my_X/rad; R_CYL(3,2) = my_Y/rad  ! Radial
         ELSE
            R_CYL(1,3) = 1.0D0; R_CYL(2,2) = 1.0D0; R_CYL(3,1) = 1.0D0
         ENDIF
      ELSE
         R_CYL(1,1) = 1.0D0; R_CYL(2,2) = 1.0D0; R_CYL(3,3) = 1.0D0
      ENDIF

      ! Rotar entradas a sistema local [FIX1: fórmulas corregidas en subrutinas]
      CALL ROT_STRESS_G2L(STRESS, STR_LOC, R_CYL, NTENS)
      CALL ROT_STRAIN_G2L(DSTRAN, DST_LOC, R_CYL, NTENS)

      ! ----------------------------------------------------------------
      ! 3. ELASTICIDAD ORTÓTROPA (COG-00-092)
      ! ----------------------------------------------------------------
      E_a = PROPS(1); E_t = PROPS(2); E_r = PROPS(3)
      nu_at = PROPS(4); nu_ar = PROPS(5); nu_tr = PROPS(6)
      G_at = PROPS(7); G_ar = PROPS(8); G_tr = PROPS(9)

      CALL GET_ELASTIC_STIFFNESS_ORTO(E_a,E_t,E_r, nu_at,nu_ar,nu_tr, &
                                      G_at,G_ar,G_tr, D_EL_LOC, NTENS)

      ! ----------------------------------------------------------------
      ! 4. PREDICTOR ELÁSTICO
      !    σ_trial = σ_old + D_el : Δε_total
      ! ----------------------------------------------------------------
      DO I = 1, NTENS
         SIG_TRIAL(I) = STR_LOC(I)
         DO J = 1, NTENS
            SIG_TRIAL(I) = SIG_TRIAL(I) + D_EL_LOC(I,J)*DST_LOC(J)
         ENDDO
      ENDDO

      ! Sin creep o sin flujo: solución elástica pura
      IF (PROPS(45) < 0.5D0 .OR. dt_h < 1.0D-9 .OR. flux_SI < 1.0D-12) THEN
         STR_LOC = SIG_TRIAL
         CALL ROT_STRESS_L2G(STR_LOC, STRESS, R_CYL, NTENS)
         CALL ROT_STIFFNESS_ENG_L2G(D_EL_LOC, DDSDDE, R_CYL, NTENS)
         ! [FIX3] fluencia correcta: flux[n/m²/s] * dt[h] * 3600[s/h]
         STATEV(1) = fluence + flux_SI * dt_h * S3600
         RETURN
      ENDIF

      ! ----------------------------------------------------------------
      ! 5. GROWTH (strain-driven, tratado en el predictor)
      !    Δε_growth = ε̇_growth(σ_trial) · Δt  [solo para predictor inicial]
      ! ----------------------------------------------------------------
      CALL CALC_RATES(STR_LOC, PROPS, NPROPS, x_m, xbar, flux_SI, temp_K, fluence, &
                      RATE_CREEP, RATE_GROWTH, NTENS, .FALSE.)

      DO I = 1, NTENS
         DST_GROWTH(I) = RATE_GROWTH(I) * dt_h
      ENDDO

      ! σ_trial_eff = σ_trial - D_el : Δε_growth
      DO I = 1, NTENS
         DO J = 1, NTENS
            SIG_TRIAL(I) = SIG_TRIAL(I) - D_EL_LOC(I,J)*DST_GROWTH(J)
         ENDDO
      ENDDO

      ! ----------------------------------------------------------------
      ! 6. NEWTON-RAPHSON — Integración Implícita Backward Euler
      !    Residual: R(σₙ₊₁) = σₙ₊₁ - σ_trial_eff + D_el:ε̇_cr(σₙ₊₁)·Δt = 0
      ! ----------------------------------------------------------------
      SIG_CURR = SIG_TRIAL   ! Valor inicial del iterador

      ! Protección: esfuerzos muy bajos → sin creep
      IF (MAXVAL(ABS(SIG_TRIAL)) < 0.1D0) THEN
         STR_LOC = SIG_TRIAL
         CALL ROT_STRESS_L2G(STR_LOC, STRESS, R_CYL, NTENS)
         CALL ROT_STIFFNESS_ENG_L2G(D_EL_LOC, DDSDDE, R_CYL, NTENS)
         ! [FIX3]
         STATEV(1) = fluence + flux_SI * dt_h * S3600
         RETURN
      ENDIF

      DO ITER = 1, MAX_ITER

         ! A) Tasas de creep evaluadas en σₖ (solo creep, no growth)
         CALL CALC_RATES(SIG_CURR, PROPS, NPROPS, x_m, xbar, flux_SI, temp_K, fluence, &
                         RATE_CREEP, RATE_GROWTH, NTENS, .TRUE.)
         DO I = 1, NTENS
            DST_CREEP(I) = RATE_CREEP(I) * dt_h
         ENDDO

         ! B) Residual: R = σₖ - σ_trial_eff + D_el:Δε_cr
         DO I = 1, NTENS
            RESID(I) = SIG_CURR(I) - SIG_TRIAL(I)
            DO J = 1, NTENS
               RESID(I) = RESID(I) + D_EL_LOC(I,J)*DST_CREEP(J)
            ENDDO
         ENDDO

         ! Convergencia (norma infinito)
         IF (MAXVAL(ABS(RESID)) < TOL_RES) EXIT

         IF (ITER == MAX_ITER) THEN
            PNEWDT = 0.25D0   ! Newton no convergió: cortar paso
            RETURN
         ENDIF

         ! C) Jacobiano numérico por diferencias finitas hacia adelante
         !    J = I + D_el · Δt · (∂ε̇_cr/∂σ)
         JAC = 0.0D0
         DO I = 1, NTENS
            JAC(I,I) = 1.0D0
         ENDDO

         DO J = 1, NTENS
            VAL_PERT  = SIG_CURR(J)
            H_PERT    = MAX(ABS(VAL_PERT)*EPS_PERT, 1.0D-6)
            SIG_PERT  = SIG_CURR
            SIG_PERT(J) = VAL_PERT + H_PERT

            CALL CALC_RATES(SIG_PERT, PROPS, NPROPS, x_m, xbar, flux_SI, temp_K, fluence, &
                            RATE_PERT, RATE_GROWTH, NTENS, .TRUE.)

            DO I = 1, NTENS
               DO K = 1, NTENS
                  JAC(I,J) = JAC(I,J) + D_EL_LOC(I,K) * &
                             ((RATE_PERT(K) - RATE_CREEP(K)) / H_PERT) * dt_h
               ENDDO
            ENDDO
         ENDDO

         ! D) Resolver: Δσ = -J⁻¹ · R  [FIX4: pivoteo parcial en MAT_INV_6X6]
         CALL MAT_INV_6X6(JAC, INV_JAC)

         DO I = 1, NTENS
            DELTA_SIG(I) = 0.0D0
            DO J = 1, NTENS
               DELTA_SIG(I) = DELTA_SIG(I) - INV_JAC(I,J)*RESID(J)
            ENDDO
            SIG_CURR(I) = SIG_CURR(I) + DELTA_SIG(I)
         ENDDO

      ENDDO   ! Fin Newton-Raphson

      ! ----------------------------------------------------------------
      ! 7. ACTUALIZACIÓN DE VARIABLES DE ESTADO
      ! ----------------------------------------------------------------
      ! [FIX3] Fluencia en unidades correctas: flux [n/m²/s] * dt [h] * 3600 [s/h]
      STATEV(1) = fluence + flux_SI * dt_h * S3600

      ! Deformaciones inelásticas acumuladas (indices 2 a 7, si hay espacio)
      IF (NSTATV >= 7) THEN
         DO I = 1, NTENS
            STATEV(1+I) = STATEV(1+I) + (DST_CREEP(I) + DST_GROWTH(I))
         ENDDO
      ENDIF

      ! Tangente consistente local: D_cons = J⁻¹ · D_el
      DO I = 1, NTENS
         DO J = 1, NTENS
            D_CONS_LOC(I,J) = 0.0D0
            DO K = 1, NTENS
               D_CONS_LOC(I,J) = D_CONS_LOC(I,J) + INV_JAC(I,K)*D_EL_LOC(K,J)
            ENDDO
         ENDDO
      ENDDO

      ! ----------------------------------------------------------------
      ! 8. ROTACIÓN DE VUELTA A GLOBAL [FIX1: fórmulas corregidas]
      ! ----------------------------------------------------------------
      CALL ROT_STRESS_L2G(SIG_CURR, STRESS, R_CYL, NTENS)
      CALL ROT_STIFFNESS_ENG_L2G(D_CONS_LOC, DDSDDE, R_CYL, NTENS)

      RETURN
      END SUBROUTINE UMAT


!=======================================================================
!  SUBRUTINAS FÍSICAS Y MATEMÁTICAS
!=======================================================================

      SUBROUTINE CALC_RATES(S, P, NP, xm, xb, flx, T, flu, RC, RG, NT, DO_C)
      ! Tasas de Creep (RC) y Growth (RG) según CRNL-4003 Rev.1
      ! [FIX2] Añadidas tasas de cizallamiento RC(4:6) — antes eran cero.
      !        Derivada Hill para cizallamiento (s_eff cancela):
      !        ε̇_τ_n1  = K1·2·L·τ·exp(-Q1/T)
      !        ε̇_τ_n2  = K2·2·L·τ·s2·exp(-Q1/T)   [s2 no cancela en n=2]
      !        ε̇_τ_K3  = K3·2·L·τ·exp(-Q3/T)
      !        ε̇_τ_irr = Kc·K4x·2·L·τ·Φ·(exp(-Q4/T)+K5)
      !        RG(4:6) = 0 (growth no tiene componentes de cizallamiento)
      !
      ! Voigt local: (1=axial, 2=trans, 3=radial, 4=at-shear, 5=ar-shear, 6=tr-shear)
      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NT, NP
      DOUBLE PRECISION, INTENT(IN) :: S(NT), P(NP)
      DOUBLE PRECISION, INTENT(IN) :: xm, xb, flx, T, flu
      DOUBLE PRECISION, INTENT(OUT) :: RC(NT), RG(NT)
      LOGICAL, INTENT(IN) :: DO_C

      ! Variables internas
      DOUBLE PRECISION :: K1,K2,K3,Q1,Q3, F1,G1,H1, F2,G2,H2
      DOUBLE PRECISION :: Kc,Q4,K5, K41,K42, Fb,Ff,Gb,Gf, F4,G4,H4
      DOUBLE PRECISION :: Kg,Q6, K61,K62, Cg,Bg
      DOUBLE PRECISION :: Ga_b,Ga_f,Gt_b,Gt_f, C6a,C6t,C6r
      DOUBLE PRECISION :: sa, st, sr, K4x, K6x
      DOUBLE PRECISION :: xb_loc, xm_loc
      DOUBLE PRECISION :: s1, s2, s4, C1a,C1t,C1r, C2a,C2t,C2r, C4a,C4t,C4r
      DOUBLE PRECISION :: s2_sq, s2_safe
      DOUBLE PRECISION :: exp_Q1, exp_Q3, exp_Q4, exp_Q6

      ! [FIX2] Variables para cizallamiento
      DOUBLE PRECISION :: tau_at, tau_ar, tau_tr
      DOUBLE PRECISION :: TH1L_p, TH1M_p, TH1N_p
      DOUBLE PRECISION :: TH2L_p, TH2M_p, TH2N_p
      DOUBLE PRECISION :: CREEPLB_p, CREEPLF_p
      DOUBLE PRECISION :: CREEPMB_p, CREEPMF_p
      DOUBLE PRECISION :: CREEPNB_p, CREEPNF_p
      DOUBLE PRECISION :: CREEPL, CREEPM, CREEPN

      RC = 0.0D0; RG = 0.0D0

      ! Guardas de entrada
      IF (.NOT. IEEE_IS_FINITE(T))   RETURN
      IF (T <= 1.0D-12)              RETURN
      IF (.NOT. IEEE_IS_FINITE(flx) .OR. .NOT. IEEE_IS_FINITE(flu) .OR. &
          .NOT. IEEE_IS_FINITE(xm)  .OR. .NOT. IEEE_IS_FINITE(xb)) RETURN
      IF (ABS(flx) < 1.0D-30)       RETURN

      ! Stresses normales locales (1=axial, 2=trans, 3=radial)
      IF (.NOT. IEEE_IS_FINITE(S(1)) .OR. .NOT. IEEE_IS_FINITE(S(2)) .OR. &
          .NOT. IEEE_IS_FINITE(S(3))) RETURN
      sa = S(1); st = S(2); sr = S(3)

      ! Stresses cizallantes locales (4=at, 5=ar, 6=tr)
      tau_at = 0.0D0; tau_ar = 0.0D0; tau_tr = 0.0D0
      IF (NT >= 4) THEN
         IF (IEEE_IS_FINITE(S(4))) tau_at = S(4)
      ENDIF
      IF (NT >= 6) THEN
         IF (IEEE_IS_FINITE(S(5))) tau_ar = S(5)
         IF (IEEE_IS_FINITE(S(6))) tau_tr = S(6)
      ENDIF

      ! Clamps posición axial
      xb_loc = MAX(0.0D0, MIN(1.0D0, xb))
      xm_loc = MAX(-1.0D6, MIN(1.0D6, xm))

      ! Lectura de PROPS (índices 10-37)
      K1=P(10); K2=P(11); K3=P(12); Q1=P(13); Q3=P(14)
      F1=P(15); G1=P(16); H1=1.5D0-F1-G1
      F2=P(17); G2=P(18); H2=1.5D0-F2-G2
      Kc=P(19); Q4=P(20); K5=P(21); K41=P(22); K42=P(23)
      Fb=P(24); Ff=P(25); Gb=P(26); Gf=P(27)
      Kg=P(28); Q6=P(29); K61=P(30); K62=P(31); Cg=P(32); Bg=P(33)
      Ga_b=P(34); Ga_f=P(35); Gt_b=P(36); Gt_f=P(37)

      ! [FIX2] Parámetros de cizallamiento (PROPS 46-57 o defaults CRNL-4003 Rev.1)
      IF (NP >= 57) THEN
         TH1L_p   = P(46); TH1M_p   = P(47); TH1N_p   = P(48)
         TH2L_p   = P(49); TH2M_p   = P(50); TH2N_p   = P(51)
         CREEPLB_p= P(52); CREEPLF_p= P(53)
         CREEPMB_p= P(54); CREEPMF_p= P(55)
         CREEPNB_p= P(56); CREEPNF_p= P(57)
      ELSE
         ! Valores estándar tubos zircaloy CRNL-4003 Rev.1
         TH1L_p   = 0.71D0;  TH1M_p   = 0.72D0;  TH1N_p   = 0.75D0
         TH2L_p   = 0.48D0;  TH2M_p   = 0.47D0;  TH2N_p   = 0.52D0
         CREEPLB_p= 0.475D0; CREEPLF_p= 0.53D0
         CREEPMB_p= 1.363D0; CREEPMF_p= 1.458D0
         CREEPNB_p= 1.154D0; CREEPNF_p= 1.126D0
      ENDIF

      ! ---------------------------------------------------------------
      ! GROWTH (Ec. 3.4, 3.8, 3.9 — stress-independent, RG(4:6) = 0)
      ! ---------------------------------------------------------------
      C6a = Ga_b + (Ga_f - Ga_b)*xb_loc
      C6t = Gt_b + (Gt_f - Gt_b)*xb_loc
      C6r = -C6a - C6t   ! Iso-volumétrico: C6a + C6t + C6r = 0

      IF (ABS(Bg) > 1.0D-100) THEN
         K6x = (K61 + K62*xm_loc) * (1.0D0 + (Cg/Bg)*flu)
      ELSE
         K6x = (K61 + K62*xm_loc) * (1.0D0 + Cg*flu)
      ENDIF

      CALL SAFE_EXP(-Q6/T, exp_Q6)
      RG(1) = Kg * K6x * C6a * flx * exp_Q6
      RG(2) = Kg * K6x * C6t * flx * exp_Q6
      RG(3) = Kg * K6x * C6r * flx * exp_Q6
      ! RG(4:6) = 0 — growth no tiene cizallamiento

      IF (.NOT. DO_C) RETURN

      ! ---------------------------------------------------------------
      ! CREEP TÉRMICO (Ec. 3.2)
      ! ---------------------------------------------------------------
      CALL HILL_CALC(F1,G1,H1, sa,st,sr, s1, C1a,C1t,C1r)
      CALL HILL_CALC(F2,G2,H2, sa,st,sr, s2, C2a,C2t,C2r)

      ! Clamp s2 para evitar overflow en rc_shear_n2 (s2 <= s2_safe)
      s2_safe = s2
      IF (.NOT. IEEE_IS_FINITE(s2_safe)) s2_safe = 0.0D0
      IF (s2_safe > 1.0D50) s2_safe = 1.0D50

      s2_sq = s2*s2
      IF (s2_sq .NE. s2_sq) s2_sq = 0.0D0
      IF (s2_sq > 1.0D100)  s2_sq = 1.0D100

      CALL SAFE_EXP(-Q1/T, exp_Q1)
      CALL SAFE_EXP(-Q3/T, exp_Q3)

      ! Tasas normales
      RC(1) = (K1*C1a*s1 + K2*C2a*s2_sq)*exp_Q1 + K3*C1a*s1*exp_Q3
      RC(2) = (K1*C1t*s1 + K2*C2t*s2_sq)*exp_Q1 + K3*C1t*s1*exp_Q3
      RC(3) = (K1*C1r*s1 + K2*C2r*s2_sq)*exp_Q1 + K3*C1r*s1*exp_Q3

      ! [FIX2] Tasas de cizallamiento térmico
      !  Derivando  seff respecto al cizallamiento: ∂s/∂τ = 2L·τ/seff → seff cancela
      !  n=1:  ε̇ = K1·(2L·τ/s1)·s1·exp = K1·2L·τ·exp   (s1 cancela)
      !  n=2:  ε̇ = K2·(2L·τ/s2)·s2²·exp = K2·2L·τ·s2·exp  (una s2 queda)
      !  K3:   ε̇ = K3·(2L·τ/s1)·s1·exp = K3·2L·τ·exp   (s1 cancela)
      IF (NT >= 4) THEN
         RC(4) = (K1*2.0D0*TH1L_p*tau_at + K2*2.0D0*TH2L_p*s2_safe*tau_at)*exp_Q1 &
                 + K3*2.0D0*TH1L_p*tau_at*exp_Q3
      ENDIF
      IF (NT >= 6) THEN
         RC(5) = (K1*2.0D0*TH1M_p*tau_ar + K2*2.0D0*TH2M_p*s2_safe*tau_ar)*exp_Q1 &
                 + K3*2.0D0*TH1M_p*tau_ar*exp_Q3
         RC(6) = (K1*2.0D0*TH1N_p*tau_tr + K2*2.0D0*TH2N_p*s2_safe*tau_tr)*exp_Q1 &
                 + K3*2.0D0*TH1N_p*tau_tr*exp_Q3
      ENDIF

      ! ---------------------------------------------------------------
      ! CREEP IRRADIACIÓN (Ec. 3.3)
      ! ---------------------------------------------------------------
      F4 = Fb + (Ff - Fb)*xb_loc
      G4 = Gb + (Gf - Gb)*xb_loc
      H4 = 1.5D0 - F4 - G4
      K4x = K41 + K42*xm_loc

      CALL HILL_CALC(F4,G4,H4, sa,st,sr, s4, C4a,C4t,C4r)
      IF (s4 > 1.0D10)  s4 = 1.0D10
      IF (s4 < -1.0D10) s4 = -1.0D10

      CALL SAFE_EXP(-Q4/T, exp_Q4)

      ! Tasas normales irradiación
      RC(1) = RC(1) + Kc * K4x * C4a * s4 * flx * (exp_Q4 + K5)
      RC(2) = RC(2) + Kc * K4x * C4t * s4 * flx * (exp_Q4 + K5)
      RC(3) = RC(3) + Kc * K4x * C4r * s4 * flx * (exp_Q4 + K5)

      ! [FIX2] Tasas cizallamiento irradiación
      !  ε̇_τ_irr = Kc·K4x·(2L·τ/s4)·s4·Φ·(exp+K5) = Kc·K4x·2L·τ·Φ·(exp+K5)
      CREEPL = CREEPLB_p + (CREEPLF_p - CREEPLB_p)*xb_loc
      CREEPM = CREEPMB_p + (CREEPMF_p - CREEPMB_p)*xb_loc
      CREEPN = CREEPNB_p + (CREEPNF_p - CREEPNB_p)*xb_loc

      IF (NT >= 4) THEN
         RC(4) = RC(4) + Kc * K4x * 2.0D0*CREEPL * tau_at * flx * (exp_Q4 + K5)
      ENDIF
      IF (NT >= 6) THEN
         RC(5) = RC(5) + Kc * K4x * 2.0D0*CREEPM * tau_ar * flx * (exp_Q4 + K5)
         RC(6) = RC(6) + Kc * K4x * 2.0D0*CREEPN * tau_tr * flx * (exp_Q4 + K5)
      ENDIF

      END SUBROUTINE CALC_RATES


      SUBROUTINE HILL_CALC(F,G,H, sa,st,sr, seff, Ca,Ct,Cr)
      ! Esfuerzo efectivo de Hill (componentes normales) y gradientes de flujo
      ! Nota: los componentes de cizallamiento son tratados por separado en CALC_RATES
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: F,G,H, sa,st,sr
      DOUBLE PRECISION, INTENT(OUT) :: seff, Ca,Ct,Cr
      DOUBLE PRECISION :: val
      DOUBLE PRECISION, PARAMETER :: SEFF_MIN = 1.0D-3  ! 0.001 MPa mínimo

      val  = F*(sa-st)**2 + G*(st-sr)**2 + H*(sr-sa)**2
      seff = SQRT(MAX(val, 1.0D-20))

      IF (seff < SEFF_MIN) THEN
         ! Estado casi hidrostático: gradientes cero
         Ca = 0.0D0; Ct = 0.0D0; Cr = 0.0D0
         seff = SEFF_MIN
      ELSE
         Ca = (F*(sa-st) - H*(sr-sa)) / seff
         Ct = (G*(st-sr) - F*(sa-st)) / seff
         Cr = (H*(sr-sa) - G*(st-sr)) / seff
      ENDIF
      END SUBROUTINE HILL_CALC


      SUBROUTINE SAFE_EXP(arg, result)
      ! Exponencial protegida contra overflow/underflow
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: arg
      DOUBLE PRECISION, INTENT(OUT) :: result
      DOUBLE PRECISION, PARAMETER :: ARG_MAX =  100.0D0
      DOUBLE PRECISION, PARAMETER :: ARG_MIN = -700.0D0

      IF (arg .NE. arg) THEN
         result = 0.0D0
      ELSE IF (arg > ARG_MAX) THEN
         result = EXP(ARG_MAX)
      ELSE IF (arg < ARG_MIN) THEN
         result = 0.0D0
      ELSE
         result = EXP(arg)
      ENDIF
      END SUBROUTINE SAFE_EXP


      SUBROUTINE GET_ELASTIC_STIFFNESS_ORTO(E1,E2,E3,n12,n13,n23,G12,G13,G23, D, NT)
      ! Rigidez ortótropa completa (9 parámetros independientes, COG-00-092)
      ! Invierte la matriz de compliance para obtener la rigidez.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NT
      DOUBLE PRECISION, INTENT(IN)  :: E1,E2,E3,n12,n13,n23,G12,G13,G23
      DOUBLE PRECISION, INTENT(OUT) :: D(NT,NT)
      DOUBLE PRECISION :: S(6,6), C(6,6)
      DOUBLE PRECISION :: E1s, E2s, E3s, G12s, G13s, G23s

      E1s  = MAX(E1,  1.0D3); E2s  = MAX(E2,  1.0D3); E3s  = MAX(E3,  1.0D3)
      G12s = MAX(G12, 1.0D3); G13s = MAX(G13, 1.0D3); G23s = MAX(G23, 1.0D3)

      S = 0.0D0
      S(1,1) =  1.0D0/E1s;  S(1,2) = -n12/E1s; S(1,3) = -n13/E1s
      S(2,1) = S(1,2);       S(2,2) =  1.0D0/E2s; S(2,3) = -n23/E2s
      S(3,1) = S(1,3);       S(3,2) = S(2,3);       S(3,3) =  1.0D0/E3s
      S(4,4) = 1.0D0/G12s
      S(5,5) = 1.0D0/G13s
      S(6,6) = 1.0D0/G23s

      CALL MAT_INV_6X6(S, C)

      D = 0.0D0
      D(1:3,1:3) = C(1:3,1:3)
      D(4,4) = C(4,4)
      IF (NT >= 6) THEN
         D(5,5) = C(5,5); D(6,6) = C(6,6)
      ENDIF
      END SUBROUTINE GET_ELASTIC_STIFFNESS_ORTO


      SUBROUTINE MAT_INV_6X6(A, AINV)
      ! Inversión Gauss-Jordan con PIVOTEO PARCIAL por filas [FIX4]
      ! El pivoteo parcial mejora la estabilidad numérica cuando
      ! el elemento diagonal es pequeño frente a elementos fuera de diagonal.
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: A(6,6)
      DOUBLE PRECISION, INTENT(OUT) :: AINV(6,6)
      DOUBLE PRECISION :: M(6,12), fac, tmp(12)
      INTEGER :: i, j, k, imax

      M = 0.0D0
      DO i = 1, 6
         DO j = 1, 6
            M(i,j) = A(i,j)
         ENDDO
         M(i,i+6) = 1.0D0
      ENDDO

      DO k = 1, 6
         ! [FIX4] Búsqueda del pivote máximo en la columna k (filas k..6)
         imax = k
         DO i = k+1, 6
            IF (ABS(M(i,k)) > ABS(M(imax,k))) imax = i
         ENDDO
         ! Intercambio de filas k e imax
         IF (imax /= k) THEN
            tmp(1:12)     = M(k,1:12)
            M(k,1:12)     = M(imax,1:12)
            M(imax,1:12)  = tmp(1:12)
         ENDIF

         IF (ABS(M(k,k)) < 1.0D-12) THEN
            ! Matriz singullar: usar identidad como fallback seguro
            AINV = 0.0D0
            DO i = 1, 6; AINV(i,i) = 1.0D0; ENDDO
            RETURN
         ENDIF

         ! Normalizar fila k
         fac = 1.0D0 / M(k,k)
         DO j = 1, 12
            M(k,j) = M(k,j) * fac
         ENDDO
         ! Eliminar columna k en las demás filas
         DO i = 1, 6
            IF (i /= k) THEN
               fac = M(i,k)
               DO j = 1, 12
                  M(i,j) = M(i,j) - fac*M(k,j)
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      DO i = 1, 6
         DO j = 1, 6
            AINV(i,j) = M(i,j+6)
         ENDDO
      ENDDO
      END SUBROUTINE MAT_INV_6X6


      SUBROUTINE ROT_STRESS_G2L(VG, VL, R, NT)
      ! Rota tensor de esfuerzos: Global → Local
      !
      ! [FIX1] CORRECCIÓN: La fórmula correcta para Global→Local es
      !        T_loc = R · T_glob · Rᵀ
      ! donde las FILAS de R son los vectores base locales en el sistema global.
      ! El código original tenía TRANSPOSE(R)·T·R que es la transformación inversa (L→G).
      !
      ! Voigt:  (σ_11, σ_22, σ_33, σ_12, σ_13, σ_23)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NT
      DOUBLE PRECISION, INTENT(IN)  :: VG(NT), R(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: VL(NT)
      DOUBLE PRECISION :: T(3,3), TR(3,3)

      T = 0.0D0
      T(1,1) = VG(1); T(2,2) = VG(2); T(3,3) = VG(3)
      IF (NT >= 4) THEN; T(1,2) = VG(4); T(2,1) = VG(4); ENDIF
      IF (NT >= 6) THEN
         T(1,3) = VG(5); T(3,1) = VG(5)
         T(2,3) = VG(6); T(3,2) = VG(6)
      ENDIF

      ! [FIX1] T_local = R · T_global · Rᵀ
      TR = MATMUL(R, MATMUL(T, TRANSPOSE(R)))

      VL = 0.0D0
      VL(1) = TR(1,1); VL(2) = TR(2,2); VL(3) = TR(3,3)
      IF (NT >= 4) VL(4) = TR(1,2)
      IF (NT >= 6) THEN; VL(5) = TR(1,3); VL(6) = TR(2,3); ENDIF
      END SUBROUTINE ROT_STRESS_G2L


      SUBROUTINE ROT_STRESS_L2G(VL, VG, R, NT)
      ! Rota tensor de esfuerzos: Local → Global
      !
      ! [FIX1] CORRECCIÓN: La fórmula correcta para Local→Global es
      !        T_glob = Rᵀ · T_loc · R
      ! El código original tenía R·T·Rᵀ que es la transformación inversa (G→L).
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NT
      DOUBLE PRECISION, INTENT(IN)  :: VL(NT), R(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: VG(NT)
      DOUBLE PRECISION :: T(3,3), TR(3,3)

      T = 0.0D0
      T(1,1) = VL(1); T(2,2) = VL(2); T(3,3) = VL(3)
      IF (NT >= 4) THEN; T(1,2) = VL(4); T(2,1) = VL(4); ENDIF
      IF (NT >= 6) THEN
         T(1,3) = VL(5); T(3,1) = VL(5)
         T(2,3) = VL(6); T(3,2) = VL(6)
      ENDIF

      ! [FIX1] T_global = Rᵀ · T_local · R
      TR = MATMUL(TRANSPOSE(R), MATMUL(T, R))

      VG = 0.0D0
      VG(1) = TR(1,1); VG(2) = TR(2,2); VG(3) = TR(3,3)
      IF (NT >= 4) VG(4) = TR(1,2)
      IF (NT >= 6) THEN; VG(5) = TR(1,3); VG(6) = TR(2,3); ENDIF
      END SUBROUTINE ROT_STRESS_L2G


      SUBROUTINE ROT_STRAIN_G2L(VG, VL, R, NT)
      ! Rota tensor de deformaciones: Global → Local
      ! Maneja la convención ingeniería (γ = 2ε tensorial para los shear).
      !
      ! [FIX1] CORRECCIÓN: misma fórmula que ROT_STRESS_G2L.
      !        T_loc = R · T_tens · Rᵀ  (T_tens con shear tensorial ε=γ/2)
      !        Resultado convertido de vuelta a γ al extraer componentes.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NT
      DOUBLE PRECISION, INTENT(IN)  :: VG(NT), R(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: VL(NT)
      DOUBLE PRECISION :: T(3,3), TR(3,3)

      T = 0.0D0
      T(1,1) = VG(1); T(2,2) = VG(2); T(3,3) = VG(3)
      IF (NT >= 4) THEN
         T(1,2) = 0.5D0*VG(4); T(2,1) = 0.5D0*VG(4)   ! γ → ε=γ/2
      ENDIF
      IF (NT >= 6) THEN
         T(1,3) = 0.5D0*VG(5); T(3,1) = 0.5D0*VG(5)
         T(2,3) = 0.5D0*VG(6); T(3,2) = 0.5D0*VG(6)
      ENDIF

      ! [FIX1] T_local = R · T_global · Rᵀ
      TR = MATMUL(R, MATMUL(T, TRANSPOSE(R)))

      VL = 0.0D0
      VL(1) = TR(1,1); VL(2) = TR(2,2); VL(3) = TR(3,3)
      IF (NT >= 4) VL(4) = 2.0D0*TR(1,2)               ! ε → γ=2ε
      IF (NT >= 6) THEN
         VL(5) = 2.0D0*TR(1,3)
         VL(6) = 2.0D0*TR(2,3)
      ENDIF
      END SUBROUTINE ROT_STRAIN_G2L


      SUBROUTINE ROT_STIFFNESS_ENG_L2G(DL, DG, R, NT)
      ! Rota tensor de rigidez de 4to orden: Local → Global
      ! Procedimiento numérico: aplica el mapa columna por columna.
      ! Con ROT_STRAIN_G2L y ROT_STRESS_L2G ya corregidos [FIX1], esta
      ! subrutina calcula automáticamente la rotación correcta.
      !
      ! DG_{ij,ab} = R^T_{im} R^T_{jn} DL_{mn,pq} R_{ap} R_{bq}
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NT
      DOUBLE PRECISION, INTENT(IN)  :: DL(NT,NT), R(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: DG(NT,NT)
      DOUBLE PRECISION :: eps_g(NT), eps_l(NT), sig_l(NT), sig_g(NT)
      INTEGER :: j, k

      DO j = 1, NT
         eps_g    = 0.0D0
         eps_g(j) = 1.0D0
         CALL ROT_STRAIN_G2L(eps_g, eps_l, R, NT)   ! G→L [FIX1 aplicado]

         sig_l = 0.0D0
         DO k = 1, NT
            sig_l = sig_l + DL(:,k)*eps_l(k)
         ENDDO

         CALL ROT_STRESS_L2G(sig_l, sig_g, R, NT)   ! L→G [FIX1 aplicado]
         DG(:,j) = sig_g
      ENDDO
      END SUBROUTINE ROT_STIFFNESS_ENG_L2G
