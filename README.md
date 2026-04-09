# UMAT en Code_Aster 17.4 — Guía de referencia

## 1. Entorno de ejecución

| Componente | Detalle |
|------------|---------|
| Solver | Code_Aster v17.4.0 (`runaster` en contenedor Docker `simvia-aster:latest`) |
| Compilador Fortran | `gfortran` via `mpif90`, invocado por `make_shared()` |
| Flags implícitos | `-fdefault-integer-8 -fdefault-real-8 -fdefault-double-8` |
| Montaje | El contenedor monta **solo el CWD** como `/work` |
| Ejecución | `cd carpeta && runaster caso.export` |

> **Importante**: symlinks no funcionan dentro del contenedor. Copiar archivos reales.

---

## 2. Estructura de un caso

Cada caso vive en su carpeta con 4 archivos mínimos:

```
mi_caso/
├── Cubo.med            # Malla (copiar, no symlink)
├── mi_umat.for         # Subrutina UMAT en Fortran 77
├── caso.comm           # Archivo de comandos Code_Aster (Python)
└── caso.export         # Descriptor de ejecución
```

### 2.1 Archivo `.export`

```
P time_limit 3000
P memory_limit 1024
P ncpus 1
P mpi_nbcpu 1
P mpi_nbnoeud 1
F comm caso.comm D 1
F med Cubo.med D 20
F rmed caso.rmed R 80
F mess caso.message R 6
```

- `D` = input, `R` = output
- Las UNITE (1, 20, 80, 6) son las que el `.comm` usa por defecto

### 2.2 Malla (`Cubo.med`)

Cubo 1×1×1, HEXA8, 10 nodos, 20 elementos (1 volumen = EL20, FPG8 = 8 GPs).

| Grupo | Tipo | Uso |
|-------|------|-----|
| `Face_inf` | GROUP_MA | Cara Z=0, apoyo empotrado |
| `Face_sup` | GROUP_MA | Cara Z=1, carga aplicada |
| `nodo_sup` | GROUP_NO | Nodo de control para POST_RELEVE_T |
| `nodos_x` | GROUP_NO | Para BC de simetría (DY=0) |
| `nodos_y` | GROUP_NO | Para BC de simetría (DX=0) |

---

## 3. Interfaz UMAT de Code_Aster

### 3.1 Firma de la subrutina

```fortran
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
 RPL,DDSDDT,DRPLDE,DRPLDT,
 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
 KSPT,KSTEP,KINC)
```

### 3.2 Argumentos — qué funciona y qué no

| Argumento | Estado | Notas |
|-----------|--------|-------|
| `STRESS(NTENS)` | **Funciona** | Entrada/salida. Actualizar con tensión al final del incremento |
| `STATEV(NSTATV)` | **Funciona** | Variables internas. Se pueden inyectar desde `.comm` via `ETAT_INIT(VARI=...)` |
| `DDSDDE(NTENS,NTENS)` | **Funciona** | Tangente consistente. Escribirla siempre |
| `STRAN(NTENS)` | **Funciona** | Deformación al inicio del incremento |
| `DSTRAN(NTENS)` | **Funciona** | Incremento de deformación |
| `TIME(2)` | **Funciona** | `TIME(1)` = tiempo al inicio, `TIME(2)` = tiempo total acumulado |
| `DTIME` | **Funciona** | Incremento de tiempo |
| `PROPS(NPROPS)` | **Funciona** | De `UMAT=_F(LISTE_COEF=(...))` |
| `NOEL` | **Funciona** | Número de elemento |
| `NPT` | **Funciona** | Número de punto de Gauss |
| `KINC` | **Funciona** | Iteración de Newton |
| `PNEWDT` | **Funciona** | Poner < 1.0 para pedir paso más corto |
| `COORDS(*)` | **CERO** | Marcado "Nondefinite" en el manual |
| `DFGRD0(3,3)` | **CERO** | Nunca poblado (ni PETIT, ni GDEF_LOG) |
| `DFGRD1(3,3)` | **CERO** | Nunca poblado (ni PETIT, ni GDEF_LOG) |
| `DROT(3,3)` | NO probado | Probablemente cero |
| `TEMP`, `DTEMP` | NO probado | — |

### 3.3 Convención Voigt (3D, NTENS=6)

```
Posición:  1     2     3     4      5      6
Tensor:    11    22    33    12     13     23
Corte:                      γ₁₂    γ₁₃    γ₂₃   (γ = 2ε para corte)
```

> **Atención**: Code_Aster usa orden (12, 13, 23), Abaqus nativo usa (12, 23, 13). Pero en la interfaz UMAT de CA se mantiene la convención Abaqus en la firma, aunque el orden real depende de la versión. Verificar con log si hay dudas.

### 3.4 Deformaciones disponibles con UMAT

| `DEFORMATION` | STRAN/DSTRAN | STRESS | DDSDDE | Notas |
|---------------|-------------|--------|--------|-------|
| `'PETIT'` | Inf. strain ε | Cauchy σ | dσ/dε | Sin corrección geométrica |
| `'PETIT_REAC'` | Inf. strain ε | Cauchy σ | dσ/dε | Con actualización de geometría entre pasos |
| `'GDEF_LOG'` | Log strain (Hencky) | Cauchy σ | dσ/d(log ε) | CA maneja la cinemática. Añade **6 VARI extra** |
| ~~`'GROT_GDEP'`~~ | — | — | — | **INCOMPATIBLE** (error COMPOR1_44) |

> **GDEF_LOG y variables internas**: el campo VARI total tiene `NB_VARI + 6` componentes. Si se usa `ETAT_INIT(VARI=...)`, el campo inyectado debe tener este tamaño, no `NB_VARI`.

---

## 4. Patrones de `.comm`

### 4.1 Caso simple (sin coords, sin bucle)

```python
from run_aster.toolbox import make_shared
DEBUT()
make_shared('/work/libumat.so', '/work/mi_umat.for')

MAIL = LIRE_MAILLAGE(FORMAT='MED', UNITE=20)
MODE = AFFE_MODELE(MAILLAGE=MAIL,
    AFFE=_F(TOUT='OUI', PHENOMENE='MECANIQUE', MODELISATION='3D'))

MATE = DEFI_MATERIAU(
    ELAS=_F(E=1000.0, NU=0.30),
    UMAT=_F(LISTE_COEF=(1000.0, 0.30)),
)
MAT = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT='OUI', MATER=MATE))

RAMPA = DEFI_FONCTION(NOM_PARA='INST', VALE=(0., 0., 1., 1.),
                      PROL_DROITE='LINEAIRE')
CHAR = AFFE_CHAR_MECA(MODELE=MODE,
    DDL_IMPO=_F(GROUP_MA='Face_inf', DX=0., DY=0., DZ=0.),
    FORCE_FACE=_F(GROUP_MA='Face_sup', FX=20.))

LINST = DEFI_LIST_REEL(DEBUT=0., INTERVALLE=_F(JUSQU_A=1., NOMBRE=20))

RESU = STAT_NON_LINE(
    MODELE=MODE, CHAM_MATER=MAT,
    EXCIT=_F(CHARGE=CHAR, FONC_MULT=RAMPA),
    COMPORTEMENT=_F(
        RELATION='UMAT',
        LIBRAIRIE='/work/libumat.so',
        NOM_ROUTINE='umat',
        NB_VARI=1,
        DEFORMATION='PETIT',   # o 'GDEF_LOG'
        TOUT='OUI'),
    INCREMENT=_F(LIST_INST=LINST),
    NEWTON=_F(MATRICE='TANGENTE', REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.e-6, ITER_GLOB_MAXI=50),
    SOLVEUR=_F(METHODE='MUMPS'),
)

# --- Post ---
IMPR_RESU(FORMAT='MED', UNITE=80,
    RESU=_F(RESULTAT=RESU, NOM_CHAM=('DEPL','SIEF_ELGA'), TOUT_ORDRE='OUI'))
FIN()
```

### 4.2 Caso con inyección de coords en STATEV (bucle iterativo)

Cuando se necesitan coordenadas de GP dentro de la UMAT, hay que:

1. **Extraer geometría inicial** → `CREA_CHAMP(EXTR, GEOMETRIE)`
2. **Mapear a GPs** → `DISC` + `FORMULE` + `EVAL`
3. **Ensamblar como VARI** → `CREA_CHAMP(ASSE, ELGA_VARI_R)`
4. **Inyectar** → `ETAT_INIT(VARI=...)`
5. **Actualizar** cada paso → `GEOM0 + DEPL → DISC → EVAL → ASSE`

Pipeline de coordenadas:
```python
# --- Preparación (una sola vez) ---
GEOM0 = CREA_CHAMP(TYPE_CHAM='NOEU_GEOM_R', OPERATION='EXTR',
                    NOM_CHAM='GEOMETRIE', MAILLAGE=MAIL)

f_X = FORMULE(VALE='X', NOM_PARA=('X','Y','Z'))
f_Y = FORMULE(VALE='Y', NOM_PARA=('X','Y','Z'))
f_Z = FORMULE(VALE='Z', NOM_PARA=('X','Y','Z'))

CH_F = CREA_CHAMP(TYPE_CHAM='ELGA_NEUT_F', OPERATION='AFFE',
    MODELE=MODE, PROL_ZERO='OUI',
    AFFE=_F(TOUT='OUI', NOM_CMP=('X1','X2','X3'), VALE_F=(f_X,f_Y,f_Z)))

GEOM_GP0 = CREA_CHAMP(TYPE_CHAM='ELGA_GEOM_R', OPERATION='DISC',
                       MODELE=MODE, CHAM_GD=GEOM0)
CH_XYZ0  = CREA_CHAMP(TYPE_CHAM='ELGA_NEUT_R', OPERATION='EVAL',
                       CHAM_F=CH_F, CHAM_PARA=GEOM_GP0)

# --- En cada paso (dentro del for) ---
VINI_I = CREA_CHAMP(TYPE_CHAM='ELGA_VARI_R', OPERATION='ASSE',
    MODELE=MODE, PROL_ZERO='OUI',
    ASSE=(
        _F(TOUT='OUI', CHAM_GD=CH_XYZ_CUR, NOM_CMP='X1', NOM_CMP_RESU='V1'),
        _F(TOUT='OUI', CHAM_GD=CH_XYZ_CUR, NOM_CMP='X2', NOM_CMP_RESU='V2'),
        _F(TOUT='OUI', CHAM_GD=CH_XYZ_CUR, NOM_CMP='X3', NOM_CMP_RESU='V3'),
        _F(TOUT='OUI', CHAM_GD=CH_XYZ_OLD, NOM_CMP='X1', NOM_CMP_RESU='V4'),
        _F(TOUT='OUI', CHAM_GD=CH_XYZ_OLD, NOM_CMP='X2', NOM_CMP_RESU='V5'),
        _F(TOUT='OUI', CHAM_GD=CH_XYZ_OLD, NOM_CMP='X3', NOM_CMP_RESU='V6'),
    ))

# primera iteración:  ETAT_INIT=_F(VARI=VINI_I)
# siguientes:         ETAT_INIT=_F(EVOL_NOLI=RESU, VARI=VINI_I)

# --- Actualización de coords tras resolver ---
DEPL_I   = CREA_CHAMP(TYPE_CHAM='NOEU_DEPL_R', OPERATION='EXTR',
                       RESULTAT=RESU, NOM_CHAM='DEPL', INST=t_fin)
GEOM_UPD = CREA_CHAMP(TYPE_CHAM='NOEU_GEOM_R', OPERATION='ASSE',
    MAILLAGE=MAIL,
    ASSE=(_F(TOUT='OUI', CHAM_GD=GEOM0,  CUMUL='OUI', COEF_R=1.,
             NOM_CMP=('X','Y','Z')),
          _F(TOUT='OUI', CHAM_GD=DEPL_I, CUMUL='OUI', COEF_R=1.,
             NOM_CMP=('DX','DY','DZ'), NOM_CMP_RESU=('X','Y','Z'))))
GEOM_GP_N = CREA_CHAMP(TYPE_CHAM='ELGA_GEOM_R', OPERATION='DISC',
                        MODELE=MODE, CHAM_GD=GEOM_UPD)
CH_XYZ_N  = CREA_CHAMP(TYPE_CHAM='ELGA_NEUT_R', OPERATION='EVAL',
                        CHAM_F=CH_F, CHAM_PARA=GEOM_GP_N)
```

### 4.3 Bucle `STAT_NON_LINE` con `reuse`

```python
for istep in range(n_steps):
    LINST_I = DEFI_LIST_REEL(DEBUT=t_ini,
                             INTERVALLE=_F(JUSQU_A=t_fin, NOMBRE=1))
    if istep == 0:
        RESU = STAT_NON_LINE(..., ETAT_INIT=_F(VARI=VINI_I),
                             INCREMENT=_F(LIST_INST=LINST_I), ...)
    else:
        RESU = STAT_NON_LINE(reuse=RESU, ...,
                             ETAT_INIT=_F(EVOL_NOLI=RESU, VARI=VINI_I),
                             INCREMENT=_F(LIST_INST=LINST_I), ...)
    DETRUIRE(NOM=LINST_I)
    DETRUIRE(NOM=VINI_I)
```

---

## 5. Esqueleto de la UMAT (Fortran 77)

### 5.1 Inicialización obligatoria

```fortran
C     Siempre inicializar estas salidas:
      DO I=1,NTENS
         DDSDDT(I) = 0.D0
         DRPLDE(I) = 0.D0
      END DO
      DRPLDT = 0.D0
      SPD    = 0.D0
      SCD    = 0.D0
      RPL    = 0.D0
      SSE    = 0.D0
      PNEWDT = 1.D0
```

### 5.2 Material con `DEFORMATION='PETIT'`

Con PETIT, Code_Aster pasa:
- `STRAN` = deformación total al inicio del incremento (Voigt, γ para corte)
- `DSTRAN` = incremento de deformación

Para reconstruir F internamente (necesario para modelos hiperelásticos):
```fortran
C     F = I + sym(STRAN + DSTRAN)  (válido solo para PETIT)
      ETOT(I) = STRAN(I) + DSTRAN(I)
      F(1,1) = 1 + ETOT(1);  F(2,2) = 1 + ETOT(2);  F(3,3) = 1 + ETOT(3)
      F(1,2) = 0.5*ETOT(4);  F(2,1) = 0.5*ETOT(4)
      F(1,3) = 0.5*ETOT(5);  F(3,1) = 0.5*ETOT(5)
      F(2,3) = 0.5*ETOT(6);  F(3,2) = 0.5*ETOT(6)
```

### 5.3 Material con `DEFORMATION='GDEF_LOG'`

Con GDEF_LOG, Code_Aster pasa deformación logarítmica (Hencky). La UMAT es formalmente idéntica a una ley en pequeñas deformaciones: la cinemática de gran deformación la hace Code_Aster internamente.

```fortran
C     sigma = lambda*tr(eps)*I + 2*mu*eps   (pero eps = log strain)
```

### 5.4 Logging desde la UMAT

```fortran
C     Abrir fichero en /work/ (el CWD del contenedor)
      IF (NOEL.EQ.20 .AND. NPT.EQ.1) THEN
         OPEN(UNIT=99,FILE='/work/mi_log.log',
     1        POSITION='APPEND',STATUS='UNKNOWN')
         WRITE(99,'(A,I6,A,I3,A,F8.4)')
     1     ' EL=',NOEL,' GP=',NPT,' t=',TIME(2)
         CLOSE(99)
      END IF
```

Para copiar el log a `REPE_OUT` (que se extrae al CWD al terminar):
```python
# Al final del .comm, antes de FIN()
import shutil, os
src = '/work/mi_log.log'
if os.path.exists(src):
    shutil.copy(src, '/work/REPE_OUT/mi_log.log')
```

---

## 6. Problemas conocidos y soluciones

| Problema | Causa | Solución |
|----------|-------|----------|
| `<F>_ABNORMAL_ABORT` al final | `make_shared` o cleanup de `FIN()` | **Ignorar** si el `.rmed` se generó correctamente |
| `COMPOR1_44: GROT_GDEP incompatible` | UMAT solo soporta PETIT/PETIT_REAC/GDEF_LOG | Usar `GDEF_LOG` para grandes deformaciones |
| DFGRD0/DFGRD1 siempre cero | CA no los implementa para UMAT | Calcular F internamente (ver §5.2) |
| COORDS siempre cero | CA no los implementa para UMAT | Inyectar via STATEV (ver §4.2) |
| SEGFAULT en la UMAT | Orden de argumentos incorrecto o PROPS faltante | Verificar firma exacta (§3.1). Todos los args son obligatorios |
| `MECANONLINE5_82: VARI incohérent` | Tamaño de VARI no coincide con NB_VARI + extras | Con GDEF_LOG usar campo VARI de tamaño NB_VARI+6 |
| Divergencia con F desde fichero | F del paso anterior no es consistente con DSTRAN | Usar F autocontenida: `F = I + sym(STRAN+DSTRAN)` |
| `EXTR_COMP` no existe | Removido en v17.4 | Usar `getValuesWithDescription()` |
| Symlink no funciona | El contenedor copia, no sigue symlinks | Copiar archivos reales |
| FX=50 → det(F) < 0 | Carga excesiva para PETIT | Reducir carga o usar GDEF_LOG |

---

## 7. Casos existentes

### 7.1 Raíz (`/Elastica/`)

| Archivo | Descripción |
|---------|-------------|
| `elastica.for` | UMAT elástica lineal original (E, ν). Autor: M. S. Armoa |
| `stage1.comm` | Caso simple: compresión con PRES_REP, E=210000, ν=0.3 |
| `elastic.export` | Export para stage1 |
| `PT_test.comm` | Caso de test de inyección de coords via STATEV |

### 7.2 `estudio_cubo_coords/` — Neo-Hooke con PETIT

| Archivo | Descripción |
|---------|-------------|
| `umat_defgrad.for` | **Neo-Hooke compresible**. F=I+sym(ε). Cauchy σ = (μ/J)dev(b̄) + κ(J-1)I. Tangente analítico |
| `cubo_defgrad.comm` | Bucle iterativo 20 pasos, FX=20, con pipeline de coords |
| `cubo_defgrad.export` | Export |
| `umat_coords_iter.for` | UMAT elástica lineal con impresión de coords GP |
| `cubo_iter_coords.comm` | Caso iterativo original con coords |
| `plot_gauss_points.py` | Gráfica 3D de evolución de GPs |

### 7.3 `estudio_cubo_grot/` — Hencky con GDEF_LOG

| Archivo | Descripción |
|---------|-------------|
| `umat_grot.for` | **Hencky hiperelástico**: σ = λ·tr(ε)·I + 2μ·ε (ε = log strain). Con logging de DFGRD0/DFGRD1 |
| `cubo_grot.comm` | Caso simple (sin bucle), 20 pasos, GDEF_LOG, NB_VARI=1 |
| `cubo_grot.export` | Export |

**Resultado**: DFGRD0/DFGRD1 confirmados como **cero** también con GDEF_LOG.

---

## 8. Checklist para crear un caso nuevo

1. **Crear carpeta**: `mkdir estudio_cubo_mi_caso`
2. **Copiar malla**: `cp ../estudio_cubo_coords/Cubo.med .` (no symlink)
3. **Escribir `.for`**: copiar esqueleto de §5, modificar la ley constitutiva
4. **Escribir `.comm`**: copiar plantilla de §4.1 (simple) o §4.2 (con coords)
   - Cambiar nombre de la `.so` y `.for`
   - Cambiar `LISTE_COEF` con las PROPS necesarias
   - Ajustar `NB_VARI`, `DEFORMATION`, carga, pasos
5. **Escribir `.export`**: copiar plantilla de §2.1, cambiar nombres
6. **Ejecutar**: `cd estudio_cubo_mi_caso && runaster caso.export`
7. **Verificar**:
   - `grep RESI_GLOB caso.message` → residuos deben ser < 1e-6
   - `grep EXCEPTION caso.message` → no debe haber excepciones reales
   - `ls -la caso.rmed` → debe existir y tener tamaño > 0
   - Si hay log Fortran: `cat mi_log.log` o `cat REPE_OUT/mi_log.log`
