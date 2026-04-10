# Tubo con presión interna — Code_Aster 17.4 + UMAT con coords de Gauss actualizadas

## Descripción del caso

Tubo cilíndrico de pared gruesa con ambos extremos completamente empotrados
(DX = DY = DZ = 0) y presión interna aplicada de forma incremental.
El tubo "se hincha" radialmente.

La UMAT recibe en cada punto de Gauss las coordenadas actuales del punto,
inyectadas mediante `STATEV`. La técnica es idéntica a la usada en
`estudio_cubo_coords/cubo_iter_coords.comm`.

---

## Geometría

```
        z = L = 2.0
         ┌──────┐
         │  ro  │  ← cara libre exterior (no tiene CC de desplazamiento
         │  ri  │                         ni carga)
         │      │
     ████│██████│████  ← Face_sup: DX=DY=DZ=0 (empotrado)
         │      │
     ←── P ────→      ← Face_int: PRES=50 (presión interna)
         │      │
     ████│██████│████  ← Face_inf: DX=DY=DZ=0 (empotrado)
         │      │
        z = 0
```

| Parámetro          | Valor |
|--------------------|-------|
| Radio interior ri  | 0.5   |
| Radio exterior ro  | 1.0   |
| Longitud L         | 2.0   |
| Relación ro/ri     | 2.0   |

---

## Malla

Malla hexaédrica estructurada generada con `gen_tubo_med.py` usando
**MEDCoupling** (API dentro del contenedor Docker de Code_Aster).

| Parámetro        | Valor |
|------------------|-------|
| Sectores (θ)     | 16    |
| Capas radiales   | 2     |
| Rebanadas axiales| 8     |
| **Total nodos**  | **432** |
| **Total HE8**    | **256** |
| QU4 Face_inf     | 32    |
| QU4 Face_sup     | 32    |
| QU4 Face_int     | 128   |

### Grupos de elementos

| Nombre     | Contenido                          | Uso en .comm                    |
|------------|------------------------------------|---------------------------------|
| `VolAll`   | Todos los hexaedros HE8            | AFFE_MODELE, AFFE_MATERIAU      |
| `Face_inf` | Cara z=0 (32 QU4)                  | DDL_IMPO DX=DY=DZ=0             |
| `Face_sup` | Cara z=L (32 QU4)                  | DDL_IMPO DX=DY=DZ=0             |
| `Face_int` | Superficie interior r=ri (128 QU4) | PRES_REP PRES=50                |

### Indexado de nodos

Los nodos se indexan en orden `(ir, iθ, iz)`:

```
node_idx(ir, ith, iz) = ir × (ntheta × (nz+1)) + (ith % ntheta) × (nz+1) + iz
```

El orden local del hexaedro (Jacobiano positivo) es:

```
cara z=iz :  (ir,ith), (ir+1,ith), (ir+1,ith+1), (ir,ith+1)
cara z=iz+1: igual orden
```

---

## Material

Elasticidad lineal isótropa, implementada en una UMAT Fortran 77:

| Parámetro         | Valor  |
|-------------------|--------|
| Módulo de Young E | 1000.0 |
| Coeficiente nu    | 0.30   |

La UMAT es puramente elástica lineal (Voigt, 3D general):

$$
\boldsymbol{\sigma}_{n+1} = \boldsymbol{\sigma}_n + \mathbb{C} : \Delta\boldsymbol{\varepsilon}
$$

$$
\mathbb{C} = \begin{pmatrix}
\lambda+2\mu & \lambda & \lambda & & & \\
\lambda & \lambda+2\mu & \lambda & & & \\
\lambda & \lambda & \lambda+2\mu & & & \\
 & & & \mu & & \\
 & & & & \mu & \\
 & & & & & \mu
\end{pmatrix}
$$

con $\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}$, $\mu = \frac{E}{2(1+\nu)}$.

---

## Condiciones de contorno y carga

| Superficie  | Condición                                         |
|-------------|---------------------------------------------------|
| `Face_inf`  | Empotrado: DX = DY = DZ = 0                       |
| `Face_sup`  | Empotrado: DX = DY = DZ = 0                       |
| `Face_int`  | Presión interna: PRES = 50 × RAMPA(t)             |
| Exterior    | Libre (sin CC ni cargas)                          |

La rampa de carga es lineal: `RAMPA(t) = t` para `t ∈ [0, 1]`.

**Convención de signo de la presión en Code_Aster:**
`PRES_REP` aplica una fuerza $\mathbf{f} = -P \cdot \hat{n}_{ext}$, donde
$\hat{n}_{ext}$ es la normal exterior a la cara. La cara `Face_int` tiene
normal apuntando hacia el eje (−r), por lo que `PRES > 0` produce una fuerza
en +r → el tubo se hincha.

---

## Coordenadas de puntos de Gauss en STATEV

Code_Aster **no** pasa las coordenadas del punto de Gauss a la UMAT:

- `COORDS` → marcado "Nondefinite" en el manual (llega siempre como cero)
- `DFGRD0`, `DFGRD1` → ídem, cero tanto en `PETIT` como en `GDEF_LOG`

### Solución: inyección via STATEV

Se usa el patrón de `PT_test.comm` de Code_Aster:

```
GEOM0 (coords nodales de la malla sin deformar)
  +
DEPL_I (desplaz. nodales convergidos del paso anterior)
  = GEOM_UPD (coords nodales actualizadas)
      ↓  CREA_CHAMP DISC
  GEOM_GP_N (ELGA_GEOM_R — interpoladas a puntos de Gauss)
      ↓  CREA_CHAMP EVAL  (formulas identidad f_X='X', f_Y='Y', f_Z='Z')
  CH_XYZ_N (ELGA_NEUT_R — mismo valor, tipo neutro)
      ↓  CREA_CHAMP ASSE  (X1→V1, X2→V2, X3→V3, X1→V4, X2→V5, X3→V6)
  VINI_I (ELGA_VARI_R — listo para ETAT_INIT)
```

### Layout de STATEV (NSTATV = 6)

| Variable | Contenido                             | Gestión      |
|----------|---------------------------------------|--------------|
| V1       | Coord X (o r cos θ) del GP — actual   | Solo `.comm` |
| V2       | Coord Y (o r sin θ) del GP — actual   | Solo `.comm` |
| V3       | Coord Z del GP — actual               | Solo `.comm` |
| V4       | Coord X del GP — paso anterior        | Solo `.comm` |
| V5       | Coord Y del GP — paso anterior        | Solo `.comm` |
| V6       | Coord Z del GP — paso anterior        | Solo `.comm` |

La UMAT **solo lee** V1–V6 para registro de diagnóstico; no los modifica.
Toda la actualización se hace desde el `.comm`.

---

## Estrategia de resolución

El `.comm` ejecuta un bucle Python sobre `n_steps = 5` incrementos de tiempo
de igual tamaño (`dt = 0.2`):

```
Para cada paso istep ∈ {0..4}:
  1. Construir VINI_I con CH_XYZ_CUR (V1:V3) y CH_XYZ_OLD (V4:V6)
  2. STAT_NON_LINE para el intervalo [istep·dt, (istep+1)·dt]
       - istep=0: ETAT_INIT=_F(VARI=VINI_I)
       - istep>0: ETAT_INIT=_F(EVOL_NOLI=RESU, VARI=VINI_I)
  3. Extraer DEPL convergido del paso
  4. Calcular GEOM_UPD = GEOM0 + DEPL  (campo nodal GEOM_R)
     → DISC → EVAL → CH_XYZ_N  (coords de Gauss deformadas)
  5. Rotar referencias: OLD ← CUR ← NEW
```

---

## Resultados

Corrida en Code_Aster 17.4 con `DEFORMATION='PETIT'`.

| Cantidad                       | Valor        |
|--------------------------------|--------------|
| Desplazamiento radial máx DR   | **0.0460**   |
| Desplazamiento axial máx \|DZ\|| **0.0098**   |
| Diagnóstico final              | `<A>_ALARM`  |

Las alarmas `QUALITY1_3` son informativas: Code_Aster avisa que el campo DEPL
se lee/inicializa desde ETAT_INIT en cada paso del bucle — comportamiento
esperado y correcto.

### Comparación con solución analítica (cilindro libre)

Para un cilindro de pared gruesa **sin restricción axial** (plano de
deformación generalizado), la solución de Lamé es:

$$
u_r(r) = \frac{P r_i^2}{E(r_o^2 - r_i^2)} \left[(1-2\nu)\,r + \frac{r_o^2}{r}\right]
$$

Con los valores del caso (E=1000, ν=0.30, P=50, ri=0.5, ro=1.0):

| Radio   | u_r analítico (libre) |
|---------|-----------------------|
| r = ri  | 0.0367                |
| r = ro  | 0.0233                |

El caso simulado tiene **extremos empotrados** (DZ=0 en ambas tapas), lo que
impide la contracción axial libre que produciría Poisson. El resultado
esperado es una mayor expansión radial en la zona media del tubo respecto
al caso libre, lo que es consistente con el valor obtenido DR_max = 0.046.

---

## Archivos

| Archivo            | Descripción                                        |
|--------------------|----------------------------------------------------|
| `gen_tubo_med.py`  | Genera `Tubo.med` con MEDCoupling (ejecutar en contenedor) |
| `gen_tubo.py`      | Versión alternativa con h5py puro (solo referencia, incompatible con HDF5 1.10.7 de Code_Aster) |
| `Tubo.med`         | Malla hexaédrica del tubo (formato MED/HDF5 1.10.7)|
| `umat_tubo.for`    | UMAT Fortran 77: elasticidad lineal + registro de coords de Gauss |
| `tubo.comm`        | Control de Code_Aster: bucle iterativo de actualización de coords |
| `tubo.export`      | Fichero de configuración de la ejecución           |
| `tubo.rmed`        | Resultados (DEPL, SIEF_ELGA, VARI_ELGA, SIGM_ELNO, EPSI_ELNO) |
| `tubo.message`     | Log de ejecución de Code_Aster                     |

---

## Cómo ejecutar

### 1. Generar la malla

```bash
cd estudio_tubo
docker run --rm -v "$(pwd)":/work simvia/code_aster:17.4.0 bash -c \
  'export PYTHONPATH=/opt/spack/var/spack/environments/simvia_env/.spack-env/view/lib/python3.11/site-packages
   export LD_LIBRARY_PATH=/opt/spack/var/spack/environments/simvia_env/.spack-env/view/lib
   /opt/spack/var/spack/environments/simvia_env/.spack-env/view/bin/python3.11 /work/gen_tubo_med.py'
```

> **Nota:** `gen_tubo_med.py` debe ejecutarse dentro del contenedor porque
> usa MEDCoupling, que está ligado al HDF5 1.10.7 del contenedor. El h5py
> del sistema (HDF5 2.0.0) genera archivos incompatibles con el lector MED
> de Code_Aster 17.4.

### 2. Resolver el caso

```bash
cd estudio_tubo
runaster tubo.export
```

`runaster` es un alias de Docker definido en el entorno:
```bash
runaster() {
    docker run --rm -v "$(pwd)":/work simvia/code_aster:17.4.0 bash -c \
      'export PYTHONPATH=...
       /opt/.../run_aster /work/'"$1"
}
```

### 3. Ver resultados

El archivo `tubo.rmed` puede abrirse con ParaViS (ParaView integrado en
Salome-Meca) o con cualquier lector MED/HDF5 como h5py:

```python
import h5py, numpy as np

with h5py.File("tubo.rmed", "r") as f:
    # Último paso (t=1.0 → paso 5)
    key = "CHA/RESU____DEPL/0000000000000000000500000000000000000005" \
          "/NOE/MED_NO_PROFILE_INTERNAL/CO"
    depl = f[key][:]
    n = len(depl) // 3
    DX, DY, DZ = depl[:n], depl[n:2*n], depl[2*n:]
    DR = np.sqrt(DX**2 + DY**2)
    print(f"DR_max = {DR.max():.6f}")
```

---

## Notas técnicas relevantes

- **h5py vs MEDCoupling**: h5py 3.x usa HDF5 2.0.0 que es incompatible
  con el lector MEDfile de Code_Aster 17.4 (HDF5 1.10.7), aun especificando
  `libver='earliest'`. La única solución fiable es generar la malla con
  **MEDCoupling dentro del contenedor** Docker.

- **DEFORMATION='PETIT'**: única opción práctica para UMAT en Code_Aster
  17.4 que no cause errores de compor. `GROT_GDEP` es incompatible (error
  COMPOR1_44). `GDEF_LOG` añade 6 STATEV extra al final.

- **PRES_REP y orientación de normales**: MEDCoupling genera las normales
  de las caras según el orden de los nodos. Para `Face_int` el orden elegido
  produce una normal en −r, lo que hace que `PRES > 0` empuje hacia +r
  (expansión del tubo). Si el tubo se contrajera en lugar de expandirse,
  habría que invertir el signo de PRES o el orden de los nodos de las caras.

- **Alarma QUALITY1_3**: completamente normal en el bucle iterativo. Code_Aster
  avisa que el campo DEPL se sobreescribe en ETAT_INIT en pasos sucesivos.
  No afecta a la convergencia ni a los resultados.
