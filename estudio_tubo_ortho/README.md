# estudio_tubo_ortho — Tubo con ortotropía cilíndrica

Caso de referencia que extiende `estudio_tubo/` (isótropo) añadiendo un material
**ortótropo en ejes cilíndricos (r, θ, z)** implementado como UMAT de Fortran 77
para Code_Aster 17.4.

---

## Geometría y malla

| Parámetro | Valor |
|-----------|-------|
| Radio interior `ri` | 0.8 |
| Radio exterior `ro` | 1.0 |
| Longitud `L` | 5.0 |
| Tipo de elemento | HE8 (hexaedros de 8 nodos) |
| Num. elementos | 256 (nθ=16, nr=2, nz=8) |

Malla: `Tubo.med` (misma que `estudio_tubo/`).

---

## Condiciones de contorno y carga

- **Extremos** (`Face_inf`, `Face_sup`): desplazamientos DX=DY=DZ=0 (empotrados).
- **Carga**: presión interna P=50 en `Face_int`, rampa lineal de 0→1 en 5 pasos (dt=0.2).

---

## Material

Ortotropía cilíndrica con 9 constantes elásticas independientes:

| Propiedad | Símbolo     | Valor |
|-----------|-------------|-------|
| Módulo radial | E_r | 500 |
| Módulo circunferencial | E_θ | 1500 |
| Módulo axial | E_z | 1000 |
| Poisson r-θ | ν_rθ | 0.15 |
| Poisson r-z | ν_rz | 0.20 |
| Poisson θ-z | ν_θz | 0.25 |
| Corte r-θ | G_rθ | 300 |
| Corte r-z | G_rz | 350 |
| Corte θ-z | G_θz | 400 |

La rigidez en ejes cilíndricos se calcula a partir de la compliance S⁻¹ estándar
de un material ortótropo con simetría ortorrómbica.

---

## Estrategia de rotación del tensor de rigidez

Se sigue la misma estrategia que `UMAT4001.f90`:

1. Se construye la matriz de rotación **R_CYL** cuyas **filas** son los vectores
   de la base local expresados en coordenadas globales cartesianas:

   ```
   Fila 1 — e_r  = ( x/r,  y/r, 0 )
   Fila 2 — e_θ  = (-y/r,  x/r, 0 )
   Fila 3 — e_z  = ( 0,    0,   1 )
   ```

   con `r = √(x²+y²)` y `(x,y)` = coordenadas actuales del punto de Gauss.

2. Las rotaciones de tensores de 2.º orden (notación tensorial) se aplican como:

   - Global → Local: `T_loc = R · T_glob · Rᵀ`
   - Local → Global: `T_glob = Rᵀ · T_loc · R`

3. Las deformaciones de ingeniería de Voigt (γᵢⱼ = 2εᵢⱼ) se convierten a
   tensoriales **antes** de rotar y se reconvierten **después** (subrut.
   `ROT_STRAIN_G2L6`).

4. La matriz de rigidez cartesiana se obtiene rotando la matriz cilíndrica
   columna a columna mediante los mismos procedimientos de rotación
   (`ROT_STIFF_L2G6`).

### Subrrutinas de la UMAT

| Subrrutina | Función |
|------------|---------|
| `UMAT` (principal) | Lee STATEV(1:2), construye R_CYL, rota ε, actualiza σ, rota C_cyl → C_cart |
| `ROT_STRESS_G2L6` | Rota tensor de tensión Global → Local (T_loc = R·T·Rᵀ) |
| `ROT_STRESS_L2G6` | Rota tensor de tensión Local → Global (T_glob = Rᵀ·T·R) |
| `ROT_STRAIN_G2L6` | Igual que G2L pero con conversión γ↔ε antes/después |
| `ROT_STIFF_L2G6`  | Rota la rigidez 6×6 columna a columna |

---

## Coordenadas en puntos de Gauss

Code_Aster 17.4 **no popula** el argumento `COORDS()` en la interfaz UMAT.
Las coordenadas se inyectan vía `STATEV(1:3)` usando el pipeline del `.comm`:

```
DISC → CREA_CHAMP(FORMULE, x/y/z) → CREA_CHAMP(EVAL) → CREA_CHAMP(ASSE) → ETAT_INIT(VARI)
```

Layout de STATEV (NSTATV=6):

| STATEV | Contenido |
|--------|-----------|
| V1–V3 | Coordenadas del GP en el instante actual |
| V4–V6 | Coordenadas del GP en el instante anterior |

---

## Archivos

| Archivo | Descripción |
|---------|-------------|
| `umat_ortho.for` | UMAT Fortran 77 — ortotropía cilíndrica con rotación R_CYL |
| `tubo_ortho.comm` | Control Code_Aster (compila UMAT, inyecta coords, resuelve) |
| `tubo_ortho.export` | Descripción de ficheros para `runaster` |
| `Tubo.med` | Malla HE8 del tubo (copiada de `estudio_tubo/`) |
| `tubo_ortho.rmed` | Resultados MED (generado al ejecutar) |
| `tubo_ortho.message` | Log de Code_Aster (generado al ejecutar) |
| `umat_ortho.log` | Log de la UMAT: coords. de Gauss por paso y elemento |
| `plot_gauss_points.py` | Script Python para visualizar las trayectorias de los GPs |
| `gauss_points_3d.png` | Vista 3D de los anillos de GPs (generada al ejecutar el script) |
| `gauss_points_cylindrical.png` | Evolución de r(t) y z(t) (generada al ejecutar el script) |

---

## Ejecución

```bash
# Desde este directorio (el directorio se monta como /work en el contenedor)
runaster tubo_ortho.export
```

Resultado esperado: `DIAGNOSTIC JOB : <A>_ALARM` (solo informativo).

### Visualización

```bash
python3 plot_gauss_points.py
```

Genera `gauss_points_3d.png` y `gauss_points_cylindrical.png`.

---

## Resultado esperado

Con E_θ=1500 > E_r=500 (tubo más rígido en hoop que en dirección radial), se
espera una **expansión radial mayor** que en el caso isótropo de referencia
(E=1000 en `estudio_tubo/`), ya que la dirección radial es relativamente blanda.
Los extremos permanecen fijos por las condiciones de contorno empotradas.
