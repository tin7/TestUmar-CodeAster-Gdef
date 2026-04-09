# Coordenadas actualizadas en UMAT — Code_Aster 17.4

## Problema

Code_Aster **no** pasa coordenadas al UMAT de Abaqus:

- `COORDS` → marcado "Nondefinite" en el manual (llega como cero)
- `DFGRD0`, `DFGRD1` → también "Nondefinite" (llegan como cero, tanto en `PETIT` como en `GDEF_LOG`)

## Solución

Se inyectan las coordenadas de los puntos de Gauss mediante `STATEV`, usando
la técnica del caso de test `PT_test.comm` de Code_Aster:

```
GEOM0 + DEPL → GEOM_UPD (coords nodales actualizadas)
    → DISC (interpolar a puntos de Gauss)
    → FORMULE/EVAL (extraer X, Y, Z como NEUT_R)
    → ASSE (ensamblar en ELGA_VARI_R como V1, V2, V3, ...)
    → ETAT_INIT(VARI=...) en STAT_NON_LINE
```

El `.comm` ejecuta un bucle Python sobre los incrementos de tiempo. En cada
paso extrae el desplazamiento, recalcula las coordenadas y las reinyecta en
`STATEV` para el paso siguiente.

La UMAT **lee** las coordenadas de `STATEV(1:3)` y `STATEV(4:6)` — **no las
modifica**. Toda la gestión de coordenadas se hace desde el `.comm`.

## Layout de STATEV (NSTATV = 6)

| Variables | Contenido | Quién escribe |
|-----------|-----------|---------------|
| V1:V3     | Coordenadas Gauss **actuales** (X, Y, Z) | `.comm` |
| V4:V6     | Coordenadas Gauss **del paso anterior** | `.comm` |

## Archivos

| Archivo | Descripción |
|---------|-------------|
| `Cubo.med` | Malla hexaédrica 1×1×1 con grupos `Face_inf`, `Face_sup` |
| `cubo_iter_coords.comm` | Caso de cálculo con bucle iterativo y rampa de carga |
| `cubo_iter_coords.export` | Fichero de control para `run_aster` |
| `umat_coords_iter.for` | UMAT elástica lineal (E=1000, ν=0.30) con impresión de coords |

## Ejecución

```bash
runaster cubo_iter_coords.export
```

Donde `runaster` ejecuta el contenedor Docker `simvia/code_aster:17.4.0`.

## Resultados esperados

Con rampa de carga lineal (FX=120 en `Face_sup`, empotrada `Face_inf`):

- Las coordenadas V1:V3 crecen progresivamente en cada paso de tiempo
- V4:V6 contienen las coordenadas del paso anterior
- La UMAT escribe un log (`umat_coords.log`) con las coords para EL/GP=1:

```
 EL=    20 GP=  1 t=  0.0000 act=    0.211325    0.211325    0.788675 prev=    0.211325    0.211325    0.788675
 EL=    20 GP=  1 t=  0.2000 act=    0.296159    0.214137    0.809241 prev=    0.211325    0.211325    0.788675
 EL=    20 GP=  1 t=  0.4000 act=    0.380994    0.216949    0.829807 prev=    0.296159    0.214137    0.809241
 EL=    20 GP=  1 t=  0.6000 act=    0.465828    0.219761    0.850373 prev=    0.380994    0.216949    0.829807
 EL=    20 GP=  1 t=  0.8000 act=    0.550662    0.222574    0.870938 prev=    0.465828    0.219761    0.850373
```

- El archivo de resultados se guarda en `cubo_iter_coords.rmed`

## Lectura de resultados MED

Los campos `ELGA` en MED se almacenan en orden **component-first** (primero
todos los GPs de V1, luego todos los GPs de V2, etc.):

```python
import h5py, numpy as np

f = h5py.File('cubo_iter_coords.rmed', 'r')
grp = f['/CHA/RESU____VARI_ELGA']
nco = grp.attrs['NCO']  # 6 componentes
ngp = 8                   # puntos de Gauss por elemento (HE8)

for ts in sorted(grp.keys()):
    vals = grp[ts]['MAI.HE8']['MED_NO_PROFILE_INTERNAL']['CO'][:]
    arr = vals.reshape(nco, ngp).T  # → (ngp, nco) = (8, 6)
    pdt = grp[ts].attrs['PDT']
    print(f't={pdt:.2f}, GP0 coords: {arr[0, 0:3]}')

f.close()
```
