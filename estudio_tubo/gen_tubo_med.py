#!/usr/bin/env python3
# coding=utf-8
"""
gen_tubo_med.py — Genera Tubo.med usando MEDCoupling.
Se ejecuta DENTRO del contenedor Docker de Code_Aster:

  docker run --rm -v "$(pwd)":/work simvia/code_aster:17.4.0 bash -c \
    'export PYTHONPATH=...;
     python3.11 /work/gen_tubo_med.py'

Geometría:  ri, ro, L  (radio interior, radio exterior, longitud)

Parámetros de discretización
-----------------------------
  ntheta  divisiones en dirección circunferencial (θ).
          Controla la resolución de la sección transversal del tubo.
          Mínimo recomendado: 16. Mejor simetría con 24 o 32.

  nr      capas de elementos en dirección radial (r, de ri a ro).
          Controla la resolución del gradiente de tensiones a través
          de la pared. Con nr=2 hay 2 elementos entre ri y ro;
          para capturar bien la distribución de Lamé se recomienda nr≥4.

  nz      divisiones en dirección axial (z, de 0 a L).
          Controla la resolución a lo largo del eje del tubo.
          Los extremos empotrados generan gradientes axiales altos;
          se recomienda nz≥16 para tubos largos o extremos restringidos.

Total de HE8 = ntheta × nr × nz
Total de nodos = ntheta × (nr+1) × (nz+1)

Grupos de elementos:
  VolAll    — todos los hexaedros
  Face_inf  — cara z=0  (condición de contorno DX=DY=DZ=0)
  Face_sup  — cara z=L  (condición de contorno DX=DY=DZ=0)
  Face_int  — superficie interior r=ri (presión interna PRES=50)
"""

import numpy as np
import medcoupling as mc
import os

# ─── Parámetros ──────────────────────────────────────────────────────────────
ri     = 0.8
ro     = 1.0
L      = 5.0
ntheta = 24
nr     = 2
nz     = 16

OUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Tubo.med")

# ─── Índice de nodo ──────────────────────────────────────────────────────────
def node_idx(ir_, ith_, iz_):
    return ir_ * (ntheta * (nz + 1)) + (ith_ % ntheta) * (nz + 1) + iz_

n_nodes = (nr + 1) * ntheta * (nz + 1)

# ─── Coordenadas ─────────────────────────────────────────────────────────────
r_vals     = np.linspace(ri, ro, nr + 1)
theta_vals = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
z_vals     = np.linspace(0.0, L, nz + 1)

coords = np.empty((n_nodes, 3), dtype=np.float64)
for ir_ in range(nr + 1):
    for ith_ in range(ntheta):
        for iz_ in range(nz + 1):
            n = node_idx(ir_, ith_, iz_)
            coords[n, 0] = r_vals[ir_] * np.cos(theta_vals[ith_])
            coords[n, 1] = r_vals[ir_] * np.sin(theta_vals[ith_])
            coords[n, 2] = z_vals[iz_]

# ─── Malla 3D: hexaedros ─────────────────────────────────────────────────────
n_hexa = nr * ntheta * nz
conn_arr = np.empty(n_hexa * 9, dtype=np.int64)   # tipo + 8 nodos por hexa
k = 0
for ir_ in range(nr):
    for ith_ in range(ntheta):
        ith1 = (ith_ + 1) % ntheta
        for iz_ in range(nz):
            conn_arr[k]   = mc.NORM_HEXA8
            conn_arr[k+1] = node_idx(ir_,   ith_,  iz_  )
            conn_arr[k+2] = node_idx(ir_+1, ith_,  iz_  )
            conn_arr[k+3] = node_idx(ir_+1, ith1,  iz_  )
            conn_arr[k+4] = node_idx(ir_,   ith1,  iz_  )
            conn_arr[k+5] = node_idx(ir_,   ith_,  iz_+1)
            conn_arr[k+6] = node_idx(ir_+1, ith_,  iz_+1)
            conn_arr[k+7] = node_idx(ir_+1, ith1,  iz_+1)
            conn_arr[k+8] = node_idx(ir_,   ith1,  iz_+1)
            k += 9

conn_idx = np.arange(0, n_hexa * 9 + 1, 9, dtype=np.int64)

m3d = mc.MEDCouplingUMesh("Tubo", 3)
m3d.setCoords(mc.DataArrayDouble(coords))
m3d.setConnectivity(mc.DataArrayInt(conn_arr), mc.DataArrayInt(conn_idx))
m3d.checkConsistencyLight()

# ─── Malla 2D: quads de cara inferior, superior e interior ───────────────────
# Face_inf (z=0): ir×ith pares, normal −z
face_inf_ids = []
for ir_ in range(nr):
    for ith_ in range(ntheta):
        ith1 = (ith_ + 1) % ntheta
        face_inf_ids.append([
            node_idx(ir_,   ith_,  0),
            node_idx(ir_,   ith1,  0),
            node_idx(ir_+1, ith1,  0),
            node_idx(ir_+1, ith_,  0),
        ])
n_inf = len(face_inf_ids)

# Face_sup (z=L): normal +z
face_sup_ids = []
for ir_ in range(nr):
    for ith_ in range(ntheta):
        ith1 = (ith_ + 1) % ntheta
        face_sup_ids.append([
            node_idx(ir_,   ith_,  nz),
            node_idx(ir_+1, ith_,  nz),
            node_idx(ir_+1, ith1,  nz),
            node_idx(ir_,   ith1,  nz),
        ])
n_sup = len(face_sup_ids)

# Face_int (r=ri): normal −r
face_int_ids = []
for ith_ in range(ntheta):
    ith1 = (ith_ + 1) % ntheta
    for iz_ in range(nz):
        face_int_ids.append([
            node_idx(0, ith_,  iz_  ),
            node_idx(0, ith_,  iz_+1),
            node_idx(0, ith1,  iz_+1),
            node_idx(0, ith1,  iz_  ),
        ])
n_int = len(face_int_ids)

all_faces = face_inf_ids + face_sup_ids + face_int_ids
n_faces = len(all_faces)

conn2d_arr = np.empty(n_faces * 5, dtype=np.int64)
for i, quad in enumerate(all_faces):
    conn2d_arr[i*5]   = mc.NORM_QUAD4
    conn2d_arr[i*5+1] = quad[0]
    conn2d_arr[i*5+2] = quad[1]
    conn2d_arr[i*5+3] = quad[2]
    conn2d_arr[i*5+4] = quad[3]

conn2d_idx = np.arange(0, n_faces * 5 + 1, 5, dtype=np.int64)

m2d = mc.MEDCouplingUMesh("Tubo", 2)
m2d.setCoords(mc.DataArrayDouble(coords))
m2d.setConnectivity(mc.DataArrayInt(conn2d_arr), mc.DataArrayInt(conn2d_idx))
m2d.checkConsistencyLight()

# ─── MEDFileUMesh ─────────────────────────────────────────────────────────────
mf = mc.MEDFileUMesh()
mf.setName("Tubo")
mf.setMeshAtLevel(0, m3d)
# La malla 2D debe compartir las MISMAS coordenadas que la 3D
m2d.tryToShareSameCoords(m3d, 1e-12)
mf.setMeshAtLevel(-1, m2d)

# ─── Familia y grupos de volumen ─────────────────────────────────────────────
# Todos los hexaedros pertenecen al grupo VolAll
fam3d = mc.DataArrayInt(np.full(n_hexa, -1, dtype=np.int64))
fam3d.setName("FamilyFieldOnCells")
mf.setFamilyFieldArr(0, fam3d)
mf.addFamily("VolAll", -1)
mf.addFamilyOnGrp("VolAll", "VolAll")

# ─── Familias de caras ───────────────────────────────────────────────────────
fam2d = np.empty(n_faces, dtype=np.int64)
fam2d[:n_inf]          = -2   # Face_inf
fam2d[n_inf:n_inf+n_sup] = -3  # Face_sup
fam2d[n_inf+n_sup:]    = -4   # Face_int
fam_arr = mc.DataArrayInt(fam2d)
fam_arr.setName("FamilyFieldOnFaces")
mf.setFamilyFieldArr(-1, fam_arr)

for fid, fname in [(-2, "Face_inf"), (-3, "Face_sup"), (-4, "Face_int")]:
    mf.addFamily(fname, fid)
    mf.addFamilyOnGrp(fname, fname)

# ─── Escribir ─────────────────────────────────────────────────────────────────
mf.write(OUT_FILE, 2)  # mode=2 → overwrite
print(f"Malla escrita en: {OUT_FILE}")
print(f"  Nodos  : {n_nodes}")
print(f"  HE8    : {n_hexa}")
print(f"  QU4 inf: {n_inf}")
print(f"  QU4 sup: {n_sup}")
print(f"  QU4 int: {n_int}")

# ─── Verificación rápida ─────────────────────────────────────────────────────
mesh_r = mc.MEDFileMesh.New(OUT_FILE)
print("\nGrupos nivel 0 (vol):", mesh_r.getGroupsOnSpecifiedLev(0))
print("Grupos nivel -1 (face):", mesh_r.getGroupsOnSpecifiedLev(-1))
