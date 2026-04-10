#!/usr/bin/env python3
# coding=utf-8
"""
gen_tubo.py — Genera malla hexaédrica estructurada (HE8) de un tubo cilíndrico
y la escribe en formato MED (HDF5) listo para Code_Aster.

Geometría
---------
  ri = 0.5   radio interior
  ro = 1.0   radio exterior
  L  = 2.0   longitud axial

Discretización
--------------
  ntheta = 16   divisiones circunferenciales
  nr     = 2    capas radiales
  nz     = 8    divisiones axiales

Grupos de elementos (GROUP_MA en Code_Aster)
--------------------------------------------
  VolAll    — todos los hexaedros (volumen)
  Face_inf  — cara inferior (z=0),  condición de contorno DX=DY=DZ=0
  Face_sup  — cara superior (z=L),  condición de contorno DX=DY=DZ=0
  Face_int  — superficie interior (r=ri), donde se aplica presión interna

Almacenamiento MED
------------------
  Coordenadas NOE/COO  : component-major  [X1..Xn, Y1..Yn, Z1..Zn]
  Conectividad MAI/NOD : column-major  [n1_e1..n1_eN, n2_e1..n2_eN, ...]
  Familias              : enteros negativos en NOD/FAM, descritas en FAS/
"""

import numpy as np
import h5py
import os

# ─── Parámetros geométricos ────────────────────────────────────────────────
ri     = 0.5   # radio interior
ro     = 1.0   # radio exterior
L      = 2.0   # longitud

ntheta = 16    # sectores circunferenciales
nr     = 2     # capas radiales
nz     = 8     # rebanadas axiales

MESH_NAME = "Tubo"
OUT_FILE  = os.path.join(os.path.dirname(__file__), "Tubo.med")

# ─── Índice de nodo ────────────────────────────────────────────────────────
# ir ∈ [0, nr], ith ∈ [0, ntheta-1], iz ∈ [0, nz]
def node_idx(ir, ith, iz):
    return ir * (ntheta * (nz + 1)) + (ith % ntheta) * (nz + 1) + iz

n_nodes = (nr + 1) * ntheta * (nz + 1)

# ─── Coordenadas (component-major para MED) ───────────────────────────────
r_vals     = np.linspace(ri, ro, nr + 1)
theta_vals = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
z_vals     = np.linspace(0.0, L, nz + 1)

X = np.empty(n_nodes); Y = np.empty(n_nodes); Z = np.empty(n_nodes)
for ir_ in range(nr + 1):
    for ith_ in range(ntheta):
        for iz_ in range(nz + 1):
            n = node_idx(ir_, ith_, iz_)
            X[n] = r_vals[ir_] * np.cos(theta_vals[ith_])
            Y[n] = r_vals[ir_] * np.sin(theta_vals[ith_])
            Z[n] = z_vals[iz_]

COO = np.concatenate([X, Y, Z])   # shape (3*n_nodes,)

# ─── Conectividad HE8 ─────────────────────────────────────────────────────
# Ordenamiento estándar (Jacobiano positivo):
#   cara inferior (z=iz): n0(ir,ith), n1(ir+1,ith), n2(ir+1,ith+1), n3(ir,ith+1)
#   cara superior (z=iz+1): n4..n7  (mismo orden)
# Almacenado en column-major para MED: array (n_hexa, 8) → .T.flatten()

hexa_list = []
for ir_ in range(nr):
    for ith_ in range(ntheta):
        ith1 = (ith_ + 1) % ntheta
        for iz_ in range(nz):
            hexa_list.append([
                node_idx(ir_,   ith_,  iz_  ) + 1,   # 0  (1-based)
                node_idx(ir_+1, ith_,  iz_  ) + 1,   # 1
                node_idx(ir_+1, ith1,  iz_  ) + 1,   # 2
                node_idx(ir_,   ith1,  iz_  ) + 1,   # 3
                node_idx(ir_,   ith_,  iz_+1) + 1,   # 4
                node_idx(ir_+1, ith_,  iz_+1) + 1,   # 5
                node_idx(ir_+1, ith1,  iz_+1) + 1,   # 6
                node_idx(ir_,   ith1,  iz_+1) + 1,   # 7
            ])
hexa_arr = np.array(hexa_list, dtype=np.int64)   # (n_hexa, 8)
NOD_HE8  = hexa_arr.T.flatten()                  # column-major → (8*n_hexa,)
n_hexa   = len(hexa_list)
FAM_HE8  = np.full(n_hexa, -1, dtype=np.int64)
NUM_HE8  = np.arange(1, n_hexa + 1, dtype=np.int64)

# ─── Conectividad QU4 ─────────────────────────────────────────────────────
# Face_inf (z=0): normal en -z → order: (ir,ith,0)→(ir,ith+1,0)→(ir+1,ith+1,0)→(ir+1,ith,0)
# Face_sup (z=L): normal en +z → order: (ir,ith,nz)→(ir+1,ith,nz)→(ir+1,ith+1,nz)→(ir,ith+1,nz)
# Face_int (r=ri): normal en -r → order: (0,ith,iz)→(0,ith,iz+1)→(0,ith+1,iz+1)→(0,ith+1,iz)

face_inf_list = []
for ir_ in range(nr):
    for ith_ in range(ntheta):
        ith1 = (ith_ + 1) % ntheta
        face_inf_list.append([
            node_idx(ir_,   ith_,  0) + 1,
            node_idx(ir_,   ith1,  0) + 1,
            node_idx(ir_+1, ith1,  0) + 1,
            node_idx(ir_+1, ith_,  0) + 1,
        ])

face_sup_list = []
for ir_ in range(nr):
    for ith_ in range(ntheta):
        ith1 = (ith_ + 1) % ntheta
        face_sup_list.append([
            node_idx(ir_,   ith_,  nz) + 1,
            node_idx(ir_+1, ith_,  nz) + 1,
            node_idx(ir_+1, ith1,  nz) + 1,
            node_idx(ir_,   ith1,  nz) + 1,
        ])

face_int_list = []
for ith_ in range(ntheta):
    ith1 = (ith_ + 1) % ntheta
    for iz_ in range(nz):
        face_int_list.append([
            node_idx(0, ith_,  iz_  ) + 1,
            node_idx(0, ith_,  iz_+1) + 1,
            node_idx(0, ith1,  iz_+1) + 1,
            node_idx(0, ith1,  iz_  ) + 1,
        ])

n_face_inf = len(face_inf_list)   # nr*ntheta
n_face_sup = len(face_sup_list)   # nr*ntheta
n_face_int = len(face_int_list)   # ntheta*nz

quad_arr = np.array(face_inf_list + face_sup_list + face_int_list, dtype=np.int64)
NOD_QU4  = quad_arr.T.flatten()   # column-major
FAM_QU4  = np.concatenate([
    np.full(n_face_inf, -2, dtype=np.int64),
    np.full(n_face_sup, -3, dtype=np.int64),
    np.full(n_face_int, -4, dtype=np.int64),
])
n_quad    = len(quad_arr)
NUM_QU4   = np.arange(n_hexa + 1, n_hexa + n_quad + 1, dtype=np.int64)

# ─── Familias de nodos (todos en familia 0, sin grupo) ───────────────────
FAM_NOE = np.zeros(n_nodes, dtype=np.int64)
NUM_NOE = np.arange(1, n_nodes + 1, dtype=np.int64)

# ─── Auxiliar: codificar nombre de grupo en array int8 de 80 bytes ────────
def encode_group_name(name: str) -> np.ndarray:
    b = name.encode('ascii')[:79]
    row = np.zeros(80, dtype=np.int8)
    row[:len(b)] = np.frombuffer(b, dtype=np.int8)
    return row

# ─── Grupos de elementos: {fam_id: (fam_label, [group_names])} ───────────
elem_families = {
    -1: ("VolAll",   ["VolAll"]),
    -2: ("Face_inf", ["Face_inf"]),
    -3: ("Face_sup", ["Face_sup"]),
    -4: ("Face_int", ["Face_int"]),
}

# ─── Escribir fichero MED ─────────────────────────────────────────────────
step_key = "-0000000000000000001-0000000000000000001"

with h5py.File(OUT_FILE, "w", libver=('earliest', 'v110')) as f:

    # ── INFOS_GENERALES ────────────────────────────────────────────────
    ig = f.create_group("INFOS_GENERALES")
    ig.attrs["MAJ"] = np.int32(4)
    ig.attrs["MIN"] = np.int32(1)
    ig.attrs["REL"] = np.int32(1)

    # ── ENS_MAA / Tubo / step ─────────────────────────────────────────
    mesh_grp = f.create_group(f"ENS_MAA/{MESH_NAME}")
    mesh_grp.attrs["DES"] = np.bytes_(b"")
    mesh_grp.attrs["DIM"] = np.int32(3)
    mesh_grp.attrs["ESP"] = np.int32(3)
    mesh_grp.attrs["NOM"] = np.bytes_(b"")
    mesh_grp.attrs["NXI"] = np.int32(-1)
    mesh_grp.attrs["NXT"] = np.int32(-1)
    mesh_grp.attrs["REP"] = np.int32(0)
    mesh_grp.attrs["SRT"] = np.int32(0)
    mesh_grp.attrs["TYP"] = np.int32(0)
    mesh_grp.attrs["UNI"] = np.bytes_(b"")
    mesh_grp.attrs["UNT"] = np.bytes_(b"")

    base = f"ENS_MAA/{MESH_NAME}/{step_key}"

    # Nodos
    f.create_dataset(f"{base}/NOE/COO", data=COO)
    f.create_dataset(f"{base}/NOE/FAM", data=FAM_NOE)
    f.create_dataset(f"{base}/NOE/NUM", data=NUM_NOE)

    # Hexaedros
    f.create_dataset(f"{base}/MAI/HE8/NOD", data=NOD_HE8)
    f.create_dataset(f"{base}/MAI/HE8/FAM", data=FAM_HE8)
    f.create_dataset(f"{base}/MAI/HE8/NUM", data=NUM_HE8)

    # Quads (faces)
    f.create_dataset(f"{base}/MAI/QU4/NOD", data=NOD_QU4)
    f.create_dataset(f"{base}/MAI/QU4/FAM", data=FAM_QU4)
    f.create_dataset(f"{base}/MAI/QU4/NUM", data=NUM_QU4)

    # ── FAS / Tubo ────────────────────────────────────────────────────
    fas_base = f"FAS/{MESH_NAME}"

    # FAMILLE_ZERO (obligatoria)
    f0 = f.create_group(f"{fas_base}/FAMILLE_ZERO")
    f0.attrs["NUM"] = np.int64(0)

    # Familias de elementos
    for fid, (label, groups) in elem_families.items():
        full_label = f"FAM_{fid}_{label}"
        fam_grp = f.create_group(f"{fas_base}/ELEME/{full_label}")
        fam_grp.attrs["NUM"] = np.int64(fid)

        # GRO/NOM debe ser shape (n,) con dtype array ('int8',(80,)),
        # identico al formato que escribe MEDfile/Salome.
        nom_dtype = np.dtype(('int8', (80,)))
        nom_data = np.empty(len(groups), dtype=nom_dtype)
        for k, g in enumerate(groups):
            nom_data[k] = encode_group_name(g)
        f.create_dataset(f"{fas_base}/ELEME/{full_label}/GRO/NOM", data=nom_data)

print(f"Malla escrita en: {OUT_FILE}")
print(f"  Nodos    : {n_nodes}")
print(f"  HE8      : {n_hexa}")
print(f"  QU4 inf  : {n_face_inf}")
print(f"  QU4 sup  : {n_face_sup}")
print(f"  QU4 int  : {n_face_int}")
