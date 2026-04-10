"""
plot_local_axes.py  —  Verificacion visual de los ejes locales del material

Lee umat_axes.log (escrito en KINC=1, NPT=1) y grafica en cada punto de Gauss
los tres vectores de la base local:
  e_r  (rojo)  — debe apuntar radialmente hacia afuera
  e_th (azul)  — debe ser tangente al circulo (counterclockwise)
  e_z  (verde) — debe apuntar a lo largo del eje del tubo

Figuras generadas:
  1. axes_top_view.png   — Vista XY desde arriba, mostrando e_r y e_th
  2. axes_3d.png         — Vista 3D con los tres vectores en algunos puntos
  3. checks.txt          — Verificacion numerica: ortonormalidad y alineacion
"""

import re
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

LOG_FILE = 'umat_axes.log'

# Geometria del tubo
ri, ro, L = 0.8, 1.0, 5.0

# ── Parsear log ─────────────────────────────────────────────────────────────
pat = re.compile(
    r'Elem=\s*(\d+)\s+'
    r'P=\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+'
    r'er=\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+'
    r'eth=\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+'
    r'ez=\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)'
)

records = []
seen_elems = set()
with open(LOG_FILE) as f:
    for line in f:
        m = pat.search(line)
        if m:
            elem = int(m.group(1))
            if elem not in seen_elems:
                seen_elems.add(elem)
                vals = list(map(float, m.groups()[1:]))
                records.append({
                    'elem': elem,
                    'P':   np.array(vals[0:3]),
                    'er':  np.array(vals[3:6]),
                    'eth': np.array(vals[6:9]),
                    'ez':  np.array(vals[9:12]),
                })

print(f"Puntos de Gauss leidos: {len(records)}")

# Convertir a arrays
P   = np.array([r['P']   for r in records])
ER  = np.array([r['er']  for r in records])
ETH = np.array([r['eth'] for r in records])
EZ  = np.array([r['ez']  for r in records])

# ── Verificacion numerica ───────────────────────────────────────────────────
checks = []
checks.append("=" * 60)
checks.append("VERIFICACION DE EJES LOCALES (KINC=1, NPT=1)")
checks.append("=" * 60)

# 1. Normas
norm_er  = np.linalg.norm(ER,  axis=1)
norm_eth = np.linalg.norm(ETH, axis=1)
norm_ez  = np.linalg.norm(EZ,  axis=1)
checks.append(f"\n1. NORMAS (deben ser 1.0):")
checks.append(f"   |e_r|  : min={norm_er.min():.6f}  max={norm_er.max():.6f}")
checks.append(f"   |e_th| : min={norm_eth.min():.6f}  max={norm_eth.max():.6f}")
checks.append(f"   |e_z|  : min={norm_ez.min():.6f}  max={norm_ez.max():.6f}")

# 2. Ortogonalidad
dot_r_th = np.einsum('ij,ij->i', ER,  ETH)
dot_r_z  = np.einsum('ij,ij->i', ER,  EZ)
dot_th_z = np.einsum('ij,ij->i', ETH, EZ)
checks.append(f"\n2. ORTOGONALIDAD (productos punto, deben ser ~0):")
checks.append(f"   e_r·e_th: max|·|={np.abs(dot_r_th).max():.2e}")
checks.append(f"   e_r·e_z:  max|·|={np.abs(dot_r_z).max():.2e}")
checks.append(f"   e_th·e_z: max|·|={np.abs(dot_th_z).max():.2e}")

# 3. Sistema dextrorso: e_r x e_th debe ser e_z
cross = np.cross(ER, ETH)
diff_ez = np.linalg.norm(cross - EZ, axis=1)
checks.append(f"\n3. SISTEMA DEXTRORSO: |e_r x e_th - e_z| (debe ser ~0):")
checks.append(f"   max error = {diff_ez.max():.2e}")

# 4. Radial: e_r debe ser paralelo a P proyectado en el plano perp al eje
#    (para eje Z: la proyeccion es simplemente (x,y,0) normalizado)
Pxy = P.copy()
Pxy[:, 2] = 0.0
r_xy  = np.linalg.norm(Pxy, axis=1, keepdims=True)
er_expected = Pxy / np.where(r_xy > 1e-9, r_xy, 1.0)
diff_er = np.linalg.norm(ER - er_expected, axis=1)
checks.append(f"\n4. ALINEACION RADIAL: |e_r - (x,y,0)/r| (eje=Z, debe ser ~0):")
checks.append(f"   max error = {diff_er.max():.2e}")

# 5. Dirección e_z constante
checks.append(f"\n5. COMPONENTE CONSTANTE e_z:")
checks.append(f"   e_z (todas las entradas): mean={EZ.mean(axis=0)}, std={EZ.std(axis=0)}")

checks.append("\n" + "=" * 60)
report = "\n".join(checks)
print(report)
with open('checks.txt', 'w') as f:
    f.write(report + "\n")
print("Guardado: checks.txt")

# ── Figura 1: Vista superior XY ─────────────────────────────────────────────
# Mostrar solo una capa z para no sobrecargar
z_vals = np.unique(np.round(P[:, 2], 2))
# Elegir 3 secciones representativas
z_show = [z_vals[0], z_vals[len(z_vals)//2], z_vals[-1]]

scale = 0.06  # longitud de las flechas

fig, axes_arr = plt.subplots(1, len(z_show), figsize=(15, 5))
fig.suptitle('Vista XY — ejes locales del material\n(rojo=e_r, azul=e_θ, solo puntos GP=1 de KINC=1)', fontsize=12)

for ax, z_target in zip(axes_arr, z_show):
    mask = np.abs(P[:, 2] - z_target) < 0.05
    px, py = P[mask, 0], P[mask, 1]
    er_x, er_y   = ER[mask, 0],  ER[mask, 1]
    eth_x, eth_y = ETH[mask, 0], ETH[mask, 1]

    # Contorno del tubo (referencia)
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(ri*np.cos(theta), ri*np.sin(theta), 'k--', lw=0.8, alpha=0.4, label='ri')
    ax.plot(ro*np.cos(theta), ro*np.sin(theta), 'k-',  lw=0.8, alpha=0.4, label='ro')

    # Flechas e_r y e_th
    # scale=8: flechas de longitud 1/8 = 0.125 unidades de datos
    ax.quiver(px, py, er_x,  er_y,  color='red',  scale=8, scale_units='xy',
              angles='xy', width=0.004, headwidth=4, label='e_r',  alpha=0.9)
    ax.quiver(px, py, eth_x, eth_y, color='blue', scale=8, scale_units='xy',
              angles='xy', width=0.004, headwidth=4, label='e_θ', alpha=0.9)

    ax.scatter(px, py, s=15, color='gray', zorder=5)
    ax.set_aspect('equal')
    ax.set_title(f'z ≈ {z_target:.2f}')
    ax.set_xlabel('X'); ax.set_ylabel('Y')
    ax.grid(True, alpha=0.3)
    if ax == axes_arr[0]:
        ax.legend(fontsize=8, loc='upper right')

plt.tight_layout()
plt.savefig('axes_top_view.png', dpi=150)
print("Guardado: axes_top_view.png")

# ── Figura 2: Vista 3D (muestra un anillo completo) ─────────────────────────
# Elegir el anillo de la sección central (z más cercano al centro del tubo)
z_mid   = z_vals[len(z_vals)//2]
mask_3d = np.abs(P[:, 2] - z_mid) < 0.05
# Añadir también el extremo para comparar
mask_end = np.abs(P[:, 2] - z_vals[0]) < 0.05

fig = plt.figure(figsize=(12, 6))
ax3 = fig.add_subplot(111, projection='3d')

for mask, label_sfx, alpha in [(mask_3d, f'z≈{z_mid:.2f}', 0.9),
                                 (mask_end, f'z≈{z_vals[0]:.2f}', 0.5)]:
    px, py, pz = P[mask, 0], P[mask, 1], P[mask, 2]
    ax3.quiver(px, py, pz,
               ER[mask,0]*scale*8, ER[mask,1]*scale*8, ER[mask,2]*scale*8,
               color='red', alpha=alpha, label=f'e_r  {label_sfx}')
    ax3.quiver(px, py, pz,
               ETH[mask,0]*scale*8, ETH[mask,1]*scale*8, ETH[mask,2]*scale*8,
               color='blue', alpha=alpha, label=f'e_θ  {label_sfx}')
    ax3.quiver(px, py, pz,
               EZ[mask,0]*scale*8, EZ[mask,1]*scale*8, EZ[mask,2]*scale*8,
               color='green', alpha=alpha, label=f'e_z  {label_sfx}')
    ax3.scatter(px, py, pz, s=20, color='gray', zorder=5)

# Wireframe del tubo
theta = np.linspace(0, 2*np.pi, 32)
z_wire = np.array([0, L])
for r_w in [ri, ro]:
    for z_w in z_wire:
        ax3.plot(r_w*np.cos(theta), r_w*np.sin(theta),
                 np.full_like(theta, z_w), 'k-', alpha=0.15, lw=0.8)

ax3.set_xlabel('X'); ax3.set_ylabel('Y'); ax3.set_zlabel('Z')
ax3.set_title('Ejes locales en 3D (rojo=e_r, azul=e_θ, verde=e_z)')
# Mostrar solo la leyenda de la primera sección para no saturar
handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles[:6], labels[:6], fontsize=7, loc='upper left')

plt.tight_layout()
plt.savefig('axes_3d.png', dpi=150)
print("Guardado: axes_3d.png")
