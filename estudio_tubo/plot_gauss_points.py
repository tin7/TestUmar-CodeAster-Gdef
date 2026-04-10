"""
plot_gauss_points.py  —  Evolución de coordenadas de Gauss - Tubo con presión interna
Parsea umat_tubo.log y grafica las trayectorias 3D de puntos de Gauss representativos.

Grupos mostrados:
  Anillo interior z≈2.63 (centro): 16 elementoss, ir=0, iz=4
  Anillo interior z≈0.13 (extremo clamped): 8 elementos, ir=0, iz=0
"""

import re
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

LOG_FILE = 'umat_tubo.log'
# LOG_FILE = 'umat_tubo.log'

# Geometría nominal
ri = 0.8
ro = 1.0
L  = 5.0

# Grupos de elementos por anillo
# ir=0 (inner, r≈0.81), nz=8 → paso entre ith es 8 elementos
# iz=4 → z≈2.63 (centro), iz=0 → z≈0.13 (extremo clamped)
# Fórmula: elem_base_ir0 = 193 + ith*8 + iz
GRUPOS = {
    'Anillo interno z≈2.63 (centro)' : {
        'elems': [197 + ith * 8 for ith in range(16)],
        'color': 'tab:red',
        'lw': 2.0, 'alpha': 0.85,
    },
    'Anillo interno z≈0.13 (extremo)': {
        'elems': [193 + ith * 8 for ith in range(16)],
        'color': 'tab:blue',
        'lw': 1.5, 'alpha': 0.60,
    },
}

# Mantener ELEMS para el segundo gráfico cilíndrico (4 puntos representativos)
ELEMS = {
    193: ('Interno z≈0.13 (extremo)', 'tab:red'),
    197: ('Interno z≈2.63 (centro)',  'tab:orange'),
    321: ('Externo z≈0.13 (extremo)', 'tab:blue'),
    325: ('Externo z≈2.63 (centro)',  'tab:cyan'),
}

# ── Parsear log ─────────────────────────────────────────────────────────────
raw = {}  # {(el, gp): {t: (x, y, z)}}
with open(LOG_FILE) as f:
    for line in f:
        m = re.match(
            r'\s+Elemento=\s+(\d+)\s+GP=\s+(\d+)\s+t=\s+([\d.]+)'
            r'\s+Coord_actual=\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)',
            line
        )
        if m:
            el = int(m.group(1))
            gp = int(m.group(2))
            t  = float(m.group(3))
            x, y, z = float(m.group(4)), float(m.group(5)), float(m.group(6))
            raw.setdefault((el, gp), {})[t] = (x, y, z)


def get_traj(el, gp=1):
    """Devuelve (times_array, coords_array (N,3)) deduplicado y ordenado por t."""
    d = raw.get((el, gp), {})
    ts = sorted(d.keys())
    return np.array(ts), np.array([d[t] for t in ts])


# ── Estructura cilíndrica de fondo ──────────────────────────────────────────
def draw_tube_wireframe(ax, ri, ro, L, ntheta=32, nz=4):
    theta = np.linspace(0, 2 * math.pi, ntheta + 1)
    # Círculos a distintas alturas
    for z_level in np.linspace(0, L, nz + 1):
        for r in (ri, ro):
            xs = r * np.cos(theta)
            ys = r * np.sin(theta)
            zs = np.full_like(xs, z_level)
            ax.plot(xs, ys, zs, color='gray', lw=0.5, alpha=0.35, linestyle='--')
    # Generatrices verticales
    for th in np.linspace(0, 2 * math.pi, 8, endpoint=False):
        for r in (ri, ro):
            ax.plot([r * math.cos(th)] * 2,
                    [r * math.sin(th)] * 2,
                    [0, L],
                    color='gray', lw=0.5, alpha=0.35, linestyle='--')


# ── Figura 1: vista 3D con tubo de fondo ────────────────────────────────────
fig = plt.figure(figsize=(12, 8))
ax  = fig.add_subplot(111, projection='3d')
draw_tube_wireframe(ax, ri, ro, L)

cmap = plt.cm.plasma
times_ref = None

for grupo, cfg in GRUPOS.items():
    color = cfg['color']
    first = True
    for el in cfg['elems']:
        ts, coords = get_traj(el)
        if len(ts) == 0:
            continue
        if times_ref is None:
            times_ref = ts

        x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
        ax.plot(x, y, z, '-', color=color, linewidth=cfg['lw'],
                alpha=cfg['alpha'], zorder=6)
        # Marcar punto final con círculo vacío
        ax.scatter([x[-1]], [y[-1]], [z[-1]], s=40, facecolors='none',
                   edgecolors=color, linewidths=1.2, zorder=7)
        if first:
            # Etiqueta solo una vez por grupo
            ax.scatter([], [], [], s=40, facecolors='none',
                       edgecolors=color, linewidths=1.2, label=grupo)
            first = False

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim(-ro * 1.05, ro * 1.05)
ax.set_ylim(-ro * 1.05, ro * 1.05)
ax.set_zlim(0, L)
ax.set_title('Evolución de puntos de Gauss — Tubo con presión interna\n'
             f'ri={ri}, ro={ro}, L={L}, PRES=50, 5 pasos (Δt=0.2)')
ax.legend(loc='upper left', fontsize=9)
ax.view_init(elev=20, azim=-60)

plt.tight_layout()
fig.savefig('gauss_points_3d.png', dpi=150)
print("Guardado: gauss_points_3d.png")


# ── Figura 2: coordenadas cilíndricas r(t) y z(t) ───────────────────────────
fig2, (ax_r, ax_z) = plt.subplots(1, 2, figsize=(12, 5))

for el, (label, color) in ELEMS.items():
    ts, coords = get_traj(el)
    r = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
    z = coords[:, 2]

    ax_r.plot(ts, r, 'o-', color=color, linewidth=1.8, markersize=5, label=label)
    ax_z.plot(ts, z, 'o-', color=color, linewidth=1.8, markersize=5, label=label)

ax_r.axhline(ri, color='gray', lw=0.8, ls=':', label=f'ri={ri}')
ax_r.axhline(ro, color='gray', lw=0.8, ls='--', label=f'ro={ro}')
ax_r.set_xlabel('Tiempo')
ax_r.set_ylabel('r (radio)')
ax_r.set_title('Evolución del radio r(t)')
ax_r.legend(fontsize=8)
ax_r.grid(True, alpha=0.3)

ax_z.set_xlabel('Tiempo')
ax_z.set_ylabel('z (axial)')
ax_z.set_title('Evolución de la coordenada axial z(t)')
ax_z.legend(fontsize=8)
ax_z.grid(True, alpha=0.3)

plt.tight_layout()
fig2.savefig('gauss_points_cylindrical.png', dpi=150)
print("Guardado: gauss_points_cylindrical.png")

plt.show()
