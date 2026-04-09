import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection

# --- Datos (únicos, sin duplicados) ---
# GP 1: coords actuales en cada instante
gp1 = np.array([
    [0.211325, 0.211325, 0.788675],  # t=0.0
    [0.246673, 0.212497, 0.797244],  # t=0.2
    [0.282020, 0.213668, 0.805813],  # t=0.4
    [0.317368, 0.214840, 0.814382],  # t=0.6
    [0.352716, 0.216012, 0.822952],  # t=0.8
])

# GP 5: coords actuales en cada instante
gp5 = np.array([
    [0.211325, 0.211325, 0.211325],  # t=0.0
    [0.220796, 0.211639, 0.213621],  # t=0.2
    [0.230268, 0.211953, 0.215917],  # t=0.4
    [0.239739, 0.212267, 0.218213],  # t=0.6
    [0.249210, 0.212581, 0.220509],  # t=0.8
])

times = [0.0, 0.2, 0.4, 0.6, 0.8]

# --- Dibujar cubo unitario ---
def draw_cube(ax):
    vertices = np.array([
        [0,0,0],[1,0,0],[1,1,0],[0,1,0],
        [0,0,1],[1,0,1],[1,1,1],[0,1,1]
    ])
    edges = [
        [0,1],[1,2],[2,3],[3,0],  # base
        [4,5],[5,6],[6,7],[7,4],  # tapa
        [0,4],[1,5],[2,6],[3,7],  # laterales
    ]
    lines = [[vertices[e[0]], vertices[e[1]]] for e in edges]
    lc = Line3DCollection(lines, colors='gray', linewidths=0.8, linestyles='--', alpha=0.5)
    ax.add_collection3d(lc)

# --- Figura ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
draw_cube(ax)

# Colormap para el tiempo
cmap = plt.cm.viridis
colors = [cmap(t / 0.8) for t in times]

# Trayectoria GP 1
ax.plot(gp1[:,0], gp1[:,1], gp1[:,2], '-', color='tab:red', linewidth=1.5, alpha=0.6)
ax.scatter(gp1[:,0], gp1[:,1], gp1[:,2], s=30, facecolors='none',
           edgecolors='red', linewidths=1.0, marker='o', zorder=5, label='GP 1')

# Trayectoria GP 5
ax.plot(gp5[:,0], gp5[:,1], gp5[:,2], '-', color='tab:blue', linewidth=1.5, alpha=0.6)
ax.scatter(gp5[:,0], gp5[:,1], gp5[:,2], s=30, facecolors='none',
           edgecolors='blue', linewidths=1.0, marker='o', zorder=5, label='GP 5')

# Anotaciones de tiempo
for i, t in enumerate(times):
    ax.text(gp1[i,0]+0.01, gp1[i,1]+0.01, gp1[i,2]+0.01, f't={t}', fontsize=7, color='red')
    ax.text(gp5[i,0]+0.01, gp5[i,1]+0.01, gp5[i,2]+0.01, f't={t}', fontsize=7, color='blue')

# Ejes y leyenda
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_zlim(0, 1)
ax.set_title('Evolución de puntos de Gauss — Elemento 20')
ax.legend(loc='upper left', fontsize=8)

# Colorbar para el tiempo
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 0.8))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.1)
cbar.set_label('Tiempo')

plt.tight_layout()
plt.savefig('/media/javier/BACKUP1/TestUMAT/Elastica/estudio_cubo_coords/gauss_points_evolution.png', dpi=150)
plt.show()
print("Guardado en gauss_points_evolution.png")
