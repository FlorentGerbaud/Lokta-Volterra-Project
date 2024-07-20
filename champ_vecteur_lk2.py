import matplotlib.pyplot as plt
import numpy as np

# Définir la grille de coordonnées
a = 4/3
b = 2/3
c = 0.2
d = 0.1
x, y = np.meshgrid(np.arange(0, 2, 0.1), np.arange(0, 5, 0.1))

# Définir les composantes x et y du champ de vecteurs
u = a*x-b*x*y
v = c*x*y-d*y

# Tracer le champ de vecteurs
fig, ax = plt.subplots()
ax.quiver(x, y, u, v)

ax.scatter(d/c, a/b, c='red')
ax.scatter(0, 0, c='green')

# Afficher le graphique
plt.show()