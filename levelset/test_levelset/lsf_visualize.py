import numpy as np
import matplotlib.pyplot as plt

# Set the figure size
fig, ax = plt.subplots(figsize=(8, 6))

# Generate grid points
x = np.linspace(-2, 2, 400)
y = np.linspace(-1.5, 1.5, 300)
X, Y = np.meshgrid(x, y)

# Define the level set function phi
# For simplification, we'll define the interface as the curve y = x^2
phi = Y - X**2

# Use contour to plot the interface phi = 0
ax.contour(X, Y, phi, levels=[0], colors='black', linewidths=2)

# Use contourf to fill the areas where phi > 0 and phi < 0
c = ax.contourf(X, Y, phi, levels=[-np.inf, 0, np.inf], colors=['blue', 'red'], alpha=0.6)

# Add labels for regions and interface
ax.text(0, 1, r'$\phi > 0$', horizontalalignment='center', verticalalignment='center', fontsize=12, color='white')
ax.text(0, -1, r'$\phi < 0$', horizontalalignment='center', verticalalignment='center', fontsize=12, color='black')
ax.text(0, 0, r'interface $\phi = 0$', horizontalalignment='center', verticalalignment='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.7, edgecolor='white'))

# Enhancements for better visualization
ax.set_title('Visualization of Level Set Function')
ax.axis('off')  # Turn off the axis
plt.tight_layout()  # Adjust the layout for better appearance

plt.show()
