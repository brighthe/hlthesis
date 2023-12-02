import numpy as np
import matplotlib.pyplot as plt

# Creating a grid of points
x = np.linspace(-2, 2, 400)
y = np.linspace(-2, 2, 400)
X, Y = np.meshgrid(x, y)

# Defining the function
Z = np.sin(np.pi * X) * np.sin(np.pi * Y)

# Plotting the function
plt.figure(figsize=(6, 6))
plt.contour(X, Y, Z, levels=[0], colors='blue')
plt.title(r'$\sin(\pi x) \cdot \sin(\pi y) = 0$')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.show()
