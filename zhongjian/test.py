import numpy as np
import sympy as sp

z = sp.symbols('z', real=True)

node = np.array([0, 0.300000012])
l = node[1] - node[0]

l0 = 1 - z / l
l1 = z / l

h0 = 1 - 3 * z**2 / l**2 + 2 * z**3 / l**3
h1 = z - 2 * z**2 / l + z**3 / l**2
h2 = 3 * z**2 / l**2 - 2 * z**3 / l**3
h3 = -z**2 / l + z**3 / l**2

u0 = 0
v0 = 0
w0 = 0
alpha0 = 0
beta0 = 0
gamma0 = 0
u1 = 0
v1 =0.000127427
w1 = 0
#alpha1 = -0.000791378
alpha1 = 0
beta1 = 0
#gamma1 = 0
gamma1 = -0.000791378

u = h0*u0 + h1*beta0 + h2*u1 + h3*beta1
uz = sp.diff(u, z)
v = h0*v0 - h1*alpha0 + h2*v1 - h3*alpha1
vz = sp.diff(v, z)
print("vz:", vz)
w = l0*w0 + l1*w1
wz = sp.diff(w, z)

alpha_ = -sp.diff(h0, z)*v0 + sp.diff(h1, z)*alpha0 - sp.diff(h2, z)*v1 + sp.diff(h3, z)*alpha1
print("alpha_:", alpha_)
alphaz = sp.diff(alpha_, z)
beta_ = sp.diff(h0, z)*u0 + sp.diff(h1, z)*beta0 + sp.diff(h2, z)*u1 + sp.diff(h3, z)*beta1
betaz = sp.diff(beta_, z)
gamma_ = l0*gamma0 + l1*gamma1
gammaz = sp.diff(gamma_, z)
print("gammaz:", gammaz)

x_values = [0.1, 0.1 - 0.01, 0.1 - 0.01/2]
y_values = [0, 0, 0]
z_values = np.linspace(0, 0.3, 100)

results = []
for x_val, y_val in zip(x_values, y_values):
    ezz_current = wz - x_val * betaz + y_val * alphaz
    exz_current = uz - y_val * gammaz - beta_
    eyz_current = vz + x_val * gammaz + alpha_
    print("vz:", vz)
    print("eyz_current:", vz + x_val*gammaz)
    print("tets:", alpha_)
    
    ezz_min = min([ezz_current.subs(z, val) for val in z_values])
    ezz_max = max([ezz_current.subs(z, val) for val in z_values])
    
    exz_min = min([exz_current.subs(z, val) for val in z_values])
    exz_max = max([exz_current.subs(z, val) for val in z_values])
    
    eyz_min = min([eyz_current.subs(z, val) for val in z_values])
    eyz_max = max([eyz_current.subs(z, val) for val in z_values])
    
    results.append({
        "x": x_val,
        "y": y_val,
        "ezz": {"min": ezz_min, "max": ezz_max},
        "exz": {"min": exz_min, "max": exz_max},
        "eyz": {"min": eyz_min, "max": eyz_max},
    })
for result in results:
    print(f"x: {result['x']}, y: {result['y']}")
    print(f"  ezz: min: {result['ezz']['min']}, max: {result['ezz']['max']}")
    print(f"  exz: min: {result['exz']['min']}, max: {result['exz']['max']}")
    print(f"  eyz: min: {result['eyz']['min']}, max: {result['eyz']['max']}\n")

h00 = -6/l**2 + 12*z/l**3
h11 = -4/l + 6*z/l**2
h22 = 6/l**2 - 12*z/l**3
h33 = -2/l + 6*z/l**2
x = 0
y = 0.1
gamma_zx = -x * (h00*u0 + h11*beta0 + h22*u1 + h33*beta1)
print("gamma_zx:", gamma_zx)
gamma_zy = -y * (h00*v0 - h11*alpha0 + h22*v1 - h33*alpha1)
print("gamma_zy:", gamma_zy)

gamma_zy_03 = gamma_zy.subs(z, 0.3)
print("gamma_zy_03:", gamma_zy_03)
