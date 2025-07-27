import sympy as sym
from sympy.vector import CoordSys3D, Del, curl
C = CoordSys3D('C')

def compute(f, alpha, beta):
    x = sym.symbols("x")
    y = sym.symbols("y")

    # 构造 f
    fx = f.dot(C.i).subs({C.x:x, C.y:y})
    fy = f.dot(C.j).subs({C.x:x, C.y:y})
    print('fx:', fx)
    print('fy:', fy)

    # 构造 curl(f)
    cf = curl(f)
    cfz = cf.dot(C.k).subs({C.x:x, C.y:y})
    print('cfz:', cfz)

    # 构造 curl(curl(f))
    ccf = curl(cf)
    ccfx = ccf.dot(C.i).subs({C.x:x, C.y:y})
    ccfy = ccf.dot(C.j).subs({C.x:x, C.y:y})
    print('ccfx:', ccfx)
    print('ccfy:', ccfy)

    sour = (ccfx*alpha - fx*beta)*C.i + (ccfy*alpha - fy*beta)*C.j
    print('sour:', sour)

a = sym.symbols("a")
b = sym.symbols("b")
c = sym.symbols("c")

u0 = 2*(sym.cos(a*C.x)*sym.cos(a*C.y)*C.i + sym.sin(a*C.x)*sym.sin(a*C.y)*C.j)
u1 = (sym.cos(b*C.x)*sym.cos(b*C.y)*C.i + sym.sin(b*C.x)*sym.sin(b*C.y)*C.j)

u0 = 2*(sym.cos(a*C.x)*sym.cos(c*C.y)*C.i + sym.sin(a*C.x)*sym.sin(c*C.y)*C.j)
u1 = (sym.cos(b*C.x)*sym.cos(c*C.y)*C.i + sym.sin(b*C.x)*sym.sin(c*C.y)*C.j)

alpha = 1
beta0 = 2*a**2
beta1 = 2*b**2



compute(u0, alpha, beta0)
compute(u1, alpha, beta1)






