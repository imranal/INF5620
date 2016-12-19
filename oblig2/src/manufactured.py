from sympy import *

x = Symbol('x')
y = Symbol('y')
A = Symbol('A')
B = Symbol('B')
c = Symbol('c')
kx = Symbol('kx')
ky = Symbol('ky')
t = Symbol('t')
w = Symbol('w')
q = x
u = (A*cos(w*t) + B*sin(w*t))*exp(-c*t)*cos(kx*x)*cos(ky*y)
ut = diff(u, t)  #Derivative u_t
utt = diff(diff(u, t), t)  #Derivative u_tt
rxx = diff(q*diff(u, x), x)  #Derivative (q*u_x)_x
ryy = diff(q*diff(u, y), y)  #Derivative (q*u_y)_y
fxyt = utt + 2*c*ut - rxx - ryy#print simplify(fxyt)
"""fxyt = (-A*c**2*cos(kx*x)*cos(t*w) + A*kx**2*x*cos(kx*x)*cos(t*w) + A*kx*sin(kx*x)*cos(t*w) + A*ky**2*x*cos(kx*x)*cos(t*w) - A*w**2*cos(kx*x)*cos(t*w) - B*c**2*sin(t*w)*cos(kx*x) + B*kx**2*x*sin(t*w)*cos(kx*x) + B*kx*sin(kx*x)*sin(t*w) + B*ky**2*x*sin(t*w)*cos(kx*x) - B*w**2*sin(t*w)*cos(kx*x))*exp(-c*t)*cos(ky*y)
-A*c*cos(kx*x)*cos(ky*y) + B*w*cos(kx*x)*cos(ky*y)
"""
t = 0
# check all these computations
V = -c*(A*cos(t*w) + B*sin(t*w))*exp(-c*t)*cos(kx*x)*cos(ky*y) + (-A*w*sin(t*w) + B*w*cos(t*w))*exp(-c*t)*cos(kx*x)*cos(ky*y)  # then V(x,y) = -A*c*cos(kx*x)*cos(ky*y) + B*w*cos(kx*x)*cos(ky*y)

I =  (A*cos(w*t) + B*sin(w*t))*exp(-c*t)*cos(kx*x)*cos(ky*y)
print I

