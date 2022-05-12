import matplotlib.pyplot as plt
import numpy as np

#Experimental Setup
n1 = 1
n2 = 1.5

d1 = 0.1
d2 = 0.005
dL = 0.001
d4 = 0.1
l = d1+d2+dL+d4

offset = 0 #offset to the right - ignore this




# Vector Calculations

euclideanBasis3D = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

def iProd(vec1, vec2):
    l = 0
    for cc in range(len(vec1)):
        l += vec1[cc] * vec2[cc]
    return l


def sMult(vec1, l):
    return [vec1[cc] * l for cc in range(len(vec1))]


def vAdd(vec1, vec2):
    return [vec1[cc] + vec2[cc] for cc in range(len(vec1))]


def vAbs(v):
    return np.sqrt(iProd(v, v))


def changeBasis(v1, B1, B2):  # ORTHONORMAL!
    v2 = np.zeros(len(B2))
    for cc in range(len(v2)):
        for dd in range(len(v1)):
            v2[cc] += v1[dd] * iProd(B1[dd], B2[cc])
    return v2

def matMult2x2_2x1(mat, vec):
    return [mat[0, 0] * vec[0] + mat[0, 1] * vec[1], mat[1, 0] * vec[0] + mat[1, 1] * vec[1]]


# Propagation Functions

def linearPropagation(r, v, d):
    return vAdd(r, sMult(v, d))


def refract(v, alpha, n1, n2):
    v = v/vAbs(v)
    n = (np.sin(alpha), -np.cos(alpha),0)
    v2 = (0,0,0)
    testScalar = v[0]/n[0]
    if(testScalar*n[1]==v[1] and testScalar*n[2] == v[2]):
        v2 = v
    else:
        e1, e2 = n, sMult(vAdd(v, sMult(n,-iProd(n,v))),-1)
        e2 = e2/vAbs(e2)
        omega1, omega2 = (iProd(v,e1))/(vAbs(v)*vAbs(e1)), 0
        if(True): #Enter domain-breaking processes here
            omega2 = np.arcsin(n1/n2*np.sin(omega1))
        vprime = (-np.cos(omega2),-np.sin(omega2))
        v2 = changeBasis(vprime, [e1,e2], euclideanBasis3D)
    return v2

def crossWithCircle(r, v, dL, offset):
    deltax,deltay =0,0
    y0 = (r[1] - offset) % dL
    if(v[1]==0):
        deltay = 0
        deltax = np.sqrt(dL**2-y0**2)
    else:
        k = v[0]/v[1]
        x0 = -k*y0
        y1,y2 = (-2*k*x0+np.sqrt(4*k**2*x0**2-4*(k**2+1)*(x0**2-dL**2)))/(2*(k**2+1)), (-2*k*x0-np.sqrt(4*k**2*x0**2-4*(k**2+1)*(x0**2-dL**2)))/(2*(k**2+1))
        x1,x2 = np.sqrt(dL**2-y1**2), np.sqrt(dL**2-y2**2)
        if(y1 >= 0):
            deltax, deltay = x1, y1
        else:
            deltax, deltay = x2, y2
    return [deltax,deltay]



def raytrace(r,v,n1,n2):
    r2 = linearPropagation(r, v, d1/iProd(v,[1,0,0]))
    v2 = refract(v,0,n1,n2)
    r3 = linearPropagation(r2, v2, d2/iProd(v2,[1,0,0]))
    radial = crossWithCircle(r3, v2, dL, offset)
    r4 = linearPropagation(r3, v2, radial[0])
    v3 = refract(v2,np.arctan(radial[1]/radial[0]),n2,n1)
    r5 = linearPropagation(r4,v3,l-r4[0])
    return r5

def evaluate(r):
    if(r[1]<0.01 and r[1]>-0.01):
        return 1
    else:
        return 0


v = [np.cos(np.pi/6), np.sin(np.pi/6), 0]
r = [0,0,0]
#raytrace(r,v,n1,n2)


def plot1():
    phi_0, phi_1 = -np.pi/4, np.pi/4
    theta_0, theta_1 = np.pi/4, 3*np.pi/4
    nphi, ntheta = 1000,1000
    rawimage = np.zeros((nphi,ntheta))
    evaluatedimage = rawimage
    for cc in range(ntheta):
        print(int(cc/nphi*100), '%')
        for dd in range(nphi):
            phi = cc / nphi * (phi_1 - phi_0) + phi_0
            theta = cc / ntheta * (theta_1 - theta_0) + theta_0
            v_0 = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
            rawimage[dd,cc] = raytrace([0,0,0],v_0, n1, n2)[2]
            evaluatedimage[dd,cc] = evaluate(raytrace([0,0,0],v_0, n1, n2))
    return evaluatedimage



plt.imshow(plot1())
plt.show()
