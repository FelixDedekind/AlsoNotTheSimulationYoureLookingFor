import matplotlib.pyplot as plt
import numpy as np


#Experimental Setup
global n1, n2
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



# Propagation Functions

def linearPropagation(r, v, d):
    return vAdd(r, sMult(v, d))


def refract(v, alpha, n1, n2):
    print(n1,n2)
    v = sMult(v,1/vAbs(v))
    n = (-np.cos(alpha), np.sin(alpha), 0)
    if(vAbs(vAdd(v,n)) < 0.00001):        #perpendicular to surface
        v2 = v
    else:
        e1 = n
        e2 = sMult(vAdd(v, sMult(e1,-iProd(v,e1))),-1)
        e2 = sMult(e2, 1/vAbs(e2))
        omega1 = np.arccos(iProd(v,e1)/(vAbs(v)*vAbs(e1)))
        if(abs(n1/n2 * np.sin(omega1)) > 1):
            v2 = (0,0,0)
        else:
            print('now')
            omega2 = np.arcsin(np.abs(n1/n2 * np.sin(omega1)))
            v2 = vAdd(sMult(e1,-np.cos(omega2)), sMult(e2, -np.sin(omega2)))
            print(v2)
    return v2

def refract_kaysquared(v, alpha, n1, n2):
    print(f'refract: {n1/n2}')
    v = sMult(v, 1 / vAbs(v))
    n = (-np.cos(alpha), np.sin(alpha), 0)
    t1a = np.cross(sMult(n,-1),v)
    t1b = np.cross(n,t1a)
    t1 = sMult(t1b,n1/n2)

    t2 = sMult(n,np.sqrt(1-(n1/n2)**2*iProd(np.cross(n,v),np.cross(n,v))))
    return vAdd(t1, sMult(t2,-1))

def crossWithCircle(r, v, dL, offset):
    x = (r[1] - offset) % (2*dL) - dL
    if(v[1]==0):

    else:

    return



#raytracing

def raytrace(r0,v0,n1,n2):
    print(f'raytrace: {n1,n2}')
    v0 = sMult(v0, 1/vAbs(v0))

    l1 = d1/iProd(v0, (1,0,0))      #this distance is incorrect for a tilted lense
    r1 = linearPropagation(r0, v0, l1)

    v1 = refract_kaysquared(v0, 0, n1, n2)

    l2 = d1/iProd(v1, (1,0,0))
    r2 = linearPropagation(r1,v1,l2)
    #l2 = d2/iProd(v1, (1,0,0))
    #r2 = linearPropagation(r1, v1, l2)

    #circleCrossing = crossWithCircle(r2, v1, dL, 0)
    #l3 = circleCrossing[1]
    #r3 = linearPropagation(r2, v1, l3)

    #v2 = refract(v1, np.arctan(circleCrossing[1]/circleCrossing[0]), n2, n1)

    #l4 = l - r3[0]
    #r4 = linearPropagation(r3, v2, l4)
    #return r4
    return r2


#evaluation and plotting

def evaluate(r):
    if(abs(r[1])<0.01):
        return 1
    else:
        return 0

def plot1():
    phi_0, phi_1 = -np.pi/4, np.pi/4
    theta_0, theta_1 = np.pi/4, 3*np.pi/4
    nphi, ntheta = 300,300
    rawimage = np.zeros((nphi,ntheta))
    evaluatedimage = rawimage
    for cc in range(ntheta):
        print(int(cc/ntheta*100), '%')
        theta = cc / ntheta * (theta_1 - theta_0) + theta_0
        for dd in range(nphi):
            phi = dd / nphi * (phi_1 - phi_0) + phi_0
            v_0 = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
            rawimage[dd,cc] = raytrace([0,0,0],v_0, n1, n2)[1]
            evaluatedimage[dd,cc] = evaluate(raytrace([0,0,0],v_0, n1, n2))
    return rawimage

#plt.imshow(plot1())
#plt.xlabel('phi')
#plt.ylabel('theta')
#plt.show()


print(raytrace((0,0,0), (1/np.sqrt(2),1/np.sqrt(2),0), n1,n2))

origin = (0,0,0)
initvec = (1/np.sqrt(2),1/np.sqrt(2),0.05)
initvec /= np.sqrt(np.dot(initvec,initvec))
#initvec = (1,0,0)

print(raytrace(origin, initvec, n1,n2))
#print(  np.tan(np.arcsin(np.sin(np.pi/4)/1.5))  *0.1 )

zList = list(np.linspace(0,0.2,10))
outvecX,outvecY,outvecZ=[],[],[]

for z in zList:
    initvec = (1 / np.sqrt(2), 1 / np.sqrt(2), z)
    initvec /= np.sqrt(np.dot(initvec, initvec))
    outvecX.append( raytrace(origin, initvec, n1,n2)[0] )
    outvecY.append(raytrace(origin, initvec, n1, n2)[1])
    outvecZ.append(raytrace(origin, initvec, n1, n2)[2])

print(outvecX)
print(outvecY)
print(outvecZ)

#plt.plot(zList, outvecX)
#plt.plot(zList, outvecZ)
#plt.show()
