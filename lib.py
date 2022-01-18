# Photoelectron_refraction
# Calculation library

import math
import numpy as np
import Config

def gauss(x, s):
    return 1.0/(math.sqrt(2*math.pi)*s)*math.exp(-x*x/(s*s*2))

def profileCube(dkx, dky, de, sigmak, sigmae, sigmaMax):
    kxLast=math.ceil(sigmak*sigmaMax/dkx)
    kyLast=math.ceil(sigmak*sigmaMax/dky)
    eLast=math.ceil(sigmae*sigmaMax/de)

    cube=np.zeros((kxLast*2+1, kyLast*2+1, eLast*2+1))
    for i in range(0, kxLast+1):
        for j in range(0, kyLast+1):
            for k in range(0, eLast+1):
                weight=gauss(i*dkx, sigmak)*gauss(j*dky, sigmak)*gauss(k*de, sigmae)
                cube[kxLast+i][kyLast+j][eLast+k]=weight
                cube[kxLast+i][kyLast+j][eLast-k]=weight
                cube[kxLast+i][kyLast-j][eLast+k]=weight
                cube[kxLast+i][kyLast-j][eLast-k]=weight
                cube[kxLast-i][kyLast+j][eLast+k]=weight
                cube[kxLast-i][kyLast+j][eLast-k]=weight
                cube[kxLast-i][kyLast-j][eLast+k]=weight
                cube[kxLast-i][kyLast-j][eLast-k]=weight

    # print(cube.sum()*dkx*dky*de)
    return cube, kxLast, kyLast, eLast


def calc1(W, V0, k0, a, V1, kFlat, kFlat_kz, kCurved_k, kxMin, kxMax, kxCount, dkx, kyMin, kyMax, kyCount, dky, eMin, eMax, eCount, de, sigmak, sigmae, dispCube1):
    profile, kxCenter, kyCenter, eCenter=profileCube(dkx, dky, de, sigmak, sigmae, Config.sigmaMax)

    k=np.zeros((3))
    for i in range(kxCount+1):
        k[0]=kxMin+dkx*i
        for j in range(kyCount+1):
            k[1]=kyMin+dky*j
            if kFlat==True:
                k[2]=kFlat_kz
            if kFlat==False:
                k[2]=math.sqrt(kCurved_k*kCurved_k-k[0]*k[0]-k[1]*k[1])
            
            kdiff=k-k0

            ek=np.inner(kdiff, kdiff)*a/2.0-W+V1
            ek_index=round((ek-eMin)/de)

            for i2 in range(-kxCenter, kxCenter+1):
                i3=i+i2
                for j2 in range(-kyCenter, kyCenter+1):
                    j3=j+j2
                    for k2 in range(-eCenter, eCenter+1):
                        k3=ek_index+k2
                        if 0<=i3 and i3<kxCount and 0<=j3 and j3<kyCount and 0<=k3 and k3<eCount:
                            dispCube1[i3][j3][k3]+=profile[i2+kxCenter][j2+kyCenter][k2+eCenter]








            

