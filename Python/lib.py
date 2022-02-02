# Photoelectron_refraction
# Calculation library

import math
import numpy as np
import Config
import random

# gauss function
def gauss(x, s):
    return 1.0/(math.sqrt(2*math.pi)*s)*math.exp(-x*x/(s*s*2))

# broadening profile
def profileCube(dkx, dky, de, sigmak, sigmae, sigmaMax):
    kxLast=math.ceil(sigmak*sigmaMax/dkx)
    kyLast=math.ceil(sigmak*sigmaMax/dky)
    eLast=math.ceil(sigmae*sigmaMax/de)

    cube=np.zeros((kxLast*2+1, kyLast*2+1, eLast*2+1))
    for i in range(0, kxLast+1):
        for j in range(0, kyLast+1):
            for k in range(0, eLast+1):
                weight=gauss(i*dkx, sigmak)*gauss(j*dky, sigmak)*gauss(k*de, sigmae)*dkx*dky*de
                cube[kxLast+i][kyLast+j][eLast+k]=weight
                cube[kxLast+i][kyLast+j][eLast-k]=weight
                cube[kxLast+i][kyLast-j][eLast+k]=weight
                cube[kxLast+i][kyLast-j][eLast-k]=weight
                cube[kxLast-i][kyLast+j][eLast+k]=weight
                cube[kxLast-i][kyLast+j][eLast-k]=weight
                cube[kxLast-i][kyLast-j][eLast+k]=weight
                cube[kxLast-i][kyLast-j][eLast-k]=weight

    # print(cube.sum()*dkx*dky*de) # should be close to 1
    return cube, kxLast, kyLast, eLast

# Original dispersion
def calc1(W, V0, k0, a, V1, kFlat, kFlat_kz, kCurved_k, kxMin, kxMax, kxCount, dkx, kyMin, kyMax, kyCount, dky, eMin, eMax, eCount, de, sigmak, sigmae, dispCube1):
    profile, kxCenter, kyCenter, eCenter=profileCube(dkx, dky, de, sigmak, sigmae, Config.sigmaMax)

    k=np.zeros((3))
    for i in range(kxCount):
        k[0]=kxMin+dkx*i
        for j in range(kyCount):
            k[1]=kyMin+dky*j
            if kFlat==True:
                k[2]=kFlat_kz
            if kFlat==False:
                k[2]=math.sqrt(kCurved_k*kCurved_k-k[0]*k[0]-k[1]*k[1])
            
            kdiff=k-k0

            esk=np.inner(kdiff, kdiff)*a/2.0-W+V1
            esk_index=round((esk-eMin)/de)

            for i2 in range(-kxCenter, kxCenter+1):
                i3=i+i2
                for j2 in range(-kyCenter, kyCenter+1):
                    j3=j+j2
                    for k2 in range(-eCenter, eCenter+1):
                        k3=esk_index+k2
                        if 0<=i3 and i3<kxCount and 0<=j3 and j3<kyCount and 0<=k3 and k3<eCount:
                            dispCube1[i3][j3][k3]+=profile[i2+kxCenter][j2+kyCenter][k2+eCenter]

# Refracted dispersion
def calc2(W, V0, k0, a, V1, kFlat, kFlat_kz, kCurved_k, surfaceConst, surfaceConst_theta, surfaceConst_phi, surfaceRandom_samples, kxMin, kxMax, kxCount, dkx, kyMin, kyMax, kyCount, dky, eMin, eMax, eCount, de, sigmak, sigmae, dispCube2):
    profile, kxCenter, kyCenter, eCenter=profileCube(dkx, dky, de, sigmak, sigmae, Config.sigmaMax)

    k=np.zeros((3))
    n=np.zeros((3))
    if surfaceConst==True:
        n[0]=math.sin(surfaceConst_theta)*math.cos(surfaceConst_phi)
        n[1]=math.sin(surfaceConst_theta)*math.sin(surfaceConst_phi)
        n[2]=math.cos(surfaceConst_theta)

    # generate random surfaces
    if surfaceConst==False:
        nList=np.zeros((surfaceRandom_samples, 3))
        for i in range(surfaceRandom_samples):
            nList[i]=genSurface()
        
        # print(nList)

    for i in range(kxCount):
        print(i)
        k[0]=kxMin+dkx*i
        for j in range(kyCount):
            k[1]=kyMin+dky*j
            if kFlat==True:
                k[2]=kFlat_kz
            if kFlat==False:
                k[2]=math.sqrt(kCurved_k*kCurved_k-k[0]*k[0]-k[1]*k[1])

            kdiff=k-k0
            esk=np.inner(kdiff, kdiff)*a/2.0-W+V1
            epk=np.inner(k, k)/2.0-V0
            eK=round((esk-eMin)/de)

            # Constant surface 
            if surfaceConst==True:
                K=calcK(k, epk, n)
                # Exclude full reflection
                if K is None:
                    continue

                iK=round((K[0]-kxMin)/dkx)
                jK=round((K[1]-kyMin)/dky)

                # Append dispersion
                for i2 in range(-kxCenter, kxCenter+1):
                    i3=iK+i2
                    for j2 in range(-kyCenter, kyCenter+1):
                        j3=jK+j2
                        for k2 in range(-eCenter, eCenter+1):
                            k3=eK+k2
                            if 0<=i3 and i3<kxCount and 0<=j3 and j3<kyCount and 0<=k3 and k3<eCount:
                                dispCube2[i3][j3][k3]+=profile[i2+kxCenter][j2+kyCenter][k2+eCenter]

            # Random surfaces
            else:
                for t in range(surfaceRandom_samples):
                    n=nList[t]
                    K=calcK(k, epk, n)
                    # Exclude full reflection
                    if K is None:
                        # print("!!Full reflection")
                        continue

                    iK=round((K[0]-kxMin)/dkx)
                    jK=round((K[1]-kyMin)/dky)

                    # Append dispersion
                    for i2 in range(-kxCenter, kxCenter+1):
                        i3=iK+i2
                        for j2 in range(-kyCenter, kyCenter+1):
                            j3=jK+j2
                            for k2 in range(-eCenter, eCenter+1):
                                k3=eK+k2
                                if 0<=i3 and i3<kxCount and 0<=j3 and j3<kyCount and 0<=k3 and k3<eCount:
                                    dispCube2[i3][j3][k3]+=profile[i2+kxCenter][j2+kyCenter][k2+eCenter]/surfaceRandom_samples

# Calculate refracted wavevector K
def calcK(k, epk, n):
    KLength=math.sqrt(2*epk)

    kPerpLength=np.inner(k, n)
    kPara=k-kPerpLength*n

    
    KPerpLength2=KLength*KLength-np.inner(kPara, kPara)
    if KPerpLength2<0:
        return None

    KPerpLength=math.sqrt(KPerpLength2)
    if kPerpLength<0:
        # KPerpLength*=-1
        return None

    K=kPara+KPerpLength*n

    return K

# Generate random surface
def genSurface():
    nt=np.zeros((3))
    while True:
        nt[0]=random.uniform(-1, 1)
        nt[1]=random.uniform(-1, 1)
        nt[2]=random.uniform(0, 1)

        nt_length=math.sqrt(np.inner(nt, nt))
        if 0.1 < nt_length and nt_length < 1:
            n=nt/nt_length
            return n
        