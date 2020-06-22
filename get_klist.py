#!/usr/bin/env python
import numpy as np
def get_kpath():
    kstic = np.loadtxt('./KPATH')
    kbase = np.loadtxt('./KBASE')
    dim = kbase.shape[0]
    nkstic = len(kstic)
    nkpath = 20
    
    kdis = []
    klist = []

    k = np.arange(dim)
    for ik in np.arange(nkstic):
        if (ik != nkstic-1):
            kstart = kstic[ik]
            kend = kstic[ik+1]
            kstep = (kend - kstart) / (nkpath - 1)

            # the start k distance 
            if (ik == 0):
                kdis.append(0)
            else:
                kdis.append(kdis[-1])

            for ip in np.arange(nkpath):
                kdirect = kstart + ip*kstep
                kcart = np.matmul(kbase.transpose(), kdirect)
                klist.append(kcart)
                if (ip != 0):
                    kdis.append(kdis[-1] + np.linalg.norm(klist[-1]-klist[-2]))
            # include the end points
            # kcart = np.matmul(kbase.transpose(), kdirect)
            # klist.append(kcart)
            # kdis.append( np.linalg.norm(kdis[-1] + np.linalg.norm(kcart[-1]-klist[-2])) )
    return klist, kdis



if __name__ == '__main__':
    klist, kdis = get_kpath()
    klist = np.array(klist)
    kdis = np.array(kdis)
    np.savetxt('./klist.dat', klist, fmt="%18.9f")
    np.savetxt('./kdis.dat', kdis, "%18.9f")

