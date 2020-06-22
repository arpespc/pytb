#!/usr/bin/env python
import numpy as np
I = 1j
class tbsystem():
    def __init__(self, num_wann):
        self.rbase = np.loadtxt('./RBASE')
        self.kbase = np.loadtxt('./KBASE')
        self.kstic = np.loadtxt('./KPATH')
        self.dim = self.rbase.shape[0]   # demension of space
        self.num_wann = num_wann
        self.tbdata = None # just declare here
        self.klist = None
        self.kdis = None

    def get_tbdata(self):
        """
        get the tight binding data
        """
        hopping_dates = np.loadtxt('./hopping_dates')
        hopping_dates_cart = []
        for hopping_term in hopping_dates:
            site1_cart = np.matmul(self.rbase.transpose(), hopping_term[0:self.dim])
            site2_cart = np.matmul(self.rbase.transpose(), hopping_term[self.dim:2*self.dim])
            hopping_term[0:self.dim] = site1_cart
            hopping_term[self.dim:2*self.dim] = site2_cart
            hopping_dates_cart.append(hopping_term)
        hopping_dates_cart = np.array(hopping_dates_cart)
        self.tbdata = hopping_dates_cart

    def get_klist(self):
        """
        return the klist in cart
        and the disstance list
        """
        kstic = self.kstic
        kbase = self.kbase
        dim = self.dim
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
        self.klist = klist
        self.kdis = kdis
        return klist, kdis
    def ham_k(self, k):
        """
        construct k space tight binding hamitonian
        k in cartisein coordinates
        """

        num_wann = self.num_wann
        hopping_dates_cart = self.tbdata
        dim = self.dim
        Ham = np.zeros((num_wann, num_wann), dtype=np.complex)
        for hopping in hopping_dates_cart:
            site1_cart = hopping[0:dim]
            site2_cart = hopping[dim:2*dim]
            delta_r = site2_cart - site1_cart
            i = int(hopping[-2])
            j = int(hopping[-1])
            phase_factor = np.exp(-I*np.dot(k, delta_r))
            t = hopping[2*dim]
            Ham[j, i] += t*phase_factor
        return Ham + Ham.conj().transpose() # np.conj(Ham.transpose())

    def bands(self):
        """
        calculate the eigenvals on a list of k points
        and output them to filei
        """
        klist = self.klist
        kdis = self.kdis
        band_data = []
        for ik , k in enumerate(klist):
            ham = self.ham_k(k)
            eigenvals = np.real(np.linalg.eigvals(ham))
            line = [kdis[ik]]
            line += sorted(eigenvals.tolist())
            band_data.append(line)
        
        np.savetxt('band.dat', band_data, "%18.6f")

if __name__ == "__main__":
    tb_model = tbsystem(3) ### dimension of hamitinian
    tb_model.get_tbdata()
    tb_model.get_klist()
    tb_model.bands()



