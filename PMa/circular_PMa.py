import numpy as np
import json, os, sys
import matplotlib.pyplot as plt

from PMa.Database import Database
from PMa.plotting import Plotting

class PMa(Database, Plotting):
    ### units
    Msun=1.989e+30 # kg
    au=1.496e+11   # m
    year=365.24 # days
    G=6.67408e-11 # mks
    Mjup=1.898e27 # kg

    # epochs
    epochs_dt = {"DR2":668., "eDR3":1038., "Hipparcos":1227.} # days
    epochs_dt_HG = {"DR2":24.25, "eDR3":24.75} # years

    def __init__(self, discdata, gaia, iparam, odir=None):
        self.odir = odir
        self.gaia = gaia
        Database.__init__(self, discdata, gaia=gaia, iparam=iparam)
        Plotting.__init__(self)
        print(self.data)
        print('-------')

    def gamma(self, Pp):
        """
        Define the function that corrects for the window smearing effect (eq. 13 in section 3.6.2 of Kervella+2019)
        """
        return Pp/(np.sqrt(2.)*np.pi) * np.sqrt(1. - np.cos(2.*np.pi/Pp))
    
    def period(self, r, m):
        """
        orbital period as a function of radial separation and mass of the primary (if mass of secondary << mass of primary)
        r in au
        m in Msun
        """
        return np.sqrt(r**3 / m)

    def eta(self, PMa_orb, inc, PA):
        """
        Deprojection (section 3.6.1)
        function of the 3D orbital proper motion anomaly vector, the inclination and the parallactic angle
        """
        # projections of pma normalised to unity
        pma_ra = np.sin(PMa_orb*np.pi/180.)
        pma_dec = np.cos(PMa_orb*np.pi/180.)

        # projections of pma normalised to unity and rotated such that x is in direction of PA
        pma_x = pma_dec * np.cos(PA*np.pi/180.) + pma_ra * np.sin(PA*np.pi/180.)
        pma_y = pma_dec * np.sin(PA*np.pi/180.) - pma_ra * np.cos(PA*np.pi/180.)

        # deprojected pma with x along PA of system and y perpendicular 
        pma_x_dep=pma_x
        pma_y_dep=pma_y/np.cos(inc*np.pi/180.)
        
        eta = 1.0 / np.sqrt( pma_x_dep**2. + pma_y_dep**2. )
        return eta
    
    def zeta(self, period, Nrandom, epoch, gaia):
        """
        When inclination and PA is known
        correcting for the orbital period of the system (see eq. 14 and section 3.6.3 in Kervella+2019)
        """

        # Randomly pick position of the star in orbit Nrandom times (Mean Anomaly)
        M1 = np.random.uniform(0.0, 2.*np.pi, Nrandom) # mean anomaly at Hipparcos
        M2 = M1 + 2.0*np.pi * self.epochs_dt_HG[gaia]/period # mean anomaly at Gaia epoch (after epochs_dt_HG)

        # Assume circular orbit and get tangential velocity vector at different epochs
        if epoch == 'DR2' or epoch == 'eDR3':
            vx = -np.sin(M2)
            vy = np.cos(M2)
        elif epoch == 'Hipparcos':
            vx = -np.sin(M1)
            vy = np.cos(M1)
        else:
            print('wrong epoch')
            return -1.
        
        # calculate difference in position between position in Hipparcos and Gaia
        deltar_x = np.cos(M2) - np.cos(M1)
        deltar_y = np.sin(M2) - np.sin(M1)

        C = period/(2.*np.pi)
        zeta =  np.sqrt( (vx - C * deltar_x/self.epochs_dt_HG[gaia])**2. +  (vy - C*deltar_y/self.epochs_dt_HG[gaia])**2. )

        return zeta
    
    def zeta_v2(self, period, Nrandom, cosi, epoch, gaia): # this incorporates uncertainty through cosi. 
        """
        When inclination and PA is unknown
        cosi is one or a set of random inclinations distributed as uniform if no information or as something else
        e.g. cosi=np.random.uniform(0.0, 1.0, Nrandom)
        """
        
        M1 = np.random.uniform(0.0, 2.*np.pi, Nrandom) # mean anomaly at Hipparcos
        M2 = M1 + 2.0*np.pi * self.epochs_dt_HG[gaia]/period # mean anomaly at Gaia epoch (after epochs_dt_HG)

        # circular orbit
        if epoch == 'DR2' or epoch == 'eDR3':
            vx = -np.sin(M2)
            vy = np.cos(M2)
        elif epoch == 'Hipparcos':
            vx = -np.sin(M1)
            vy = np.cos(M1)
        else:
            print('wrong epoch')
            return -1.

        deltar_x = np.cos(M2) - np.cos(M1)
        deltar_y = (np.sin(M2) - np.sin(M1)) * cosi
        vy = vy * cosi
        C = period / (2.*np.pi)

        # zeta_x=vx - C*deltar_x/delta_t_HG
        # zeta_y=vy - C*deltar_y/delta_t_HG
        
        zeta =  np.sqrt( (vx - C*deltar_x/self.epochs_dt_HG[gaia])**2. +  (vy - C*deltar_y/self.epochs_dt_HG[gaia])**2. )/np.sqrt(vx**2+vy**2)
        return zeta

    def calculate_zeta(self, P, epoch, known_inc, inc, Nrandom):
        zeta = np.zeros((len(P), Nrandom))

        # for i, period in enumerate(P):
        #     if known_inc: # when inc and PA are known
        #         zeta[i,:] = self.zeta(period, Nrandom=Nrandom, epoch=epoch, gaia=self.gaia)
        #     else: # when inc and PA are unknown
        #         zeta[i,:] = self.zeta_v2(period, Nrandom=Nrandom, cosi=np.cos(inc), epoch=epoch, gaia=self.gaia)

        for i, period in enumerate(P):
            zeta[i,:] = self.zeta(period, Nrandom=Nrandom, epoch=epoch, gaia=self.gaia)

        return zeta

    def m2(self, rs, v, PA_PMa, epoch, inc, PA, know_inc, Nrandom):
        """
        links the mass of the secondary (m2), the radial distance (rs), the mass of the primary (mstar) and the v (the norm of the
        tangential PMa vector converted to linear velocity using the Gaia parallax, "dVt" in the Vizier catalogue)
        """

        eta = self.eta(PA_PMa, inc, PA)
        periods = self.period(rs, m=self.data['params']['mstar'][0]) # years        
        zeta = self.calculate_zeta(P=periods, epoch=epoch, known_inc=know_inc, inc=inc, Nrandom=Nrandom)

        const = np.sqrt(self.data['params']['mstar'][0] * self.Msun / self.G)
        r_term = const * np.sqrt(rs * self.au) / self.gamma(periods/(self.epochs_dt[epoch] / self.year))
        corr_vel = v / eta

        mass = np.matmul(np.array([r_term]).T, np.array([corr_vel])) / zeta # Mstar
        return mass / self.Mjup # masses in Mjup
    
    def get_mass_MC(self, data, rs, epoch, Nrandom):
        """
        Calculates the mass of companion causing the PMa
        rs [au], PMa [mas], parallax [mas], inc [deg], PA [deg], mstar [Msun] (default is 1.0), default gaia is 'eDR3'
        returns ±1,2,3 sigma and mean of m2 at the given rs
        """
        PMa_ras = np.random.normal(data[epoch]['PMa_ra'][0], data[epoch]['PMa_ra_err'][0], Nrandom)
        PMa_decs=np.random.normal(data[epoch]['PMa_dec'][0], data[epoch]['PMa_dec_err'][0], Nrandom)

        if data['params']['inc'][0] is None or data['params']['inc_err'][0] is None:
            print('unknown system inc')
            know_inc = False
            incs=np.arccos(np.random.uniform(0.0, 1.0, Nrandom))*180/np.pi 
        else:
            incs=np.random.normal(data['params']['inc'][0], data['params']['inc_err'][0], Nrandom)
            know_inc = True

        if data['params']['PA'][0] is None or data['params']['PA_err'][0] is None:
            print('unknown system PA')
            PAs=np.random.uniform(0.0, 180.0, Nrandom)
        else:
            PAs=np.random.normal(data['params']['PA'][0], data['params']['PA_err'][0], Nrandom)

        msMC = np.zeros((len(rs), Nrandom))
        ms = np.zeros((len(rs), 7)) # ±1,2,3 sigma + mean

        PAs_PMa = np.zeros(Nrandom)
        ###### MC

        v_ = np.sqrt(PMa_ras**2. + PMa_decs**2.) / data['params']['parallax'][0] *4740.470 # m/s        
        PAs_PMa[:] = np.arctan2(PMa_ras, PMa_decs)*180./np.pi # deg

        msMC = self.m2(rs=rs, v=v_, PA_PMa=PAs_PMa, epoch=epoch, inc=incs, PA=PAs, know_inc=know_inc, Nrandom=Nrandom)

        PAs_PMa[PAs_PMa<0.0]=PAs_PMa[PAs_PMa<0.0]+360.
        print('PA mean and std = %1.2f +- %1.2f'%(np.mean(PAs_PMa), np.std(PAs_PMa)))

        ### Retrieve percentiles
        for ir in range(len(rs)):
            ms[ir,:]=np.percentile(msMC[ir,:], [0.135, 2.28, 15.9, 50., 84.2, 97.7, 99.865 ]) # ±1,2,3 sigma

        return ms
    
    def mass_retrieval(self, aps, Nrandom, epoch, subdir=None):

        try:
            dir = self.odir + subdir + '/'
            if self.discdata['PA'] is None:
                if self.discdata['inc'] is None:
                    file_path = 'PMa_limits_{}_{}_noPA_noInc.txt'.format(epoch, self.target)
                else:
                    file_path = 'PMa_limits_{}_{}_noPA.txt'.format(epoch, self.target)
            else:
                if self.discdata['inc'] is None:
                    file_path = 'PMa_limits_{}_{}_noInc.txt'.format(epoch, self.target)
                else:
                    file_path = 'PMa_limits_{}_{}.txt'.format(epoch, self.target)

            if os.path.exists(dir+file_path):
                answer = False
                while answer == False:
                    recalculate = input('Results already ran for this star, do you want to recalculate? [y/n]')
                    if recalculate == 'y':
                        answer = True
                        ms = self.get_mass_MC(data=self.data, rs=aps, epoch=epoch, Nrandom=Nrandom)  
                        saved_ms = np.vstack([aps, ms.T]).T
                        print('Saving ms for {}'.format(epoch))
                        np.savetxt(dir+file_path, saved_ms, 
                    header='Percentiles from Monte Carlo simulation\n semi-major axis\t  0.135\t 2.28\t 15.9\t 50\t 84.2\t 97.7\t 99.865')
                    elif recalculate == 'n':
                        answer = True
                        saved_ms = np.loadtxt(dir+file_path)
                    else:
                        print('Answer again...')
                        answer = False
            else:
                if not os.path.exists(dir):
                    os.makedirs(dir)
                ms = self.get_mass_MC(data=self.data, rs=aps, epoch=epoch, Nrandom=Nrandom)
                saved_ms = np.vstack([aps, ms.T]).T
                print('Saving ms for {}'.format(epoch))
                np.savetxt(dir+file_path, saved_ms, 
                    header='Percentiles from Monte Carlo simulation\n semi-major axis\t  0.135\t 2.28\t 15.9\t 50\t 84.2\t 97.7\t 99.865')

        except: # no saving
            ms = self.get_mass_MC(data=self.data, rs=aps, epoch=epoch, Nrandom=Nrandom)
            saved_ms = np.vstack([aps, ms.T]).T

        return saved_ms
