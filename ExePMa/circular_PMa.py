"""
Companion mass constraints from the Hipparcos-Gaia proper motion anomaly (PMa).

Implements the circular-orbit method of Kervella et al. (2019, A&A 623, A72) using the
PMa values of the Kervella et al. (2022, A&A 657, A7) catalogue. The main entry point is
the PMa class, which combines catalogue access (Database)
and plotting (Plotting).
"""
import numpy as np
import os

from ExePMa.Database import Database
from ExePMa.plotting import Plotting

class PMa(Database, Plotting):
    """
    Companion mass versus separation from the PMa of a single star.

    On initialisation the star's PMa and stellar parameters are retrieved (from Vizier or a
    local catalogue copy) and stored in self.data. Masses are then obtained with
    mass_retrieval() and plotted with plot_pma().

    Parameters
    ----------
    **kwargs
        Passed to Database; see that class for the full list
        (star, hip, gaia, geometry, mstar, local_catalog).

    Attributes
    ----------
    data : pandas.DataFrame
        One row (indexed by HIP identifier) with column groups 'Hipparcos',
        'eDR3' and 'params'. See get_data().

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from ExePMa import circular_PMa
    >>> star = circular_PMa.PMa(star='HD 107146')
    >>> fig, ax = plt.subplots()
    >>> ax, aps = star.plot_pma(ax, epoch='eDR3', pma_params={'savelog': False})
    >>> ax = star.add_axis_info(ax)
    """
    ### units
    M_SUN=1.989e+30 # kg
    AU=1.496e+11   # m
    YEAR=365.24 # days
    G=6.67408e-11 # mks
    M_JUP=1.898e27 # kg

    # epochs
    EPOCHS_DT = {"DR2":668., "eDR3":1038., "Hipparcos":1227.} # days
    EPOCHS_DT_HG = {"DR2":24.25, "eDR3":24.75} # years

    def __init__(self, **kwargs):
        """Retrieve the PMa data of the star (prints self.data). See class docstring."""

        ### Get all information from Vizier catalog
        Database.__init__(self, **kwargs) # creating data DataFrame
        print(self.data)
        print('-----------------')

        ### Initialise Plotting
        Plotting.__init__(self)

    def gamma(self, Pp):
        """
        Correction for smearing of the PMa over the observing window.

        Eq. 13 (Sect. 3.6.2) of Kervella et al. (2019).

        Parameters
        ----------
        Pp : float or ndarray
            Orbital period in units of the observing window length
            (P / EPOCHS_DT[epoch], dimensionless).

        Returns
        -------
        float or ndarray
            Smearing factor gamma (dimensionless, -> 1 for long periods).
        """
        return Pp/(np.sqrt(2.)*np.pi) * np.sqrt(1. - np.cos(2.*np.pi/Pp))
    
    def period(self, r, m):
        """
        Keplerian orbital period, assuming the companion mass is negligible.

        Parameters
        ----------
        r : float or ndarray
            Orbital radius [au].
        m : float
            Mass of the primary [Msun].

        Returns
        -------
        float or ndarray
            Orbital period [yr].
        """
        return np.sqrt(r**3 / m)

    def eta(self, PMa_orb, inc, PA):
        """
        Deprojection factor of the PMa vector onto the orbital plane.

        Sect. 3.6.1 of Kervella et al. (2019). The unit PMa vector is rotated into the
        frame of the system's major axis (PA) and the minor-axis component is
        deprojected by cos(inc).

        Parameters
        ----------
        PMa_orb : float or ndarray
            Position angle of the PMa vector [deg, east of north].
        inc : float or ndarray
            Inclination of the orbital plane [deg].
        PA : float or ndarray
            Position angle of the orbital plane's major axis [deg, east of north].

        Returns
        -------
        float or ndarray
            eta = |projected| / |deprojected| PMa (dimensionless, <= 1).
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
        Efficiency of the PMa as a proxy for the instantaneous orbital velocity.

        Eq. 14 (Sect. 3.6.3) of Kervella et al. (2019) for a circular orbit: the PMa is
        measured relative to the long-term Hipparcos-Gaia proper motion, which partly
        absorbs the orbital motion when the period is comparable to the baseline. The
        orbital phase at the Hipparcos epoch is drawn uniformly Nrandom times.

        Parameters
        ----------
        period : float
            Orbital period [yr].
        Nrandom : int
            Number of random orbital phases.
        epoch : {'Hipparcos', 'DR2', 'eDR3'}
            Epoch at which the PMa is measured.
        gaia : {'DR2', 'eDR3'}
            Gaia release defining the Hipparcos-Gaia baseline (EPOCHS_DT_HG).

        Returns
        -------
        ndarray, shape (Nrandom,)
            zeta (dimensionless).
        """

        # Randomly pick position of the star in orbit Nrandom times (Mean Anomaly)
        M1 = np.random.uniform(0.0, 2.*np.pi, Nrandom) # mean anomaly at Hipparcos
        M2 = M1 + 2.0*np.pi * self.EPOCHS_DT_HG[gaia]/period # mean anomaly at Gaia epoch (after epochs_dt_HG)

        # Assume circular orbit and get tangential velocity vector at different epochs
        if epoch == 'DR2' or epoch == 'eDR3':
            vx = -np.sin(M2)
            vy = np.cos(M2)
        elif epoch == 'Hipparcos':
            vx = -np.sin(M1)
            vy = np.cos(M1)
        else:
            return NotImplementedError('Wrong epoch chosen')
        
        # calculate difference in position between position in Hipparcos and Gaia
        deltar_x = np.cos(M2) - np.cos(M1)
        deltar_y = np.sin(M2) - np.sin(M1)

        C = period/(2.*np.pi)
        zeta =  np.sqrt( (vx - C * deltar_x/self.EPOCHS_DT_HG[gaia])**2. +  (vy - C*deltar_y/self.EPOCHS_DT_HG[gaia])**2. )

        return zeta
    
    def zeta_v2(self, period, Nrandom, cosi, epoch, gaia): # this incorporates uncertainty through cosi. 
        """
        Variant of zeta() that folds the inclination into zeta (currently unused).

        Intended for when the inclination and PA are unknown; see the TODO in
        calculate_zeta().

        Parameters
        ----------
        period : float
            Orbital period [yr].
        Nrandom : int
            Number of random orbital phases.
        cosi : float or ndarray, shape (Nrandom,)
            cos(inclination), e.g. np.random.uniform(0, 1, Nrandom) for isotropic orbits.
        epoch : {'Hipparcos', 'DR2', 'eDR3'}
            Epoch at which the PMa is measured.
        gaia : {'DR2', 'eDR3'}
            Gaia release defining the Hipparcos-Gaia baseline.

        Returns
        -------
        ndarray, shape (Nrandom,)
            zeta (dimensionless).
        """
        
        M1 = np.random.uniform(0.0, 2.*np.pi, Nrandom) # mean anomaly at Hipparcos
        M2 = M1 + 2.0*np.pi * self.EPOCHS_DT_HG[gaia]/period # mean anomaly at Gaia epoch (after epochs_dt_HG)

        # circular orbit
        if epoch == 'DR2' or epoch == 'eDR3':
            vx = -np.sin(M2)
            vy = np.cos(M2)
        elif epoch == 'Hipparcos':
            vx = -np.sin(M1)
            vy = np.cos(M1)
        else:
            return NotImplementedError('Wrong epoch chosen')

        deltar_x = np.cos(M2) - np.cos(M1)
        deltar_y = (np.sin(M2) - np.sin(M1)) * cosi
        vy = vy * cosi
        C = period / (2.*np.pi)
        
        zeta =  np.sqrt( (vx - C*deltar_x/self.EPOCHS_DT_HG[gaia])**2. +  (vy - C*deltar_y/self.EPOCHS_DT_HG[gaia])**2. )/np.sqrt(vx**2+vy**2)
        return zeta

    def calculate_zeta(self, P, epoch, known_inc, inc, Nrandom):
        """
        Evaluate zeta() for each orbital period.

        Parameters
        ----------
        P : ndarray, shape (Na,)
            Orbital periods [yr].
        epoch : {'Hipparcos', 'DR2', 'eDR3'}
            Epoch at which the PMa is measured.
        known_inc : bool
            Whether the inclination is known. Currently unused (see TODO).
        inc : ndarray
            Inclinations [deg]. Currently unused (see TODO).
        Nrandom : int
            Number of random orbital phases per period.

        Returns
        -------
        ndarray, shape (Na, Nrandom)
            zeta for each period and random phase.
        """
        zeta = np.zeros((len(P), Nrandom))
        # [TODO: check what to do about zeta_v2]
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
        Companion mass needed to produce a given tangential velocity anomaly.

        m2 = sqrt(mstar / G) * sqrt(r) * (v / eta) / (gamma * zeta), following
        Kervella et al. (2019). mstar is taken from self.data.

        Parameters
        ----------
        rs : ndarray, shape (Na,)
            Orbital radii [au].
        v : ndarray, shape (Nrandom,)
            Norm of the tangential PMa as a linear velocity [m/s] ("dVt" in the catalogue).
        PA_PMa : ndarray, shape (Nrandom,)
            Position angle of the PMa vector [deg].
        epoch : {'Hipparcos', 'eDR3'}
            Epoch at which the PMa is measured.
        inc : ndarray, shape (Nrandom,)
            Inclination samples [deg].
        PA : ndarray, shape (Nrandom,)
            Major-axis position angle samples [deg].
        know_inc : bool
            Whether the inclination is known (passed to calculate_zeta()).
        Nrandom : int
            Number of Monte Carlo samples.

        Returns
        -------
        ndarray, shape (Na, Nrandom)
            Companion masses [Mjup].
        """

        eta = self.eta(PA_PMa, inc, PA)
        periods = self.period(rs, m=self.data['params']['mstar'].iloc[0]) # years        
        zeta = self.calculate_zeta(P=periods, epoch=epoch, known_inc=know_inc, inc=inc, Nrandom=Nrandom)

        const = np.sqrt(self.data['params']['mstar'].iloc[0] * self.M_SUN / self.G)
        r_term = const * np.sqrt(rs * self.AU) / self.gamma(periods/(self.EPOCHS_DT[epoch] / self.YEAR))
        corr_vel = v / eta

        mass = np.matmul(np.array([r_term]).T, np.array([corr_vel])) / zeta # Mstar
        return mass / self.M_JUP # masses in Mjup
    
    def get_mass_MC(self, data, rs, epoch, Nrandom):
        """
        Monte Carlo distribution of companion masses at each orbital radius.

        Draws Nrandom PMa vectors from the catalogue values and uncertainties, plus
        inclinations and PAs (Gaussian if given in geometry, otherwise isotropic
        inclinations and uniform PA in [0, 180) deg), and computes m2() for each.

        Parameters
        ----------
        data : pandas.DataFrame
            Star data, normally self.data.
        rs : ndarray, shape (Na,)
            Orbital radii [au].
        epoch : {'Hipparcos', 'eDR3'}
            Epoch at which the PMa is measured.
        Nrandom : int
            Number of Monte Carlo samples.

        Returns
        -------
        ndarray, shape (Na, 7)
            Mass percentiles [Mjup] at each radius, columns
            0.135, 2.28, 15.9, 50 (median), 84.2, 97.7, 99.865
            (i.e. -3, -2, -1, 0, +1, +2, +3 sigma).
        """
        PMa_ras = np.random.normal(data[epoch]['PMa_ra'].iloc[0], data[epoch]['PMa_ra_err'].iloc[0], Nrandom)
        PMa_decs=np.random.normal(data[epoch]['PMa_dec'].iloc[0], data[epoch]['PMa_dec_err'].iloc[0], Nrandom)

        if data['params']['inc'].iloc[0] is None:
            print('unknown system inc')
            know_inc = False
            incs=np.arccos(np.random.uniform(0.0, 1.0, Nrandom))*180/np.pi 
        else:
            incs=np.random.normal(data['params']['inc'].iloc[0], data['params']['inc_err'].iloc[0], Nrandom)
            know_inc = True

        if data['params']['PA'].iloc[0] is None:
            print('unknown system PA')
            PAs=np.random.uniform(0.0, 180.0, Nrandom)
        else:
            PAs=np.random.normal(data['params']['PA'].iloc[0], data['params']['PA_err'].iloc[0], Nrandom)

        msMC = np.zeros((len(rs), Nrandom))
        ms = np.zeros((len(rs), 7)) # ±1,2,3 sigma + mean

        PAs_PMa = np.zeros(Nrandom)

        ### MC
        v_ = np.sqrt(PMa_ras**2. + PMa_decs**2.) / data['params']['parallax'].iloc[0] *4740.470 # m/s        
        PAs_PMa[:] = np.arctan2(PMa_ras, PMa_decs)*180./np.pi # deg

        msMC = self.m2(rs=rs, v=v_, PA_PMa=PAs_PMa, epoch=epoch, inc=incs, PA=PAs, know_inc=know_inc, Nrandom=Nrandom)

        PAs_PMa[PAs_PMa<0.0]=PAs_PMa[PAs_PMa<0.0]+360.
        print('PA mean and std = %1.2f +- %1.2f'%(np.mean(PAs_PMa), np.std(PAs_PMa)))

        ### Retrieve percentiles
        for ir in range(len(rs)):
            ms[ir,:]=np.percentile(msMC[ir,:], [0.135, 2.28, 15.9, 50., 84.2, 97.7, 99.865 ]) # ±1,2,3 sigma

        return ms
    
    def mass_retrieval(self, epoch, **kwargs):
        """
        Companion mass (percentiles) needed to explain the PMa, versus separation.

        Method of Kervella et al. (2019), assuming a circular orbit. Optionally saves
        the result to <odir><subdir>/PMa_limits_<epoch>_<star>[_noPA][_noInc].txt.

        Parameters
        ----------
        epoch : {'Hipparcos', 'eDR3'}
            Epoch at which the PMa is measured.
        pma_params : dict, optional
            Keyword argument with any of:

            - Na (int, default 300): number of separations to sample.
            - aps (ndarray, default np.logspace(-0.5, 2.5, Na)): separations [au].
            - Nrandom (int, default 10000): Monte Carlo samples per separation.
            - savelog (bool, default True): save the result to a text file.
            - subdir (str, default 'log'): directory name for the saved file.
            - odir (str, default './'): parent path of subdir (needs trailing /).

        Returns
        -------
        ndarray, shape (Na, 8)
            Column 0: separations [au]; columns 1-7: mass percentiles [Mjup] as returned
            by get_mass_MC().
        """
        # read in any kwargs params
        pma_params = kwargs.get('pma_params', {})
        # default values if not in pma_params
        Na = pma_params.get('Na', 300)
        aps = pma_params.get('aps', np.logspace(-0.5, 2.5, Na))
        Nrandom = pma_params.get('Nrandom', 10000)
        savelog = pma_params.get('savelog', True)
        subdir = pma_params.get('subdir', 'log')
        odir = pma_params.get('odir', './')

        if savelog:
            PMadir = odir + subdir + '/'
            if not os.path.exists(PMadir):
                os.makedirs(PMadir)

            ### create filepath for saving PMa info
            if self.data['params']['PA'].iloc[0] is None:
                if self.data['params']['inc'].iloc[0] is None:
                    file_path = f'PMa_limits_{epoch}_{self.star}_noPA_noInc.txt'
                else:
                    file_path = f'PMa_limits_{epoch}_{self.star}_noPA.txt'
            else:
                if self.data['params']['inc'].iloc[0] is None:
                    file_path = f'PMa_limits_{epoch}_{self.star}_noInc.txt'
                else:
                    file_path = f'PMa_limits_{epoch}_{self.star}.txt'

        ms = self.get_mass_MC(data=self.data, rs=aps, epoch=epoch, Nrandom=Nrandom)  
        saved_ms = np.vstack([aps, ms.T]).T

        if savelog:
            print(f'Saving ms for {epoch}')
            np.savetxt(PMadir+file_path, saved_ms, 
            header='Percentiles from Monte Carlo simulation\n semi-major axis\t  0.135\t 2.28\t 15.9\t 50\t 84.2\t 97.7\t 99.865')

        return saved_ms
