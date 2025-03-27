import pandas as pd
from astroquery.vizier import Vizier
import astropy.units as u
import numpy as np

class Database():
    V_EDR3 = Vizier(columns=["**"], catalog="J/A+A/657/A7") # Kervella+22 Vizier Catalog

    def __init__(self, star, **kwargs):
        """ Potential kwargs:
            - gaia [default:'eDR3']: gaia epoch required from the Vizier catalogs
            - geometry [default: {}]: geometry info about the system. """
        self.star = star
        self.geometry = self.build_geometryinfo(kwargs.get("geometry", {}))
        self.gaia = kwargs.get("gaia", 'eDR3')

        self.data = self.get_data(self.gaia)

    def build_geometryinfo(self, discinfo):
        """ Build dictionary with geometry info 
            inc and PA are given in deg.
            Returns: Dictionary with geometry information about the system. Values depend on the geometry disctionary given to the function """
        return {"inc":discinfo.get('inc', None), # default None if no inc is given
                "inc_err":discinfo.get('inc_err', 1.0), # default 1.0 if no inc_err is given [TODO: check default]
                "PA":discinfo.get('PA', None), # default None if no PA is given
                "PA_err":discinfo.get('PA_err', 1.0), # default 1.0 if no PA_err is given [TODO: check default]
                }

    def vizier_query(self, gaia):
        """ Query Vizier catalog for Kervella PMa information. 
            Get info from both gaia and hipparcos. 
            Returns: info for hip, info for gaia, info on star """
        
        # Initialise the dictionaries
        dhip = {'star':[self.star]}
        dgaia = {'star':[self.star]}
        dparams = {'star':[self.star]}

        if gaia == 'eDR3':
            # Query catalog
            results = self.V_EDR3.query_object(self.star, radius=1*u.arcmin)
            query_okay = True
            if len(results) == 0:
                query_okay = False
            elif len(results[0]) == 0:
                query_okay = False
            
            if query_okay == False:
                results = self.V_EDR3.query_object(self.star, radius=1*u.deg)
                
            # Get Hip data
            hip_name = int(self.star.replace('HIP ',''))
            query_hip_names = np.array(results[0]['HIP'])
            
            w, = np.where(query_hip_names == hip_name)
            if len(w) == 0:
                raise ValueError(f"The star {self.star} could not be found in Kervella's catalog. Please do a manual query to understand why.")
            w = w[0]
            
            dhip['S_N'] = results[0]['snrPMaH2EG3a'].data[w] # SNR hip - HipeDR3
            dhip['PMa_ra'] = results[0]['PMaRAH2EG3a'].data[w] # PMa Right Ascension hip - HipeDR3
            dhip['PMa_ra_err'] = results[0]['e_PMaRAH2EG3a'].data[w] # PMa Right Ascension err hip - HipeDR3
            dhip['PMa_dec'] = results[0]['PMaDEH2EG3a'].data[w] # PMa Dec hip - HipeDR3
            dhip['PMa_dec_err'] = results[0]['e_PMaDEH2EG3a'].data[w] # PMa Dec err hip - HipeDR3
            
            # Get eDR3 data
            dgaia['S_N'] = results[0]['snrPMaH2EG3b'].data[w] # SNR eDR3 - HipeDR3
            dgaia['PMa_ra'] = results[0]['PMaRAH2EG3b'].data[w] # PMa Right Ascension eDR3 - HipeDR3
            dgaia['PMa_ra_err'] = results[0]['e_PMaRAH2EG3b'].data[w] # PMa Right Ascension err eDR3 - HipeDR3
            dgaia['PMa_dec'] = results[0]['PMaDEH2EG3b'].data[w] # PMa Dec eDR3 - HipeDR3
            dgaia['PMa_dec_err'] = results[0]['e_PMaDEH2EG3b'].data[w] # PMa Dec err eDR3 - HipeDR3

            # Get star info
            dparams['parallax'] = results[0]['PlxG3'].data[w] # mas
            dparams['dpc'] = 1./(dparams['parallax']*1.0e-3)
            dparams['mstar'] = results[0]['M1'].data[w]

            # Add geometry info to star info
            dparams['PA'] = self.geometry['PA']
            dparams['PA_err'] = self.geometry['PA_err']
            dparams['inc'] = self.geometry['inc']
            dparams['inc_err'] = self.geometry['inc_err']

        elif gaia == 'DR2':
            return NotImplementedError('DR2 epoch still needs to be implemented from Kervella+2019')
        else:
            return NotImplementedError('Wrong choice of gaia epoch.')

        return dhip, dgaia, dparams
    
    def get_data(self, gaia):
        d = {}
        
        # Query data from Kervella vizier catalog
        dhip, dgaia, dparams = self.vizier_query(gaia)

        d['Hipparcos'] = pd.DataFrame(data=dhip).set_index('star')
        d[gaia] = pd.DataFrame(data=dgaia).set_index('star')
        d['params'] = pd.DataFrame(data=dparams).set_index('star')

        return pd.concat(d, axis=1)

    def create_single_star_table(self, local_catalog):

        dhip = {'star':[self.star]}
        dgaia = {'star':[self.star]}
        dparams = {'star':[self.star]}
        
        hip_names = np.array([f'HIP {i}' for i in local_catalog['HIP']])
        w, = np.where(hip_names == self.star)
        if len(w) == 0:
            raise ValueError(f'Star {self.star} was not found!')
        w = w[0]

        dhip['S_N'] = float(local_catalog['snrPMaH2EG3a'].data[w]) # SNR hip - HipeDR3
        dhip['PMa_ra'] = float(local_catalog['PMaRAH2EG3a'].data[w]) # PMa Right Ascension hip - HipeDR3
        dhip['PMa_ra_err'] = float(local_catalog['e_PMaRAH2EG3a'].data[w]) # PMa Right Ascension err hip - HipeDR3
        dhip['PMa_dec'] = float(local_catalog['PMaDEH2EG3a'].data[w]) # PMa Dec hip - HipeDR3
        dhip['PMa_dec_err'] = float(local_catalog['e_PMaDEH2EG3a'].data[w]) # PMa Dec err hip - HipeDR3

        # Get eDR3 data
        dgaia['S_N'] = float(local_catalog['snrPMaH2EG3b'].data[w]) # SNR eDR3 - HipeDR3
        dgaia['PMa_ra'] = float(local_catalog['PMaRAH2EG3b'].data[w]) # PMa Right Ascension eDR3 - HipeDR3
        dgaia['PMa_ra_err'] = float(local_catalog['e_PMaRAH2EG3b'].data[w]) # PMa Right Ascension err eDR3 - HipeDR3
        dgaia['PMa_dec'] = float(local_catalog['PMaDEH2EG3b'].data[w]) # PMa Dec eDR3 - HipeDR3
        dgaia['PMa_dec_err'] = float(local_catalog['e_PMaDEH2EG3b'].data[w]) # PMa Dec err eDR3 - HipeDR3

        # Get star info
        dparams['parallax'] = float(local_catalog['PlxG3'].data[w]) # mas
        dparams['dpc'] = 1./(dparams['parallax']*1.0e-3)
        dparams['mstar'] = float(local_catalog['M1'].data[w])

        # Add geometry info to star info
        dparams['PA'] = self.geometry['PA']
        dparams['PA_err'] = self.geometry['PA_err']
        dparams['inc'] = self.geometry['inc']
        dparams['inc_err'] = self.geometry['inc_err']

        return dhip, dgaia, dparams
    
