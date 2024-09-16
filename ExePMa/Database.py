import pandas as pd
from astroquery.vizier import Vizier
import astropy.units as u

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
            # Get Hip data
            dhip['S_N'] = results[0]['snrPMaH2EG3a'].data[0] # SNR hip - HipeDR3
            dhip['PMa_ra'] = results[0]['PMaRAH2EG3a'].data[0] # PMa Right Ascension hip - HipeDR3
            dhip['PMa_ra_err'] = results[0]['e_PMaRAH2EG3a'].data[0] # PMa Right Ascension err hip - HipeDR3
            dhip['PMa_dec'] = results[0]['PMaDEH2EG3a'].data[0] # PMa Dec hip - HipeDR3
            dhip['PMa_dec_err'] = results[0]['e_PMaDEH2EG3a'].data[0] # PMa Dec err hip - HipeDR3
            
            # Get eDR3 data
            dgaia['S_N'] = results[0]['snrPMaH2EG3b'].data[0] # SNR eDR3 - HipeDR3
            dgaia['PMa_ra'] = results[0]['PMaRAH2EG3b'].data[0] # PMa Right Ascension eDR3 - HipeDR3
            dgaia['PMa_ra_err'] = results[0]['e_PMaRAH2EG3b'].data[0] # PMa Right Ascension err eDR3 - HipeDR3
            dgaia['PMa_dec'] = results[0]['PMaDEH2EG3b'].data[0] # PMa Dec eDR3 - HipeDR3
            dgaia['PMa_dec_err'] = results[0]['e_PMaDEH2EG3b'].data[0] # PMa Dec err eDR3 - HipeDR3

            # Get star info
            dparams['parallax'] = results[0]['PlxG3'].data[0] # mas
            dparams['dpc'] = 1./(dparams['parallax']*1.0e-3)
            dparams['mstar'] = results[0]['M1'].data[0]

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
        d= {}
        
        # Query data from Kervella vizier catalog
        dhip, dgaia, dparams = self.vizier_query(gaia)

        d['Hipparcos'] = pd.DataFrame(data=dhip).set_index('star')
        d[gaia] = pd.DataFrame(data=dgaia).set_index('star')
        d['params'] = pd.DataFrame(data=dparams).set_index('star')

        return pd.concat(d, axis=1)
