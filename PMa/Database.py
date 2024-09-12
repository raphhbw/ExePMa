import pandas as pd
from astroquery.vizier import Vizier
import astropy.units as u
import sys

class Database():
    v_eDR3 = Vizier(columns=["**"], catalog="J/A+A/657/A7") # Kervella+22 Vizier Catalog

    def __init__(self, discdata, gaia, iparam):
        
        self.discdata = discdata
        self.target = discdata['name']
        self.data = self.get_data(gaia, iparam)


    def vizier_query(self, gaia):
        dhip = {'star':[self.target]}
        dgaia = {'star':[self.target]}
        dparams = {'star':[self.target]}
        if gaia == 'eDR3':
            results = self.v_eDR3.query_object(self.target, radius=1*u.arcmin)
            # Get Hip data
            dhip['S_N'] = results[0]['snrPMaH2EG3a'].data[0]    
            dhip['PMa_ra'] = results[0]['PMaRAH2EG3a'].data[0]
            dhip['PMa_ra_err'] = results[0]['e_PMaRAH2EG3a'].data[0]
            dhip['PMa_dec'] = results[0]['PMaDEH2EG3a'].data[0]
            dhip['PMa_dec_err'] = results[0]['e_PMaDEH2EG3a'].data[0]
            # Get eDR3 data
            dgaia['S_N'] = results[0]['snrPMaH2EG3b'].data[0]
            dgaia['PMa_ra'] = results[0]['PMaRAH2EG3b'].data[0]
            dgaia['PMa_ra_err'] = results[0]['e_PMaRAH2EG3b'].data[0]
            dgaia['PMa_dec'] = results[0]['PMaDEH2EG3b'].data[0]
            dgaia['PMa_dec_err'] = results[0]['e_PMaDEH2EG3b'].data[0]

            dparams['parallax'] = results[0]['PlxG3'].data[0] # mas
            dparams['dpc'] = 1./(dparams['parallax']*1.0e-3)
            
            dparams['mstar'] = results[0]['M1'].data[0]

        elif gaia == 'DR2':
            # TODO update with Kervella+19 (if needed)
            print('Need to update')
            sys.exit()
            # # Get Hip data
            # d['S_N_hip'] = results[0]['snrPMaH2G2'].data[0]    
            # d['PMa_hip_ra'] = results[0]['PMaRAH2G2'].data[0]
            # d['PMa_hip_ra_err'] = results[0]['e_PMaRAH2G2'].data[0]
            # d['PMa_hip_dec'] = results[0]['PMaDEH2EG3a'].data[0]
            # d['PMa_hip_dec_err'] = results[0]['e_PMaDEH2EG3a'].data[0]
            # # Get eDR3 data
            # d['S_N_dr3'] = results[0]['snrPMaH2EG3b'].data[0]
            # d['PMa_dr3_ra'] = results[0]['PMaRAH2EG3b'].data[0]
            # d['PMa_dr3_ra_err'] = results[0]['e_PMaRAH2EG3b'].data[0]
            # d['PMa_dr3_dec'] = results[0]['PMaDEH2EG3b'].data[0]
            # d['PMa_dr3_dec_err'] = results[0]['e_PMaDEH2EG3b'].data[0]

        else:
            print('wrong gaia epoch')
            sys.exit()

        return dhip, dgaia, dparams
    
    def read_param(self, d, iparam):
        if iparam:
            d['PA'] = self.discdata['PA']['median'] # deg
            d['PA_err'] = self.discdata['PA']['sigma'] # deg

            if 'l' in self.discdata['inc'].keys():
                # lower limit
                print('inc lower limit of ', self.discdata['inc']['l'])
                d['inc'] = (self.discdata['inc']['l']+90.)/2.
                d['inc_err'] = (90.-self.discdata['inc']['l'])/3.
                
            elif 'u' in self.discdata['inc'].keys():
                # upper limit
                print('inc upper limit of ', self.discdata['inc']['u']) 
                d['inc'] = self.discdata['inc']['u']/2.
                d['inc_err'] = d['inc']/3.
                
            else: 
                d['inc'] = self.discdata['inc']['median'] # deg
                d['inc_err'] = self.discdata['inc']['sigma'] # deg
        else:
            d['PA'] = self.discdata['PA']
            d['PA_err'] = self.discdata['PA_err']
            d['inc'] = self.discdata['inc']
            d['inc_err'] = self.discdata['inc_err']

        return d
    
    def get_data(self, gaia, iparam):
        d= {}
        # dparams = {'star':[self.target]}
        dhip, dgaia, dparams = self.vizier_query(gaia)  
        dparams = self.read_param(dparams, iparam)

        d['Hipparcos'] = pd.DataFrame(data=dhip).set_index('star')
        d[gaia] = pd.DataFrame(data=dgaia).set_index('star')
        d['params'] = pd.DataFrame(data=dparams).set_index('star')

        return pd.concat(d, axis=1)
