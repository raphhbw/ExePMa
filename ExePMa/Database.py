"""
Retrieval of Hipparcos-Gaia PMa data and stellar parameters for a single star.

Data come from the Kervella et al. (2022, A&A 657, A7) catalogue, either queried on
Vizier (J/A+A/657/A7) or read from a local copy.
"""
import pandas as pd
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import astropy.units as u
import numpy as np

class Database():
    """
    PMa data and stellar parameters of one star, stored in self.data.

    Normally used through PMa rather than directly.

    Parameters
    ----------
    star : str, optional
        Star name, resolved to a HIP number through SIMBAD (e.g. 'HD 107146').
    hip : str, optional
        Hipparcos identifier, e.g. 'HIP 60074'; skips the SIMBAD query.
        One of star or hip is required.
    gaia : {'eDR3'}, default 'eDR3'
        Gaia release used for the PMa. Only eDR3 is implemented.
    geometry : dict, optional
        Orientation of the system (e.g. of its debris disc) used to deproject the PMa.
        Keys: inc, inc_err, PA, PA_err [deg]. inc/PA default to
        None (unknown, marginalised over), their errors to 1 deg.
    mstar : float, optional
        Stellar mass [Msun] overriding the catalogue value M1.
    local_catalog : astropy.table.Table, optional
        Local copy of the Kervella+22 table (same column names as Vizier). If given,
        data are read from it instead of querying Vizier.

    Attributes
    ----------
    star : str or None
        Star name as given.
    hip : str
        Hipparcos identifier, 'HIP <number>'.
    geometry : dict
        Output of build_geometryinfo().
    gaia : str
        Gaia release.
    mstar : float or None
        Stellar mass override [Msun].
    data : pandas.DataFrame
        See get_data().
    """
    V_EDR3 = Vizier(columns=["**"], catalog="J/A+A/657/A7") # Kervella+22 Vizier Catalog

    def __init__(self, **kwargs):
        """Resolve the HIP identifier and retrieve the data. See class docstring for kwargs."""
        self.star = kwargs.get("star", None)
        self.hip = kwargs.get("hip", None)

        if self.star is None and self.hip is None:
            raise ValueError("Please provide either a star name (star) or a Hipparcos number (hip).")
        elif self.hip is None:
            self.hip = self.get_hip_number(self.star)

        self.geometry = self.build_geometryinfo(kwargs.get("geometry", {}))

        self.gaia = kwargs.get("gaia", 'eDR3')
        self.mstar = kwargs.get("mstar", None)

        self.data = self.get_data(self.gaia, kwargs.get("local_catalog", None))

    def build_geometryinfo(self, discinfo):
        """
        Fill in defaults for the system geometry.

        Parameters
        ----------
        discinfo : dict
            May contain inc, inc_err, PA, PA_err [deg].

        Returns
        -------
        dict
            Keys inc and PA (default None = unknown) and inc_err and
            PA_err (default 1.0) [deg].
        """
        return {"inc":discinfo.get('inc', None), # default None if no inc is given
                "inc_err":discinfo.get('inc_err', 1.0), # default 1.0 if no inc_err is given [TODO: check default]
                "PA":discinfo.get('PA', None), # default None if no PA is given
                "PA_err":discinfo.get('PA_err', 1.0), # default 1.0 if no PA_err is given [TODO: check default]
                }
    
    def get_hip_number(self, star_name):
        """
        Find the Hipparcos identifier of a star through SIMBAD.

        Parameters
        ----------
        star_name : str
            Any identifier SIMBAD can resolve.

        Returns
        -------
        str
            Identifier of the form 'HIP <number>'.

        Raises
        ------
        ValueError
            If the star is not found, is ambiguous, or has no HIP identifier.
        """
        # Customize Simbad to include the HIP identifier in the results
        custom_simbad = Simbad()
        custom_simbad.TIMEOUT = 10  # Set a timeout for the query
        custom_simbad.add_votable_fields('ids')  # Include identifiers in the query

        # Query SIMBAD for the star
        result = custom_simbad.query_object(star_name)
        
        if result is None:
            raise ValueError(f"Star '{star_name}' not found in SIMBAD.")
        elif len(result)>1:
            raise ValueError(f"Multiple entries found for star '{star_name}'. Please specify further.")

        # Extract the HIP number from the identifiers
        identifiers = result['IDS'][0]
        for identifier in identifiers.split('|'):
            if identifier.strip().startswith('HIP'):
                print(f"The HIP number of {star_name} is: {identifier.strip()}")
                return identifier.strip()
        
        raise ValueError(f"HIP number not found for star '{star_name}'.")

    def vizier_query(self, gaia):
        """
        Query the Kervella+22 catalogue on Vizier for self.hip.

        Searches within 1 arcmin of the star (1 deg if nothing is found) and selects
        the row matching the HIP number.

        Parameters
        ----------
        gaia : {'eDR3'}
            Gaia release. Only eDR3 is implemented.

        Returns
        -------
        dhip, dgaia, dparams : dict
            Hipparcos-epoch PMa, Gaia-epoch PMa and stellar parameters; see
            get_data() for the keys and units.
        """
        
        # Initialise the dictionaries
        dhip = {'star':[self.hip]}
        dgaia = {'star':[self.hip]}
        dparams = {'star':[self.hip]}

        if gaia == 'eDR3':
            # Query catalog
            results = self.V_EDR3.query_object(self.hip, radius=1*u.arcmin)
            query_okay = True
            if len(results) == 0:
                query_okay = False
            elif len(results[0]) == 0:
                query_okay = False
            
            if query_okay == False:
                results = self.V_EDR3.query_object(self.hip, radius=1*u.deg)
                
            # Get Hip data
            hip_name = int(self.hip.replace('HIP ',''))
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
            dparams['ruwe'] = results[0]['RUWE'].data[w] # RUWE
            if self.mstar is None:
                dparams['mstar'] = results[0]['M1'].data[w]
            else:
                dparams['mstar'] = self.mstar

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
    
    def get_data(self, gaia, local_catalog=None):
        """
        Build the data table of the star.

        Parameters
        ----------
        gaia : {'eDR3'}
            Gaia release.
        local_catalog : astropy.table.Table, optional
            If given, read from it with create_single_star_table();
            otherwise query Vizier with vizier_query().

        Returns
        -------
        pandas.DataFrame
            One row indexed by HIP identifier, with two-level columns:

            - 'Hipparcos' and gaia (e.g. 'eDR3'): PMa at that epoch, with
              S_N (signal-to-noise of the PMa), PMa_ra, PMa_ra_err,
              PMa_dec, PMa_dec_err [mas/yr].
            - 'params': parallax [mas], dpc [pc], ruwe, mstar [Msun],
              PA, PA_err, inc, inc_err [deg].

            Access a value with e.g. data['params']['dpc'].iloc[0].
        """
        d = {}
        
        # Query data from Kervella vizier catalog, or read it from a local copy
        if local_catalog is None:
            dhip, dgaia, dparams = self.vizier_query(gaia)
        else:
            dhip, dgaia, dparams = self.create_single_star_table(local_catalog)

        d['Hipparcos'] = pd.DataFrame(data=dhip).set_index('star')
        d[gaia] = pd.DataFrame(data=dgaia).set_index('star')
        d['params'] = pd.DataFrame(data=dparams).set_index('star')

        return pd.concat(d, axis=1)

    def create_single_star_table(self, local_catalog):
        """
        Same output as vizier_query(), read from a local copy of the catalogue.

        Parameters
        ----------
        local_catalog : astropy.table.Table
            Kervella+22 table with Vizier column names (HIP, snrPMaH2EG3a, ...).
            ruwe is set to NaN if there is no RUWE column. eDR3 only.

        Returns
        -------
        dhip, dgaia, dparams : dict
            As in vizier_query().

        Raises
        ------
        ValueError
            If self.hip is not in the table.
        """

        dhip = {'star':[self.hip]}
        dgaia = {'star':[self.hip]}
        dparams = {'star':[self.hip]}
        
        hip_names = np.array([f'HIP {i}' for i in local_catalog['HIP']])
        w, = np.where(hip_names == self.hip)
        if len(w) == 0:
            raise ValueError(f'Star {self.hip} was not found in the local catalog!')
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
        dparams['ruwe'] = float(local_catalog['RUWE'].data[w]) if 'RUWE' in local_catalog.colnames else np.nan # RUWE
        dparams['mstar'] = float(local_catalog['M1'].data[w]) if self.mstar is None else self.mstar

        # Add geometry info to star info
        dparams['PA'] = self.geometry['PA']
        dparams['PA_err'] = self.geometry['PA_err']
        dparams['inc'] = self.geometry['inc']
        dparams['inc_err'] = self.geometry['inc_err']

        return dhip, dgaia, dparams
    
