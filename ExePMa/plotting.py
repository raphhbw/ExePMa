import numpy as np

class Plotting():
    ### units
    M_SUN=1.989e+30 # kg
    AU=1.496e+11   # m
    YEAR=365.24 # days
    G=6.67408e-11 # mks
    M_JUP=1.898e27 # kg
    
    def __init__(self):
        pass

    def au2arc(self, x):
        return x / self.data['params']['dpc'][0]
    
    def arc2au(self, x):
        return x * self.data['params']['dpc'][0]

    def plot_pma(self, ax, epoch, **kwargs):
        """ Plotting PMa curve function.
            Choose your epoch [eDR3 or Hipparcos] for your curve. If snr > 3 [default] curve will be plotted. 
            Kwargs:
                - plotting_params [dict]:
                    - snr [default: 3.0]: snr determines the strength of the PMa
                    - upper_limit [default: False]: plots upper limit when PMa not significant enough
                    - min_x, min_y [default None]: mask for aps
                    - color [default: [C1, C0]]: colour for Hipparcos, Gaia curves
                    - alpha [default: 0.1]
                - pma_params [dict] given to mass_retrieval function
        """
        # read in any plotting kwargs params
        plotting_params = kwargs.get('plotting_params', {})
        # default values if not in plotting_params
        snr = plotting_params.get('snr', 3.0)
        upper_limit = plotting_params.get('upper_limit', False)
        color = plotting_params.get('color', ['C1', 'C0'])
        alpha = plotting_params.get('alpha', 0.1)
        min_x = plotting_params.get('min_x', None)
        max_x = plotting_params.get('max_x', None)

        pma_data = self.mass_retrieval(epoch=epoch, **kwargs)

        aps = pma_data[:,0]
        ms = pma_data[:,1:]

        # Conditions on aps based on min_x and max_x
        if min_x is not None and max_x is not None: #min_x and max_x have values
            mask = (aps >= min_x) & (aps <= max_x)
            aps = aps[mask]
            ms = ms[mask]
        elif min_x is not None and max_x is None: #min_x has value, max_x is None
            mask = aps >= min_x
            aps = aps[mask]
            ms = ms[mask]
        elif min_x is None and max_x is not None: #min_x is None, max_x has value
            mask = aps <= max_x
            aps = aps[mask]
            ms = ms[mask]
        else: #min_x and max_x are None
            pass
        # Plot the PMa curve
        if self.data[epoch].iloc[0]['S_N'] >= snr:
            ax.plot(aps, ms[:,3], color=color[0] if epoch=='Hipparcos' else color[1], label=f'{epoch} PMa') # mean ms value
            ax.fill_between(aps, ms[:,2], ms[:,4], alpha=alpha, color=color[0] if epoch=='Hipparcos' else color[1]) # +/- 1sigma
            ax.fill_between(aps, ms[:,1], ms[:,5], alpha=alpha, color=color[0] if epoch=='Hipparcos' else color[1]) # +/- 2sigma
            ax.fill_between(aps, ms[:,0], ms[:,6], alpha=alpha, color=color[0] if epoch=='Hipparcos' else color[1]) # +/- 3sigma
        else: # upper limits when significance < snr
            if upper_limit:
                ax.plot(aps, ms[:,6], color=color[0] if epoch=='Hipparcos' else color[1], label=f'Upper limit {epoch} PMa')

        return ax, aps
    
    def add_axis_info(self, ax, **kwargs):
        """ Add legends + axis information for plot.
        Kwargs:
            - plotting_params [dict]:
                - xlim [default [0.3, 300]]: xlimit with [xmin, xmax]
                - ylim [default [3.0e-1, 1.0e2]]: ylimit with [ymin, ymax]
        """
        # read in any plotting kwargs params
        plotting_params = kwargs.get('plotting_params', {})
        # default values if not in plotting_params
        xlim = plotting_params.get('xlim', [0.3, 300])
        ylim = plotting_params.get('ylim', [3.0e-1, 1.0e2])

        ax.set_ylabel(r'$M_{\rm p}$ [$M_\mathrm{jup}$]')
        ax.set_xlabel(r'$a_{\rm p}$ [au]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.legend(frameon=False, loc=3)

        ax.tick_params(which='both', top=False)
        secax = ax.secondary_xaxis('top', functions=(self.au2arc, self.arc2au))
        secax.set_xlabel(r'$a_{\rm p}$ [arcsec]')
        return ax
    
    def plot_disc_extent(self, ax, aps, discdata, mstar, **kwargs):
        """ Plot the disc extent + rhill argument for edge truncation.
            Kwargs options:
                - plotting_params [dict]:
                    - disc_color [default: C1]: colour given to the disc
                    - disc_alpha [default 0.4]
                - discdata [dict]:
                    - r_in: inner edge location in au
                    - r_out: outer edge location in au
                    - NRhill [default 3.0]: number of hill radii 
        """
        # read in any plotting kwargs params
        plotting_params = kwargs.get('plotting_params', {})
        # default values if not in plotting_params
        disc_color = plotting_params.get('disc_color', 'C1')
        disc_alpha = plotting_params.get('disc_alpha', 0.4)

        r_in = discdata.get('r_in')
        r_out = discdata.get('r_out')
        NRhill = discdata.get('NRhill', 3.0)

        # Plot disc edges
        ax.axvline(r_in, color='grey', ls='--', lw=1.)
        ax.axvline(r_out, color='grey', ls='--', lw=1.)
        
        Mpldisc=np.ones_like(aps)*1.0e-2
        Mpldisc[aps<r_in]=np.maximum(3/NRhill**3 *mstar*self.M_SUN/self.M_JUP * ( (r_in/aps[aps<r_in] -1.) )**(3.), Mpldisc[aps<r_in])
        Mpldisc[aps>r_out]=np.maximum(3/NRhill**3 *mstar*self.M_SUN/self.M_JUP * ( ( 1.-r_out/aps[aps>r_out]))**(3.),Mpldisc[aps>r_out])
        ax.fill_between(aps, Mpldisc, np.ones(len(aps))*1.0e3, color=disc_color, alpha=disc_alpha, hatch='//', label='Disc')
        return ax