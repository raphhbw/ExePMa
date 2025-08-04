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
        pma_label = plotting_params.get('label', f'{epoch} PMa')

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
    
    def plot_disc_extent(self, ax, aps, discdata, mstar, plot_lines=False, **kwargs):
        """ Plot the disc extent + rhill argument for edge truncation.
            Kwargs options:
                - plotting_params [dict]:
                    - disc_color [default: C1]: colour given to the disc
                    - disc_alpha [default 0.3]
                - discdata [dict]:
                    - r_in: inner edge location in au
                    - r_out: outer edge location in au
                    - NRhill [default 3.0]: number of hill radii 
                    - gaps [default 0]: number of gaps in the disc
                    - gaps_in [default None]: inner gaps in the disc (as a list), if gaps > 0
                    - gaps_out [default None]: outer gaps in the disc (as a list), if gaps > 0
        """
        # read in any plotting kwargs params
        plotting_params = kwargs.get('plotting_params', {})
        # default values if not in plotting_params
        disc_color = plotting_params.get('disc_color', 'C1')
        disc_alpha = plotting_params.get('disc_alpha', 0.3)

        r_in = discdata.get('r_in')
        r_out = discdata.get('r_out')
        NRhill = discdata.get('NRhill', 3.0)
        gaps = discdata.get('gaps', 0)
        gaps_in = discdata.get('gaps_in', None)
        gaps_out = discdata.get('gaps_out', None)

        disc_extent = [r_in, r_out]
        
        if gaps_in is None:
            if gaps > 0:
                raise ValueError('gaps_in values for the disc not provided')
        else:
            if len(gaps_in) != gaps:
                raise ValueError('Review gaps_in values for the disc, not equal to gaps')
            else:
                disc_extent = np.concatenate((disc_extent, gaps_in))
        
        if gaps_out is None:
            if gaps > 0:
                raise ValueError('gaps_out values for the disc not provided')
        else:
            if len(gaps_out) != gaps:
                raise ValueError('Review gaps_out values for the disc, not equal to gaps')
            else:
                disc_extent = np.concatenate((disc_extent, gaps_out))

        # reorder disc extent
        disc_extent = np.sort(disc_extent)

        Mpldisc = np.ones((gaps+1, len(aps))) * 1.0e-2

        for i in range(gaps+1):
            j = i+i
            rin_i = disc_extent[j]
            rout_i = disc_extent[j+1]

            # Plot edges
            ax.axvline(rin_i, color='grey', ls='--', lw=1., zorder=1)
            ax.axvline(rout_i, color='grey', ls='--', lw=1., zorder=1)

            Mpldisc_i = Mpldisc[i,:]
            Mpldisc_i[aps<rin_i]=np.maximum(3/(NRhill**3) *mstar*self.M_SUN/self.M_JUP * ( (rin_i/aps[aps<rin_i] -1.) )**(3.), Mpldisc_i[aps<rin_i])
            Mpldisc_i[aps>rout_i]=np.maximum(3/(NRhill**3) *mstar*self.M_SUN/self.M_JUP * ( ( 1.-rout_i/aps[aps>rout_i]))**(3.),Mpldisc_i[aps>rout_i])
            if plot_lines:
                ax.plot(aps, 3/(NRhill**3) *mstar*self.M_SUN/self.M_JUP * ( (rin_i/aps -1.) )**(3.))
                ax.plot(aps, 3/(NRhill**3) *mstar*self.M_SUN/self.M_JUP * ( ( 1.-rout_i/aps))**(3.))

            Mpldisc[i,:] = Mpldisc_i

        if gaps == 0:
            ax.fill_between(aps, Mpldisc[0], np.ones(len(aps))*1.0e3, color=disc_color, alpha=disc_alpha, hatch='//', label='Disc 3R$_\mathrm{Hill}$', zorder=1)
        else:
            ax.fill_between(aps, np.min(Mpldisc, axis=0), np.ones(len(aps))*1.0e3, color=disc_color, alpha=disc_alpha, hatch='//', label='Disc 3R$_\mathrm{Hill}$', zorder=1)

        return ax
    
    def ruwe_cutoff(self, ax, dpc, mstar, ruwe):
        ''' Plot ruwe cutoff based on Limbach+2024 and Kiefer+2024 (paperI) '''

        if ruwe < 1.4:
            ap_gaia = ((1038*24*3600/(2*np.pi))**2 * self.G * mstar*self.M_SUN)**(1/3) / self.AU #1038 days from Gaia into au
            mp_gaia_cutoff = 2*np.sqrt(2*np.log(2))*0.19e-3 * dpc * mstar*self.M_SUN/self.M_JUP/ap_gaia

            ruwe_aps = np.linspace(0.1, 20, 300)

            inner_mps = 2*np.sqrt(2*np.log(2))*0.19e-3 * dpc * mstar*self.M_SUN/self.M_JUP/ruwe_aps
            b = np.log10(mp_gaia_cutoff)-2*np.log10(ap_gaia)

            outer_mps = 10**(2*np.log10(ruwe_aps) + b)

            ax.fill_between(ruwe_aps, np.maximum(inner_mps, outer_mps), np.ones(len(ruwe_aps))*1.0e3, color='grey', alpha=0.2, hatch='..', zorder=1, linestyle='--', edgecolor='k', label='RUWE cut')
            return ax
        else:
            print(f'RUWE > 1.4 ({ruwe}), no cut applied')
            return ax