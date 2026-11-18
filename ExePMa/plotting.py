"""
Plotting helpers for PMa companion-mass curves and complementary constraints
(disc extent / Hill stability, inner-edge sculpting, diffusion, Gaia RUWE).

All functions draw on a matplotlib Axes with separation [au] on x and planet
mass [Mjup] on y, and return the Axes.
"""
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
        """Convert separation [au] to angular separation [arcsec] at the star's distance."""
        return x / self.data['params'].iloc[0]['dpc']
    
    def arc2au(self, x):
        """Convert angular separation [arcsec] to separation [au] at the star's distance."""
        return x * self.data['params'].iloc[0]['dpc']

    def plot_pma(self, ax, epoch, **kwargs):
        """
        Plot the companion mass needed to explain the PMa versus separation.

        If the PMa signal-to-noise S_N >= snr, the median mass is drawn with shaded
        1/2/3-sigma bands. Otherwise nothing is drawn, unless upper_limit is set, in
        which case the +3sigma mass is drawn as an upper limit with the region above it
        hatched.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to draw on.
        epoch : {'eDR3', 'Hipparcos'}
            Epoch of the PMa.
        plot_sigmas : list of int, default [1, 2, 3]
            Sigma bands to shade: [1], [1, 2] or [1, 2, 3] (1 sigma is always drawn).
        plotting_params : dict, optional
            Any of:
            - snr (float, default 3.0): minimum PMa S/N for a detection curve.
            - upper_limit (bool, default False): draw upper limit when S_N < snr.
            - min_x, max_x (float, default None): only plot separations in this
              range (applied after x_scaling).
            - color (list, default ['C1', 'C0']): colours for Hipparcos, Gaia.
            - alpha (float, default 0.1): opacity of each sigma band.
            - label (str, default '<epoch> PMa'): legend label.
            - x_scaling (float, default 1): multiplies the separations [au] before
              plotting, e.g. 1/r_disc to plot a_p/r_disc. Only the PMa curve is
              rescaled: disc extent, RUWE cut and the default axis labels stay in au,
              and the top arcsec axis is wrong when x_scaling != 1 (hide it with
              upper_label=False, upper_ticks=False and set xlabel in
              add_axis_info()).
        pma_params : dict, optional
            Passed to mass_retrieval().

        Returns
        -------
        ax : matplotlib.axes.Axes
        aps : ndarray
            Plotted separations (scaled and masked) [au, times x_scaling].
        """
        pma_levels = kwargs.get('plot_sigmas', [1, 2, 3])
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
        x_scaling = plotting_params.get('x_scaling', 1)

        pma_data = self.mass_retrieval(epoch=epoch, **kwargs)

        aps = x_scaling * pma_data[:,0]
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
            ax.plot(aps, ms[:,3], color=color[0] if epoch=='Hipparcos' else color[1], label=pma_label) # median ms value
            ax.fill_between(aps, ms[:,2], ms[:,4], alpha=alpha, color=color[0] if epoch=='Hipparcos' else color[1]) # +/- 1sigma
            if 2 in pma_levels:
                ax.fill_between(aps, ms[:,1], ms[:,5], alpha=alpha, color=color[0] if epoch=='Hipparcos' else color[1]) # +/- 2sigma
            if 3 in pma_levels:
                ax.fill_between(aps, ms[:,0], ms[:,6], alpha=alpha, color=color[0] if epoch=='Hipparcos' else color[1]) # +/- 3sigma
        else: # upper limits when significance < snr
            if upper_limit:
                ax.plot(aps, ms[:,6], color=color[0] if epoch=='Hipparcos' else color[1], label=f'Upper limit {epoch} PMa')
                ax.fill_between(aps, ms[:,6], np.ones_like(aps)*1.0e3, color=color[0] if epoch=='Hipparcos' else color[1], alpha=0.3, hatch='//')

        return ax, aps
    
    def add_axis_info(self, ax, **kwargs):
        r"""
        Set log axes, limits, labels and legend, and add a top axis in arcsec.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to format.
        plotting_params : dict, optional
            Any of:
            - xlim (list, default [0.3, 300]): x limits [au].
            - ylim (list, default [0.3, 100]): y limits [Mjup].
            - xlabel (str, default r'$a_{\rm p}$ [au]'): bottom x-axis label.
            - ylabel (str, default r'$M_{\rm p}$ [$M_\mathrm{jup}$]'): y-axis label.
            - upper_xlabel (str, default r'$a_{\rm p}$ [arcsec]'): top-axis label.
            - upper_label (bool, default True): show upper_xlabel.
            - upper_ticks (bool, default True): show ticks on the top axis.

        Returns
        -------
        matplotlib.axes.Axes
        """
        # read in any plotting kwargs params
        plotting_params = kwargs.get('plotting_params', {})
        # default values if not in plotting_params
        xlim = plotting_params.get('xlim', [0.3, 300])
        ylim = plotting_params.get('ylim', [3.0e-1, 1.0e2])
        xlabel = plotting_params.get('xlabel', r'$a_{\rm p}$ [au]')
        ylabel = plotting_params.get('ylabel', r'$M_{\rm p}$ [$M_\mathrm{jup}$]')
        upper_xlabel = plotting_params.get('upper_xlabel', r'$a_{\rm p}$ [arcsec]')
        show_uplabel = plotting_params.get('upper_label', True)
        show_upticks = plotting_params.get('upper_ticks', True)

        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.legend(frameon=False, loc=3)

        ax.tick_params(which='both', top=False)
        secax = ax.secondary_xaxis('top', functions=(self.au2arc, self.arc2au))
        if show_uplabel:
            secax.set_xlabel(upper_xlabel)
        if show_upticks == False:
            secax.set_xticks([])
        
        return ax
    
    def plot_disc_extent(self, ax, aps, discdata, mstar, plot_lines=False, **kwargs):
        """
        Plot the disc edges and the region excluded by disc stability.

        A planet of mass mp at a clears NRhill Hill radii,
        R_H = a (mp / 3 mstar)^(1/3). Planets whose cleared zone would overlap the
        disc (or any disc ring between gaps) are excluded; this region is hatched.
        Dashed vertical lines mark the disc and gap edges.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to draw on.
        aps : ndarray
            Separations [au] at which to evaluate the constraint.
        discdata : dict
            - r_in, r_out (float): disc inner and outer edges [au].
            - NRhill (float, default 3.0): number of Hill radii cleared.
            - gaps (int, default 0): number of gaps in the disc.
            - gaps_in, gaps_out (list of float): inner and outer edges of each
              gap [au]; required if gaps > 0.
        mstar : float
            Stellar mass [Msun].
        plot_lines : bool, default False
            Also draw the individual inner/outer-edge curves.
        plotting_params : dict, optional
            disc_color (default 'C1'), disc_alpha (default 0.3),
            disc_zorder (default 1).

        Returns
        -------
        matplotlib.axes.Axes

        Raises
        ------
        ValueError
            If gaps > 0 and gaps_in/gaps_out are missing or of wrong length.
        """
        # read in any plotting kwargs params
        plotting_params = kwargs.get('plotting_params', {})

        # default values if not in plotting_params
        disc_color = plotting_params.get('disc_color', 'C1')
        disc_alpha = plotting_params.get('disc_alpha', 0.3)
        disc_zorder = plotting_params.get('disc_zorder', 1)

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
            ax.axvline(rin_i, color='grey', ls='--', lw=1., zorder=disc_zorder)
            ax.axvline(rout_i, color='grey', ls='--', lw=1., zorder=disc_zorder)

            Mpldisc_i = Mpldisc[i,:]
            Mpldisc_i[aps<rin_i]=np.maximum(3/(NRhill**3) *mstar*self.M_SUN/self.M_JUP * ( (rin_i/aps[aps<rin_i] -1.) )**(3.), Mpldisc_i[aps<rin_i])
            Mpldisc_i[aps>rout_i]=np.maximum(3/(NRhill**3) *mstar*self.M_SUN/self.M_JUP * ( ( 1.-rout_i/aps[aps>rout_i]))**(3.),Mpldisc_i[aps>rout_i])
            if plot_lines:
                ax.plot(aps, 3/(NRhill**3) *mstar*self.M_SUN/self.M_JUP * ( (rin_i/aps -1.) )**(3.))
                ax.plot(aps, 3/(NRhill**3) *mstar*self.M_SUN/self.M_JUP * ( ( 1.-rout_i/aps))**(3.))

            Mpldisc[i,:] = Mpldisc_i

        if gaps == 0:
            ax.fill_between(aps, Mpldisc[0], np.ones(len(aps))*1.0e3, color=disc_color, alpha=disc_alpha, hatch='//',
                            label=rf'Disc {NRhill:g}R$_\mathrm{{Hill}}$',
                            zorder=disc_zorder)
        else:
            # for i in range(gaps):
            #     gap_rin = disc_extent[2*i+1]
            #     gap_rout = disc_extent[2*i+2]
            #     mask = (aps >= gap_rin) & (aps <= gap_rout)
            #     Mpldisc[:, mask] = 1e3

            ax.fill_between(aps, np.min(Mpldisc, axis=0), np.ones(len(aps))*1.0e3, color=disc_color, alpha=disc_alpha, hatch='//',
                            label=rf'Disc {NRhill:g}R$_\mathrm{{Hill}}$',
                            zorder=disc_zorder)

        return ax

    def single_planet_inner_edge_sculptor(self, ax, aps, discdata, mstar, min_mp, **kwargs):
        """
        Plot the mass-separation relation of a single planet sculpting the disc inner edge.

        The planet's NRhill Hill radii reach exactly the inner edge:
        mp = 3 mstar / NRhill^3 * (r_in / a - 1)^3.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to draw on.
        aps : ndarray
            Planet separations [au].
        discdata : dict
            r_in (disc inner edge [au]) and optionally NRhill (default 3.0).
        mstar : float
            Stellar mass [Msun].
        min_mp : float or None
            Only draw the relation where mp >= min_mp [Mjup] (e.g. the minimum mass
            able to sculpt the edge); None to skip the curve.
        ap_err : float, optional
            If given, draw a horizontal band of half-width ap_err [au] instead of a line.
        min_mp_multiple : float, optional
            Minimum planet mass [Mjup] in the multi-planet scenario; drawn as a point on the
            relation at the corresponding separation. Works with or without min_mp.
        mp_err : float, default 1
            Vertical error bar on the min_mp_multiple point [Mjup].
        plotting_params : dict, optional
            line_color (default 'r'), multi_mp_color (default
            'darkorange'), alpha (default 0.5), label
            (default 'Inner edge sculptor').

        Returns
        -------
        matplotlib.axes.Axes
        """
        plotting_params = kwargs.get('plotting_params', {})
        line_color = plotting_params.get('line_color', 'r')
        multi_mp_color = plotting_params.get('multi_mp_color', 'darkorange')
        alpha = plotting_params.get('alpha', 0.5)
        label = plotting_params.get('label', 'Inner edge sculptor')

        r_in = discdata.get('r_in')
        NRhill = discdata.get('NRhill', 3.0)

        ap_err = kwargs.get('ap_err', None)
        min_mp_multiple = kwargs.get('min_mp_multiple', None)
        mp_err = kwargs.get('mp_err', 1)

        mstar_mjup = mstar*self.M_SUN/self.M_JUP # stellar mass [Mjup]

        if min_mp is not None:
            mps = 3/(NRhill**3) * mstar_mjup * ( (r_in/aps -1.) )**(3.)
            cond = (mps >= min_mp)
            if ap_err is None:
                ax.plot(aps[cond], mps[cond], color=line_color, alpha=alpha, ls='-', lw=2., zorder=1, label=label)
            else:
                ax.fill_betweenx(mps[cond], aps[cond] - ap_err, aps[cond] + ap_err, alpha=alpha, color=line_color, zorder=1, label=label)

        if min_mp_multiple is not None:
            # separation at which a planet of mass min_mp_multiple sculpts r_in (inverse of mps relation above)
            ap_multiple = r_in / (1. + (np.asarray(min_mp_multiple) * NRhill**3 / (3. * mstar_mjup))**(1./3.))
            ax.errorbar(ap_multiple, min_mp_multiple, xerr=ap_err, yerr=mp_err, fmt='o', markersize=5, mec='black', color=multi_mp_color, label=r'min m$_{plt}$ (multiple planets)')
        
        return ax

    def plot_diffusion(self, ax, aps, discdata, mstar, age):
        """
        Plot the minimum planet mass to sculpt the disc inner edge within the system age.

        Diffusion-timescale argument, Eq. 7 of Pearce et al. (2022).

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to draw on.
        aps : ndarray
            Planet separations [au].
        discdata : dict
            r_in: disc inner edge [au].
        mstar : float
            Stellar mass [Msun].
        age : float
            System age [Myr].

        Returns
        -------
        matplotlib.axes.Axes
        """
        r_in = discdata.get('r_in')

        diff = 0.331*aps*r_in**(-1/4)*age**(-1/2)*mstar**(3/4) # Eq 7 Pearce2022

        ax.plot(aps, diff, color='r', ls='--', lw=1., zorder=1, label=r'Diffusion timescale')
        return ax
    
    def ruwe_cutoff(self, ax, dpc, mstar, ruwe):
        """
        Plot the region excluded by a low Gaia RUWE (Limbach+2024; Kiefer+2024, paper I).

        If ruwe < 1.4 the star shows no significant astrometric excess noise, which
        rules out companions producing a large photocentre wobble. For periods shorter
        than the Gaia eDR3 baseline (1038 d) the limit scales as mp ∝ dpc*mstar / a
        (fixed astrometric signal); for longer periods it rises as a^2. The
        excluded region (0.1-20 au) is hatched. Nothing is drawn if ruwe >= 1.4.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to draw on.
        dpc : float
            Distance [pc].
        mstar : float
            Stellar mass [Msun].
        ruwe : float
            Gaia RUWE of the star (data['params']['ruwe'] or manual input).

        Returns
        -------
        matplotlib.axes.Axes
        """

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
