import matplotlib.pyplot as plt
import numpy as np
import os

class Plotting():
    def __init__(self):
        pass

    def au2arc(self, x):
        return x / self.data['params']['dpc'][0]
    
    def arc2au(self, x):
        return x * self.data['params']['dpc'][0]
    
    def saving_path(self, subdir):
        
        dir = self.odir + subdir + '/'
        if self.discdata['PA'] is None:
            if self.discdata['inc'] is None:
                file_path = 'Figure_PMa_{}_noPA_noInc.pdf'.format(self.target)
            else:
                file_path = 'Figure_PMa_{}_noPA.pdf'.format(self.target)
        else:
            if self.discdata['inc'] is None:
                file_path = 'Figure_PMa_{}_noInc.pdf'.format(self.target)
            else:
                file_path = 'Figure_PMa_{}.pdf'.format(self.target)

        if not os.path.exists(dir):
            os.makedirs(dir)
        
        return dir+file_path
        
    def plot(self, data_hip, data_gaia, figsize, only_save = False, subdir=None, snr=3.0, disc=False, mstar=False):
        # Retrieve radius and percentiles
        aps = data_gaia[:,0]
        ms_gaia = data_gaia[:,1:]
        ms_hip = data_hip[:,1:]

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        if self.data['Hipparcos']['S_N'][0] >= snr:
            plt.plot(aps, ms_hip[:,3], color='C1', label='Hipp PMa')
            plt.fill_between(aps, ms_hip[:,2], ms_hip[:,4], alpha=0.1, color='C1')
            plt.fill_between(aps, ms_hip[:,1], ms_hip[:,5], alpha=0.1, color='C1')
            plt.fill_between(aps, ms_hip[:,0], ms_hip[:,6], alpha=0.1, color='C1')
        else:
            plt.plot(aps, ms_hip[:,6], color='C1', label='Upper limit Hipp PMa')

        if self.data[self.gaia]['S_N'][0] >= snr:
            ax.plot(aps, ms_gaia[:,3], color='C0', label='eDR3 PMa')
            ax.fill_between(aps, ms_gaia[:,2], ms_gaia[:,4], alpha=0.2, color='C0')
            ax.fill_between(aps, ms_gaia[:,1], ms_gaia[:,5], alpha=0.2, color='C0')
            ax.fill_between(aps, ms_gaia[:,0], ms_gaia[:,6], alpha=0.2, color='C0')
        else:
            ax.plot(aps, ms_gaia[:,6], color='C0', label='Upper limit eDR3 PMa')

        if disc != False:
            rin, rout = disc
            ax.axvline(rin, color='grey', ls='--', lw=1.)
            ax.axvline(rout, color='grey', ls='--', lw=1.)

            Mpldisc=np.ones_like(aps)*1.0e-2
            NRhill=3.

            Mpldisc[aps<rin]=np.maximum(3/NRhill**3 *mstar * ( (rin/aps[aps<rin] -1.) )**(3.), Mpldisc[aps<rin])
            Mpldisc[aps>rout]=np.maximum(3/NRhill**3 *mstar * ( ( 1.-rout/aps[aps>rout]))**(3.),Mpldisc[aps>rout])
            ax.fill_between(aps, Mpldisc, np.ones(len(aps))*1.0e3, color='C1', alpha=0.4, hatch='//')
        
        ax.set_ylabel(r'$M_{\rm p}$ [$M_\mathrm{jup}$]')
        ax.set_xlabel(r'$a_{\rm p}$ [au]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(aps[0], aps[-1])
        ax.set_ylim(3.0e-1, 1.0e3)
        ax.legend(frameon=False, loc=3, fontsize=8)

        ax.tick_params(which='both', top=False)

        secax = ax.secondary_xaxis('top', functions=(self.au2arc, self.arc2au))
        secax.set_xlabel(r'$a_{\rm p}$ [arcsec]')

        plt.tight_layout()
        if subdir is not None:
            save_path = self.saving_path(subdir=subdir)
            if only_save:
                print('Saving figure...')
                plt.savefig(save_path)
            else:
                plt.savefig(save_path)
                print('Saving figure...')
                plt.show()
        else:
            plt.show()
