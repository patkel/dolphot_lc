from glob import glob

fs = glob('driz_cr/*15936*/diffs/singlestar')

import subprocess
print(fs)

import numpy

dates = []
fluxes = []
flux_errs = []
visits = []
chis = []


summary_dates = []
summary_fluxes = []
summary_flux_errs = []



for f in fs:
    import numpy


    try:

        statf = f.replace('singlestar','stats.txt')            

        sf = numpy.loadtxt(statf)



        lf = numpy.loadtxt(f)[0,:]                                                                                                                                  
                                                                                                                                                                    
        fl = open(f + '.columns').readlines()

        im_index = 0
                                                                                                                                                                    
        summary_added = False
                                                                                                                                                                    
        visit_dates = []
                                                                                                                                                                    
        for l in fl:
                                                                                                                                                                    
                                                                                                                                                                    
            if l.find('Normalized count rate') != -1 and l.find('_flc') == -1 and l.find('uncertainty') == -1:
                                                                                                                                                                    
                                                                                                                                                                    
                spl = l.split(' ')
                ind = int(float(spl[0][:-1])) - 1
                                                                                                                                                                    
                summary_flux = lf[ind]                                                                               
                summary_flux_err = lf[ind + 1]
                                                                                                                                                                    
            elif l.find('Normalized count rate') != -1 and l.find('_flc') != -1 and l.find('uncertainty') == -1:
                print(l)

                fake_err = sf[im_index]
                                                                                                                                                                    
                spl = l.split(' ')
                ind = int(float(spl[0][:-1])) - 1
                                                                                                                                                                    
                x_ref = lf[2]
                y_ref = lf[3]
                                                                                                                                                                    
                                                                                                                                                                    
                #big_x, big_y = w.wcs_world2pix(small_ra, small_dec, 1,
                #                           ra_dec_order=True)
                                                                                                                                                                    
                from astropy.wcs import WCS
                from astropy.io import fits
                                                                                                                                                                    
                                                                                                                                                                    
                if True:
                    w = WCS(fits.open('driz_cr/15936_42/diffs/coadd_sky_UVIS_tweak_SNAP_F200LP.chip1.fits')[0])                  
                    ra, dec = w.wcs_pix2world( x_ref + 0.5, y_ref + 0.5, 0 )
                                                                                                                      
                    print(ra,dec)
                                                                                                                      
                    direct = f.replace('singlestar', '')
                                                                                                                      
                    file = direct + spl[4] + '.fits'
                                                                                                                      
                    print(file)
                                                                                                                                                                    
                                                                                                                                                                    
                    forig = file.replace('diffs','Images').replace('.chip1','').replace('.chip2','')
                                                                                                                                                                    
                    print(forig)
                                                                                                                                                                    
                    if False:
                                                                                                                                                                    
                        if file.find('chip1') != -1:                                                                       
                            wf = WCS(fits.open(forig)[2])
                                                                                                                           
                        else:
                                                                                                                           
                            wf = WCS(fits.open(forig)[5])
                                                                                                                           
                                                                                                                          
                                                                                                                          
                                                                                                                          
                        big_x, big_y = wf.wcs_world2pix(ra, dec, 1,
                                                   ra_dec_order=True)
                                                                                                                          
                                                                                                                          
                        data = fits.open(file)
                                                                                                                          
                        val = data[0].data[int(big_y - 0.5), int(big_x - 0.5)]
                                                                                                                      
                                                                                                                      
                    flag = lf[ind + 10]
                                                                                                                      
                    print('flag', flag)
                                                                                                                      
                    chi = lf[ind + 5]

                    im_index += 1
                                                                                                                      
                    if flag <= 3: # and chi < 3.: # and val > 0: # manual says flags 1--3 are usuable
                                                                                                                      
                        flux = lf[ind]                                                                               
                        flux_err = lf[ind + 1]
                                                                                                                     
                                                                                                                      
                                                                                                                     
                                                                                                                      
                                                                                                                      
                                                                                                                      
                                                                                                                     
                        output = subprocess.run( [f'gethead',f'{file}', 'MJD-OBS'], capture_output=True, text=True )
                        mjdobs = float(output.stdout[:-1])
                                                                                                                     
                                                                                                                     
                        print(file, mjdobs, flux, flux_err, chi)
                                                                                                                                                                    
                        
                                                                                                                     
                        visit_dates.append(mjdobs)                                                                                                                 
                        dates.append(mjdobs)
                        fluxes.append(flux)
                        flux_errs.append(fake_err) #flux_err)
                        visits.append( file.split('/')[0].split('_')[1] + ' ' + spl[4] )
                                                                                                                                                                    
        if summary_added == False:
            summary_dates.append( numpy.mean(visit_dates) )
            summary_fluxes.append(summary_flux)
            summary_flux_errs.append(summary_flux_err)
            summary_added = True

    except:
        print('failed')

a = list(zip(dates, fluxes, flux_errs))

a.sort()

for b in a:
    print(b)

import pylab
from scipy.optimize import minimize


# ── Fold caustic light curve model ──────────────────────────────────────────────

def fold_caustic_flux(t, t_c, t_E, t_star, F0, F_caustic):
    """
    Flux from a uniform-disk source crossing a fold caustic.

    Parameters
    ----------
    t : array-like
        Times (same units as t_c / t_E / t_star).
    t_c : float
        Caustic crossing time (source centre on caustic).
    t_E : float
        Einstein ring crossing timescale (> 0).
    t_star : float
        Source angular-radius crossing time (> 0, controls peak width).
    F0 : float
        Baseline flux (outside caustic).
    F_caustic : float
        Caustic magnification amplitude (> 0).

    Returns
    -------
    numpy.ndarray
    """
    t = numpy.atleast_1d(numpy.asarray(t, dtype=float))

    # Gauss-Legendre quadrature over the uniform-disk surface
    n_quad = 64
    s_pts, s_wts = numpy.polynomial.legendre.leggauss(n_quad)
    disk_wts = s_wts * numpy.sqrt(numpy.maximum(1.0 - s_pts**2, 0.0))
    disk_norm = disk_wts.sum()   # normalisation ≈ π/2

    # Effective caustic-crossing time for each disk element  (n_quad, N_t)
    t_c_eff = t_c + s_pts[:, numpy.newaxis] * t_star
    u = (t[numpy.newaxis, :] - t_c_eff) / t_E

    # Point-source fold magnification: 1/sqrt(−u) inside (u<0), 0 outside
    eps = 1e-8
    mag_point = numpy.where(u < 0, 1.0 / numpy.sqrt(numpy.maximum(-u, eps)), 0.0)

    mag_disk = (disk_wts[:, numpy.newaxis] * mag_point).sum(axis=0) / disk_norm
    return F0 + F_caustic * mag_disk


def _chi2_caustic(params, t_data, F_data, F_err):
    t_c, t_E, t_star, F0, F_caustic = params
    if t_E <= 0 or t_star <= 0 or F_caustic <= 0:
        return 1e30
    F_model = fold_caustic_flux(t_data, t_c, t_E, t_star, F0, F_caustic)
    return float(numpy.sum(((F_data - F_model) / F_err) ** 2))


# ── Fit ─────────────────────────────────────────────────────────────────────────

t_data = numpy.array(dates) * 24        # MJD × 24 h  (matches existing x-axis)
F_data = numpy.array(fluxes)
F_err  = numpy.array(flux_errs)

# Shift to relative time to keep parameter values numerically well-conditioned
t0_ref = numpy.median(t_data)
t_rel  = t_data - t0_ref

idx_peak   = int(numpy.argmax(F_data))
t_c0       = t_rel[idx_peak]
F0_0       = numpy.median(F_data)
F_caustic0 = max(F_data[idx_peak] - F0_0, 1e-8)
t_span     = t_rel.max() - t_rel.min() if len(t_rel) > 1 else 1.0
t_E0       = max(t_span / 10.0, 0.1)
t_star0    = t_E0 / 4.0

p0 = [t_c0, t_E0, t_star0, F0_0, F_caustic0]
print(f"Initial params: t_c={t_c0:.3f} h  t_E={t_E0:.3f} h  t_*={t_star0:.3f} h  "
      f"F0={F0_0:.5f}  F_caustic={F_caustic0:.5f}")

result = minimize(
    _chi2_caustic,
    p0,
    args=(t_rel, F_data, F_err),
    method='Nelder-Mead',
    options={'xatol': 1e-8, 'fatol': 1e-8, 'maxiter': 100000, 'maxfev': 200000},
)

t_c_fit, t_E_fit, t_star_fit, F0_fit, F_caustic_fit = result.x
chi2_dof = result.fun / max(len(t_data) - 5, 1)

print("Best-fit caustic transit parameters:")
print(f"  t_c      = {t_c_fit + t0_ref:.4f}  (relative: {t_c_fit:+.4f} h from median)")
print(f"  t_E      = {t_E_fit:.4f} h")
print(f"  t_star   = {t_star_fit:.4f} h")
print(f"  F0       = {F0_fit:.6f}")
print(f"  F_caustic= {F_caustic_fit:.6f}")
print(f"  chi2/dof = {chi2_dof:.3f}")

# ── Plot ─────────────────────────────────────────────────────────────────────────

t_pad  = 0.25 * t_span
t_fine = numpy.linspace(t_rel.min() - t_pad, t_rel.max() + t_pad, 1000)
F_fine = fold_caustic_flux(t_fine, t_c_fit, t_E_fit, t_star_fit, F0_fit, F_caustic_fit)

pylab.figure(figsize=(10, 5))

# Individual-image data points
pylab.errorbar(t_rel, F_data, yerr=F_err,
               fmt='o', color='red', capsize=3, zorder=3, label='per-image data')

# Summary (visit-averaged) points
if summary_dates:
    t_sum = numpy.array(summary_dates) * 24 - t0_ref
    pylab.errorbar(t_sum, summary_fluxes, yerr=summary_flux_errs,
                   fmt='s', color='black', capsize=3, zorder=4, label='visit average')

# Best-fit model
pylab.plot(t_fine, F_fine, '-', color='royalblue', lw=2,
           label=(f'fold caustic fit\n'
                  f'$t_c$={t_c_fit:+.2f} h  '
                  f'$t_E$={t_E_fit:.2f} h  '
                  f'$t_*$={t_star_fit:.2f} h\n'
                  f'$F_0$={F0_fit:.4f}  '
                  f'$\chi^2$/dof={chi2_dof:.2f}'))

pylab.axhline(F0_fit, color='gray', linestyle='--', alpha=0.5, label=f'baseline $F_0$')
pylab.axhline(0, color='k', linestyle='-', alpha=0.2, lw=0.8)

pylab.xlabel(f'Time − {t0_ref:.1f}  (MJD × 24 h)')
pylab.ylabel('Normalized count rate')
pylab.legend(fontsize=8)
pylab.title('Stellar caustic transit')
pylab.tight_layout()
pylab.show()
