#Assignment 2
import matplotlib.pyplot as pt
import numpy as np
import astropy.units as u
from astropy.table import Table
from dust_extinction.parameter_averages import CCM89, F99
from synphot import units, config
from synphot import SourceSpectrum,SpectralElement,Observation,ExtinctionModel1D
from synphot.models import BlackBodyNorm1D
from synphot.spectrum import BaseUnitlessSpectrum
from synphot.reddening import ExtinctionCurve
from astroquery.simbad import Simbad
from astroquery.mast import Observations
from astropy import visualization

# ================= accuaring spectrum data =======================================================
# Data downloaded to folder /home/jibak/mastDownload/IUE

# reading HDE 251204 data (UV spectrum)
t_lwr11609 = Table.read('/home/jibak/mastDownload/IUE/lwr11609/lwr11609mxlo_vo.fits') #1800A - 3000 A
t_swp15081 = Table.read('/home/jibak/mastDownload/IUE/swp15081/swp15081mxlo_vo.fits') # 1200 - 2000 A
wav_UV_2 = t_lwr11609['WAVE'][0,].quantity
UVflux_2 = t_lwr11609['FLUX'][0,].quantity #reading table as quantity vector
wav_UV_1 = t_swp15081['WAVE'][0,].quantity
UVflux_1 = t_swp15081['FLUX'][0,].quantity

# reading HD 63922 data (UV spectrum)
t_lwr08237 = Table.read('/home/jibak/mastDownload/IUE/lwr08237/lwr08237mxlo_vo.fits') # 1800 - 3000 A
t_swp09511 = Table.read('/home/jibak/mastDownload/IUE/swp09511/swp09511mxlo_vo.fits') # 1200 - 2000 A
wav_UV_4 = t_lwr08237['WAVE'][0,].quantity
UVflux_4 = t_lwr08237['FLUX'][0,].quantity
wav_UV_3 = t_swp09511['WAVE'][0,].quantity
UVflux_3 = t_swp09511['FLUX'][0,].quantity

# U,V,B band data from simbad
custom_query = Simbad()
custom_query.add_votable_fields('fluxdata(U)','fluxdata(B)','fluxdata(V)')
#--------------------------------
phot_table=custom_query.query_object('HD 63922')
Umag_HD63922=phot_table['FLUX_U']
Bmag_HD63922=phot_table['FLUX_B']
Vmag_HD63922=phot_table['FLUX_V']
#---------------------------------
phot_table=custom_query.query_object('HDE 251204')
Umag_HDE251204=phot_table['FLUX_U']
Bmag_HDE251204=phot_table['FLUX_B']
Vmag_HDE251204=phot_table['FLUX_V']
#------------------------------------------------
#zero flux pass band
wav_U = 0.3660 * u.micron
zeroflux_U_nu = 1.81E-23 * u.Watt/(u.m*u.m*u.Hz)
wav_B = 0.4400 * u.micron
zeroflux_B_nu = 4.26E-23 * u.Watt/(u.m*u.m*u.Hz)
wav_V = 0.5530 * u.micron
zeroflux_V_nu = 3.64E-23 * u.Watt/(u.m*u.m*u.Hz)
# converting to flux per wavelength unit
zeroflux_U = zeroflux_U_nu.to(u.erg/u.AA/u.cm/u.cm/u.s,equivalencies=u.spectral_density(wav_U))
zeroflux_B = zeroflux_B_nu.to(u.erg/u.AA/u.cm/u.cm/u.s,equivalencies=u.spectral_density(wav_B))
zeroflux_V = zeroflux_V_nu.to(u.erg/u.AA/u.cm/u.cm/u.s,equivalencies=u.spectral_density(wav_V))

#--------------------------------------------------
#Flux of U,V,B band
Uflux_HD63922 = zeroflux_U * 10.**(-0.4*Umag_HD63922)
Bflux_HD63922 = zeroflux_B * 10.**(-0.4*Bmag_HD63922)
Vflux_HD63922 = zeroflux_V * 10.**(-0.4*Vmag_HD63922)
Uflux_HDE251204 = zeroflux_U * 10.**(-0.4*Umag_HDE251204)
Bflux_HDE251204 = zeroflux_B * 10.**(-0.4*Bmag_HDE251204)
Vflux_HDE251204 = zeroflux_V * 10.**(-0.4*Vmag_HDE251204)
#--------------------------------------------------
#======================================================================================================

#11111111111111111111        Extinction curve      111111111111111111111111111111111111111111111
Evb=0.71 # fromthe paper
Rv=3.1 #guessed
E_swp=-2.5*(np.log(UVflux_1.value/UVflux_3.value))/Evb/np.log(10)-(Vmag_HDE251204-Vmag_HD63922)/Evb
E_lwr=-2.5*np.log(UVflux_2.value[:562]/UVflux_4.value[:562])/Evb/np.log(10)-(Vmag_HDE251204-Vmag_HD63922)/Evb
pt.plot(2E4/(wav_UV_2.value[:562]+wav_UV_4.value[:562]),E_lwr,'g-',label="LWR")
pt.plot(2E4/(wav_UV_1.value+wav_UV_3.value),E_swp,'y-',label="SWP")
wav = np.arange(0.1, 3.0, 0.001)*u.micron
ext=CCM89(Rv=3.1)
pt.plot(1/wav, (ext(wav)-1)*Rv,'r--',label="CCM89")
pt.xlabel("$\lambda^{-1} (\mu^{-1})$")
pt.ylabel("$\dfrac{E(\lambda - V)}{E(B-V)}$")
pt.title("Extinction curve of HDE 251204")
pt.legend()
pt.show()
#________________________________________________________________________________________________

# 22222222222222222222    Blackbody spectrum  2222222222222222222222222222222222222222222222222222
# effective temperature of B0III star = 29200 K
sp = SourceSpectrum(BlackBodyNorm1D, temperature=29200)
v_band = SpectralElement.from_filter('johnson_v') 
vega = SourceSpectrum.from_vega() # For unit conversion  
sp_norm = sp.normalize(4.11 * units.VEGAMAG, v_band, vegaspec=vega) # crating blackbody spectrum with V magnitude = HD63922
wave=np.arange(1100,3300,10)  # Angstrom
pt.plot(wave,sp(wave,flux_unit=units.FLAM),'r',label="Blackbody 29200K")
pt.xlabel("$\lambda (\AA)$")
pt.ylabel("Flux $erg s^{−1} cm^{−2} \AA ^{-1}$")
pt.title("Blackbody spectrum in UV range of B0III type")
pt.legend()
pt.show()
#_________________________________________________________________________________________________
#3333333333333333333333 Redenning Blackbody curve 333333333333333333333333333333333333
Rd=sp_norm(wave,flux_unit=units.FLAM)*ext.extinguish(wave*u.AA,Ebv=0.71) # Reddened Blackbody curve
pt.subplot(1,2,1)
pt.semilogy(wave*u.AA,Rd,'r-',label="Reddened spectra")
pt.semilogy(wave*u.AA,sp_norm(wave,flux_unit=units.FLAM),'b',label="Unextinguished spectra")
pt.title("Theoretical plot")
pt.xlabel("$\lambda (\AA)$")
pt.legend()
pt.ylabel("Flux $erg  s^{−1} cm^{−2} \AA ^{-1}$")
pt.subplot(1,2,2)
pt.semilogy(wav_UV_1,UVflux_1,'m-',wav_UV_2,UVflux_2,'c-',label="HDE 251204")
pt.semilogy(wav_UV_3,UVflux_3,'y-',wav_UV_4,UVflux_4,'g-',label="HD 63922")
pt.xlabel("$\lambda (\AA)$")
pt.ylabel("Flux $erg s^{−1} cm^{−2} \AA ^{-1}$")
pt.title("Observed spectrum")
pt.tight_layout()
pt.legend()
pt.show()
#======= plotting data =====================


