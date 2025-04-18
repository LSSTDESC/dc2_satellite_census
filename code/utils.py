#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner, Kabelo Tsiane"

import os
from os.path import join, exists
from collections import OrderedDict as odict

import pylab as plt
import numpy as np
import healpy as hp
import fitsio
import pandas as pd

from healpy import read_map
from healsparse import HealSparseMap
from IPython.core.debugger import set_trace

def setdefaults(kwargs,defaults):
    for k,v in defaults.items():
        kwargs.setdefault(k,v)
    return kwargs


def get_dirname(dirname_in):
    # First see if the path exists
    basedir = '/project/shared/data/satsim/'
    dirname = os.path.join(basedir,dirname_in)
    if os.path.exists(dirname):
        return dirname

    raise Exception("Could not find %s"%os.path.basename(dirname))


def get_datadir():
    return get_dirname('lsst_dc2_v6')


def get_plotdir():
    return get_dirname('plots')


def get_tabledir():
    return get_dirname('tables')


def get_texdir():
    return get_dirname('tex')


BITS = odict([
    ('GOOD',   0b000000000), # No flags
    ('DWARF4', 0b000000001), # near known dwarf (type2 >= 4)
    ('DWARF3', 0b000000010), # near known dwarf (type2 == 3)
    ('DWARF2', 0b000000100), # near known dwarf (type2 == 2)
    ('ASSOC',  0b000001000), # near object in catalogs (excluding DWARFs)
    ('STAR',   0b000010000), # near bright star
    ('EBV',    0b000100000), # E(B-V) > 0.2
    ('FOOT',   0b001000000), # outside of footprint
    ('FAIL',   0b010000000), # location of failure in ugali PS1 processing
    ('ART',    0b100000000), # location of artifact in PS1 footprint
])
BADBITS = (BITS['ASSOC'] | BITS['STAR'] | BITS['EBV'] | BITS['FOOT'] | \
           BITS['FAIL'] | BITS['ART'])
BADBITS2 = BADBITS | BITS['DWARF2']

# Satellite classification
TYPE = odict([
    (0, 'not included'),
    (1, 'not real'),
    (2, 'likely cluster'),
    (3, 'likely dwarf'),
    (4, 'kinematic but strange'),
    (5, 'kinematic dwarf')
])

DATADIR  = get_datadir()
DES_SIMS = join(DATADIR,'sim_population_v12.2.0_des_results_0000001-0100000.fits')
PS1_SIMS = join(DATADIR,'sim_population_v13.0.1_ps1_results_0000001-1000000.fits')
# wondering what to do about lsst_sims...
LSST_SIMS = '/home/kb/software/simple_adl/notebooks/simsv6_results.csv'
DWARFS   = join(DATADIR,'dwarf_literature_values.txt')
EBV      = join(DATADIR,'ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz')

MASK_VERSION  = '6.1'
DES_MASK = join(DATADIR,'healpix_mask_des_v{}.fits.gz'.format(MASK_VERSION))
PS1_MASK = join(DATADIR,'healpix_mask_ps1_v{}.fits.gz'.format(MASK_VERSION))
LSST_MASK = '/home/kb/lsst_wfd_toy.fits'

DES_CLASS = join(DATADIR,'classifier_des.yaml')
PS1_CLASS = join(DATADIR,'classifier_ps1.yaml')

DES_DENSITY = join(DATADIR,'des_stellar_density_r22_equ_n128_v0.fits.gz')
#DES_DENSITY = join(DATADIR,'gaia_stellar_density_G21_equ_n128_v0.fits')
PS1_DENSITY = join(DATADIR,'ps1_stellar_density_r22_equ_n128_v0.fits.gz')
#PS1_DENSITY = join(DATADIR,'gaia_stellar_density_G21_equ_n128_v0.fits')

LVDB = join(DATADIR,'lvdb_v3.1.fits')
MASTER = join(DATADIR, 'mw_sats_master.csv')


def load_des_sims(filename=None):
    return load_sims('des',filename=filename)


def load_ps1_sims(filename=None):
    return load_sims('ps1',filename=filename)


def load_sims(survey,filename=None, lsst_version=6):
    if filename is not None: pass
    elif survey['name'] == 'des': filename = DES_SIMS
    elif survey['name'] == 'ps1': filename = PS1_SIMS
    elif survey['name'] == 'lsst':
        filename = LSST_SIMS
    else: raise Exception("Unrecognized survey: %s"%survey)

    print("Loading %s..."%filename)
    if '.csv' in filename:
        sims = read_csv(filename)
    else:
        sims = fitsio.read(filename)

    # Rename the stellar density column
    # names = list(sims.dtype.names)
    # if 'STELLAR_DENSITY' not in names:
    #     names[names.index('DENSITY')] = 'STELLAR_DENSITY'
    #     sims.dtype.names = names
    return sims


def load_hpxmap(filename,nest=False):
    """ Load the mask file """
    print("Loading %s..."%filename)
    if '.hs' in filename:
        return HealSparseMap.read(filename)
    else:
        return read_map(filename,nest)


#load_mask = load_hpxmap


def load_des_mask(nest=False):
    #return fitsio.read(DES_MASK)['T'].ravel()
    return load_hpxmap(DES_MASK,nest)


def load_ps1_mask(nest=False):
    #return fitsio.read(PS1_MASK)['T'].ravel()
    return load_hpxmap(PS1_MASK,nest)

def load_lsst_mask(nest=False):
    return load_hpxmap(LSST_MASK, nest)

def load_mask(survey,nest=False):
    if survey['name'] == 'ps1': return load_ps1_mask()
    elif survey['name'] == 'des': return load_des_mask()
    elif survey['name'] == 'lsst': return load_lsst_mask()
    else: raise Exception("Unrecognized survey: %s"%survey)


### Stellar density

def load_density(survey,nest=False):
    """ Load stellar density [stars/arcmin^2] """
    if 'ps1' in survey: filename = PS1_DENSITY
    elif 'des' in survey: filename = DES_DENSITY
    else: raise Exception("Unrecognized survey: %s"%survey)

    print("Loading %s..."%filename)
    return hp.read_map(filename)


def load_des_density():
    return load_density('des')


def load_ps1_density():
    return load_density('ps1')


### fracdet
def load_des_fracdet(nest=False):
    #https://cdcvs.fnal.gov/redmine/projects/des-y3/wiki/Systematic_maps
    #filename = get_datadir()+'/y3a2_footprint_griz_1exp_v2.0.fits.gz'
    filename = get_datadir()+'/y3a2_griz_o.4096_t.32768_coverfoot_EQU.fits.gz'
    return load_hpxmap(filename,nest)


def load_des_maps():
    mask = load_des_mask().astype('i2')
    frac = load_des_fracdet()

    foot = (mask & BITS['FOOT'] == 0)
    good = (mask & BADBITS == 0)

    return mask, frac, foot, good


def load_ps1_maps():
    mask = load_ps1_mask().astype('i2')

    foot = (mask & BITS['FOOT'] == 0)
    good = (mask & BADBITS == 0)

    frac = foot

    return mask, frac, foot, good


def load_maps(survey):
    if   survey == 'ps1': return load_ps1_maps()
    elif survey == 'des': return load_des_maps()
    else: raise Exception("Unrecognized survey: %s"%survey)


def load_dwarfs(upper=None):
    dwarfs = np.genfromtxt(DWARFS,names=True,dtype=None)
    if upper == True:
        dwarfs.dtype.names = list(map(str.upper,dwarfs.dtype.names))
    elif upper == False:
        dwarfs.dtype.names = list(map(str.lower,dwarfs.dtype.names))
    return dwarfs


def load_ebv():
    return read_map(EBV,nest=True)


def load_lvdb():
    return fitsio.read(LVDB)


def load_master_csv():
    return np.recfromcsv(MASTER)


def load_objects_from_master(survey):
    master = load_master_csv()
    cut = []
    for sat in master:
        sel = sat['type2'] > 0
        sel &= survey in sat['survey']
        cut.append(sel)
    cut = np.array(cut)
    return master[cut]


def load_dwarf_table_csv(filename,upper=None):
    """ Parse physical parameters of dwarfs out of csv tables. """
    df = pd.read_csv(filename)
    df.rename(columns={'TS':'ts','SIG':'sig','mod_actual':'mod','ah':'a_h',
                       'ellipticity':'ell','r12':'r_physical',
                       'M_V':'abs_mag','m_v':'abs_mag',
                       'abbreviation':'abbr',
                       'distance_kpc':'distance'},
              inplace=True)
    df['r_physical'] /= 1e3 # pc -> kpc
    return df.to_records(index=False)


def load_dwarf_table_tex(filename,upper=None):
    """ Parse physical parameters of dwarfs out of tex tables. """

    with open(filename,'r') as f:
        dtypes = [('name','S30'),('ts',float),('sig',float),('pdet',float),
                  ('ra',float),('dec',float),
                  ('mod',float),('a_h',float),('ell',float),
                  ('distance',float),('abs_mag',float),('r_physical',float)]

        if upper == True:
            dtypes = [(n.upper(),t) for n,t in dtypes]
        elif upper == False:
            dtypes = [(n.lower(),t) for n,t in dtypes]

        names = [k[0] for k in dtypes]
        values = []
        for line in f.read().splitlines():
            if not line: continue
            if line.startswith(('\\',' \\','%','{')): continue
            line = line.strip('\\')
            line = line.replace(r'\ldots','nan')
            items = line.split('&')
            #print(items)

            name = items[0].split('\\')[0].strip()
            ts   = float(items[1])**2
            sig  = float(items[2])
            pdet = float(items[3])
            ra = float(items[4]) # deg
            dec = float(items[5]) # deg
            mod = float(items[6])
            a_h = float(items[7])
            ell = items[8].strip().replace('...','0')
            ell = float(ell) if not ell.startswith('$<') else 0
            distance = float(items[9]) # kpc
            abs_mag = float(items[10]) # mag
            r_physical = float(items[11])/1000. # kpc
            values.append((name,ts,sig,pdet,ra,dec,mod,a_h,ell,distance,abs_mag,r_physical))

        return np.rec.fromrecords(values,dtype=dtypes)


def load_des_dwarf_table(tex=False):
    """ Load list of dwarfs in DES """
    if tex: return load_dwarf_table_tex(get_texdir()+'/table_known_des.tex')
    return load_dwarf_table_csv(get_datadir()+'/recovered_satellites_des.csv')


def load_ps1_dwarf_table(tex=False):
    """ Load list of dwarfs in PS1 """
    if tex: return load_dwarf_table_tex(get_texdir()+'/table_known_ps1.tex')
    return load_dwarf_table_csv(get_datadir()+'/recovered_satellites_ps1.csv')


def load_all_dwarfs():
    """ Load list of all Milky Way dwarfs """
    #df = pd.read_csv(get_datadir()+'/dsphs_v2.csv',comment='#')
    #data = df.to_records(index=False)
    #data.dtype.names = [str(n.lower()) for n in data.dtype.names]
    #return data

    data = load_dwarf_table_csv(MASTER)
    # select only dwarfs
    return data[(data['type2'] > 2) & (data['distance'] < 400)]


def load_datasets(survey='ps1'):
    from selection_function import SurveySelectionFunction

    if survey == 'ps1':
        dwarfs = load_ps1_dwarf_table()
        ssf = SurveySelectionFunction(PS1_CLASS)
    elif survey == 'des':
        dwarfs = load_des_dwarf_table()
        ssf = SurveySelectionFunction(DES_CLASS)
    else:
        raise Exception()

    return ssf,dwarfs


def select_mask(ra,dec,mask):
    """ Select satellites in good regions of the mask

    Parameters
    ----------
    ra   : right ascension (deg)
    dec  : declination (deg)
    mask : healpix mask

    Returns
    -------
    sel  : boolean selection array
    """
    # Mask selection
    pix = hp.ang2pix(hp.get_nside(mask), ra, dec, lonlat=True)
    mp = np.array(mask[pix], dtype=int)
    return ( (mp & BADBITS) == 0 )


def select(data,mask):
    """ Select a good subset of simulations.

    Parameters
    ----------
    data : the simulated data
    mask : the healpix mask

    Returns
    -------
    sel  : boolean selection array
    """

    # try:
    #     data['TS'][data['DIFFICULTY'] == 1] = np.nan
    #     data['TS'][data['DIFFICULTY'] == 2] = 99**2
    # except ValueError: pass
    try:
        data['SIG'][data['DIFFICULTY'] == 1] = np.nan
        data['SIG'][data['DIFFICULTY'] == 2] = 99
    except ValueError: pass

    # Mask selection
    sel = select_mask(data['RA'],data['DEC'],mask)
    # Other possible selections
    # sel &= ((data['FLAG']==0) | (data['FLAG']==8))
    # #sel &= (np.abs(data['GLAT']) > 15)
    # sel &= (data['FRACDET_HALF'] > 0.5)

    # # Other failed runs
    # sel &= (np.nan_to_num(data['SIG']) >= 0)
    # sel &= (np.nan_to_num(data['TS'])  >= 0)

    # Select everything
    #sel = np.ones(len(data),dtype=bool)

    return sel


def detect_simple(data):
    SIG = 6.0
    try:
        return (data['SIG'] > SIG)
    except ValueError:
        return (data['sig'] > SIG)


def detect_ugali(data):
    TS = 36.0
    try:
        return (data['TS'] > TS)
    except ValueError:
        return (data['ts'] > TS)


def detect(data):
    return detect_simple(data) #& detect_ugali(data)

