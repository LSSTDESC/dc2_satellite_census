from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import os

'''
Clone the throughputs repository from https://github.com/lsst/throughputs

Checkout tagged versions (e.g., git checkout DC2production or git checkout 1.1)

Copy the baseline throughputs to a local directory (e.g. 'throughputs_1p1/baseline')

'''

dir1p1 = 'throughputs_1p1/baseline'
# dirdc2 = 'throughputs_1p4/baseline'
dirdc2 = 'throughputs_dc2_1p4/baseline'
#dirdc2 = 'throughputs_1p7/baseline'

# From https://github.com/DarkEnergySurvey/descolors
from collections import OrderedDict as odict

BAND_COLORS = odict([
    ('u','#56b4e9'),
    ('g','#008060'),
    ('r','#ff4000'),
    ('i','#850000'),
    ('z','#6600cc'),
    ('y','#000000'),
])


# Read in the total ('final') throughput curves in each filter, and the atmosphere curve.
# filterdir = os.getenv('LSST_THROUGHPUTS_BASELINE')
lsst = {}
filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
# filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'orange', 'z':'r', 'y':'m'}
#for f in filterlist:
#    lsst[f] = Bandpass()
#    lsst[f].readThroughput(os.path.join(filterdir, 'total_'+f+'.dat'))
#atmos = Bandpass()
#atmos.readThroughput(os.path.join(filterdir, 'atmos_std.dat'))

#%%
# Read in each component separately.
allcommon_components = ['detector', 'lens1', 'lens2', 'lens3', 'm1', 'm2', 'm3', 'atmos_std']

allcommon_1p1 = {}
for c in allcommon_components:
    allcommon_1p1[c] = ascii.read(os.path.join(dir1p1, c+'.dat'))
    # allcommon[c].readThroughput(os.path.join(filterdir, c +'.dat'))


# We probably won't want to show each component separately, as the plot will be quite messy.
#  So let's batch up "mirrors" and "lenses" and the other components we will want to plot.

common_1p1 = {}
common_1p1['detector'] = allcommon_1p1['detector']
common_1p1['atmosphere'] = allcommon_1p1['atmos_std']
common_1p1['mirrors'] = allcommon_1p1['m1'].copy()
common_1p1['mirrors']['col2'] = allcommon_1p1['m1']['col2'] * allcommon_1p1['m2']['col2'] * allcommon_1p1['m3']['col2']
common_1p1['lenses'] = allcommon_1p1['m1'].copy()
common_1p1['lenses']['col2'] = allcommon_1p1['lens1']['col2'] * allcommon_1p1['lens2']['col2'] * allcommon_1p1['lens3']['col2']
#common['mirrors'] = copy.deepcopy(allcommon['m1'])
#common['mirrors'].sb = allcommon['m1'].sb * allcommon['m2'].sb * allcommon['m3'].sb
#common['lenses'] = copy.deepcopy(allcommon['lens1'])
#common['lenses'].sb = allcommon['lens1'].sb * allcommon['lens2'].sb * allcommon['lens3'].sb

# Note that we could have also combined these throughput curves using 'Bandpass.readThroughputList'
lsst_filters_1p1 = {}
total_1p1 = {}
for f in filterlist:
    lsst_filters_1p1[f] = ascii.read(os.path.join(dir1p1, 'filter_'+f+'.dat'))
    total_1p1[f] = ascii.read(os.path.join(dir1p1, 'total_'+f+'.dat'))
    # lsst_filters[f].readThroughput(os.path.join(filterdir, 'filter_'+f+'.dat'))

allcommon_dc2 = {}
for c in allcommon_components:
    allcommon_dc2[c] = ascii.read(os.path.join(dirdc2, c+'.dat'))

# We probably won't want to show each component separately, as the plot will be quite messy.
#  So let's batch up "mirrors" and "lenses" and the other components we will want to plot.

common_dc2 = {}
common_dc2['detector'] = allcommon_dc2['detector']
common_dc2['atmosphere'] = allcommon_dc2['atmos_std']
common_dc2['mirrors'] = allcommon_dc2['m1'].copy()
common_dc2['mirrors']['col2'] = allcommon_dc2['m1']['col2'] * allcommon_dc2['m2']['col2'] * allcommon_dc2['m3']['col2']
common_dc2['lenses'] = allcommon_dc2['lens1'].copy()
common_dc2['lenses']['col2'] = allcommon_dc2['lens1']['col2'] * allcommon_dc2['lens2']['col2'] * allcommon_dc2['lens3']['col2']

# Note that we could have also combined these throughput curves using 'Bandpass.readThroughputList'
lsst_filters_dc2 = {}
total_dc2 = {}
for f in filterlist:
    lsst_filters_dc2[f] = ascii.read(os.path.join(dirdc2, 'filter_'+f+'.dat'))
    total_dc2[f] = ascii.read(os.path.join(dirdc2, 'total_'+f+'.dat'))


params = {
   'axes.labelsize': 24,
   'font.size': 24,
   'legend.fontsize': 20,
#   'xtick.labelsize': 16,
   'xtick.major.width': 3,
   'xtick.minor.width': 2,
   'xtick.major.size': 8,
   'xtick.minor.size': 5,
   'xtick.direction': 'in',
   'xtick.top': True,
   'lines.linewidth':3,
   'axes.linewidth':3,
   'axes.labelweight':3,
   'axes.titleweight':3,
   'ytick.major.width':3,
   'ytick.minor.width':2,
   'ytick.major.size': 8,
   'ytick.minor.size': 5,
   'ytick.direction': 'in',
   'ytick.right': True,
#   'ytick.labelsize': 20,
   'text.usetex': True,
   'text.latex.preamble': r'\boldmath',
   'figure.figsize': [9, 7]
   }

for comp in allcommon_components:

    plt.rcParams.update(params)
    fig=plt.figure()

    plt.plot(allcommon_1p1[comp]['col1'], allcommon_1p1[comp]['col2'], 'k', label='v1.1')
    plt.plot(allcommon_dc2[comp]['col1'], allcommon_dc2[comp]['col2'], 'k', linestyle='--', label='v1.4')
    plt.legend()
    plt.minorticks_on()
    if comp in 'atmos_std':
        comptitle = 'atmos'
    else:
        comptitle = comp
    plt.title(comptitle)
    plt.xlabel('wavelength (nm)')
    plt.ylabel('throughput')
    plt.savefig(comp+'_compare.png')
    plt.show()

for comp in common_dc2.keys():

    plt.rcParams.update(params)
    fig=plt.figure()

    plt.plot(common_1p1[comp]['col1'], common_1p1[comp]['col2'], 'k', label='v1.1')
    plt.plot(common_dc2[comp]['col1'], common_dc2[comp]['col2'], 'k', linestyle='--', label='v1.4')
    plt.legend()
    plt.minorticks_on()
    if comp in 'atmos_std':
        comptitle = 'atmos'
    else:
        comptitle = comp
    plt.title(comptitle)
    plt.xlabel('wavelength (nm)')
    plt.ylabel('throughput')
    plt.savefig(comptitle+'_compare.png')
    plt.show()


plt.rcParams.update(params)
fig=plt.figure()

for filt in filterlist:

    plt.plot(lsst_filters_1p1[filt]['col1'], lsst_filters_1p1[filt]['col2'], color=BAND_COLORS[filt], label='__no label__')
    plt.plot(lsst_filters_dc2[filt]['col1'], lsst_filters_dc2[filt]['col2'], color=BAND_COLORS[filt], linestyle='--', label='__no label__')

plt.plot(lsst_filters_1p1[filt]['col1'], lsst_filters_1p1[filt]['col2'], color=BAND_COLORS[filt], label='v1.1')
plt.plot(lsst_filters_dc2[filt]['col1'], lsst_filters_dc2[filt]['col2'], color=BAND_COLORS[filt], linestyle='--', label='v1.4')
plt.legend()
plt.minorticks_on()
plt.title('filters')
plt.xlabel('wavelength (nm)')
plt.ylabel('throughput')
plt.savefig('filters_compare.png')
plt.show()

plt.rcParams.update(params)
fig=plt.figure()

for filt in filterlist:

    plt.plot(total_1p1[filt]['col1'], total_1p1[filt]['col2'], color=BAND_COLORS[filt], label='__no label__')
    plt.plot(total_dc2[filt]['col1'], total_dc2[filt]['col2'], color=BAND_COLORS[filt], linestyle='--', label='__no label__')

plt.plot(total_1p1[filt]['col1'], total_1p1[filt]['col2'], color=BAND_COLORS[filt], label='v1.1')
plt.plot(total_dc2[filt]['col1'], total_dc2[filt]['col2'], color=BAND_COLORS[filt], linestyle='--', label='v1.4')
plt.legend()
plt.minorticks_on()
plt.title('total throughput')
plt.xlabel('wavelength (nm)')
plt.ylabel('throughput')
plt.savefig('total_throughput_compare.png')
plt.show()
