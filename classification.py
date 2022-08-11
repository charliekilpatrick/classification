from astropy.table import Table, Column, vstack, hstack
import os
import glob
import re
import sys
import dustmaps.sfd
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

# From Table 6 of Schlafly & Finkbeiner for R_filter values given Rv=3.1
sf11rv = {
    'WFPC2_F218W':7.760,
    'WFPC2_F300W':4.902,
    'WFPC2_F336W':4.424,
    'WFPC2_F487N':3.175,
    'WFPC2_F450W':3.410,
    'WFPC2_F547M':2.735,
    'WFPC2_F555W':2.755,
    'WFPC2_F606W':2.415,
    'WFPC2_F675W':2.151,
    'WFPC2_F702W':1.948,
    'WFPC2_F814W':1.549,
    'WFC3_IR_F105W':0.969,
    'WFC3_IR_F110W':0.881,
    'WFC3_IR_F125W':0.726,
    'WFC3_IR_F140W':0.613,
    'WFC3_IR_F160W':0.512,
    'WFC3_UVIS_F200LP':2.958,
    'WFC3_UVIS_F218W':7.760,
    'WFC3_UVIS_F225W':6.989,
    'WFC3_UVIS_F275W':5.487,
    'WFC3_UVIS_F300X':5.228,
    'WFC3_UVIS_F336W':4.453,
    'WFC3_UVIS_F350LP':2.624,
    'WFC3_UVIS_F390W':3.896,
    'WFC3_UVIS_F438W':3.623,
    'WFC3_UVIS_F475W':3.248,
    'WFC3_UVIS_F475X':3.116,
    'WFC3_UVIS_F547M':2.756,
    'WFC3_UVIS_F555W':2.855,
    'WFC3_UVIS_F600LP':1.781,
    'WFC3_UVIS_F606W':2.488,
    'WFC3_UVIS_F625W':2.259,
    'WFC3_UVIS_F775W':1.643,
    'WFC3_UVIS_F814W':1.536,
    'WFC3_UVIS_F850LP':1.208,
    'ACS_CLEAR':2.436,
    'ACS_WFC_F330W':4.472,
    'ACS_WFC_F435W':3.610,
    'ACS_WFC_F475W':3.268,
    'ACS_WFC_F550M':2.620,
    'ACS_WFC_F555W':2.792,
    'ACS_WFC_F606W':2.471,
    'ACS_WFC_F625W':2.219,
    'ACS_WFC_F658N':2.224,
    'ACS_WFC_F775W':1.629,
    'ACS_WFC_F814W':1.526,
    'ACS_WFC_F850LP':1.243,
}

for pole in ['ngp', 'sgp']:
    datadir=dustmaps.sfd.data_dir()
    file = os.path.join(datadir, 'sfd',
        'SFD_dust_4096_{}.fits'.format(pole))
    if not os.path.exists(file):
        dustmaps.sfd.fetch()
        break

sfd=dustmaps.sfd.SFDQuery()

def generate_inst_filt_colname(inst_filt, colnames, keytype='mag'):

    inst_filt = inst_filt.upper()

    matches = [i for i in colnames if inst_filt in i]
    ver = len(matches)
    if keytype=='magerr':
        ver = ver-1

    colname = inst_filt+'_'+str(ver).zfill(2)

    if keytype=='magerr':
        colname = colname+'_ERR'

    return(colname)

def parse_inst_filt(var):

    if 'wfc3' in var.lower():
        inst = 'WFC3'
    elif 'acs' in var.lower():
        inst = 'ACS'
    elif 'wfpc2' in var.lower():
        inst = 'WFPC2'
    else:
        raise Exception(f'COULD NOT PARSE INST FROM {var}')

    reg='\.(f.*?)\.'
    m = re.findall(reg, var)
    if len(m)!=1:
        raise Exception(f'COULD NOT PARSE FILT FROM {var}')
    else:
        filt = m[0]

    val = inst+'_'+m[0].upper()

    return(val)


# For an input list of instruments and filters, load all relevant mist data 
# tables into a single table object
def load_mist_models(mistdir, inst_filts):

    mist_dir = os.path.join(mistdir, 'FEH_0000')
    inst_dirs = {'ACS':'ACS_WFC','WFPC2':'WFPC2','WFC3':'WFC3'}

    all_tables = {}

    inst_filt_data = {'ACS':[], 'WFPC2':[],'WFC3':[]}
    for inst_filt in inst_filts:
        inst, filt = inst_filt.split('_')
        inst_filt_data[inst].append(filt)

    for inst in ['ACS','WFPC2','WFC3']:

        id_base = inst_dirs[inst]
        subdir = os.path.join(mist_dir, id_base)

        insttable = None

        if all([inst not in i for i in inst_filts]):
            continue

        for cmd in glob.glob(os.path.join(subdir, '*.track.eep.cmd')):

            mass = os.path.basename(cmd).replace('M.track.eep.cmd','')
            mass = float(mass) * 1.0e-4

            table = Table.read(cmd, format='ascii', header_start=14)
            table.add_column(Column([mass]*len(table), name='mass'))

            use_cols = ['star_age','log_Teff','log_L', 'mass', 'phase']
            use_cols.extend(inst_filts)

            for col in table.keys():
                newcol=''
                if 'ACS_WFC' in col:
                    newcol = col.replace('ACS_WFC','ACS')
                elif 'WFC3_UVIS' in col:
                    newcol = col.replace('WFC3_UVIS','WFC3')
                elif 'WFC3_IR' in col:
                    newcol = col.replace('WFC3_IR','WFC3')
                if newcol not in use_cols and col not in use_cols:
                    table.remove_column(col)
                elif col!=newcol and newcol:
                    table.rename_column(col, newcol)

            if insttable is None:
                insttable = table
            else:
                insttable = vstack([insttable, table])

        insttable.sort('mass','star_age')

        all_tables[inst]=insttable

    return(all_tables)

def get_a_lambda(inst_filt, ebv):

    inst, filt = inst_filt.split('_')

    use_col = None
    for key in sf11rv.keys():
        if inst in key and filt in key:
            use_col = key

    if use_col is None:
        raise Exception(f'COULD NOT GET R_V for {inst_filt}')

    return(sf11rv[use_col])

def get_zpt(zpt_data, name):

    for row in zpt_data:
        if row['imagename'].replace('.fits','').replace('.drz','') in name:
            return(row['zpt'])

    raise Exception(f'COULD NOT GET ZPT FOR {name}')

def load_and_parse_table(table_name, dist=None, ebv=None, zpt_data=None):

    table = Table.read(table_name, format='ascii')

    if dist:
        dm = 5 * np.log10(float(dist)*1.0e6) - 5
    else:
        dm = 0.

    colnames = []
    for key in table.keys():
        if key=='MAG_AUTO':
            inst_filt=parse_inst_filt(table_name)
            name=generate_inst_filt_colname(inst_filt, colnames, keytype='mag')

            if zpt_data:
                zpt = get_zpt(zpt_data, table_name)
            else:
                zpt = 0.

            if ebv:
                a_lambda = get_a_lambda(inst_filt, ebv)
            else:
                a_lambda = 0.

            table.rename_column(key, name)

            table[name] = table[name] + zpt - dm - a_lambda

            colnames.append(name)

        elif key=='MAGERR_AUTO':
            name=parse_inst_filt(table_name)
            name=generate_inst_filt_colname(name, colnames, keytype='magerr')

            table.rename_column(key, name)

        elif key.endswith('mag'):
            name=parse_inst_filt(key)
            name=generate_inst_filt_colname(name, colnames, keytype='mag')

            if zpt_data:
                zpt = get_zpt(zpt_data, table_name)
            else:
                zpt = 0.

            if ebv:
                a_lambda = get_a_lambda(inst_filt, ebv)
            else:
                a_lambda = 0.

            table.rename_column(key, name)

            table[name] = table[name] + zpt - dm - a_lambda

            colnames.append(name)

        elif key.endswith('magerr'):
            name=parse_inst_filt(key)
            name=generate_inst_filt_colname(name, colnames, keytype='magerr')

            table.rename_column(key, name)

        elif key=='ALPHA_J2000':
            table.rename_column(key, 'RA')
        elif key=='DELTA_J2000':
            table.rename_column(key, 'DEC')
        else:
            table.remove_column(key)

    return(table)

def get_unique_inst_filt(table):

    inst_filts = []
    for key in table.keys():

        if key.endswith('ERR') or key=='RA' or key=='DEC':
            continue

        parts = key.split('_')

        inst_filt = parts[0]+'_'+parts[1]

        if inst_filt not in inst_filts:
            inst_filts.append(inst_filt)

    return(inst_filts)

def do_classification(filename, distance, ra, dec):

    coord = SkyCoord(ra, dec, unit=(u.hour, u.deg))
    ebv = sfd(coord) * 0.86

    zpt_data = Table.read('zeropoints.dat', format='ascii')

    table = load_and_parse_table(os.path.join('data', filename), dist=distance,
        ebv=ebv, zpt_data=zpt_data)
    inst_filts = get_unique_inst_filt(table)

    mist_dir = '/Users/ckilpatrick/scripts/python/progenitors/sed/data/mist'

    models = load_mist_models(mist_dir, inst_filts)

    classification_table = Table([[0.],[0.],[0.],[0.]], names=('mass','phase',
        'ra','dec'))

    for row in table:

        total_chi2 = np.zeros(70334)
        for key in table.keys():
            if (('WFC3' in key or 'ACS' in key or 'WFPC2' in key) and
                not key.endswith('ERR')):

                inst = key.split('_')[0]
                filt = key.split('_')[1]

                mag = row[key]
                if np.isnan(float(mag)):
                    continue
                check_key = inst+'_'+filt
                if check_key not in list(models[inst].keys()):
                    continue

                magerr = row[key+'_ERR']

                chi2 = (models[inst][check_key].data-mag)**2 / (magerr**2)

                total_chi2 += chi2

        idx = np.argmin(total_chi2)

        if 'WFC3' in models.keys():
            data=models['WFC3'][idx]
        elif 'ACS' in models.keys():
            data=models['ACS'][idx]
        elif 'WFPC2' in models.keys():
            data=models['WFPC2'][idx]

        classification_table.add_row([data['mass'],data['phase'],row['RA'],row['DEC']])

    mask = (classification_table['mass'] > 18.0) & (classification_table['phase']==0)
    if len(classification_table[mask])>0:
        coords = SkyCoord(classification_table[mask]['ra'],
            classification_table[mask]['dec'], unit=(u.deg, u.deg))
        seps = coord.separation(coords)

        idx = np.argmin(seps)
        sep_physical = seps[idx].radian * distance * 1.0e6
        print(f'Separation is {sep_physical}')
        return(sep_physical)
    else:
        print('ERROR: no stars classified as main sequence O-type')
        return(None)
