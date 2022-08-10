from astropy.table import Table
import os
import glob
import re

def parse_inst_filt(var):

    if 'wfc3' in var.lower():
        isnt = 'WFC3'
    elif 'acs' in var.lower():
        inst = 'ACS'
    elif 'wfpc2' in var.lower():
        inst = 'WFPC2'
    else:
        raise Exception(f'COULD NOT PARSE INST FROM {var}')

    reg='\.(f.+)\.'
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

    all_tables = None

    # Get 
    inst_filt_data = {'ACS':[], 'WFPC2':[],'WFC3':[]}
    for inst_filt in inst_filts:
        inst, filt = inst_filt.split('_')
        inst_filt_data[inst].append(filt)

    for inst in ['ACS','WFPC2','WFC3']:

        id_base = inst_dirs[inst]
        subdir = os.path.join(mist_dir, id_base)

        insttable = None

        for cmd in glob.glob(os.path.join(subdir, '*.track.eep.cmd')):

            mass = os.path.basename(cmd).replace('M.track.eep.cmd','')
            mass = float(mass) * 1.0e-4

            table = Table.read(cmd, format='ascii')
            table.add_column(Column([mass]*len(table), name='mass'))

            use_cols = ['star_age','log_Teff','log_L']

            for filt in inst_filt_data[inst]:
                colname = inst_filt_data[inst]+'_'+filt
                use_cols.append(colname)

            for col in table.keys():
                if col not in use_cols:
                    table.remove_column(colname)

            if insttable is None:
                insttable = table
            else:
                insttable = vstack([insttable, table])

        if all_tables is None:
            all_tables = insttable
        else:
            all_tables = hstack([all_tables, insttable])


    return(all_tables)


def load_and_parse_table(table_name):

    table = Table.read(table_name, format='ascii')

    inst_filts = []
    for key in table.keys():
        if key=='MAG_AUTO':
            name=parse_inst_filt(table_name)

            colname = inst_
    
