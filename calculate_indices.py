'''
Calculate various climate extreme indices from PRIMAVERA HighResMIP simulations
Use ICCLIM package
Write output in nice directory structure
'''


import icclim
import datetime, glob, os, subprocess
from netCDF4 import Dataset

file_format = '/group_workspaces/jasmin2/primavera5/stream1/CMIP6/HighResMIP/{}/{}/highresSST-present/{}/{}/{}/{}/*/{}_{}*_{}_{}*.nc'
file_prim_format = '/group_workspaces/jasmin2/primavera4/stream1/PRIMAVERA/HighResMIP/{}/{}/highresSST-present/{}/{}/{}/{}/*/{}_{}*_{}_{}*.nc'

file_out = '/group_workspaces/jasmin2/primavera1/model_derived_data/climate_extreme_indices/{}'
fileout_format = '/group_workspaces/jasmin2/primavera1/model_derived_data/climate_extreme_indices/{}/{}/highresSST-present/{}/{}/{}/{}/{}/{}'

years = range(1950, 2015)

mohc_LM = {'institute':'MOHC', 'model':'HadGEM3-GC31-LM', 'label':'r1i1p1f1', 'grid':'gn'}
mohc_LM_2 = {'institute':'MOHC', 'model':'HadGEM3-GC31-LM', 'label':'r1i2p1f1', 'grid':'gn'}
mohc_LM_3 = {'institute':'MOHC', 'model':'HadGEM3-GC31-LM', 'label':'r1i3p1f1', 'grid':'gn'}
mohc_HM = {'institute':'MOHC', 'model':'HadGEM3-GC31-HM', 'label':'r1i1p1f1', 'grid':'gn'}
mohc_HM_2 = {'institute':'MOHC', 'model':'HadGEM3-GC31-HM', 'label':'r1i2p1f1', 'grid':'gn'}
mohc_HM_3 = {'institute':'MOHC', 'model':'HadGEM3-GC31-HM', 'label':'r1i3p1f1', 'grid':'gn'}

ecmwf_LR = {'institute':'ECMWF', 'model':'ECMWF-IFS-LR', 'label':'r1i1p1f1', 'grid':'gr'}
ecmwf_HR = {'institute':'ECMWF', 'model':'ECMWF-IFS-HR', 'label':'r1i1p1f1', 'grid':'gr'}

ecearth_LR = {'institute':'EC-Earth-Consortium', 'model':'EC-Earth3', 'label':'r1i1p1f1', 'grid':'gr'}
ecearth_HR = {'institute':'EC-Earth-Consortium', 'model':'EC-Earth3-HR', 'label':'r1i1p1f1', 'grid':'gr'}

mpim_LR = {'institute':'MPI-M', 'model':'MPIESM-1-2-HR', 'label':'r1i1p1f1', 'grid':'gn'}
mpim_HR = {'institute':'MPI-M', 'model':'MPIESM-1-2-XR', 'label':'r1i1p1f1', 'grid':'gn'}

cmcc_LR = {'institute':'CMCC', 'model':'CMCC-CM2-HR4', 'label':'r1i1p1f1', 'grid':'gn'}
cmcc_HR = {'institute':'CMCC', 'model':'CMCC-CM2-VHR4', 'label':'r1i1p1f1', 'grid':'gn'}

# CNRM model f2 has surf temp error over sea-ice
cm61_LR = {'institute':'CNRM-CERFACS', 'model':'CNRM-CM6-1', 'label':'r21i1p1f2', 'grid':'gr'}
cm61_HR = {'institute':'CNRM-CERFACS', 'model':'CNRM-CM6-1-HR', 'label':'r2i1p1f2', 'grid':'gr'}

def define_indices():
    '''
    Define climate extreme indices with names consistent with ccclim
    '''
    global index1, index2, index3, index4, index5, index6, index7, index8, index9, index16, index17, index18, index19, index20, index21, index23, index24, index27, index31, index31_params, index32_3hr, index31_params_3hr, index32_6hr, index31_params_6hr, index33, index33_params, index17_6hr, index17_params_6hr, index18_6hr, index18_params_6hr

    index1 = {'name':'FD', 'var':'tasmin', 'period':'year', 'var_period':'day'}
    index2 = {'name':'SU', 'var':'tasmax', 'period':'year', 'var_period':'day'}
    index3 = {'name':'ID', 'var':'tasmax', 'period':'year', 'var_period':'day'}
    index4 = {'name':'TR', 'var':'tasmin', 'period':'year', 'var_period':'day'}
    index5 = {'name':'GSL', 'var':'tas', 'period':'year', 'var_period':'day'}
    index6 = {'name':'TXx', 'var':'tasmax', 'period':'month', 'var_period':'day'}
    index7 = {'name':'TNx', 'var':'tasmin', 'period':'month', 'var_period':'day'}
    index8 = {'name':'TXn', 'var':'tasmax', 'period':'month', 'var_period':'day'}
    index9 = {'name':'TNn', 'var':'tasmin', 'period':'month', 'var_period':'day'}
    index16 = {'name':'DTR', 'var':['TX','TN'], 'period':'month', 'var_period':'day'}
    index17 = {'name':'RX1day', 'var':'pr', 'period':'month', 'var_period':'day'}
    index18 = {'name':'RX5day', 'var':'pr', 'period':'month', 'var_period':'day'}
    index19 = {'name':'SDII', 'var':'pr', 'period':'year', 'var_period':'day'}
    index20 = {'name':'R10mm', 'var':'pr', 'period':'year', 'var_period':'day'}
    index21 = {'name':'R20mm', 'var':'pr', 'period':'year', 'var_period':'day'}
    index23 = {'name':'CDD', 'var':'pr', 'period':'year', 'var_period':'day'}
    index24 = {'name':'CWD', 'var':'pr', 'period':'year', 'var_period':'day'}
    index27 = {'name':'PRCPTOT', 'var':'pr', 'period':'year', 'var_period':'day'} # annual

# Michael's idea: TX3x (monthly maximum of running 3 day average of tasmax)
    index31_params = {'indice_name': 'TX3x',
                    'calc_operation': 'run_mean',
                    'extreme_mode': 'max',
                    'window_width': 3}
    index31 = {'name':'TX3x', 'var': 'tasmax', 'period':'month', 'var_period':'day', 'long_name':'Maximum 3-day running mean maximum Near-Surface Air Temperature', 'params':index31_params}

    # Xuebin's higher freq precip
    # this one for CNRM-CERFACS, EC-Earth-Consortium, MOHC
    # use 2x3hr periods
    index32_params_3hr = {'indice_name': 'RX6hour',
                    'calc_operation': 'run_sum',
                    'extreme_mode': 'max',
                    'window_width': 2}
    index32_3hr = {'name':'RX6hour', 'var': 'pr', 'period':'month', 'var_period':'3hr', 'long_name':'Maximum 6-hour running mean maximum precipitation (from 2x3hr periods)', 'params': index32_params_3hr}

    # this one for MPI-M, ECMWF, CMCC
    # use 1x6hr periods
    index32_params_6hr = {'indice_name': 'RX6hour',
                    'calc_operation': 'max',
                    'extreme_mode': 'max'}
    index32_6hr = {'name':'RX6hour', 'var': 'pr', 'period':'month', 'var_period':'Prim6hr', 'long_name':'Maximum 6-hour running mean maximum precipitation (from 1x6hr period)', 'params': index32_params_6hr}

    index33_params = {'indice_name': 'RX3hour',
                    'calc_operation': 'max',
                    'extreme_mode': 'max'}
    index33 = {'name':'RX3hour', 'var': 'pr', 'period':'month', 'var_period':'3hr', 'long_name':'Maximum 3-hour mean precipitation', 'params':index33_params}

    # RX1day for CMCC (who do not have daily precip, so use running mean over 4 x 6 hr
    index17_params_6hr = {'indice_name': 'RX1day',
                    'calc_operation': 'run_sum',
                    'extreme_mode': 'max',
                    'window_width': 4}
    index17_6hr = {'name':'RX1day', 'var': 'pr', 'period':'month', 'var_period':'Prim6hr', 'long_name':'Highest 1-day precipitation amount (derived from running mean of 4x6hr mean periods)', 'params':index17_params_6hr}

    # RX5day for CMCC (who do not have daily precip, so use running mean over 20 x 6 hr
    index18_params_6hr = {'indice_name': 'RX5day',
                    'calc_operation': 'run_sum',
                    'extreme_mode': 'max',
                    'window_width': 20}
    index18_6hr = {'name':'RX5day', 'var': 'pr', 'period':'month', 'var_period':'Prim6hr', 'long_name':'Highest 5-day precipitation amount (derived from running mean of 20x6hr mean periods)', 'params':index18_params_6hr}

    # CDD for CMCC (who do not have daily precip, only 6 hourly
    index23_params_6hr = {'indice_name': 'CDD',
                    'calc_operation': 'max_nb_consecutive_events',
                    'logical_operation': 'lt',
                    'thresh':0.0000116,
                    }
    index23_6hr = {'name':'CDD', 'var': 'pr', 'period':'year', 'var_period':'Prim6hr', 'long_name':'Maximum number of consecutive dry days (precipitation < 1 mm) (derived from running mean of 20x6hr mean periods)', 'params':index23_params_6hr}


    indices = [index6, index7, index8, index9, index17, index18, index23, index31]
#    indices = [index6, index7, index8, index9, index17_6hr, index18_6hr, index23, index31]
    #indices = [index32_3hr, index32_6hr]
    #indices = [index17_6hr, index18_6hr]
    #indices = [index17_6hr, index18_6hr, index32_3hr, index32_6hr]
    #indices = [index33, index32_3hr, index32_6hr]
    #indices = [index33]
    indices = [index17, index18]

    return indices

def metadata_from_file(fname, model_info):
    '''
    Extract metadata from file fname
    '''
    ncatted_string = []
    with Dataset(fname,'r') as x:
        for i in x.ncattrs():
            if i != 'history':
                #print i, x.getncattr(i)
                ncatted_string.append('-a '+i+',global,a,c,"'+str(x.getncattr(i))+'"')
                model_info[i] = x.getncattr(i)

        # get the calendar
        time = x.variables['time']
        calendar = time.calendar
        calendar_units = time.units
    return ncatted_string, calendar, calendar_units

def add_model_metadata(f_ref, fname):
    '''
    Add metadata from source model file to output ccclim file
    '''
    model_info = {}
    # now extract the global attributes from the netcdf file
    ncatted_string = []
    file_size = os.path.getsize(f_ref)
    print 'file_size ',f_ref, file_size
    if os.path.exists(f_ref) and file_size > 0:
        ncatted_string, calendar, calendar_units = metadata_from_file(f_ref, model_info)

        cmd = 'ncatted -h -O '+' '.join(ncatted_string)+' '+fname
        #print cmd
        subprocess.call(cmd, shell=True)

if __name__ == '__main__':

    #runs = [mohc_LM, mohc_HM, mohc_LM_2, mohc_HM_2, mohc_LM_3, mohc_HM_3 ]
    runs = [mohc_HM, mohc_HM_2, mohc_HM_3 ]
    #runs = [mpim_LR, mpim_HR]
    #runs = [ecearth_LR, ecearth_HR]
    runs = [ecearth_HR]
    #runs = [cmcc_LR, cmcc_HR]
    #runs = [cmcc_HR]
    #runs = [cm61_LR, cm61_HR]
    #runs = [mohc_LM]
    #runs = [ecmwf_LR, ecmwf_HR]
    print runs

    # for these models, index32_3hr: CNRM-CERFACS, EC-Earth-Consortium, MOHC
    # for these models, index32_6hr: MPI-M, ECMWF, CMCC
    hr3_data = ['MOHC', 'EC-Earth-Consortium', 'CNRM-CERFACS']
    hr6_data = ['MPI-M', 'CMCC', 'ECMWF']

    indices = define_indices()

    for run in runs:
        for index in indices:
            index_name = index['name']
            index_var = index['var']
            index_period = index['period']
            index_var_period = index['var_period']
            if index_name == 'RX6hour':
                if run['institute'] in hr3_data and '3hr' not in index['var_period']:
                    print '3hr data but not right model ',run['institute'], index['var_period']
                    continue
                elif run['institute'] in hr6_data and '6hr' not in index['var_period']:
                    print '6hr data but not right model ',run['institute'], index['var_period']
                    continue
                else:
                    pass
            
            for year in years:
                # CMCC does not have daily precip, so need to look elsewhere
                if 'Prim' in index_var_period:
                    fname = file_prim_format.format(run['institute'], run['model'], run['label'], index['var_period'], index_var, run['grid'], index_var, index['var_period'], run['grid'], str(year))

                else:
                    fname = file_format.format(run['institute'], run['model'], run['label'], index['var_period'], index_var, run['grid'], index_var, index['var_period'], run['grid'], str(year))
                files = sorted(glob.glob(fname))

                if len(files) == 0:
                    print 'searching, no files ',fname, run['model'], index_var
                    continue
                print files[0]
                fname_out = os.path.basename(files[0])
                if len(files) > 1:
                    fname_out_last = os.path.basename(files[-1])
                    first_period = fname_out.split('-')[-1][:-3]
                    last_period = fname_out_last.split('-')[-1][:-3]
                    fname_out = fname_out.replace(first_period, last_period)
                    
                fname_out = fname_out.replace(index['var_period'], index_period)
                fname_out = fname_out.replace(index_var, index_name, 1)
                f_out = fileout_format.format(run['institute'], run['model'], run['label'], index_period, index_name, run['grid'], 'latest', fname_out)

                # test if directory path exists, if not create
                if not os.path.exists(os.path.dirname(f_out)):
                    os.makedirs(os.path.dirname(f_out))

                # these are user-defined indices
                if 'params' in index.keys():
                    icclim.indice(user_indice=index['params'], in_files=files, var_name=index_var, slice_mode = index_period, out_file=f_out)
                # need to edit the long_name to make it consistent with this user variable
                    cmd = 'ncatted -O -h -a long_name,'+index_name+',o,c,"'+index['long_name']+'" '+f_out
                    os.system(cmd)
                    if index_name == 'CDD' and run['institute'] == 'CMCC':
                        # need to scale 6 hour period back to days
                        with Dataset(f_out,'r+') as fh:
                            cdd = fh.variables[index_name]
                            cdd_rescale = cdd[:]
                            cdd_rescale  = cdd[:] // 4
                            cdd[:] = cdd_rescale
		    elif run['institute'] == 'CMCC' and (index_name == 'RX1day' or index_name == 'RX5day'):
                        # need to scale 6 hour period back to days (devide by 4, seems ICCLIM assumes day period)
                        with Dataset(f_out,'r+') as fh:
                            rx = fh.variables[index_name]
                            rx_rescale = rx[:]
                            rx_rescale  = rx[:] / 4.0
                            rx[:] = rx_rescale
			# change units to mm (the user defined index does not do this)
			cmd = 'ncatted -O -a units,'+index_name+',m,c,"mm" '+f_out
			os.system(cmd)
                    elif index_name == 'RX6hour':
                        if run['institute'] in hr3_data:
                            with Dataset(f_out,'r+') as fh:
                                rx = fh.variables[index_name]
                                rx_rescale = rx[:]
                                rx_rescale  = rx[:] / 8.0
                                rx[:] = rx_rescale
                        elif run['institute'] in hr6_data:
                            with Dataset(f_out,'r+') as fh:
                                rx = fh.variables[index_name]
                                rx_rescale = rx[:]
                                rx_rescale  = rx[:] / 4.0
                                rx[:] = rx_rescale
                        cmd = 'ncatted -O -a units,'+index_name+',m,c,"mm" '+f_out
                        os.system(cmd)
                    elif index_name == 'RX3hour':
                        with Dataset(f_out,'r+') as fh:
                            rx = fh.variables[index_name]
                            rx_rescale = rx[:]
                            rx_rescale  = rx[:] / 8.0
                            rx[:] = rx_rescale

			# change units to mm (the user defined index does not do this)
                        cmd = 'ncatted -O -a units,'+index_name+',m,c,"mm" '+f_out
                        os.system(cmd)

                else:
                    # these are standard indices
                    icclim.indice(indice_name=index_name, in_files=files, var_name=index_var, slice_mode = index_period, out_file=f_out)

                add_model_metadata(files[0], f_out)
                cmd = 'nccopy -k 4 -d 3 '+f_out+' '+f_out+'.tmp.nc4'
                os.system(cmd)
                os.rename(f_out+'.tmp.nc4', f_out)
