'''
Calculate various climate extreme indices from PRIMAVERA HighResMIP simulations
Use ICCLIM package
Write output in nice directory structure
'''


import icclim
import datetime, glob, os, subprocess
from netCDF4 import Dataset

file_format = '/group_workspaces/jasmin2/primavera5/stream1/CMIP6/HighResMIP/{}/{}/highresSST-present/{}/{}/{}/{}/*/{}_{}*_{}_{}*.nc'

file_out = '/group_workspaces/jasmin2/primavera1/model_derived_data/climate_extreme_indices/{}'
fileout_format = '/group_workspaces/jasmin2/primavera1/model_derived_data/climate_extreme_indices/{}/{}/highresSST-present/{}/{}/{}/{}/{}/{}'

years = range(1950, 2015)

mohc_LM = {'institute':'MOHC', 'model':'HadGEM3-GC31-LM', 'label':'r1i3p1f1', 'grid':'gn'}
mohc_HM = {'institute':'MOHC', 'model':'HadGEM3-GC31-HM', 'label':'r1i3p1f1', 'grid':'gn'}
cm61_LR = {'institute':'CNRM-CERFACS', 'model':'CNRM-CM6-1', 'label':'r21i1p1f2', 'grid':'gr'}
ecmwf_LR = {'institute':'ECMWF', 'model':'ECMWF-IFS-LR', 'label':'r1i1p1f1', 'grid':'gr'}

def define_indices():
    '''
    Define climate extreme indices with names consistent with ccclim
    '''
    global index1, index2, index3, index4, index5, index6, index7, index8, index9, index16, index17, index18, index19, index20, index21, index23, index24, index27, index31, index31_params, index32_3hr, index31_params_3hr, index32_6hr, index31_params_6hr
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
    index31 = {'name':'TX3x', 'var': 'tasmax', 'period':'month', 'var_period':'day', 'long_name':'Maximum 3-day running mean maximum Near-Surface Air Temperature'}

    # this one for CNRM-CERFACS, EC-Earth-Consortium, MOHC
    # use 2x3hr periods
    index32_params_3hr = {'indice_name': 'RX6hour',
                    'calc_operation': 'run_mean',
                    'extreme_mode': 'max',
                    'window_width': 2}
    index32_3hr = {'name':'RX6hour', 'var': 'pr', 'period':'month', 'var_period':'3hr', 'long_name':'Maximum 6-hour running mean maximum precipitation (from 2x3hr periods)'}

    # this one for MPI-M, ECMWF, CMCC
    # use 1x6hr periods
    index32_params_6hr = {'indice_name': 'RX6hour',
                    'extreme_mode': 'max'}
    index32_6hr = {'name':'RX6hour', 'var': 'pr', 'period':'month', 'var_period':'Prim6hr', 'long_name':'Maximum 6-hour running mean maximum precipitation (from 1x6hr period)'}

    indices = [index17, index31]

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

    runs = [mohc_LM, mohc_HM]
    indices = define_indices()

    for run in runs:
        for index in indices:
            index_name = index['name']
            index_var = index['var']
            index_period = index['period']
            for year in years:
                fname = file_format.format(run['institute'], run['model'], run['label'], index['var_period'], index_var, run['grid'], index_var, index['var_period'], run['grid'], str(year))
                files = sorted(glob.glob(fname))
                print files[0]
                fname_out = os.path.basename(files[0])
                fname_out = fname_out.replace(index['var_period'], index_period)
                fname_out = fname_out.replace(index_var, index_name)
                f_out = fileout_format.format(run['institute'], run['model'], run['label'], index_period, index_name, run['grid'], 'latest', fname_out)
                if not os.path.exists(os.path.dirname(f_out)):
                    os.makedirs(os.path.dirname(f_out))
                if not index == index31:
                    icclim.indice(indice_name=index_name, in_files=files, var_name=index_var, slice_mode = index_period, out_file=f_out)
                else:
                    icclim.indice(user_indice=index31_params, in_files=files, var_name=index_var, slice_mode = index_period, out_file=f_out)
                # need to edit the long_name to make it consistent with this user variable
                    cmd = 'ncatted -O -h -a long_name,'+index_name+',o,c,"'+index['long_name']+'" '+f_out
                    os.system(cmd)

                add_model_metadata(files[0], f_out)
                cmd = 'nccopy -k 4 -d 3 '+f_out+' '+f_out+'.tmp.nc4'
                os.system(cmd)
                os.rename(f_out+'.tmp.nc4', f_out)
