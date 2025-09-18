"""
Development of bio uptake module
 
Magne Simonsen
MET Norway
2022 May
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

from datetime import datetime, timedelta

from bio_uptake_tools import compute_Cb_dynamic
from bio_uptake_tools import compute_Cb_instant
from bio_uptake_tools import get_Cw_from_opendrift_conc
from bio_uptake_tools import get_Cw_from_trajan_file
from bio_uptake_tools import get_water_conc_ana
from bio_uptake_tools import get_water_salinity
from bio_uptake_tools import plot_bio_map
from bio_uptake_tools import plot_bio_concentration_timeseries

import yaml
import glob
import os



# Get command line arguments
if __name__ == "__main__":

    parser  = argparse.ArgumentParser()
#    parser.add_argument('--tstart', type=str, help='First time step', required=False, default='2020-09-10')
    parser.add_argument('-f', '--configfile', type=str, help='Path to the config file', required=False, default='bio_input.yaml')   
    parser.add_argument('--verbose', action='store_true', help='Verbose output', required=False, default=False)

    args = parser.parse_args()

    verbose = args.verbose
    print('verbose:',verbose)

    if args.configfile is not None:
        input_yaml = args.configfile
        if not os.path.exists(input_yaml):
            raise FileNotFoundError(f"Configuration file not found: {input_yaml}")










print('*******************')
print('Run bio uptake model')
print('')


# read input config from yaml file
print('Reading input configuration from: {}'.format(input_yaml))

with open(input_yaml, 'r') as f:
    conf = yaml.safe_load(f)


tstart = datetime.strptime(conf.get('tstart','2000-01-01'), '%Y-%m-%d %H')
tend   = datetime.strptime(conf.get('tend','2000-02-01'), '%Y-%m-%d %H')


exp_nm = conf.get('exp_nm', '')  # default to empty string if not specified


stations = conf.get('stations_to_eval', None)  # list of stations to evaluate


positions = [conf['stations_coords'][item] for item in  stations ]
Npos = len(positions)


water_conc = conf.get('water_conc', 'model')  # default to 'model' if not specified
imp_radius = conf.get('impact_radius', 0)  # default to 0 m if not specified
if water_conc == 'model':
    transport_model_file = '{}/{}'.format( conf.get('model_folder',''),  conf.get('model_fn','').replace('EXPNM', exp_nm))
    print('Using water concentration from model output:', transport_model_file)
    fill_conc = conf.get('fill_concentration', 0.0)  # default to 0 Bq/m^3 if not specified
    print(f'Use fill concentration: {fill_conc} Bq/m^3')

    if not os.path.exists(transport_model_file):
        raise FileNotFoundError(f"Transport model file not found: {transport_model_file}")
elif water_conc == 'constant':
    print('Using constant water concentration.')
elif water_conc == 'decay':
    print('Using decay water concentration.')
elif water_conc == 'sinus':
    print('Using sinus function water concentration.')
elif water_conc == 'randomspikes':
    print('Using randomspikes function water concentration.')
elif water_conc == 'setzero':
    print('Using setzero function water concentration.')
else:
    raise ValueError(f"Unknown water concentration option: {water_conc}. Choose 'model', 'decay', 'sinus', 'randomspikes', 'setzero', or 'constant'.")

bio_uptake_model = conf.get('bio_uptake_model', 'dynamic')  # default to 'dynamic' if not specified


print('\nNumber of positions:', Npos)
for ii, pos in enumerate(stations):
    print(f'{pos}: {positions[ii]}')



saveplots     = conf.get('saveplots', False)

if saveplots:
    output_dir = conf.get('output_folder', None)
    print('Output directory: {}'.format( output_dir) )

    if output_dir is not None:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f'Created output directory: {output_dir}')



if conf.get('plot_locations_on_map', False):
    print('Plotting locations on map...')

    ofn = f'{output_dir}/bio_uptake_map_{exp_nm}.png' if saveplots else None 
    plot_bio_map(positions, transport_model_file, ofn=ofn, title=exp_nm, labels=stations, verbose=verbose)





# Initial bio concentration
ini_Cb = conf.get('initial_bio_conc', 0.)  # default to 0. Bq/kg if not specified







print('\n # ####################')
print('Bio uptake module parameters:')
# Uptake and elimination coefficients (s-1)


salt_CF            =  conf.get('salt_CF',False)  # If true, use salinity-dependent CF values, 
                                                 # otherwise, use constant CF value
salt_intervals     =  conf.get('salt_intervals_CF', None)


if salt_CF:



    salt_default  = conf.get('salt_default',0.)
    bio_thalf = conf['bio_half_life']  # days
    kel = np.log(2.) / (bio_thalf * 24 * 3600)  # convert days to seconds
    upper_limits= []
    lower_limits= []
    kup_salts   = []
    cf_salts    = []
    kup=0.
    concentration_factor=0.

    for s_int in salt_intervals:
        lower_limits.append(s_int[0])
        upper_limits.append(s_int[1])
        kup_salts.append( s_int[2] * kel )  # uptake coefficient is a factor of depuration coefficient
        cf_salts.append( s_int[2])
        

    print(f'Use salinity-dependent uptake coefficients ')
    print(f'Lower limits: {lower_limits}')
    print(f'Upper limits: {upper_limits}')
    print(f'CF values   : {cf_salts} ')
    print(f'kups        : {kup_salts}')


else:

    try: 
        kel  =  conf['depuration_coefficient']
        kup  =  conf['uptake_coefficient']
        bio_thalf = np.log(2.) / (kel*24*3600)  # days
        concentration_factor = kup / kel  # dimensionless
        print('Using provided uptake and depuration coefficients:')

    except KeyError:
        bio_thalf = conf['bio_half_life']  # days
        concentration_factor = conf['concentration_factor']  # dimensionless
        kel = np.log(2.) / (bio_thalf * 24 * 3600)  # convert days to seconds
        kup = concentration_factor * kel  # uptake coefficient is a factor of depuration coefficient
        print('Using bio half life and concentration factor to compute uptake and depuration coefficients:')



print('Bio half life:                {:4.2f} days'.format(bio_thalf))
print('Concentration factor:         {:4.2f} L/kg '.format(concentration_factor))
print('Uptake coefficient (kup):     {:4.2e} s-1'.format(kup))
print('Depuration coefficient (kel): {:4.2e} s-1'.format(kel))
print('')
















print ('\n # ####################')
print ('Get water concentration and compute bio uptake')

if water_conc == 'model':
    print('Using water concentration from model output:', exp_nm)
    print('Impact radius: {} m'.format(imp_radius) )

else: 
    print('Using analytical water concentration function: {}'.format(water_conc))

    ini_Cw                = conf.get('ana_initial_water_conc', 1.)    # Initial water concentration
    ana_decaycoeff        = conf.get('ana_decay_coeff', 1e-2)         # default to 1e-2 if not specified
    nspikes               = conf.get('ana_nspikes', 8)                # default to 8 spikes if not specified
    nwaves                = conf.get('ana_nwaves', 1)                 # default to 1 wave if not specified



# Store results in arrays
time_arr_all                   = []
seawater_concentration_all     = []
bio_concentration_all          = []
salt_time_all                  = []
seawater_salinity_all          = []
cf_ts_all                      = []

print('')
for ii, pos in enumerate(positions):
    print(f'{ii+1} of {Npos} - Processing {stations[ii]} position: {pos}')


    if water_conc == 'model':
        # Get water concentration from model output
        #[time_arr, seawater_concentration] = get_Cw_from_opendrift_conc(transport_model_file, tstart, tend, pos=pos, imp_radius=imp_radius, verbose=verbose)
        [time_arr, seawater_concentration] = get_Cw_from_trajan_file(transport_model_file, tstart, tend, pos=pos, verbose=verbose)
        dt = (time_arr[1] - time_arr[0]).total_seconds()  # Time step in seconds
        fill_time = np.array([tstart + timedelta(seconds=item) for item in range(0, int( (tend-tstart).total_seconds() )+1, int(dt))])
        if verbose:
            print('Time range from model: ',time_arr[0],time_arr[-1], len(time_arr))
            print('Total time range:      ',fill_time[0],fill_time[-1], len(fill_time))
        fill_conc = np.ones(len(fill_time)) * fill_conc  # Fill concentration in Bq/m^3
        match_ft  = (fill_time >= time_arr[0]) & (fill_time <= time_arr[-1])  # 
        fill_conc[match_ft] = seawater_concentration  # Replace fill concentration with model concentration at the correct time steps

        time_arr                 = fill_time  # Use fill time as time array
        seawater_concentration   = fill_conc  # Use fill concentration as seawater concentration



    elif water_conc == 'sinus' or water_conc == 'decay' or water_conc == 'randomspikes' or water_conc == 'setzero':
        [time_arr, seawater_concentration] = get_water_conc_ana(tstart, tend, C0=ini_Cw, func=water_conc, decaycoeff=ana_decaycoeff,nwaves=nwaves, nspikes=nspikes, verbose=verbose)
    
    
    elif water_conc == 'constant':
        print(f'Using constant water concentration: {ini_Cw} Bq/m^3')
        time_arr = np.array( [tstart+timedelta(hours=item) for item in range(int((tend-tstart).total_seconds()/3600)+1)] )
        seawater_concentration = np.ones(len(time_arr)) * ini_Cw  # 
    
    else:
        raise ValueError(f"Unknown water concentration option: {water_conc}. Choose 'model', 'analytical', or 'constant'.")

    if verbose:
        print('Convert sea water concentration from Bq/m^3 to Bq/L')
    seawater_concentration = seawater_concentration / 1000.  # Convert Bq/m^3 to Bq/L

    print('Done setting water concentration for position:', stations[ii])
    print( 'First timestep: ',time_arr[0], 'Last timestep: ', time_arr[-1], 'len(time)',len(time_arr), 'len(conc)',len(seawater_concentration) )
    
    
    

    if salt_CF:
        print('****')
        print(' Get salinity from ocean model') 
        datadir   =  conf.get('salt_datadir','')
        filenames = conf.get('salt_filenames','')
        if verbose:
            print(f'Look into files: {datadir}/{filenames}')
            print(f'Time range: {tstart} to {tend} Position: {pos}')
        ocean_model_files = np.sort( glob.glob(f'{datadir}/{filenames}') ).ravel()
        [salt_time, seawater_salinity]   = get_water_salinity(ocean_model_files, tstart, tend, pos=pos, imp_radius=imp_radius, level=[-3,-1] , verbose=verbose)


        dt = (salt_time[1] - salt_time[0]).total_seconds()  # Time step in seconds
        fill_time = np.array([tstart + timedelta(seconds=item) for item in range(0, int( (tend-tstart).total_seconds() )+1, int(dt))])

        fill_salt = np.ones(len(fill_time)) * salt_default  # Fill concentration in Bq/m^3
        match_ft  = (fill_time >= salt_time[0]) & (fill_time <= salt_time[-1])  # 
        fill_salt[match_ft] = seawater_salinity  # Replace fill concentration with model concentration at the correct time steps

        salt_time                 = fill_time  # Use fill time as time array
        seawater_salinity   = fill_salt  # Use fill concentration as seawater concentration

        
        print( f'Salt time, len, range: {len(salt_time)} {salt_time[0]} {salt_time[-1]}')
        print( f'Salt len, range :      {len(seawater_salinity)} {np.min(seawater_salinity):4.2f} {np.max(seawater_salinity):4.2f}')






    
    # Compute uptake 
    if bio_uptake_model == 'instant':
        if verbose:
            print('Using instant bio uptake model.')
        bio_concentration = compute_Cb_instant(seawater_concentration, concentration_factor, time=time_arr)

    elif bio_uptake_model == 'dynamic':
        if verbose:
            print('Using dynamic bio uptake model.')
        if salt_CF:
            bio_concentration, cf_ts = compute_Cb_dynamic(kup,kel, seawater_concentration, time=time_arr, Cb0=ini_Cb, salt_ts=seawater_salinity, salt_koeffs=[lower_limits,upper_limits,kup_salts])
        else:
            bio_concentration, cf_ts = compute_Cb_dynamic(kup,kel, seawater_concentration, time=time_arr, Cb0=ini_Cb)
    else:
        raise ValueError(f"Unknown bio uptake model: {bio_uptake_model}. Choose 'instant' or 'dynamic'.")

    time_arr_all.append(time_arr)
    seawater_concentration_all.append(seawater_concentration)
    bio_concentration_all.append(bio_concentration)
    if salt_CF:
        salt_time_all.append(salt_time)
        seawater_salinity_all.append(seawater_salinity)
        cf_ts_all.append(cf_ts)
        print(f'Mean salt: {np.mean(seawater_salinity):4.2f} Mean CF: {np.mean(cf_ts):4.2f}')













print('\n # ####################')
print('Plotting results')


to_plotting01 = {
   'time'  : time_arr_all,
   'ydata' : seawater_concentration_all,
   'ylabel': 'Sea water concentration (Bq/L)', 
   'labels': stations,
   }

to_plotting02 = {
   'time'  : time_arr_all,
   'ydata' : bio_concentration_all,
   'ylabel': 'Bio concentration (Bq/kg dw)', 
   'labels': stations,
    }

apparent_concentration_factor = []
for ii in range(Npos):
    concfact = bio_concentration_all[ii] / seawater_concentration_all[ii]
    concfact[np.isnan(concfact)] = 0.  # Replace NaN with 0
    apparent_concentration_factor.append( concfact )

to_plotting03 = {
   'time'  : time_arr_all,
   'ydata' : apparent_concentration_factor,
   'ylabel': 'Apparent concentration factor (dimensionless)', 
   'labels': stations,
    }

to_plotting = [to_plotting01, to_plotting02, to_plotting03]

if salt_CF:
    to_plotting04 = {
        'time'  : salt_time_all,
        'ydata' : seawater_salinity_all,
        'ylabel': 'Seawater salinity', 
        'labels': stations,
        }
    to_plotting05 = {
        'time'  : salt_time_all,
        'ydata' : cf_ts_all,
        'ylabel': 'Concentration Factor ($C_F$)', 
        'labels': stations,
        }
    to_plotting.append(to_plotting04)
    to_plotting.append(to_plotting05)

thalfstr = '$t_{1/2}$'
title = '{}\nWater concentration: {} \nBio uptake model: {} \n$C_F={}$,  {}={} '.format(exp_nm, water_conc, bio_uptake_model, concentration_factor,thalfstr, bio_thalf)


outfilename = f'{output_dir}/bio_model_{exp_nm}02.png' if saveplots else None
plot_bio_concentration_timeseries( to_plotting , suptitle=title, verbose=verbose, outfile=outfilename)



if not saveplots:
    plt.show()



print('FiNiSH')