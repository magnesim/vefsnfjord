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

from bio_uptake_tools import computeCfarr_dynamic
from bio_uptake_tools import get_Cw_from_opendrift_conc
from bio_uptake_tools import plot_bio_map

import yaml
import os



# Get command line arguments
if __name__ == "__main__":

    parser  = argparse.ArgumentParser()
#    parser.add_argument('--tstart', type=str, help='First time step', required=False, default='2020-09-10')
    
    args = parser.parse_args()
#    tstartS     = args.tstart







# read input config from yaml file
input_yaml = 'bio_input.yaml'
with open(input_yaml, 'r') as f:
    conf = yaml.safe_load(f)

tstart = datetime.strptime(conf['tstart'], '%Y-%m-%d %H')
tend   = datetime.strptime(conf['tend'], '%Y-%m-%d %H')

fn = '{}/{}'.format( conf['model_folder'],  conf['model_fn'].replace('EXPNM', conf['exp_nm']))


stations = conf['stations_to_eval']

print('Stations to evaluate:', stations)

positions = [conf['stations_coords'][item] for item in  conf['stations_to_eval'] ]
Npos = len(positions)


water_conc = conf.get('water_conc', 'model')  # default to 'model' if not specified
if water_conc == 'model':
    print('Using water concentration from model output:', fn)
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
    raise ValueError(f"Unknown water concentration option: {water_conc}. Choose 'model', 'analytical', or 'constant'.")

bio_uptake_model = conf.get('bio_uptake_model', 'dynamic')  # default to 'dynamic' if not specified


print('Number of positions:', Npos)
print('Positions:', positions)
imp_radius = conf.get('impact_radius', 0)  # default to 0 m if not specified
print('Impact radius:', imp_radius, 'm')


output_dir = conf.get('output_folder', None)
print('Output directory:', output_dir)

if False:
#if output_dir is not None:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f'Created output directory: {output_dir}')



if conf.get('plot_locations_on_map', False):
    print('Plotting locations on map...',output_dir)
    if output_dir is not None:
        ofn = f'{output_dir}/bio_uptake_map.png'
    else:
        ofn = None

    plot_bio_map(positions, fn, ofn, verbose=True)








# Initial bio concentration
Cf0 = conf.get('ana_Cf0', 0.)  # default to 0. Bq/kg if not specified

# Inintial water concentration
Cw0 = conf.get('ana_Cw0', 1.) 


# Uptake and elimination coefficients (s-1)
try: 
    kel  =  conf['depuration_coefficient']
    kup  =  conf['uptake_coefficient']
    bio_thalf = np.log(2.) / (kel*24*3600)  # days
    concentration_factor = kup / kel  # dimensionless

except KeyError:
    bio_thalf = conf['bio_half_life']  # days
    concentration_factor = conf['concentration_factor']  # dimensionless
    kel = np.log(2.) / (bio_thalf * 24 * 3600)  # convert days to seconds
    kup = concentration_factor * kel  # uptake coefficient is a factor of depuration coefficient



print('\n # ####################')
print('Bio uptake module parameters:')
print('Bio half life:                {:4.2f} days'.format(bio_thalf))
print('Concentration factor:         {:4.2f} '.format(concentration_factor))
print('Uptake coefficient (kup):     {:4.5e} '.format(kup))
print('Depuration coefficient (kel): {:4.5e} '.format(kel))
print('')










# Store results in arrays
time_arr_all   = []
Cwarr_all      = []
cfarr_all      = []

for ii, pos in enumerate(positions):
    print(f'{ii} of {Npos} - Processing {stations[ii]} position: {pos}')


    if water_conc == 'model':
        # Get water concentration from model output
        print('Using water concentration from model output:', fn)

        [time_arr, Cwarr] = get_Cw_from_opendrift_conc(fn, tstart, tend, pos=pos, imp_radius=imp_radius)
    
    
    elif water_conc == 'sinus' or water_conc == 'decay' or water_conc == 'randomspikes' or water_conc == 'setzero':
        from bio_uptake_tools import get_water_conc_ana
        # Get analytical water concentration
        print('Using analytical water concentration function.')
        nspikes = conf.get('ana_nspikes', 8)  # default to 8 spikes if not specified
        [time_arr, Cwarr] = get_water_conc_ana(tstart, tend, C0=Cw0, func=water_conc, decaycoeff=1e-2, nspikes=nspikes)
    
    
    elif water_conc == 'constant':
        # Use constant water concentration
        print(f'Using constant water concentration: {Cw0}')
        time_arr = np.array( [tstart+timedelta(hours=item) for item in range(int((tend-tstart).total_seconds()/3600))] )
        Cwarr = np.ones(len(time_arr)) * Cw0  # Example constant value
    
    else:
        raise ValueError(f"Unknown water concentration option: {water_conc}. Choose 'model', 'analytical', or 'constant'.")

    print( 'First timestep: ',time_arr[0], 'Last timestep: ', time_arr[-1])
    
    
    






    
    # Compute uptake 
    if bio_uptake_model == 'instant':
        from bio_uptake_tools import computeCfarr_instant
        print('Using instant bio uptake model.')
        cfarr = computeCfarr_instant(Cwarr, concentration_factor, time=time_arr)

    elif bio_uptake_model == 'dynamic':
        print('Using dynamic bio uptake model.')
        cfarr = computeCfarr_dynamic(kup,kel, Cwarr, time=time_arr, Cf0=Cf0)
    else:
        raise ValueError(f"Unknown bio uptake model: {bio_uptake_model}. Choose 'instant' or 'dynamic'.")

    time_arr_all.append(time_arr)
    Cwarr_all.append(Cwarr)
    cfarr_all.append(cfarr)
    














ncol=1
nrow=2
fig = plt.figure(figsize=[10,6])
ax1=plt.subplot(nrow,ncol,1)
for ii in range(Npos):
    ax1.plot(time_arr_all[ii], Cwarr_all[ii])
ax1.set_ylabel('Sea water concentration (Bq/m$^3$)')
ax1.grid()

ax2=plt.subplot(nrow,ncol,2)
for ii in range(Npos):
    ax2.plot(time_arr_all[ii], cfarr_all[ii], label = stations[ii])
ax2.grid()
ax2.set_ylabel('Bio concentration (Bq/kg)')
ax2.legend()

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.gcf().autofmt_xdate()
plt.tight_layout()

plt.show()



print('FiNiSH')