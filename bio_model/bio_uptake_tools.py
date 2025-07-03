import numpy as np
from datetime import timedelta, datetime


def get_water_conc_ana(t0,t1, C0=1, func='', decaycoeff = 0.0011, nspikes=2,verbose=False, nwaves=1):
    """
    Generate analytical water concentration time series.
    Parameters
    ----------
    t0 : datetime
        Start time.
    t1 : datetime
        End time.
    C0 : float
        Initial water concentration.
    func : str
        Function to use for generating the water concentration time series.
    decaycoeff : float
        Decay coefficient for the water concentration.
    nspikes : int
        Number of random spikes to add to the water concentration.
    Returns
    -------
    time_arr : np.ndarray
        Array of datetime objects representing the time steps.
    water_conc : np.ndarray
        Array of water concentration values at each time step.
    """
    
    dt=1. 
    
    time_arr = np.array( [t0 + timedelta(hours=float(dt*i)) for i in range(int((t1 - t0).total_seconds()/3600)) ] )
    Nt=len(time_arr)
    tt=np.arange(Nt)

    water_conc = np.ones(Nt) * C0 
    
    if verbose:
        print('Use analytical functions: ', func)
    
    if 'sinus' in func:        
        water_conc = water_conc + water_conc*np.random.rand() + np.sin(2*nwaves*np.pi*tt/Nt) 
    if 'decay' in func:
        if verbose:
            print('Decay coefficient:', decaycoeff)
        water_conc = water_conc*np.random.rand()  + water_conc * np.exp(-tt*decaycoeff)
    if 'setzero' in func:
        water_conc = water_conc + water_conc*np.random.rand() 
        half = int(Nt/2)
        water_conc[half:] = 0.    
    if 'randomspikes' in func:
        ran1 = [int(item) for item in np.random.rand(nspikes)*Nt]
        ran2 = np.random.rand(nspikes)*C0*5 + C0 
        water_conc[ran1] = water_conc[ran1] + ran2
        
    
    water_conc[water_conc<0.] = 0.
    
    return time_arr, water_conc





def get_Cw_from_opendrift_conc(filename='', t0=None, t1=None, pos=None, species=['Total'], imp_radius=0.):

    """
    Extract water concentration from an OpenDrift netCDF file.
    Parameters
    ----------
    filename : str
        Path to the OpenDrift netCDF file.
    t0 : datetime
        Start time for the extraction.
    t1 : datetime
        End time for the extraction.
    pos : tuple
        Latitude and longitude of the position to extract data for.
    species : list
        List of species to extract data for. Default is ['Total'].
    Returns
    -------
    time_arr : np.ndarray
        Array of datetime objects representing the time steps.
    conc : np.ndarray
        Array of water concentration values at each time step for the specified species.
    """




    from netCDF4 import Dataset, date2index, num2date 
    
    
    nc = Dataset(filename,'r')
    nct = nc.variables['time']
    lat        = nc.variables['lat'][:]
    lon        = nc.variables['lon'][:]
    t0i = date2index(t0, nct, calendar=None, select='nearest')
    t1i = date2index(t1, nct, calendar=None, select='nearest')
    time = num2date(nct,units=nct.units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)
    time_arr = time[t0i:t1i]
    
    masktmp = nc.variables['land'][:]
    mask = np.zeros(masktmp.shape, dtype=np.int8)
    mask[masktmp==1] = 0
    mask[masktmp==0] = 1
    

    stations_ji = nearest_latlon([pos[0]], [pos[1]], lat, lon, mask=mask,verbose=True, radius=imp_radius)
    print('Nearest station index:', stations_ji)
    
    if len(stations_ji) == 0:
        print('No stations found within the specified radius.')
        return time_arr, np.zeros(len(time_arr))
    
    stay=stations_ji[0]
    stax=stations_ji[1]
    

    for sp in species:
        if sp=='Total':
            if len(stax) == 1 and len(stay) == 1:
                water_conc = nc.variables['concentration_smooth'][t0i:t1i, :, -1, stay,stax].squeeze()
                water_conc = np.sum(water_conc, axis=1)
            elif len(stax) > 1 and len(stay) > 1:
                water_conc = nc.variables['concentration_smooth'][t0i:t1i, :, -1, stay,stax].sum(axis=1).mean(axis=(1,2))
    print('conc shape:', water_conc.shape)
    
    nc.close()

    
    return time_arr, water_conc
    







def compute_Cb_dynamic(kup, kel, C1, time=None, Cb0=0.):

    """
    Compute times series of biological concentration with a dynamic model.
    Use sea water concentration and uptake and depuration rate constants.
    Parameters
    ----------
    kup : float
        Uptake rate constant.
    kel : float
        Depuration rate constant.
    C1 : np.ndarray
        Array of water concentration values at each time step.
    time : np.ndarray, optional
        Array of datetime objects representing the time steps. If None, it will be computed from C1.
    Cb0 : float, optional
        Initial concentration factor. Default is 0.
    Returns
    -------
    Cfarr : np.ndarray
        Array of concentration factor values at each time step.
    """


    Nt=len(time)
        
    dts = (time[1]-time[0]).total_seconds()
    
    Cfarr = np.zeros(Nt)
    Cfp0=Cb0


    for ii in range(Nt):    
        Cfp1  = Cfp0 + C1[ii] * kup*dts  - Cfp0 *kel*dts 
        Cfarr[ii] = Cfp1
        Cfp0 = Cfp1 


    return Cfarr 




def compute_Cb_instant(C1, cr,  time=None):

    """
    Compute the water concentration time series based on an instantaneous model.
    Parameters
    ----------
    C1 : np.ndarray
        Array of water concentration values at each time step.
    cr : float
        Concentration rate constant.
    time : np.ndarray, optional
        Array of datetime objects representing the time steps. If None, it will be computed from C1.
    Returns
    -------
    Cfarr : np.ndarray
        Array of concentration factor values at each time step.
    """

    Nt=len(time)        
    dts = (time[1]-time[0]).total_seconds()
    Cfarr = C1 * cr


    return Cfarr 




def plot_bio_map(pos, filename, ofn=None, ext = [12.4, 13.3, 65.8, 66.09], verbose=False):
    """
    Plot a map with the specified positions and topography.
    Parameters
    ----------
    pos : list of tuples
        List of tuples containing latitude and longitude of positions to plot.
    filename : str
        Path to the netCDF file containing topography data.
    ofn : str
        Output filename for the plot.
    ext : list, optional
        List specifying the extent of the map in the format [lon_min, lon_max, lat_min, lat_max]. Default is [12.4, 13.3, 65.8, 66.09].
    Returns
    -------
    None
    """


    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    if verbose:
        print('Plotting map with positions:', pos)
        print('Using topography data from:', filename)

    nc = Dataset(filename,'r')
    #nct = nc.variables['time']
    lat        = nc.variables['lat'][:]
    lon        = nc.variables['lon'][:]
    topo        = nc.variables['topo'][:]
    topo[topo.mask==1] = 0
    topo.mask=None
    nc.close()
    if verbose:
        print('done')

    
    proj_pp=ccrs.PlateCarree()
    
    fig=plt.figure()
    ax=plt.subplot(projection=ccrs.Orthographic(13,66))
    ax.set_extent(ext,proj_pp)
#    ax.coastlines(resolution='10m', zorder=5, color='k')
#    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=4)
#    ax.add_feature(cfeature.OCEAN, facecolor='lightblue', zorder=3)

    c2 = ax.contour(lon,lat, topo, [0.,1.], colors='k', linewidths=1.6, 
                    zorder=6, transform=proj_pp)
    for ii in range(len(pos)):
        s2 = ax.scatter(pos[ii][1], pos[ii][0], transform=proj_pp)
    if not ofn==None:
        plt.savefig(ofn, dpi=300,bbox_inches='tight')
        plt.close()
    if verbose:
        print('Plot saved to:', ofn)




def nearest_latlon(latp1,lonp1,lata,lona,mask=None, verbose=False,radius=0.):
    """
    Find the nearest latitude and longitude grid points to given points.
    Parameters
    ----------
    latp1 : list or np.ndarray
        List or array of latitudes of points to be evaluated.
    lonp1 : list or np.ndarray
        List or array of longitudes of points to be evaluated.
    lata : np.ndarray
        Array of latitudes of the grid points.
    lona : np.ndarray
        Array of longitudes of the grid points.
    mask : np.ndarray, optional
        Mask array to avoid points on land. If None, no masking is applied.
    verbose : bool, optional
        If True, print additional information during processing. Default is False.
    radius : float, optional
        Radius in meters to include for the nearest grid points. If 0, only the nearest grid point is returned. Default is 0.
    Returns
    -------
    mltidx : np.ndarray
        Array of indices of the nearest grid points for each input point.   
    """
    
    import numpy as np
    


    Re= 6371000.  # Earth radius m
    
#    if mask.all() is not None:
    if mask is not None:
        # avoid points on land
        latam=lata*mask
        lonam=lona*mask
    
    mltidx=[]
    for pp in range(len(latp1)):
        latp=latp1[pp]
        lonp=lonp1[pp]
        # Diff between point and gridcells
        dlat=latam-latp
        dlon=lonam-lonp
        
        # Middle latitude
        mlat = 0.5*(lata+latp)

        # Distance in meters    
        dr = np.sqrt( (dlat*2*np.pi* Re / 360.)**2 + (dlon*2*np.pi*Re*np.cos(np.deg2rad(mlat))/360.)**2  )
    
        # Find index of nearest gridpoint
        
        minval = np.min(dr)
        
        if radius == 0:
            idx=np.argmin(dr)
            nearest = np.unravel_index(idx, dr.shape)
            mltidx.append([[nearest[0]], [nearest[1]]])  # Append as a list of lists for consistency
            if verbose:
                print(f"Point {pp}: lat={latp}, lon={lonp}, nearest point at index={nearest} with distance={minval:.2f} m")

        elif radius > 0.:
            # Find all grid points within the specified radius
            within_radius = np.where(dr <= radius)
            if verbose:
                print(f"Point {pp}: lat={latp}, lon={lonp}, points within radius({radius}m): {len(within_radius[0])}, min distance={minval:.2f} m")
            if len(within_radius[0]) > 0:
                # If there are points within the radius, append them to mltidx
                if verbose:
                    print(f"Point {pp}: lat={latp}, lon={lonp}, nearest points within radius at indices={within_radius}")
                mltidx.append(within_radius)
            else:
                idx=np.argmin(dr)  
                nearest = np.unravel_index(idx, dr.shape)
                print("No points found within the specified radius. Use nearest",nearest)
                mltidx.append([[nearest[0]], [nearest[1]]])


    # Convert mltidx to a numpy array for consistency
    if isinstance(mltidx, list):
        mltidx = [np.array(item) for item in mltidx]
#    else:
#        mltidx = np.array(mltidx)
#    if isinstance(mltidx, list) and mltidx == []:
    
    if verbose:
        print("Nearest grid points indices:", mltidx)
    
    if isinstance(mltidx, list) and len(mltidx) == 1:
        mltidx = np.array(mltidx[0])  # If only one point, convert to a single array
    elif isinstance(mltidx, list) and len(mltidx) > 1:
        mltidx = np.array(mltidx)  # Convert list of arrays to a 2D array


    if verbose:
        print("Final nearest grid points indices:", mltidx)
    if len(mltidx) == 0:
        return mltidx  # If no points found, return empty mltidx

    plotting = True
    if plotting:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        try:
            npoints = len(mltidx[0])
        except TypeError:
            npoints = 1  # If mltidx is a single point, set npoints to 1
        for ii in range(npoints):
            idx = np.array([mltidx[0][ii], mltidx[1][ii]])
            ax.plot(lona[idx[0],idx[1]], lata[idx[0],idx[1]], color='red', ls='none', marker='o', label='Nearest Point')
        ax.scatter(lonp1, latp1, color='green', label='Input Points')
        ax.contour(lona, lata, mask, levels=[0.01,0.5,0.99], colors='black', linewidths=0.5, linestyles='dashed')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_title('Nearest Grid Points')
#        ax.legend()

    return mltidx





# #PLOTTING FUNCTIONS

def plot_bio_concentration_timeseries(data, suptitle='', verbose=False):

    """
    Plot bio concentration time series.
    Parameters
    ----------
    data : list of dicts
        List containing dictionaries with 'time', 'ydata', and 'ylabel' keys.
    suptitle : str, optional
        Super title for the plot. Default is an empty string.
    Returns
    -------
    None
    """
    
    import matplotlib.pyplot as plt

    ncol = 1
    nrow = len(data)
    
    fig = plt.figure(figsize=[10, 3*nrow])
    fig.suptitle(suptitle)

    for ii, d in enumerate(data):
        ax = plt.subplot(nrow, ncol, ii + 1)
        nstations = len(data[ii]['ydata'])
        if verbose:
            print('Plotting data for:', d['ylabel'], 'with', nstations, 'stations')
        for jj in range(nstations):
            ax.plot(d['time'][jj], d['ydata'][jj], label=d['labels'][jj])
        if ii == 0:
            ax.legend()
        ax.set_ylabel(d['ylabel'])
        ax.grid()
        ax.set_xlim([d['time'][0][0] , d['time'][0][-1]])

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    #plt.show()