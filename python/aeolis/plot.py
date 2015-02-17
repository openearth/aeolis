import os
import numpy as np
import matplotlib.pyplot as plt

import filesys

def plot_distribution(mass, locs=None, max_locs=5,
                      colormap='coolwarm', figsize=(15,15)):
    '''Plot sediment distribution in model domain

    Parameters
    ----------
    mass : pandas.DataFrame
        Dataframe with mass contents of grid cells
    locs : iterable, optional
        Indices of spatial grid cells to include in the plot. If not given max_locs
        locations are chosen at a constant distance form each other
    max_locs : int, optional
        Maximum number of spatial grid cell locations to include
    colormap : string or matplotlib colormap, optional
        Colormap for coloring time evolution of sediment distributions
    figsize : 2-tuple, optional
        Figure dimensions

    Returns
    -------
    matplotlib.Figure
        Figure object containing axes
    numpy.ndarray of matplotlib.Axes objects
        Array with axes objects
    '''
    
    if locs is None:
        locs = [int(x) for x in np.round(np.linspace(0, mass.shape[1]-1, max_locs))]

    nlocs = len(locs)
    nlayers = mass.shape[2]

    fig, axs = plt.subplots(nlayers, nlocs, figsize=figsize)
    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):

            m = mass.iloc[:,locs[j],i,:]
            (m / m.sum()).plot(legend=False, colormap=colormap, ax=axs[i,j])
            
            if i == 0:
                axs[i,j].set_title('location: %s' % mass.items[locs[j]])
            if j == 0:
                axs[i,j].set_ylabel('fraction')
            if i == axs.shape[0]-1:
                axs[i,j].set_xlabel('grain fraction')

    return fig, axs


def plot_profiles(fpath, start=0, stop=np.inf, step=1, height_markers=[-1, 1],
                  max_subplots=10, clim=.3e-3, colormap='coolwarm', figsize=(20,30)):
    '''Plot profile evolution including sediment sorting and armoring

    Parameters
    ----------
    fpath : string
        Path to model output
    start : int, optional
        Start position used for data read
    stop : int, optional
        Stop position used for data read
    step : int, optional
        Step or stride used for data read
    height_markers : iterable, optional
        Mark heights in list, for example tidal range
    clim : float
        Maximum grain size in color scale
    max_subplots : int
        Maximum number of subplots to be rendered. Above this value
        will result in an error.
    colormap : string or matplotlib colormap, optional
        Colormap for coloring grain size distribution
    figsize : 2-tuple, optional
        Figure dimensions
    '''

    # read dimensions
    d = filesys.load_dimensions(fpath)
    nx = d['nx']
    nl = d['nl']
    x = d['ax_x']

    # read data files
    z = filesys.load_dataframe(os.path.join(fpath, 'z.out'),
                               start=start, stop=stop, step=step).as_matrix()
    d50 = filesys.load_dataframe(os.path.join(fpath, 'd50.out'),
                                 start=start, stop=stop, step=step).as_matrix()

    # determine starting distribution
    n = d50.shape[0]
    d50_0 = np.mean(d50[0,:,:])
    d50[d50 == 1e-3] = d50_0 # FIXME ?

    if n > max_subplots:
        raise ValueError('Maximum number of %d subplots exceeded with %d, limit the data' % (max_subplots, n-max_subplots))

    fig, axs = plt.subplots(n, 2, figsize=figsize)
    for i in range(n):

        # construct pcolormesh matrices
        X = np.repeat(x, nl).reshape((-1,nl))
        Z1 = np.repeat(z[i,:], nl).reshape((-1,nl))
        Z2 = np.repeat(z[i,:] - z[0,:], nl).reshape((-1,nl))
        L = .25 / nl * np.repeat(np.arange(nl), nx).reshape((-1,nx)).T
        C = d50[i,:,:] - d50_0

        zlim = [np.max(np.abs(Z1)),
                np.max(np.abs(Z2))]

        # plot profile
        axs[i,0].pcolormesh(X, Z1-np.ptp(Z1)*L, C,
                            vmin=-clim, vmax=clim, cmap=colormap)
        axs[i,1].pcolormesh(X, Z2-np.ptp(Z2)*L, C,
                            vmin=-clim, vmax=clim, cmap=colormap)

        # plots markers
        markers = [np.interp(m, z[i,:], x) for m in height_markers]
        for j, ax in enumerate(axs[i,:]):
            if i == 0 and j == 1:
                continue
            
            for m in markers:
                ax.plot([m, m], [-zlim[j], zlim[j]], '--b', linewidth=1)
            
            ax.set_xlim((0, np.max(x)))
            ax.set_ylim((-zlim[j], zlim[j]))
            ax.grid()

        axs[i,0].set_ylabel('height [$m$]')
        axs[i,1].set_ylabel('bed level change [$m$]')

    axs[0,0].set_title('profile')
    axs[0,1].set_title('bed level change')

    axs[-1,0].set_xlabel('x [$m$]')
    axs[-1,1].set_xlabel('x [$m$]')

    axs[0,1].set_axis_off()

    return fig, axs
            
                                                                                
def plot_sediment_balance(fpath, rhom=1650., figsize=(10,5)):
    '''Compute and plot sediment balance of model

    Parameters
    ----------
    fpath : string
        Path to model output
    rhom : float, optional
        Density of soil including pores
    figsize : 2-tuple, optional
        Figure dimensions

    Returns
    -------
    matplotlib.Figure
        Figure object containing axes
    numpy.ndarray of matplotlib.Axes objects
        Array with axes objects
    '''
    
    # read dimensions
    d = filesys.load_dimensions(fpath)
    dx = d['dx']
    dt = d['dt_out']

    # read data files
    z = filesys.load_dataframe(os.path.join(fpath, 'z.out'))
    mass = filesys.load_dataframe(os.path.join(fpath, 'mass.out'))
    supply = filesys.load_dataframe(os.path.join(fpath, 'supply.sum.out'))
    Ct = filesys.load_dataframe(os.path.join(fpath, 'Ct.avg.out'))
    u = filesys.load_dataframe(os.path.join(fpath, 'u.avg.out'))
    
    # compute erosion/deposition
    v1 = np.sum(np.maximum(0., z.iloc[-1,:] - z.iloc[0,:]) * dx) # deposition
    v2 = np.sum(np.maximum(0., z.iloc[0,:] - z.iloc[-1,:]) * dx) # erosion

    v3 = -sum(z.iloc[-1,:] - z.iloc[0,:]) * dx
    v4 = -sum(mass.sum(axis=3).sum(axis=2).iloc[:,-1] - \
              mass.sum(axis=3).sum(axis=2).iloc[:,0]) * dx / rhom

    v5 = sum(supply.sum(axis=2).sum(axis=1)) * dx / rhom

    v6 = Ct.sum(axis=2).iloc[-1,:].multiply(u[0]).sum() * dt / rhom
    v7 = Ct.sum(axis=2).iloc[:,-1].sum() * dx / rhom

    # print output statistics
    print 'total erosion:                %f' % v2
    print 'total deposition:             %f' % v1
    print 'net loss from profile (z):    %f' % v3
    print 'net loss from profile (mass): %f' % v4
    print 'net loss over domain border:  %f' % (v6+v7)
    print 'net supply to wind:           %f' % v5

    # compute time series
    #s1 = supply.sum(axis=2).sum(axis=0).cumsum() * dx / rhom
    s2 = Ct.sum(axis=2).iloc[-1,:].multiply(u[0]).cumsum() * dt / rhom + \
         Ct.sum(axis=2).sum(axis=0) * dx / rhom
    s3 = -z.iloc[:,:].subtract(z.iloc[0,:]).sum(axis=1) * dx

    # plot time series
    fig, axs = plt.subplots(1, 2, figsize=figsize)
    
    s3.plot(legend=False, label='loss from bed', ax=axs[0])
    #s1.plot(legend=False, label='supply to wind', ax=axs[0])
    s2.plot(legend=False, label='loss over domain border +\nsaltation', ax=axs[0])
    axs[0].set_xlabel('time')
    axs[0].set_ylabel('supply / loss [$m^3$]')
    axs[0].set_title('supply vs. loss')
    axs[0].legend(loc='lower right')
    
    (s3-s2).plot(legend=False, ax=axs[1])
    axs[1].set_xlabel('time')
    axs[0].set_ylabel('loss [$m^3$]')
    axs[1].set_title('difference')

    return fig, axs
