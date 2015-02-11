import matplotlib as plt

def plot_distribution(mass, locs=None, max_locs=5,
                      colormap='coolwarm', figsize=(15,15)):
    
    if locs is None:
        locs = np.round(np.linspace(0, mass.shape[1], max_locs))

    nlocs = len(locs)
    nlayers = mass.shape[2]

    fig, axs = subplots(nlayers, nlocs, figsize=figsize)
    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            m = mass.iloc[:,idxs[j],i,:]
            (m / sum(m)).plot(legend=False, colormap=colormap, ax=axs[i,j])
            if i == 0:
                axs[i,j].set_title('location: %s' % mass.columns[idxs[j]])

    return axs, fig
