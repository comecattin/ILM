import pyemma
from dimension_reduction import *
from load_feat import *

def implied_time_scale(cluster,lags,nits):

    its = pyemma.msm.its(cluster.dtrajs,lags=lags,nits=nits,errors='bayes')

    return its

def plot_its(its,data, cluster, save=False, display=False,outdir=''):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    data_concatenated = np.concatenate(data.get_output())
    pyemma.plots.plot_feature_histograms(data_concatenated, feature_labels=['Feat 1', 'Feat 2'], ax=axes[0])
    pyemma.plots.plot_density(*data_concatenated.T, ax=axes[1], cbar=False, alpha=0.1)
    axes[1].scatter(*cluster.clustercenters.T, s=15, c='C1')
    axes[1].set_xlabel('Feat 1')
    axes[1].set_ylabel('Feat 2')
    pyemma.plots.plot_implied_timescales(its, ax=axes[2], units='ps')
    
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/its_validation.pdf', dpi=300, bbox_inches="tight")
    if display:
        plt.show()



if __name__ == '__main__':

     # Path
    pdb = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb'
    traj = ['/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS01_md_all_fitBB_protonly.xtc','/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS02_md_all_fitBB_protonly.xtc']
    outdir='/home/ccattin/Documents/Code/outputs'
    # Feat
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    # Parameters
    save = False
    display = False
    T = 300
    lag=1000
    nits = 4
    lags=[1, 2, 5, 10, 20, 50]

    pair_indices = tools.create_pairIndices_from_pairNames(pdb,pairNames)
    feat = create_feat(pdb,pair_indices)
    data = load_data(traj,feat)

    plot_feat_hist(data,feat,
                   display=display,
                   save=save,
                   outdir=outdir)
    plot_density_energy(data=data,
                        T=T,
                        pairNames=pairNames,
                        display=display,
                        save=save,
                        outdir=outdir
                        )
    pca = pca_reduction(data=data,
                        T=T,
                        save=save,
                        display=display,
                        outdir=outdir)
    tica = tica_reduction(data=data,
                          lag=lag,
                          T=T,
                          save=save,
                          display=display,
                          outdir=outdir)
    cluster = clustering(reduction=tica,
                         method='kmeans',
                         k=200,
                         stride=1,
                         save=save,
                         display=display,
                         outdir=outdir)
    
    its = implied_time_scale(cluster=cluster,
                             lags=lags,
                             nits=nits)
    

    plot_its(its=its,
             data=tica,
             cluster=cluster,
             save=save,
             display=display,
             outdir=outdir)