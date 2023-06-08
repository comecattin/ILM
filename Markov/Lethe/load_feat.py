import pyemma
import tools
import numpy as np
import matplotlib.pyplot as plt


def create_feat(pdb):
    """Create a PyEmma featurizer

    Parameters
    ----------
    pdb : str
        Path to a reference topology .pdb file
    
    Returns
    -------
    feat : PyEmma object
        PyEmma featurizer
    """
    # Create the feat
    feat = pyemma.coordinates.featurizer(pdb)
    return feat

def offset(feat):
    """Get the offset for residue

    Parameters
    ----------
    feat : PyEmma object
        PyEmma featurizer

    Returns
    -------
    offset : int
        offset to subtract
    """
    with open(feat.topologyfile) as f:
        lines = f.readlines()
        for line in lines:
            if 'ATOM' in line:
                return int(line.split()[5])

def feat_atom_distances(feat,pair_indices):
    """Add atom distances to the featurizer

    Parameters
    ----------
    feat : PyEmma object
        PyEmma featurizer
    pair_indices : list
        List containing the indices of the pair to compute distances

    Returns
    -------
    feat : PyEmma object
        PyEmma featurizer
    """
    # Add the right coordinates to the feat
    feat.add_distances(indices=pair_indices, periodic=True, indices2=None)
    print(f"PyEmma feat description:\n{feat.describe()}")
    return feat

def feat_residue_midist(feat,pair_indices):
    """Add the minimum distance between residues

    Parameters
    ----------
    feat : PyEmma object
        PyEmma featurizer
    pair_indices : list
        List containing the indices of residue pairs to compute distances

    Returns
    -------
    feat : PyEmma object
        PyEmma featurizer
    """
    # Remove the offset that PyEmma put
    pair_indices = pair_indices - offset(feat)
    # Add the correct residues
    feat.add_residue_mindist(residue_pairs=pair_indices)
    print(f"PyEmma feat description:\n{feat.describe()}")
    return feat


def load_data(traj, feat, stride, ram):
    """Load all the given trajectories in PyEmma

    Parameters
    ----------
    traj : list
        List that contain the path of all the trajectories to analyze
    feat : PyEmma object
        PyEmma featurizer
    stride : int
        Number of stride to consider. Only read every stride'th frame.
    ram : bool
        Load the data in the ram of not

    Returns
    -------
    data : List
        List of all the value of the feat for each snapshot
    """
    if not ram:
        data = pyemma.coordinates.source(traj, features=feat, stride=stride)
        out = data.get_output()
        print("Lengths (number of trajectories ):", len(out))
        print(
            "Shape of elements (for each trajectories, number of timestep, number of features):",
            out[0].shape,
        )

    if ram:
        data = pyemma.coordinates.load(traj, features=feat, stride=stride)

        print("Lengths (number of trajectories ):", len(data))
        print(
            "Shape of elements (for each trajectories, number of timestep, number of features):",
            data[0].shape,
        )

    return data


def plot_feat_hist(data, feat, display=False, save=False, outdir=""):
    """Make a histogram plot of the feat values

    Parameters
    ----------
    data : List
        List of all the value of the feat for each snapshot
    feat : PyEmma object
        PyEmma featurizer
    display : bool, optional
        Display or not the plots, by default False
    save : bool, optional
        Save or not the plots, by default False
    outdir : str, optional
        Path where to save the files, by default ''
    """

    if type(data) == pyemma.coordinates.data.feature_reader.FeatureReader:
        data = data.get_output()

    data_concatenated = np.concatenate(data)
    fig, ax = pyemma.plots.plot_feature_histograms(
        data_concatenated, feature_labels=feat, ignore_dim_warning=True
    )

    if save:
        if outdir == "":
            raise Exception("Please provide a directory to save the file")
        else:
            plt.savefig(f"{outdir}/feat_hist.pdf", dpi=300, bbox_inches="tight")
    if display:
        plt.show()


def plot_density_energy(
        data,
        T,
        pairNames,
        save=False,
        display=False,
        outdir="",
        ij=(0,1)):
    """Plot the density map and the energy map

    Parameters
    ----------
    data : List
        List of all the value of the feat for each snapshot
    T : float
        Temperature
    pairNames : list
        List containing the name of the pairs
    save : bool, optional
        Save or not the plots, by default False
    display : bool, optional
        Display or not the plots, by default False
    outdir : str, optional
        Path where to save the files, by default ''
    ij : tuple, optional
        Index to project the representation, by default (0,1)
    """
    if type(data) == pyemma.coordinates.data.feature_reader.FeatureReader:
        data = data.get_output()
    data_concatenated = np.concatenate(data)
    
    # Index for the projection
    i,j=ij
    
    # Density plot
    fig, axes = plt.subplots(1, 1, sharex=True, sharey=True)
    
    pyemma.plots.plot_density(*data_concatenated.T[[i,j]], ax=axes)

    axes.set_xlabel(f" {pairNames[i]} (nm)")
    axes.set_ylabel(f" {pairNames[j]} (nm)")

    if save:
        if outdir == "":
            raise Exception("Please provide a directory to save the file")
        else:
            plt.savefig(f"{outdir}/data_density.pdf", dpi=300, bbox_inches="tight")
    if display:
        plt.show()

    # Energy plot
    fig, axes = plt.subplots(1, 1, sharex=True, sharey=True)
    kT = tools.get_kT(T)
    fig, axes = pyemma.plots.plot_free_energy(
        *data_concatenated.T[[i,j]], ax=axes, kT=kT, cbar_label="free energy / kJ.mol-1"
    )

    axes.set_xlabel(f" {pairNames[i]} (nm)")
    axes.set_ylabel(f" {pairNames[j]} (nm)")

    if save:
        if outdir == "":
            raise Exception("Please provide a directory to save the file")
        else:
            plt.savefig(
                f"{outdir}/data_free_energy_direct_from_density.pdf",
                dpi=300,
                bbox_inches="tight",
            )
    if display:
        plt.show()


def vamp_score(data, dim):
    """Compute the vamp score of the feat

    Parameters
    ----------
    data : List
        List of all the value of the feat for each snapshot
    dim : int
        Dimension of the VAMP reduction

    Returns
    -------
    score : float
        VAMP2 score
    """

    score = pyemma.coordinates.vamp(data[:-1], dim=dim).score(
        test_data=data[-1], score_method="VAMP2"
    )

    print("VAMP2 score: {:.2f}".format(score))

    return score

def score_cv(data, dim, lag, number_of_splits=10, validation_fraction=0.5):
    """Compute a cross-validated VAMP2 score.
    
    We randomly split the list of independent trajectories into
    a training and a validation set, compute the VAMP2 score,
    and repeat this process several times.
    
    Parameters
    ----------
    data : list of numpy.ndarrays
        The input data.
    dim : int
        Number of processes to score; equivalent to the dimension
        after projecting the data with VAMP2.
    lag : int
        Lag time for the VAMP2 scoring.
    number_of_splits : int, optional, default=10
        How often do we repeat the splitting and score calculation.
    validation_fraction : int, optional, default=0.5
        Fraction of trajectories which should go into the validation
        set during a split.
    """
    nval = int(len(data) * validation_fraction)
    scores = np.zeros(number_of_splits)
    for n in range(number_of_splits):
        ival = np.random.choice(len(data), size=nval, replace=False)
        vamp = pyemma.coordinates.vamp(
            [d for i, d in enumerate(data) if i not in ival], lag=lag, dim=dim)
        scores[n] = vamp.score([d for i, d in enumerate(data) if i in ival])
    return scores


if __name__ == "__main__":
    # Path
    pdb = "/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb"
    traj = [
        "/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS01_md_all_fitBB_protonly.xtc",
        "/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS02_md_all_fitBB_protonly.xtc",
    ]
    outdir = "/home/ccattin/Documents/Code/outputs"
    # Feat
    pairNames = ["64_CA-130_CA", "119_CA-24_CA","123_CA-24_CA","119_CA-27_CA"]
    # Parameters
    save = True
    display = True
    T = 300
    lag = 1000
    ij = (1,3)

    pair_indices = tools.create_pairIndices_from_pairNames(pdb, pairNames)
    feat = create_feat(pdb)
    feat = feat_atom_distances(feat,pair_indices)
    data = load_data(traj, feat, stride=5, ram=True)

    plot_feat_hist(data, feat, display=display, save=save, outdir=outdir)
    plot_density_energy(
        data=data, T=T, pairNames=pairNames, display=display, save=save, outdir=outdir,ij=ij
    )

    print(vamp_score(data=data, dim=2))
    print(offset(pair_indices,feat))
