#! /usr/bin/env python3
"""LETHE's parser"""

import argparse


def parsing():
    """Parsing of the option given

    Returns
    -------
    parser : argparse.parser
        Parser used
    args : argparse.args
        Argument given by the user
    """

    # Init
    parser = argparse.ArgumentParser(
        "CLI LETHE program",
        description="Automatize MSM building using PyEmma.",
    )

    # List of files
    parser.add_argument(
        "-f",
        "--files",
        required=True,
        nargs="+",
        help="List of files"
    )
    # Distances to add to the feat
    parser.add_argument(
        "-d",
        "--distances",
        nargs="+",
        help="List of the pair names",
    )
    parser.add_argument(
        "--residue",
        nargs="+",
        help="Shortest distance between residue indices to consider. Pairs of residue are separated by a space and '-' separate the two residues inside a pair"
    )
    # Compute the VAMP2 score
    parser.add_argument(
        "--vamp-score",
        action="store_true",
        help="Compute the VAMP score of the features",
    )
    # Topology .pdb file
    parser.add_argument("-t", "--topology", required=True, help="Reference .pdb file")
    # Load into RAM
    parser.add_argument(
        "--ram",
        action="store_true",
        help="Store the data into the RAM. More efficient but maybe not possible for large data set",
    )
    # What properties to plot
    parser.add_argument(
        "-p",
        "--plot",
        nargs="+",
        help="Plot wanted (feat_hist, density_energy, vamp_lag_dim, pca, tica, vamp, its, cluster, cktest, stationary, eigenvectors, metastable_membership, mfpt, committor)",
        choices=[
            "feat_hist",
            "density_energy",
            "vamp_lag_dim",
            "pca",
            "tica",
            "vamp",
            "its",
            "cluster",
            "cktest",
            "stationary",
            "eigenvectors",
            "metastable_membership",
            "mfpt",
            "committor",
        ],
        required=True,
    )
    # Do not display the plot
    parser.add_argument("--no-plot", action="store_true", help="Do not display plots")
    # Output directory to save plots
    parser.add_argument("-o", "--outdir", help="Path to save the plots")
    # Temperature in K of the system
    parser.add_argument(
        "--T", "--temperature", type=float, help="Temperature of the system"
    )
    # Get the lag time and the dimension for dimension reduction
    # Lag time
    parser.add_argument(
        "--vamp-lags",
        nargs="+",
        type=int,
        help="Dimension reduction lags to test"
    )
    # Dimension
    parser.add_argument(
        "--vamp-dim",
        type=int,
        help="Number of dimension to test the dimension reduction"
    )
    # Dimension reduction
    parser.add_argument(
        "--reduction",
        type=str,
        choices=["pca", "tica", "vamp", "none"],
        help="Do a PCA, TICA or VAMP dimension reduction. Select PCA, TICA, VAMP or none",
        default="none",
    )
    # Select the TICA lag time
    parser.add_argument(
        "--tica-lag", type=int, help="Lag used for the TICA dimension reduction"
    )
    # Select the VAMP lag time
    parser.add_argument(
        "--vamp-lag", type=int, help="Lag used for the VAMP dimension reduction"
    )
    # Select the number of dimension for the reduction
    parser.add_argument(
        "--dim", type=int, help="Dimension to project for the dimension reduction"
    )
    # Lag time
    parser.add_argument("--lag", help="Lag time for the MSM", type=int, required=True)
    # Clustering method
    parser.add_argument(
        "--cluster",
        choices=["kmeans", "regspace"],
        help="Clustering method (kmeans or regspace)",
        type=str,
    )
    # VAMP2 score as a function of the number of cluster
    parser.add_argument(
        "--vamp-cluster",
        nargs="+",
        type=int,
        help="Compute the VAMP2 score as a function of the number of cluster center. Put the number of cluster center after"
    )
    # Number of cluster
    parser.add_argument(
        "-k", "--cluster-number", help="Number of cluster wanted", type=int
    )
    # Number of stride
    parser.add_argument("--stride", help="Number of stride", type=int, default=1)
    # Perform ITS analysis
    parser.add_argument(
        "--its",
        nargs="+",
        type=int,
        help="Perform ITS analysis on the following lags",
    )
    # Number of iteration for ITS analysis
    parser.add_argument(
        "--nits", type=int, help="Number of iteration for the ITS validation"
    )
    # ITS analysis on different cluster size
    parser.add_argument(
        "--its-cluster",
        nargs="+",
        type=int,
        help="Perform a cluster size analysis using the ITS analysis",
    )
    # CK Test
    parser.add_argument(
        "--cktest",
        type=bool,
        choices=[True, False],
        help="Perform a CK test on the MSM builded given the number of metastable states. True for having the error, False for not",
    )
    # Meta-stable states
    parser.add_argument(
        "--state", type=int, help="Number of state to consider in the MSM"
    )
    # 95% confidence interval
    parser.add_argument(
        "--confidence",
        action="store_true",
        help="Create the MSM with a Bayesian estimation of the error. This is presented in a 95 confidence interval",
    )
    # Save the model
    parser.add_argument(
        "--save",
        nargs=2,
        type=str,
        help="File name and model name to save the clustering and MSM",
    )
    # Load the model
    parser.add_argument(
        "--load",
        nargs=2,
        type=str,
        help="Load previous model. File name (.pyemma) and then model name.",
    )
    # PCCA analysis, display the stationary proba
    parser.add_argument(
        "--pcca",
        action="store_true",
        help="PCCA and TPT analysis. Display the stationary probabilities",
    )
    # TPT between two states
    parser.add_argument(
        "--state-path",
        nargs=2,
        type=int,
        help="Get the TPT between the two given states",
    )
    # Save .pdb files of the meta stable states
    parser.add_argument(
        '--pdb',
        type=int,
        help="Generate .pdb files of the sampled meta stable states. The number of this option correspond to the number of frame to write in the .pdb files."
    )

    # Get the arguments
    args = parser.parse_args()

    return parser, args


def LETHE_handle_error(parser, args):
    """Get the standard error in LETHE

    Parameters
    ----------
    parser : argparse.parser
        Parser used
    args : argparse.args
        Argument given by the user
    """

    # No file given
    if not args.files:
        parser.error("No trajectories file given")

    # No distance given
    if not args.distances and not args.residue:
        parser.error("No distances given")

    # Forgot to give the temperature
    if "density_energy" in args.plot:
        if not args.T:
            parser.error("Please provide a temperature")

    # Forgot to give dimension
    if args.reduction and not args.dim:
        parser.error("Please provide a dimension for dimension reduction")
    if args.vamp_score and not args.dim:
        parser.error("Please provide a dimension for VAMP dimension reduction")

    if 'vamp_lag_dim' in args.plot and (not args.vamp_lags or not args.vamp_dim):
        parser.error("Please provide a correct maximum dimension and list of lags for VAMP2 score analysis")

    # Forgot to give the lag time in reduction
    if args.reduction == "tica" and not args.tica_lag:
        parser.error("Please provide a TICA lag")
    if args.reduction == "vamp" and not args.vamp_lag:
        parser.error("Please provide a VAMP lag")

    if args.cluster and not args.cluster_number:
        parser.error("Please provide a number of cluster")

    # ITS analysis error
    if args.its:
        if not args.nits:
            parser.error("Please provide a number of iteration")

    if not args.state:
        if args.cktest or args.pcca:
            parser.error("Please provide a number of meta-stable state")

    if args.state <= max(args.state_path):
        parser.error("This path does not exist. The state counting begin at 0. Please check your path")
