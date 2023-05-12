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

    #Init
    parser = argparse.ArgumentParser(
        "CLI LETHE program",
        description="Automatize MSM building using PyEmma.",
    )


    # List of files
    parser.add_argument(
        '-f',
        '--files',
        required=True,
        nargs='+',
        help='List of files'
        )
    # Distances to add to the feat
    parser.add_argument(
        '-d',
        '--distances',
        nargs='+',
        help='List of the pair names'
    )
    # Topology .pdb file
    parser.add_argument(
        '-t',
        '--topology',
        required=True,
        help='Reference .pdb file'
    )
    # What properties to plot
    parser.add_argument(
        '-p',
        '--plot',
        nargs='+',
        help='Plot wanted (feat_hist, density_energy, its, cluster, cktest, stationary, eigenvectors, metastable_membership, mfpt, committor)',
        choices=[
             'feat_hist','density_energy','its','cluster','cktest','stationary','eigenvectors','metastable_membership','mfpt','committor'
             ],
        required=True
    )
    # Do not display the plot 
    parser.add_argument(
        '--no-plot',
        action='store_true',
        help='Do not display plots'
    )
    # Output directory to save plots
    parser.add_argument(
        '-o',
        '--outdir',
        help='Path to save the plots'
    )
    # Temperature in K of the system
    parser.add_argument(
        '--T',
        '--temperature',
        type=float,
        help='Temperature of the system'
    )
    # Dimension reduction
    parser.add_argument(
        '--reduction',
        type=str,
        choices=['pca','tica','none'],
        help='Do a PCA or TICA dimension reduction. Select PCA, TICA or none',
        default='none'
    )
    parser.add_argument(
         '--tica-lag',
         type=int,
         help='Lag used for the TICA dimension reduction'
    )
    # Lag time
    parser.add_argument(
        '--lag',
        help='Lag time for the MSM',
        type=int,
        required=True
    )
    # Clustering method
    parser.add_argument(
        '--cluster',
        choices=['kmeans', 'regspace'],
        help='Clustering method (kmeans or regspace)',
        type=str
    )
    # Number of cluster
    parser.add_argument(
        '-k',
        '--cluster-number',
        help='Number of cluster wanted',
        type=int
    )
    # Number of stride
    parser.add_argument(
        '--stride',
        help='Number of stride',
        type=int
    )
    # Perform ITS analysis
    parser.add_argument(
        '--its',
        nargs='+',
        type=int,
        help='Perform ITS analysis on the following lags',
    )
    # Number of iteration for ITS analysis
    parser.add_argument(
        '--nits',
        type=int,
        help='Number of iteration for the ITS validation'
    )
    # ITS analysis on different cluster size
    parser.add_argument(
        '--its-cluster',
        nargs='+',
        type=int,
        help='Perform a cluster size analysis using the ITS analysis'
    )
    # CK Test
    parser.add_argument(
        '--cktest',
        type=bool,
        choices=[True,False],
        help='Perform a CK test on the MSM builded given the number of metastable states. True for having the error, False for not'
    )
    # Meta-stable states
    parser.add_argument(
        '--state',
        type=int,
        help='Number of state to consider in the MSM'
    )
    # 95% confidence interval
    parser.add_argument(
        '--confidence',
        action='store_true',
        help='Create the MSM with a Bayesian estimation of the error. This is presented in a 95 confidence interval'
    )
    # Save the model
    parser.add_argument(
        '--save',
        nargs=2,
        type=str,
        help='File name and model name to save the clustering and MSM'
    )
    # Load the model
    parser.add_argument(
        '--load',
        nargs=2,
        type=str,
        help='Load previous model. File name (.pyemma) and then model name.'
    )
    # PCCA analysis, display the stationary proba
    parser.add_argument(
        '--pcca',
        action='store_true',
        help='PCCA and TPT analysis. Display the stationary probabilities'
    )
    parser.add_argument(
        '--state-path',
        nargs=2,
        type=int,
        help='Get the TPT between the two given states'
    )

    #Get the arguments
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
    if not args.distances:
        parser.error("No distances given")

    # Forgot to give the temperature
    if 'density_energy' in args.plot:
            if not args.T:
                parser.error('Please provide a temperature')
    
    if args.reduction == 'tica' and not args.tica_lag:
         parser.error('Please provide a TICA lag')

    if args.cluster and not args.cluster_number:
         parser.error('Please provide a number of cluster')
    
    # ITS analysis error
    if args.its:
        if not args.nits:
            parser.error("Please provide a number of iteration")
    
    if not args.state:
         if args.cktest or args.pcca:
              parser.error('Please provide a number of meta-stable state')