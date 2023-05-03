#! /usr/bin/env python3
"""LETHE automatizes the MSM analysis"""

import argparse

def main():
    """
    CLI main function
    """
    
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
    parser.add_argument(
        '-d',
        '--distances',
        nargs='+',
        help='List of the pair names'
    )
    
    args = parser.parse_args()

    #Handles error
    if not args.files:
        parser.error("No trajectories file given")
    if not args.distances:
        parser.error("No distances given")

    file_list = args.files
    pairNames = args.distances
    print(file_list)
    print(pairNames)

if __name__ == '__main__':
    main()