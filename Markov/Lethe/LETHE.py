"""LETHE automatizes the MSM analysis"""

import argparse

def main():
    """
    CLI main function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--files', nargs='+', help='List of files')
    args = parser.parse_args()

    file_list = args.files
    print(file_list)

if __name__ == '__main__':
    main()