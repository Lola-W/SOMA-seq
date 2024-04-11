"""
generateBarcode.py

Description:
This script processes single-cell RNA-seq data to calculate the total counts per barcode
from a matrix file from STARSolo and outputs a CSV file with these counts. This helper script is
is not necessarily required for CellRanger alignment.

Usage:
python utils/generateBarcode.py <directory_path>

Arguments:
<directory_path>: The directory containing STARSolo alignment output, it is where 'matrix.mtx' and 'barcodes.tsv' are located,
                and 'barcode_counts.csv' will be saved.

Example:
python generateBarcode.py /path/to/alignment_output

This will read the matrix and barcodes from the specified STARSolo alignment output directory, calculate the counts per barcode,
and save the results in the same directory as 'barcode_counts.csv'.
"""

import sys
from scipy.io import mmread
import pandas as pd

def generate_barcode_counts(matrix_path, barcodes_path, output_path):
    """
    Generates a CSV file containing the sum of counts for each barcode.

    Parameters:
    - matrix_path (str): Path to the .mtx file containing the matrix.
    - barcodes_path (str): Path to the .tsv file containing the barcodes.
    - output_path (str): Path where the barcode counts CSV will be saved.
    """
    # Load the matrix and barcodes
    matrix = mmread(matrix_path).tocsc()
    barcodes = pd.read_csv(barcodes_path, header=None)[0]

    # Sum the counts for each barcode
    counts_per_barcode = matrix.sum(axis=0).A1  # Convert to 1D numpy array

    # Create a DataFrame
    barcode_counts = pd.DataFrame({
        'cell': barcodes,
        'id': counts_per_barcode
    })

    # Save to file
    barcode_counts.to_csv(output_path, index=False)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generateBarcode.py <directory_path>")
        sys.exit(1)

    directory_path = sys.argv[1]
    matrix_file = f'{directory_path}/matrix.mtx'
    barcodes_file = f'{directory_path}/barcodes.tsv'
    output_file = f'{directory_path}/barcode_counts.csv'

    generate_barcode_counts(matrix_file, barcodes_file, output_file)
    print(f'Barcode counts file saved to: {output_file}')
