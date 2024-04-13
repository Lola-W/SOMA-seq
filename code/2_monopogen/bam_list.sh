#!/bin/bash

# Usage: 
# bash bam_list.sh <source_directory> <destination_directory>
#
# <source_directory>: Directory containing subdirectories for BAM files.
# <destination_directory>: Directory where .bam.lst files will be saved.
#
# Example:
# bash bam_list.sh /home/wengjiam/scratch/soma_seq/Monopogen/AIBS/human-var-gru /home/wengjiam/scratch/soma_seq/Monopogen/bamlsts

if [ "$#" -ne 2 ]; then
    echo "Incorrect usage!"
    echo "Correct usage: bash bam_list.sh <source_directory> <destination_directory>"
    exit 1
fi

source_dir="$1"
dest_dir="$2"

# Ensure the destination directory exists
mkdir -p "$dest_dir"

# Process each directory within the source directory
for dir in "$source_dir"/*; do
    if [ -d "$dir" ]; then
        # Extract the directory name
        dir_name=$(basename "$dir")

        # Create and overwrite the .bam.lst file with the specified content
        echo "${dir_name},$source_dir/STARsolo/STAR_results/${dir_name}/Aligned.sortedByCoord.out.bam" > "${dest_dir}/${dir_name}.bam.lst"
    fi
done

echo "Generated .bam.lst files in: $dest_dir"
