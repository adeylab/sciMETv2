import h5py
import os
import numpy as np
import argparse

# to run, type in the command line:
# python /home/groups/ravnica/src/sciMET/sciMET_cellCalls2h5.py [input folder path] [output prefix]

# Define and parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_folder")
parser.add_argument("output_prefix")
args = parser.parse_args()

# Input folder and output file paths
folder_path = args.input_folder
output_file = args.output_prefix + ".h5"

with h5py.File(output_file, 'w') as hdf5_file:
    cg_group = hdf5_file.create_group("CG")
    ch_group = hdf5_file.create_group("CH")

    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)

        if file_name.endswith(".CG.cov") or file_name.endswith(".CH.cov"):
            cell_id = file_name.replace('.CG.cov', '').replace('.CH.cov', '')

            cov = np.genfromtxt(file_path, delimiter='\t', dtype=[('chr', 'S10'), ('pos', int), ('pct', float), ('t', int), ('c', int)])
            cov = np.sort(cov, order=['chr', 'pos'])

            if file_name.endswith(".CG.cov"):
                cg_group.create_dataset(cell_id, data=cov, compression='gzip', compression_opts=9)

            elif file_name.endswith(".CH.cov"):
                ch_group.create_dataset(cell_id, data=cov, compression='gzip', compression_opts=9)

print("Files combined and written to HDF5.")






