import sys
import h5py

def merge_hdf5_files(file_list, output_file):
    with h5py.File(output_file, 'w') as merged_h5:
        for file_name in file_list:
            with h5py.File(file_name, 'r') as file:
                for group_name in file.keys():
                    group = file[group_name]
                    if group_name not in merged_h5:
                        merged_group = merged_h5.create_group(group_name)
                    else:
                        merged_group = merged_h5[group_name]
                    # Copy datasets from source file to the merged file
                    for dataset_name, dataset in group.items():
                        if dataset_name not in merged_group:
                            merged_group.create_dataset(dataset_name, data=dataset[...],
                                                        compression="gzip", compression_opts=9)
                        else:
                            print(f"Dataset '{dataset_name}' already exists in '{group_name}' group.")
                            # Handle if the dataset already exists in the destination file

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python merge_hdf5.py output_file file1.h5 file2.h5 ...")
        sys.exit(1)
    
    output_file = sys.argv[1]
    file_list = sys.argv[2:]
    
    merge_hdf5_files(file_list, output_file)

