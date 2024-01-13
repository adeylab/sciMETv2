import h5py
import sys

def rename_datasets(input_file, output_file, prefix):
    with h5py.File(input_file, 'r') as file_in, h5py.File(output_file, 'w') as file_out:
        def recursive_copy(source, destination):
            for name, item in source.items():
                if isinstance(item, h5py.Group):
                    group = destination.create_group(name)
                    recursive_copy(item, group)
                elif isinstance(item, h5py.Dataset):
                    new_name = prefix + '_' + name
                    destination.copy(item, new_name)

        recursive_copy(file_in, file_out)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input_file output_file prefix")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    prefix = sys.argv[3]

    rename_datasets(input_file, output_file, prefix)
