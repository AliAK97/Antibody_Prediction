## part 1: separate .cif and .json files

def separate_files(pdb_id):
    counter_cif = 0
    counter_json = 0
    base_path = f'/Users/ali/Desktop/Reincke_Lab/Antigen_Prediction/AF3/HPC_files/output_files/raw_files/new_predictions/{pdb_id}'
    dest_path = '/Users/ali/Desktop/Reincke_Lab/Antigen_Prediction/AF3/HPC_files/output_files/sorted_output/new_sorted'
    for first_level_name in os.listdir(base_path):
        first_level_path = os.path.join(base_path, first_level_name)
        if os.path.isdir(first_level_path):
            #print(f"Processing directory: {first_level_name}")
            for file in os.listdir(first_level_path):
                if file.endswith('.cif'):
                    src_file = os.path.join(first_level_path, file)
                    dst_file = os.path.join(f'{dest_path}/{pdb_id}_cif', file)
                    shutil.copy2(src_file, dst_file)
                    counter_cif += 1
                elif file.endswith('summary_confidences.json'):
                    src_file = os.path.join(first_level_path, file)
                    dst_file = os.path.join(f'{dest_path}/{pdb_id}_json', file)
                    shutil.copy2(src_file, dst_file)
                    counter_json += 1
    print(f"Total .cif files copied: {counter_cif}")
    print(f"Total .json files copied: {counter_json}")


pdb_id = '8j7e'
separate_files(pdb_id)

## part 2: rename files in each folder
def rename_files(pdb_id):
    cif_path = f'/Users/ali/Desktop/Reincke_Lab/Antigen_Prediction/AF3/HPC_files/output_files/sorted_output/new_sorted/{pdb_id}_cif'
    json_path = f'/Users/ali/Desktop/Reincke_Lab/Antigen_Prediction/AF3/HPC_files/output_files/sorted_output/new_sorted/{pdb_id}_json'

    for i, file in enumerate(os.listdir(cif_path)):
        if file.endswith('.cif'):
            new_name = f'{pdb_id}_{i}.cif'
            os.rename(os.path.join(cif_path, file), os.path.join(cif_path, new_name))

    for i, file in enumerate(os.listdir(json_path)):
        if file.endswith('.json'):
            new_name = f'{pdb_id}_{i}.json'
            os.rename(os.path.join(json_path, file), os.path.join(json_path, new_name))


rename_files(pdb_id)