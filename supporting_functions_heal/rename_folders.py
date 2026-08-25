import os
import pydicom
import shutil
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_folder', type=str, help='Path to the input directory', required=True)
    args = parser.parse_args()
    input_folder = args.input_folder

    dicom_folder = os.path.join(input_folder, "sorted_dcms")
    if not os.path.exists(dicom_folder):
        os.makedirs(dicom_folder)

    # Move all subfolders into sorted_dcms (except sorted_dcms itself)
    for item in os.listdir(input_folder):
        src = os.path.join(input_folder, item)
        dst = os.path.join(dicom_folder, item)
        # CHANGED FOR WRAPPER: only move DICOM series folders. Moving everything
        # swallowed pipeline folders (logs, nifti, derivatives, maps) on a re-run
        if os.path.isdir(src) and src != dicom_folder and "SER" in item:
            if not os.path.exists(dst):
                shutil.move(src, dicom_folder)

    subfolders = os.listdir(dicom_folder)
    subfolders = [s for s in subfolders if "SER" in s]
    subfolders.sort()
    
    for folder in subfolders:
        current_folder = os.path.join(dicom_folder, folder)
        dicoms = [d for d in os.listdir(current_folder) if not d.startswith('.')]
        if not dicoms:
            continue

        # Read first DICOM in folder
        temp = pydicom.dcmread(os.path.join(current_folder,dicoms[0]))
        series = temp[0x008,0x103e].value

        # CHANGED FOR WRAPPER: skip folders that already carry the series description,
        # otherwise re-running the pipeline prefixes the name a second time
        if folder.startswith(series+'_'):
            print('already renamed, skipping: '+folder)
            continue

        new_folder_name = os.path.join(dicom_folder,series+'_'+folder)
        print(new_folder_name)
        shutil.move(current_folder,new_folder_name)

if __name__ == "__main__":
    main()