#!/bin/bash

# Input is a folder of sorted dicoms
# CHANGED FOR WRAPPER: parent_path / segmentation_filename / center come from the
# arguments when called by run_heal_pipeline.sh, otherwise the defaults below are used
parent_path="${1:-/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial2075}"
segmentation_filename="${2:-WCMyofascial2075.nii.gz}"
center="${3:-cornell}"
rescan_rpbm="n"
rescan_nonrpbm="n"
python_bin="${python_bin:-python3}" # CHANGED FOR WRAPPER: run_heal_pipeline.sh exports this

### ------------------------------------------------

# # Before this use dicom sorting application e.g DICOM sort
dcm_path=$parent_path/sorted_dcms

# # Cornell- need to rename folders so they have series description
# CHANGED FOR WRAPPER: rename_folders.py moves EVERY sub-folder of parent_path into
# sorted_dcms, so it must not run again once nifti/ or derivatives/ exist
if [ "$center" = "cornell" ]; then
    if [ -d "$parent_path/nifti" ] || [ -d "$parent_path/derivatives" ]; then
        echo "Skipping rename_folders.py - $parent_path already has nifti/ or derivatives/"
    else
        $python_bin ./supporting_functions_heal/rename_folders.py -input_folder $parent_path
    fi
fi

# # Remove superflous folders and files
find "$dcm_path" -depth -path '*scout*' -delete
find "$dcm_path" -depth -path '*_ADC_*' -delete
find "$dcm_path" -depth -path '*_TENSOR*' -delete
find "$dcm_path" -depth -path '*TRACEW_*' -delete
find "$dcm_path" -depth -path '*TRACEW_*' -delete
find "$dcm_path" -depth -path '*FA_*' -delete
rm -f $parent_path/DICOMDIR
rm -f $parent_path/LOG.txt
rm -f $parent_path/README.txt

# move retro recon to another folder
mkdir -p "$parent_path/retro_recon"
mkdir -p "$parent_path/retro_recon/sorted_dcms"
find "$dcm_path" -type d -name '*_RR_*' -exec mv {} "$parent_path/retro_recon/sorted_dcms" \;
# move rescan to another folder
mkdir -p "$parent_path/rescan"
mkdir -p "$parent_path/rescan/sorted_dcms"
find "$dcm_path" -type d -name '*_rescan_*' -exec mv {} "$parent_path/rescan/sorted_dcms" \;

# # # Convert to nifti
nifti_path="$parent_path/nifti"
mkdir -p "$nifti_path"
dcm2niix -f "%e_%p_%d" -m y -px y -z y -o "$nifti_path" "$dcm_path"

# # B matrix and bvec extraction
$python_bin bmatrix_export.py -input_folder $parent_path -sequence STEAM -center $center
$python_bin bmatrix_export.py -input_folder $parent_path -sequence PGSE -center $center
$python_bin bmatrix_export.py -input_folder $parent_path -sequence OGSE -center $center

# CHANGED FOR WRAPPER: designer_heal.sh and resample_roi.py are now run as their own
# steps by run_heal_pipeline.sh (step 2 and step 3) so they can be re-run on their own.
# Uncomment the block below to go back to running everything from this script.
# # # Run heal designer script
# chmod +x designer_heal.sh # so always have permissions
# ./designer_heal.sh "$parent_path"
#
# python3 ./resample_roi.py \
#     -parent_folder "$parent_path" \
#     -segmentation_filename "$segmentation_filename"

# # Run T2 Mapping script?
# # Find and convert CPMG sequences from dcm to nifti --> seperate by echo files for t2 mapping
# for cpmg_folder in "$dcm_path"/*CPMG*; do
#     if [ -d "$cpmg_folder" ]; then
#         echo "Converting CPMG sequence: $(basename "$cpmg_folder")"
#         dcm2niix -z y -f "%p_%s" -o "$nifti_path" "$cpmg_folder"
#     fi
# done

# python3 ./t2_mapping.py "$parent_path"


# ---- RESCAN -----
# Convert to nifti
# parent_path_rescan="$parent_path/rescan"
# dcm_path_rescan="$parent_path_rescan/sorted_dcms"
# nifti_path_rescan="$parent_path_rescan/nifti"
# mkdir -p "$nifti_path_rescan"
# dcm2niix -f "%e_%p_%d" -m y -px y -z y -o "$nifti_path_rescan" "$dcm_path_rescan"

# if [ "$rescan_rpbm" = "y" ]; then
#     # # # B matrix and bvec extraction
#     $python_bin bmatrix_export.py -input_folder $parent_path_rescan -sequence STEAM -center $center
#     $python_bin bmatrix_export.py -input_folder $parent_path_rescan -sequence PGSE -center $center
#     $python_bin bmatrix_export.py -input_folder $parent_path_rescan -sequence OGSE -center $center

#     # run designer processing on rescan
#     # chmod +x designer_heal_rescan.sh # so always have permissions
#     # ./designer_heal_rescan.sh "$parent_path_rescan"
# fi


