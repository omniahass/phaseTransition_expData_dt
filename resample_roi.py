#!/usr/bin/env python3
# Resample ROI from MPRAGE to DWI- first 

import os
import argparse
import SimpleITK as sitk


def main():
    parser = argparse.ArgumentParser(
        description="Resample a segmentation to match the DWI image grid."
    )

    parser.add_argument(
        "-parent_folder",
        required=True,
        help="Path to the subject folder (e.g. /path/to/WCMyofascial2075)",
    )
    parser.add_argument(
        "-segmentation_filename",
        required=True,
        help="Segmentation filename (e.g. WCMyofascial2075.nii.gz)",
    )

    args = parser.parse_args()

    parent_folder = args.parent_folder.rstrip("/")

    # Accept either a filename inside the subject folder or a full path to the
    # segmentation, so ROIs kept in a shared folder can be used without copying
    if os.path.isabs(args.segmentation_filename):
        seg_filepath = args.segmentation_filename
    else:
        seg_filepath = os.path.join(parent_folder, args.segmentation_filename)
    if not os.path.exists(seg_filepath):
        raise SystemExit(f"Segmentation not found: {seg_filepath}")
    roi = sitk.ReadImage(seg_filepath)

    # Load DWI
    nifti_image = os.path.join(
        parent_folder,
        "derivatives",
        "all",
        "0022ms",
        "dwiec.nii",
    )
    if not os.path.exists(nifti_image):
        raise SystemExit(f"DWI not found, run the designer step first: {nifti_image}")
    ref_img = sitk.ReadImage(nifti_image)

    # If DWI is 4D, extract the first volume
    if ref_img.GetDimension() == 4:
        ref_img = ref_img[:, :, :, 0]

    # Resample ROI to reference image
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(ref_img)
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    resampler.SetTransform(sitk.Transform())  # Identity transform
    resampler.SetDefaultPixelValue(0)

    roi_resampled = resampler.Execute(roi)

    # Save output - the name is built from the subject folder, because that is how
    # rpbm_bayesian_singlecase.m and read_roi_values_singlecase.ipynb look it up
    subject = os.path.basename(parent_folder)
    output_filename = subject + "_reformatdwi.nii.gz"
    output_file = os.path.join(parent_folder, output_filename)
    sitk.WriteImage(roi_resampled, output_file)

    print(f"Saved resampled ROI to:\n{output_file}")


if __name__ == "__main__":
    main()