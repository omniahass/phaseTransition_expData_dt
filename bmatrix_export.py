# b matrix and b vector extraction script
#
# For a given parent folder containing sub-folders of sorted dicoms
# - extracts the b matrix
# - extracts the scanner provided b vectors
# - performs eigenvalue decomposition for the b matrix and saves vectors 
#
#
# Gabrielle Baxter 09/18/2025 gabrielle.baxter@nyulangone.org
#
# Required inputs:
# -parent_folder
# -sequence_type
# -center (default is cornell)

import os
import pydicom
import numpy as np
import re
import argparse

def main():
        parser = argparse.ArgumentParser()

        # Add a flag-style argument
        parser.add_argument('-input_folder', type=str, help='Path to the input directory', required=True)
        parser.add_argument('-sequence',type=str,help='Type of sequence',default='STEAM')
        parser.add_argument('-center',type=str,help='cornell or nyu',default='cornell')
        args = parser.parse_args()

        parent_folder = args.input_folder
        parent_dcm_folder = parent_folder + '/sorted_dcms'

        dcm_folders = os.listdir(parent_dcm_folder)
        dcm_folders.sort()

        # Get just the standard images for different diffusion times
        sequence_type = args.sequence
        print(sequence_type)
        dt_folders = [s for s in dcm_folders if sequence_type in s] # get only STEAM sequences
        dt_folders = [s for s in dt_folders if not 'RR' in s] # remove all retro recon phase folders
        dt_folders = [s for s in dt_folders if not 'b0' in s] # remove all reverse phase encoding folders
        # dt_folders = [s for s in dt_folders if not '_rescan' in s] # remove all reverse phase encoding folders
        dt_folders.sort()
        print(dt_folders)

        center = args.center

        # Create folders for files
        bmat_path = os.path.join(parent_folder,'bmatrix')
        if not os.path.exists(bmat_path):
                os.makedirs(bmat_path)

        # eigen_path = os.path.join(parent_folder,'vectors_eigendecomp')
        # if not os.path.exists(eigen_path):
        #         os.makedirs(eigen_path)

        # header_path = os.path.join(parent_folder,'vectors_header')
        # if not os.path.exists(header_path):
        #         os.makedirs(header_path)

        # Iterate through diffusion times
        for dt_folder in dt_folders:
                print(dt_folder)
                if sequence_type == 'OGSE':
                        match = re.search(r'(.{4})Hz', dt_folder)
                else:
                        match = re.search(r'(.{4})ms', dt_folder) # Get the diffusion time from the folder name

                diffusion_time = match.group(1)
                dt_folder_path = os.path.join(parent_dcm_folder,dt_folder)
                dcm_files = os.listdir(dt_folder_path) # Get dicom files corresponding to diffusion time
                dcm_files.sort()

                if dcm_files[0] == '1.dcm':
                        dcm_files = os.listdir(dt_folder_path) 
                        dcm_files.sort(key=lambda x: int(os.path.splitext(x)[0])) 

                bmat_array = []
                eig_dirs = np.zeros([len(dcm_files),3])
                header_dirs = np.zeros([len(dcm_files),3])
                # New bmat and vector set for each dicom (each direction)
                for i in range(0,len(dcm_files)):
                        dcm_info = pydicom.dcmread(os.path.join(dt_folder_path,dcm_files[i])) # read dicom file

                        if center == 'nyu':
                        # Get bmatrix 
                                matrix = dcm_info[0x5200, 0x9230][1][0x0018,0x9117][0][0x0018, 0x9601]
                                xx = matrix[0][0x018,0x9602].value
                                xy = matrix[0][0x018,0x9603].value
                                xz = matrix[0][0x018,0x9604].value
                                yy = matrix[0][0x018,0x9605].value
                                yz = matrix[0][0x018,0x9606].value
                                zz = matrix[0][0x018,0x9607].value

                                # Get directions from header
                                directions = dcm_info[0x5200, 0x9230][1][0x0018,0x9117][0][0x0018, 0x9076]
                                directions = directions[0][0x0018,0x9089].value
                                header_dirs[i] = directions

                        if center == 'cornell':
                                # Get bmatrix
                                matrix = dcm_info[0x19, 0x1027].value
                                xx = matrix[0]
                                xy = matrix[1]
                                xz = matrix[2]
                                yy = matrix[3]
                                yz = matrix[4]
                                zz = matrix[5]

                                # Get directions from header
                                directions = dcm_info[0x0019,0x100e]
                                directions = directions.value
                                header_dirs[i] = directions

                        bmat = np.array([xx,xy,xz,yy,yz,zz]) # This is the standard format 
                        bmat_array.append(bmat)

                        # Reconstruct 3x3 symmetric B matrix
                        B = np.array([[xx, xy, xz],
                                        [xy, yy, yz],
                                        [xz, yz, zz]], dtype=float)
                        
                        # Eigen-decomposition
                        eigvals, eigvecs = np.linalg.eigh(B)
                        idx = np.argmax(eigvals) # Get index of max eigenvalue
                        bval = eigvals[idx] # Corresponding b val
                        bvec = eigvecs[:, idx] # Corresponding b vecs
                        bvec = bvec / np.linalg.norm(bvec)
                        eig_dirs[i] = bvec

                # Choose suffix based on sequence type
                if sequence_type == "OGSE":
                        unit_suffix = "hz"
                else:
                        unit_suffix = "ms"

                # Save bmat
                bmat_array = np.array(bmat_array)
                save_filepath = os.path.join(bmat_path, f"{diffusion_time}{unit_suffix}.txt")
                np.savetxt(save_filepath, bmat_array, fmt="%.1f")

if __name__ == "__main__":
    main()
