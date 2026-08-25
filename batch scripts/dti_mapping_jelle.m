% DTI mapping using Jelle DTI.m class

clear all
close all
clc

addpath('./supporting_functions_dti/');

parent_directory = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/';

subjects = dir(parent_directory);
subjects = {subjects.name};
% subjects = subjects(contains(subjects, 'WCMyofascial20'));
subjects = subjects(contains(subjects, 'C2_'));

diffusion_times = ["0022ms","0042ms","0081ms","0156ms","0300ms"];

%%

for k = 24:length(subjects)
    subject = subjects{k};
    % subject='WCMyofascial2059';
    disp(subject)
    
    % parent_folder = fullfile(parent_directory, subject,'/rescan/');
    parent_folder = fullfile(parent_directory, subject,'/');
    for i = 1:length(diffusion_times)
        diffusion_time = diffusion_times(i)
        nifti_folder = strcat(parent_folder,'derivatives/all/',diffusion_time,'/');

        if isfolder(nifti_folder)
            
            % Read image
            img = niftiread(strcat(nifti_folder,'dwiec.nii')); % should this be load_untouch??
            info = niftiinfo(strcat(nifti_folder,'dwiec.nii'));
            info.ImageSize = info.ImageSize(1:3);
            info.PixelDimensions = info.PixelDimensions(1:3);
            info.Datatype = 'double';
            
            % Read b matrix and bvecs and bvals
            bmat_folder = strcat(parent_folder,'bmatrix/');
            bmatrix = load(strcat(bmat_folder,diffusion_time,'.txt'));
        
            % Remove the b0 at the beginning and end
            img = img(:,:,:,2:end-1);
            bmatrix = bmatrix(2:end-1,:);
            
            % mapping using b matrix
            dti = DTI(bmatrix, 'wlls');
            outliers = dti.outliers(img);
            dt = dti.dt(img);
            eig = dti.eig(dt);
            ad = dti.ad(eig);
            rd = dti.rd(eig);
            l2 = dti.l2(eig);
            l3 = dti.l3(eig);
            fa = dti.fa(eig);
        
            maps_folder = strcat(parent_folder,"maps");
            if ~exist(maps_folder,"dir") % make map folder if it doesn't exist
                mkdir(maps_folder)
            end
            niftiwrite(clip(1000*abs(rd),0,3.5),strcat(maps_folder, "/", diffusion_time, "_RD.nii.gz"),info);
            niftiwrite(clip(1000*abs(ad),0,3.5),strcat(maps_folder, "/", diffusion_time, "_AD.nii.gz"),info);
            niftiwrite(clip(1000*abs(l2),0,3.5),strcat(maps_folder, "/", diffusion_time, "_L2.nii.gz"),info);
            niftiwrite(clip(1000*abs(l3),0,3.5),strcat(maps_folder, "/", diffusion_time, "_L3.nii.gz"),info);
            niftiwrite(fa,strcat(maps_folder, "/", diffusion_time, "_FA.nii.gz"),info);
        end
    end
end