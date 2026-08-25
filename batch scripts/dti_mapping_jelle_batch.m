% DTI mapping using Jelle DTI.m class

clear all
close all
clc

addpath('./supporting_functions_dti/');

parent_directory = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/nyu/';
% parent_directory = '/Volumes/labspace/Projects/HEAL/cornell_rescan';
nifti_folder = fullfile(parent_directory,'dwi_smoothed');
bmat_folder = fullfile(parent_directory,'bmatrix');
maps_folder = fullfile(parent_directory,'maps_smoothed');

subjects = dir(nifti_folder);
subjects = {subjects.name};
% subjects = subjects(contains(subjects, 'WCMyofascial20'));
subjects = subjects(contains(subjects, 'C2_'));
subjects = subjects(~contains(subjects, '_reg'));
temp = subjects';


diffusion_times = ["0022ms","0042ms","0081ms","0156ms","0300ms"];

%%

for k = 1:length(subjects)
    subject = subjects{k};
    subject_code = extractBefore(subject, strlength(subject)-3);
    disp(subject_code)
    disp(fullfile(nifti_folder,subject))

    % Read image
    img = niftiread(fullfile(nifti_folder,subject)); 
    info = niftiinfo(fullfile(nifti_folder,subject)); 
    info.Datatype = 'double';

    % Read b matrix and bvecs and bvals
    bmatrix = load(fullfile(bmat_folder, [subject_code '.txt']));

    % Remove the b0 at the beginning and end
    img = img(:,:,:,2:end-1);
    bmatrix = bmatrix(2:end-1,:);

    % mapping using b matrix
    dti = DTI(bmatrix, 'wlls');
    %outliers = dti.outliers(img);
    dt = dti.dt(img);
    eig = dti.eig(dt);
    % ad = dti.ad(eig);
    % rd = dti.rd(eig);
    % l2 = dti.l2(eig);
    % l3 = dti.l3(eig);
    % fa = dti.fa(eig);
    md = dti.md(eig);
    
    % % save 4D
    % info.ImageSize = [120,120,36,6];
    % niftiwrite(dt,strcat(maps_folder, "/", subject_code, "_dt.nii.gz"),info);
    % info.ImageSize = [120,120,36,3];
    % niftiwrite(clip(1000*abs(eig),0,3.5),strcat(maps_folder, "/", subject_code, "_eig.nii.gz"),info);
    
    % edit info
    info.ImageSize = info.ImageSize(1:3);
    info.PixelDimensions = info.PixelDimensions(1:3);
    % save 3D
    niftiwrite(clip(1000*abs(md),0,3.5),strcat(maps_folder, "/", subject_code, "_MD.nii.gz"),info);
    % niftiwrite(clip(1000*abs(rd),0,3.5),strcat(maps_folder, "/", subject_code, "_RD.nii.gz"),info);
    % niftiwrite(clip(1000*abs(ad),0,3.5),strcat(maps_folder, "/", subject_code, "_AD.nii.gz"),info);
    % niftiwrite(clip(1000*abs(l2),0,3.5),strcat(maps_folder, "/", subject_code, "_L2.nii.gz"),info);
    % niftiwrite(clip(1000*abs(l3),0,3.5),strcat(maps_folder, "/", subject_code, "_L3.nii.gz"),info);
    % niftiwrite(fa,strcat(maps_folder, "/", subject_code, "_FA.nii.gz"),info);

end