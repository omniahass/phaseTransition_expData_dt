% 
clear all
dvs_file = '/Users/gabriellebaxter/Downloads/DiffusionVectors_HEAL_20250820.dvs';
[dirs,mods,Ndirs] = Read_all_dvs_directions(dvs_file);
dirs = dirs{1}; % vectors but with scaling
mods = mods{1}; % factor of scaling (

unscaled_dirs = dirs./mods; 

[~,name,~] = fileparts(dvs_file);
bvec_filename = ['/Users/gabriellebaxter/Documents/code/diffusionGradients/dvsfiles/' name '.bvec'];

vecs = unscaled_dirs.'; % transpose to write

% Save to text file with space-separated columns
fid = fopen(bvec_filename,'w');
fprintf(fid, '%f ', vecs(1,:));
fprintf(fid, '\n');
fprintf(fid, '%f ', vecs(2,:));
fprintf(fid, '\n');
fprintf(fid, '%f ', vecs(3,:));
fprintf(fid, '\n');
fclose(fid);

function [dirs,mods,Ndirs] = Read_all_dvs_directions(inputDVSfile,UniqueNumber)
% [dirs,mods,Ndirs] = Read_all_dvs_directions(inputDVSfile,UniqueNumber)
%
% It is still unstable, needs further tweaking but worked when I needed it
%
%
% ONLY WORKS WITH .dvs files in the format of the product diffusion sequence
%
% By: Santiago Coelho (17/06/2021)
fileID=fopen(inputDVSfile,'r');
formatSpec = '%s';
A = textscan(fileID,formatSpec);
fclose(fileID);

% Filter directions
NdirSets=0;AllDirs_numbers=[];AllDirs_coord=[];
for ii=1:length(A{1})
    current_line=cell2mat(A{1}(ii));
    if (length(current_line)>11)&&strcmp(current_line(2:11),'directions')
        NdirSets=NdirSets+1;
        AllDirs_coord(NdirSets)=ii;
        AllDirs_numbers(NdirSets)=str2num(current_line(13:end-1));
    end
end

for ii=1:NdirSets
    clearvars current_dirs
    for jj=1:AllDirs_numbers(ii)
        x=str2num(regexprep(cell2mat(A{1}(AllDirs_coord(ii)+10+7*(jj-1))),',',''));
        y=str2num(regexprep(cell2mat(A{1}(AllDirs_coord(ii)+11+7*(jj-1))),',',''));
        z=str2num(regexprep(cell2mat(A{1}(AllDirs_coord(ii)+12+7*(jj-1))),',',''));
        current_dirs(jj,:)=[x y z];
    end
    dirs_all{ii}=current_dirs;
    mods_all{ii}=sqrt(sum(current_dirs.^2,2));
    if size(current_dirs,1)~=AllDirs_numbers(ii)
        error('Something went wrong with the number of directions')
    end
end

if ~exist('UniqueNumber', 'var') || isempty(UniqueNumber)
    dirs=dirs_all;
    mods=mods_all;
else
    if ~any(AllDirs_numbers==UniqueNumber)
        error('Requested set of directions is not in the dvs file')
    end
    dirs=dirs_all{AllDirs_numbers==UniqueNumber};
    mods=mods_all{AllDirs_numbers==UniqueNumber};
end
Ndirs=AllDirs_numbers;
end