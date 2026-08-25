#!/bin/bash

# CHANGED FOR WRAPPER: parent folder comes from the first argument when called by
# run_heal_pipeline.sh, otherwise it falls back to the hard-coded path below
dataDir="${1:-/Volumes/labspace/Projects/HEAL/C2_001}"
deltas=(0022ms 0042ms 0081ms 0156ms 0300ms)

# Initialize an empty string for concatenating paths
all_img_paths=""
dpath=$dataDir/derivatives/all
dwi_path=$dpath/dwi_designer.nii

# # # --------------------------------------------- DENOISING -------------------------------------------------
### Loop through each delta value to concatenate paths
for delta in "${deltas[@]}"; do
    echo "Processing delta = ${delta} ms..."
    img_file=$(find "$dataDir/nifti/" -type f -name "*${delta}*_32dir.nii*" | sort -V | head -n 1 | xargs basename)
    img_dir="/data/nifti/$img_file"
    echo "Image : $img_dir"
    all_img_paths+="$img_dir,"

done
# Remove trailing comma from all_img_paths
all_img_paths="${all_img_paths%,}"

# Run magdn for all deltas
echo "  > Magnitude denoising"
rm -rf $dpath/magdn
docker run --platform linux/amd64 -it \
    -v "$dataDir:/data" \
    nyudiffusionmri/designer2:v2.0.16 \
    designer -denoise -shrinkage frob -extent 5,5,5 -nthreads 24 -nocleanup \
    -scratch /data/derivatives/all/magdn \
    "$all_img_paths" \
    /data/derivatives/all/dwi_designer.nii

# # ## ---------------------------------------- Split up files --------------------------------------------
size=(`mrinfo -size $dwi_path`)
echo "$dwi_path"
echo "$size"
nvols=${size[3]}
ndeltas=${#deltas[@]}
nvol_per_dmri=$((nvols/ndeltas))
echo "$nvol_per_dmri"

c=0
for i in ${deltas[@]}; do
    (( c++ ))
    echo "${deltas[@]}"
    start=$(((c-1)*nvol_per_dmri))
    end=$(((nvol_per_dmri*c)-1))
    echo "$c $start $end"
    folder_path=$dpath/${i} # put everything into individual folders
    mkdir -p $folder_path
    img_path=$folder_path/${i}.nii
    mrconvert -force -coord 3 $start:$end $dwi_path $img_path
done
rm -f $dpath/dwi_designer.nii
rm -f $dpath/dwi_designer.bval
rm -f $dpath/dwi_designer.bvec
# Clean up files except sigma
magdn_path=$dpath/magdn
find "$magdn_path" -maxdepth 1 -type f ! -name "sigma.nii" -delete
# # ------------------------------------------------ DEGIBBS ----------------------------------------------
for delta in "${deltas[@]}"; do
    echo "Processing delta = ${delta} ms..."
    echo "Degibbs"

    img_folder=$dpath/${delta}
    img_path=$dpath/${delta}/${delta}.nii
    bvec_path=$dataDir/bvec_bval/$delta.bvec
    
    orig_bvec=(`ls -d $dataDir/nifti/*${delta}*_32dir.bvec`)
    cp $orig_bvec $img_folder/${delta}.bvec
    orig_bval=(`ls -d $dataDir/nifti/*${delta}*_32dir.bval`)
    cp $orig_bval $img_folder/${delta}.bval
    orig_json=(`ls -d $dataDir/nifti/*${delta}*_32dir.json`)
    cp $orig_json $img_folder/${delta}.json

    echo "$img_folder"
    echo "$orig_bvec"
    echo "$orig_bval"
    echo "$bvec_path"

# # # #     # -------------------------------run degibbs-------------------------------------
    echo "Degibbs"
    name=${delta}_dg
    processing_name=degibbs
    rm -rf $img_folder/${processing_name}
    docker run --rm --platform linux/amd64 -it \
        -v "$dataDir:/data" \
        nyudiffusionmri/designer2:main  \
        designer -degibbs -pf 6/8 -pe_dir k \
        -nocleanup \
        -scratch /data/derivatives/all/${delta}/${processing_name} \
        /data/derivatives/all/${delta}/${delta}.nii /data/derivatives/all/${delta}/${delta}_dg.nii

    # -----------------------------------------Convert to mif---------------------------------------------
    cp $img_folder/${delta}.bval $img_folder/${name}.bval
    bval_path=$img_folder/${delta}.bval
    cp $img_folder/${delta}.json $img_folder/${delta}_dg.json
    dmri_path=$img_folder/${name}.nii
    cp $img_folder/${delta}.bvec $img_folder/${name}.bvec
    bvec_path=$img_folder/${name}.bvec
    dmri_bids=$img_folder/${name}.json

    mrconvert $dmri_path -force -fslgrad $bvec_path $bval_path -json_import $dmri_bids $img_folder/${processing_name}/${name}.mif
done

####### -------------------- topup and eddy PGSE -----------------------

# PGSE
deltas=(0022ms)
for delta in "${deltas[@]}"; do
    echo "$delta"
    # set up topup paths
    topup_result_path=$dpath/$delta/topup
    mkdir -p $topup_result_path
    # get reverse PE image and json file and convert to mif
    reverse_pe_delta="0022ms"
    rpe_file=$(find "$dataDir/nifti/" -type f -name "*${reverse_pe_delta}*_b0_FH.nii*" | sort -V | head -n 1 | xargs basename) # find the PGSE reverse PE image
    echo "rpe file $rpe_file"
    rpe_path="$dataDir/nifti/$rpe_file"
    rpe_bids_file=$(find "$dataDir/nifti/" -type f -name "*${reverse_pe_delta}*_b0_FH.json*" | sort -V | head -n 1 | xargs basename) # find the PGSE reverse PE json
    rpe_bids="$dataDir/nifti/$rpe_bids_file"
    mrconvert $rpe_path -force -coord 3 0 -json_import $rpe_bids $topup_result_path/b0rpe.mif # take just the first b0 

    pe_path="$dpath/$delta/degibbs/${delta}_dg.mif"
    mrconvert -force -coord 3 0 $pe_path $topup_result_path/b0pe.mif

    cd $topup_result_path
    mrinfo b0pe.mif -force -export_pe_eddy topup_config_1.txt topup_indicies_1.txt # info for normal PE
    mrinfo b0rpe.mif -force -export_pe_eddy topup_config_2.txt topup_indicies_2.txt # info for reverse PE
    cat topup_config_1.txt topup_config_2.txt > topup_acqp.txt # concatenate them for topup info

    # concatenate pair
    mrconvert -force b0rpe.mif b0rpe.nii
    mrconvert -force -coord 3 0 b0pe.mif b0pe.nii
    mrcat -force -axis 3 b0pe.nii b0rpe.nii b0_pair_topup.nii

    echo "Running topup"
    topup --imain=b0_pair_topup.nii --datain=topup_acqp.txt --config=b02b0.cnf --scale=1 --out=topup_results --iout=topup_results.nii.gz
    mrmath -force topup_results.nii.gz mean topup_corrected_mean.nii -axis 3
    # CHANGED FOR WRAPPER: "-t 1" (thread count) is not accepted by FreeSurfer 7.4.1
    # and made synthstrip exit without writing the mask, which broke every eddy call
    mri_synthstrip -i topup_corrected_mean.nii -o topup_corrected_brain.nii.gz -m topup_corrected_brain_mask.nii.gz

    echo "Running eddy $delta"
    name=${delta}_dg 
    img_folder=$dpath/${delta}
    img_path=$img_folder/${name}.nii
    processing_name=degibbs
    
    #create eddy folder
    mkdir -p $img_folder/eddy
    cd $img_folder/eddy
    mrinfo $img_folder/${processing_name}/${name}.mif -force -export_pe_eddy config.txt indicies.txt

    topup=$topup_result_path/topup_results
    data=$img_path
    dmri_bvec=$img_folder/${delta}.bvec
    dmri_bval=$img_folder/${delta}.bval
    echo "$dmri_bvec"
    echo "$dmri_bval"

    # make a fake bval 
    dmri_bval_fake=$img_folder/fake.bval
    cp $dmri_bval "$dmri_bval_fake"
    awk '{for(i=1; i<=NF; i++) if ($i < 100) $i = 0; print}' "$dmri_bval" > "$dmri_bval_fake" # change first value to 0 so eddy works

    mask=$topup_result_path/topup_corrected_brain_mask.nii.gz
    acqp=config.txt
    index=indicies.txt
    out2=eddy_processing/dwiec.nii.gz
    mkdir -p eddy_processing
    
    # # ----------------------------------run eddy-------------------------------------
    eddy --imain=$data --mask=$mask --acqp=$acqp --index=$index \
    --bvecs=$dmri_bvec --bvals=$dmri_bval_fake --repol --data_is_shelled \
    --topup=$topup --out=$out2 --verbose

# -----------------------------------------------------------------------------------
    dmri_bval=$img_folder/${delta}.bval
    mrconvert -force -fslgrad $dmri_bvec $dmri_bval $out2 dwiec.mif
    mrconvert -force -export_grad_fsl $img_folder/dwiec.bvec $img_folder/dwiec.bval dwiec.mif $img_folder/dwiec.nii
done

####### -------------------- topup and eddy STEAM -----------------------
## run topup once
deltas=(0042ms)
for delta in "${deltas[@]}"; do
    echo "$delta"
    topup_result_path=$dpath/$delta/topup
    mkdir -p $topup_result_path
    reverse_pe_delta="0042ms"
    rpe_file=$(find "$dataDir/nifti/" -type f -name "*${reverse_pe_delta}*_b0_FH.nii*" | sort -V | head -n 1 | xargs basename) # find the PGSE reverse PE image
    echo "$rpe_file"
    rpe_path="$dataDir/nifti/$rpe_file"
    rpe_bids_file=$(find "$dataDir/nifti/" -type f -name "*${reverse_pe_delta}*_b0_FH.json*" | sort -V | head -n 1 | xargs basename) # find the PGSE reverse PE json
    rpe_bids="$dataDir/nifti/$rpe_bids_file"
    mrconvert $rpe_path -force -coord 3 0 -json_import $rpe_bids $topup_result_path/b0rpe.mif # take just the first b0 

    pe_path="$dpath/$delta/degibbs/${delta}_dg.mif"
    mrconvert -force -coord 3 0 $pe_path $topup_result_path/b0pe.mif

    cd $topup_result_path
    mrinfo b0pe.mif -force -export_pe_eddy topup_config_1.txt topup_indicies_1.txt # info for normal PE
    mrinfo b0rpe.mif -force -export_pe_eddy topup_config_2.txt topup_indicies_2.txt # info for reverse PE
    cat topup_config_1.txt topup_config_2.txt > topup_acqp.txt # concatenate them for topup info

    # concatenate pair
    mrconvert -force b0rpe.mif b0rpe.nii
    mrconvert -force -coord 3 0 b0pe.mif b0pe.nii
    mrcat -force -axis 3 b0pe.nii b0rpe.nii b0_pair_topup.nii

    echo "Running topup"
    topup --imain=b0_pair_topup.nii --datain=topup_acqp.txt --config=b02b0.cnf --scale=1 --out=topup_results --iout=topup_results.nii.gz
    mrmath -force topup_results.nii.gz mean topup_corrected_mean.nii -axis 3
    # mri_synthstrip -i topup_corrected_mean.nii -o topup_corrected_brain.nii.gz -m topup_corrected_brain_mask.nii.gz
done

######### -------------------- eddy -----------------------
deltas=(0042ms 0081ms 0156ms 0300ms)

for delta in "${deltas[@]}"; do
    echo "Running eddy $delta"
    name=${delta}_dg 
    img_folder=$dpath/${delta}
    img_path=$img_folder/${name}.nii
    processing_name=degibbs
    
    #create eddy folder
    mkdir -p $img_folder/eddy
    cd $img_folder/eddy
    mrinfo $img_folder/${processing_name}/${name}.mif -force -export_pe_eddy config.txt indicies.txt

    reverse_pe_delta="0042ms"
    topup_result_path=$dpath/$reverse_pe_delta/topup
    echo "Topup path = ${topup_result_path}"
    topup=$topup_result_path/topup_results
    data=$img_path
    dmri_bvec=$img_folder/${delta}.bvec
    dmri_bval=$img_folder/${delta}.bval
    echo "$dmri_bvec"
    echo "$dmri_bval"

    # make a fake bval 
    dmri_bval_fake=$img_folder/fake.bval
    cp $dmri_bval "$dmri_bval_fake"
    awk '{for(i=1; i<=NF; i++) if ($i < 100) $i = 0; print}' "$dmri_bval" > "$dmri_bval_fake" # change first value to 0 so eddy works

    mask=$dpath/0022ms/topup/topup_corrected_brain_mask.nii.gz # use brain mask from other diffusion time
    acqp=config.txt
    index=indicies.txt
    out2=eddy_processing/dwiec.nii.gz
    mkdir -p eddy_processing
    
    # # ----------------------------------run eddy-------------------------------------
    eddy --imain=$data --mask=$mask --acqp=$acqp --index=$index \
    --bvecs=$dmri_bvec --bvals=$dmri_bval_fake --repol --data_is_shelled \
    --topup=$topup --out=$out2 --verbose

# -----------------------------------------------------------------------------------
    dmri_bval=$img_folder/${delta}.bval
    mrconvert -force -fslgrad $dmri_bvec $dmri_bval $out2 dwiec.mif
    mrconvert -force -export_grad_fsl $img_folder/dwiec.bvec $img_folder/dwiec.bval dwiec.mif $img_folder/dwiec.nii

done