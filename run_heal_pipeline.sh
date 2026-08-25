#!/bin/bash
#
# HEAL pipeline wrapper - runs the whole single case pipeline in order
#
# Usage:
#   ./run_heal_pipeline.sh                                  # uses the settings below
#   ./run_heal_pipeline.sh /path/to/case                    # override the parent folder
#   ./run_heal_pipeline.sh /path/to/case case.nii.gz        # + segmentation file
#   ./run_heal_pipeline.sh /path/to/case case.nii.gz 4567   # + only run steps 4,5,6,7
#
# Steps:
#   1  heal_pipeline.sh                          dicom clean up -> nifti + bmatrix
#   2  designer_heal.sh                          denoise, degibbs, topup, eddy
#   3  resample_roi.py                           segmentation -> dwi space ROI
#   4  dti_mapping_jelle_singlecase.m            DTI maps per diffusion time
#   5  heal_maps_registration_singlecase.ipynb   register maps within scan
#   6  rpbm_bayesian_singlecase.m                RPBM a / kappa maps
#   7  read_roi_values_singlecase.ipynb          slopes + ROI values -> csv
#
# Every step reads the outputs of the previous one from $parent_path, which is
# passed to the shell scripts as an argument and to MATLAB / the notebooks through
# the HEAL_PARENT environment variable.

# ================================ USER SETTINGS ================================
#parent_path="/Volumes/labspace/omnia/projects/HEAL/data/WCMyofascial2075"
parent_path="/mnt/fieree01lab/omnia/projects/HEAL/data/WCMyofascial2075"
segmentation_filename="WCMyofascial2075.nii.gz"
center="cornell"
steps="1234567"                 # which steps to run, e.g. "1234567", "45", "7"

python_bin="python3"
jupyter_bin="jupyter"
matlab_bin="/Applications/MATLAB_R2025b.app/bin/matlab"
# ===============================================================================

# Command line arguments override the settings above
[ -n "$1" ] && parent_path="$1"
[ -n "$2" ] && segmentation_filename="$2"
[ -n "$3" ] && steps="$3"

set -o pipefail

# Always run from the folder holding the scripts - heal_pipeline.sh, the MATLAB
# addpath calls and bmatrix_export.py all use relative paths
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir" || exit 1

parent_path="${parent_path%/}"          # no trailing slash, MATLAB adds its own
subject="$(basename "$parent_path")"
deltas=(0022ms 0042ms 0081ms 0156ms 0300ms)

export HEAL_PARENT="$parent_path"       # read by the .m files and the notebooks
export python_bin                       # read by heal_pipeline.sh

log_dir="$parent_path/logs"
log_file="$log_dir/pipeline_$(date '+%Y%m%d_%H%M%S').log"

# ---------------------------------- helpers -----------------------------------
msg()  { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$log_file"; }
warn() { echo "[$(date '+%H:%M:%S')] WARNING: $*" | tee -a "$log_file"; }
fail() { echo "[$(date '+%H:%M:%S')] ERROR: $*" | tee -a "$log_file"; exit 1; }

do_step() { [[ "$steps" == *"$1"* ]]; }

# Debugging checks - fail early with the exact path that is missing
check_file() {
    if [ -s "$1" ]; then
        msg "  ok       $1"
    else
        fail "missing or empty file: $1"
    fi
}

check_dir() {
    if [ -d "$1" ]; then
        msg "  ok       $1"
    else
        fail "missing folder: $1"
    fi
}

check_glob() { # check_glob <pattern> <expected count>
    local n
    n=$(ls -1d $1 2>/dev/null | wc -l | tr -d ' ')
    if [ "$n" -ge "$2" ]; then
        msg "  ok       $n file(s) matching $1"
    else
        fail "expected at least $2 file(s) matching $1, found $n"
    fi
}

need_cmd() {
    if command -v "$1" >/dev/null 2>&1; then
        msg "  ok       $1 -> $(command -v "$1")"
    else
        missing_cmds+=("$1")
    fi
}

need_bin() { # for full paths like matlab
    if [ -x "$1" ]; then
        msg "  ok       $1"
    else
        missing_cmds+=("$1")
    fi
}

need_pymod() {
    if "$python_bin" -c "import $1" >/dev/null 2>&1; then
        msg "  ok       python module $1"
    else
        missing_cmds+=("python module $1")
    fi
}

run_step() { # run_step <number> <name> <command...>
    local number="$1" name="$2" status; shift 2
    msg ""
    msg "=============================================================="
    msg " STEP $number - $name"
    msg "=============================================================="
    msg " running: $*"
    "$@" 2>&1 | tee -a "$log_file"
    status=$?
    [ $status -eq 0 ] || fail "step $number ($name) failed with exit code $status"
    msg " step $number done"
}

run_notebook() {
    "$jupyter_bin" nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=-1 \
        --output-dir "$log_dir" \
        --output "$(basename "${1%.ipynb}")_executed.ipynb" \
        "$1"
}

run_matlab() {
    "$matlab_bin" -sd "$script_dir" -batch "$1"
}

# DTI.m calls mex_wls and mex_dti_eig directly, with no pure MATLAB fallback, so the
# compiled binaries have to match the architecture of the installed MATLAB
check_mex() {
    local matlab_root arch have
    matlab_root="$(dirname "$(dirname "$matlab_bin")")"
    if   [ -d "$matlab_root/bin/glnxa64" ]; then arch="mexa64"      # linux cluster
    elif [ -d "$matlab_root/bin/maca64"  ]; then arch="mexmaca64"   # apple silicon
    else                                        arch="mexmaci64"   # intel mac
    fi
    for m in mex_wls mex_dti_eig; do
        if [ -f "$script_dir/supporting_functions_dti/${m}.${arch}" ]; then
            msg "  ok       supporting_functions_dti/${m}.${arch}"
        else
            have=$(ls "$script_dir"/supporting_functions_dti/${m}.mex* 2>/dev/null | xargs -n1 basename 2>/dev/null | tr '\n' ' ')
            missing_cmds+=("supporting_functions_dti/${m}.${arch}  (MATLAB is ${arch}, this folder has: ${have:-nothing})")
        fi
    done
}

# these two toolboxes are used inside the executed code path
check_toolboxes() {
    local matlab_root
    matlab_root="$(dirname "$(dirname "$matlab_bin")")"
    [ -d "$matlab_root/toolbox/images" ] || warn "Image Processing Toolbox not found - padarray/medfilt3 need it"
    [ -d "$matlab_root/toolbox/stats" ]  || warn "Statistics and Machine Learning Toolbox not found - mad/nanmedian/randsample need it"
}

# ------------------------------- start up -------------------------------------
[ -d "$parent_path" ] || { echo "ERROR: parent folder does not exist: $parent_path"; exit 1; }
mkdir -p "$log_dir" || exit 1

msg "=============================================================="
msg " HEAL pipeline"
msg "=============================================================="
msg " parent folder  : $parent_path"
msg " subject        : $subject"
msg " segmentation   : $segmentation_filename"
msg " center         : $center"
msg " steps          : $steps"
msg " script folder  : $script_dir"
msg " log file       : $log_file"

# ------------------------- preflight: tools and inputs -------------------------
msg ""
msg "--- preflight checks ---"
missing_cmds=()

do_step 1 && { need_cmd "$python_bin"; need_cmd dcm2niix; need_pymod pydicom; need_pymod numpy; }
do_step 2 && { need_cmd docker; need_cmd mrinfo; need_cmd mrconvert; need_cmd mrcat
               need_cmd mrmath; need_cmd topup; need_cmd eddy; need_cmd mri_synthstrip; }
do_step 3 && { need_cmd "$python_bin"; need_pymod SimpleITK; }
do_step 4 && { need_bin "$matlab_bin"; check_mex; check_toolboxes; }
do_step 5 && { need_cmd "$jupyter_bin"; need_pymod ants; }
do_step 6 && { need_bin "$matlab_bin"; check_toolboxes; }
do_step 7 && { need_cmd "$jupyter_bin"; need_pymod nibabel; need_pymod dipy
               need_pymod scipy; need_pymod pandas; }

if [ ${#missing_cmds[@]} -gt 0 ]; then
    for c in "${missing_cmds[@]}"; do echo "  MISSING  $c" | tee -a "$log_file"; done
    fail "${#missing_cmds[@]} required tool(s) not found - fix the paths at the top of this script or source your FSL / conda setup"
fi

# scripts this run needs
do_step 1 && check_file "$script_dir/heal_pipeline.sh"
do_step 2 && check_file "$script_dir/designer_heal.sh"
do_step 3 && check_file "$script_dir/resample_roi.py"
do_step 4 && check_file "$script_dir/dti_mapping_jelle_singlecase.m"
do_step 5 && check_file "$script_dir/heal_maps_registration_singlecase.ipynb"
do_step 6 && check_file "$script_dir/rpbm_bayesian_singlecase.m"
do_step 7 && check_file "$script_dir/read_roi_values_singlecase.ipynb"

# input data
do_step 1 && check_dir "$parent_path/sorted_dcms"
if do_step 2; then
    docker info >/dev/null 2>&1 || fail "docker is not running - designer2 needs it"
    msg "  ok       docker daemon is running"
fi

# =============================== STEP 1 =======================================
# in : $parent_path/sorted_dcms
# out: $parent_path/nifti, $parent_path/bmatrix
if do_step 1; then
    chmod +x ./heal_pipeline.sh
    run_step 1 "prepare dicoms -> nifti + bmatrix" \
        ./heal_pipeline.sh "$parent_path" "$segmentation_filename" "$center"

    msg "--- checking step 1 output ---"
    check_dir "$parent_path/nifti"
    for d in "${deltas[@]}"; do
        check_glob "$parent_path/nifti/*${d}*_32dir.nii*"  1
        check_glob "$parent_path/nifti/*${d}*_32dir.bvec"  1
        check_glob "$parent_path/nifti/*${d}*_32dir.bval"  1
        check_glob "$parent_path/nifti/*${d}*_32dir.json"  1
        check_file "$parent_path/bmatrix/${d}.txt"
    done
    # designer_heal.sh needs a reverse phase encode b0 for PGSE (0022ms) and STEAM (0042ms)
    check_glob "$parent_path/nifti/*0022ms*_b0_FH.nii*"  1
    check_glob "$parent_path/nifti/*0042ms*_b0_FH.nii*"  1
    check_glob "$parent_path/nifti/*0022ms*_b0_FH.json"  1
    check_glob "$parent_path/nifti/*0042ms*_b0_FH.json"  1
fi

# =============================== STEP 2 =======================================
# in : $parent_path/nifti
# out: $parent_path/derivatives/all/<delta>/dwiec.nii, .../magdn/sigma.nii
if do_step 2; then
    chmod +x ./designer_heal.sh
    run_step 2 "designer - denoise, degibbs, topup, eddy" \
        ./designer_heal.sh "$parent_path"

    msg "--- checking step 2 output ---"
    check_file "$parent_path/derivatives/all/magdn/sigma.nii"
    for d in "${deltas[@]}"; do
        check_file "$parent_path/derivatives/all/${d}/dwiec.nii"
        check_file "$parent_path/derivatives/all/${d}/${d}.bval"
        check_file "$parent_path/derivatives/all/${d}/${d}.bvec"
        # the DTI fit uses one b matrix row per dwi volume - catch a mismatch here
        nvols=$(mrinfo -size "$parent_path/derivatives/all/${d}/dwiec.nii" | awk '{print $4}')
        nrows=$(grep -cve '^[[:space:]]*$' "$parent_path/bmatrix/${d}.txt")
        if [ "$nvols" = "$nrows" ]; then
            msg "  ok       ${d}: $nvols volumes, $nrows b matrix rows"
        else
            warn "${d}: dwiec.nii has $nvols volumes but bmatrix/${d}.txt has $nrows rows"
        fi
    done
fi

# =============================== STEP 3 =======================================
# in : $parent_path/$segmentation_filename, $parent_path/derivatives/all/0022ms/dwiec.nii
# out: $parent_path/${subject}_reformatdwi.nii.gz
if do_step 3; then
    check_file "$parent_path/derivatives/all/0022ms/dwiec.nii"
    # segmentation_filename may be a name inside the case folder or a full path
    case "$segmentation_filename" in
        /*) seg_path="$segmentation_filename" ;;
        *)  seg_path="$parent_path/$segmentation_filename" ;;
    esac
    check_file "$seg_path"

    run_step 3 "resample segmentation into dwi space" \
        "$python_bin" ./resample_roi.py \
            -parent_folder "$parent_path" \
            -segmentation_filename "$segmentation_filename"

    msg "--- checking step 3 output ---"
    check_file "$parent_path/${subject}_reformatdwi.nii.gz"
fi

# =============================== STEP 4 =======================================
# in : $parent_path/derivatives/all/<delta>/dwiec.nii, $parent_path/bmatrix
# out: $parent_path/maps/<delta>_{RD,AD,FA,L2,L3}.nii.gz
if do_step 4; then
    for d in "${deltas[@]}"; do
        check_file "$parent_path/derivatives/all/${d}/dwiec.nii"
        check_file "$parent_path/bmatrix/${d}.txt"
    done

    run_step 4 "DTI mapping (MATLAB)" run_matlab dti_mapping_jelle_singlecase

    msg "--- checking step 4 output ---"
    for d in "${deltas[@]}"; do
        for m in RD AD FA L2 L3; do
            check_file "$parent_path/maps/${d}_${m}.nii.gz"
        done
    done
fi

# =============================== STEP 5 =======================================
# in : $parent_path/maps/<delta>_<map>.nii.gz, derivatives/all/<delta>/dwiec.nii
# out: $parent_path/maps/<delta>_<map>_withinscanreg.nii.gz
if do_step 5; then
    for d in "${deltas[@]}"; do
        check_file "$parent_path/derivatives/all/${d}/dwiec.nii"
        check_file "$parent_path/maps/${d}_RD.nii.gz"
    done

    run_step 5 "register maps to the 0022ms b0" \
        run_notebook heal_maps_registration_singlecase.ipynb

    msg "--- checking step 5 output ---"
    for d in "${deltas[@]:1}"; do       # 0022ms is the fixed image, it is not registered
        for m in RD AD FA L2 L3; do
            check_file "$parent_path/maps/${d}_${m}_withinscanreg.nii.gz"
        done
    done
fi

# =============================== STEP 6 =======================================
# in : maps/<delta>_{RD,AD}[_withinscanreg].nii.gz, magdn/sigma.nii, ROI mask
# out: $parent_path/rpbm_mapping/{a,kappa,sigma_dperp}.nii
if do_step 6; then
    check_file "$parent_path/derivatives/all/0022ms/0022ms.bval"
    check_file "$parent_path/derivatives/all/0022ms/0022ms.bvec"
    check_file "$parent_path/derivatives/all/0022ms/dwiec.nii"
    check_file "$parent_path/derivatives/all/magdn/sigma.nii"
    check_file "$parent_path/${subject}_reformatdwi.nii.gz"
    check_file "$parent_path/maps/0022ms_RD.nii.gz"
    check_file "$parent_path/maps/0022ms_AD.nii.gz"
    for d in "${deltas[@]:1}"; do
        check_file "$parent_path/maps/${d}_RD_withinscanreg.nii.gz"
        check_file "$parent_path/maps/${d}_AD_withinscanreg.nii.gz"
    done

    run_step 6 "RPBM bayesian mapping (MATLAB)" run_matlab rpbm_bayesian_singlecase

    msg "--- checking step 6 output ---"
    check_file "$parent_path/rpbm_mapping/sigma_dperp.nii"
    check_file "$parent_path/rpbm_mapping/a.nii"
    check_file "$parent_path/rpbm_mapping/kappa.nii"
fi

# =============================== STEP 7 =======================================
# in : maps/, rpbm_mapping/, ROI mask
# out: $parent_path/slopes/, $parent_path/dti_rpbm_results.csv
if do_step 7; then
    check_file "$parent_path/${subject}_reformatdwi.nii.gz"
    check_file "$parent_path/rpbm_mapping/a.nii"
    check_file "$parent_path/rpbm_mapping/kappa.nii"

    run_step 7 "slopes + ROI values -> csv" \
        run_notebook read_roi_values_singlecase.ipynb

    msg "--- checking step 7 output ---"
    check_file "$parent_path/slopes/RDslope.nii.gz"
    check_file "$parent_path/slopes/RDr2.nii.gz"
    check_file "$parent_path/slopes/RDrho.nii.gz"
    check_file "$parent_path/dti_rpbm_results.csv"
    msg ""
    msg "--- $parent_path/dti_rpbm_results.csv ---"
    cut -d, -f1-6 "$parent_path/dti_rpbm_results.csv" | tee -a "$log_file"
fi

msg ""
msg "=============================================================="
msg " pipeline finished - steps $steps for $subject"
msg " log: $log_file"
msg "=============================================================="
