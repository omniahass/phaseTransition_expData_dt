**Integrating the HEAL processing pipeline onto the Server**

**Starting Point**

*Six scripts that had to be run manually one at a time, on one particular Mac, across three languages (bash, MATLAB and Jupyter). each script carried a hard-coded absolute path to a subject folder. Processing a case meant editing six files by hand and remembering the order. There were no checks between stages, so a step that failed halfway through would leave partial outputs, and the next step would consume them without complaint. None of it had been run outside macOS.*

**Changes Made**

A single wrapper, run\_heal\_pipeline.sh, with one settings block at the top: the case folder, the segmentation file, and which steps to run. Everything downstream derives its paths from that one value.

| \# | Script | What it does |
| :---- | :---- | :---- |
| 1 | heal\_pipeline.sh | DICOM cleanup, NIfTI conversion, b-matrix extraction |
| 2 | designer\_heal.sh | Denoising, Gibbs correction, topup and eddy across five diffusion times |
| 3 | resample\_roi.py | Segmentation resampled into DWI space |
| 4 | dti\_mapping\_singlecase.py  (rewritten) | Diffusion tensor fit: AD, RD, FA, L2, L3 per diffusion time |
| 5 | heal\_maps\_registration\_singlecase.ipynb | Maps registered to the shortest diffusion time |
| 6 | rpbm\_bayesian\_singlecase.m | RPBM fit: fibre diameter a and permeability kappa |
| 7 | read\_roi\_values\_singlecase.ipynb | Diffusivity slopes and per-ROI values, written to CSV |

\- Each step can be run on its own, so a failure at step 6 does not mean redoing the four hours before it. 

\- Before anything executes, a preflight check verifies every tool, Python package and input file the selected steps need, and reports all the problems at once rather than dying on the first. 

\- After each step it confirms the specific files the next step will open. Everything is logged with timestamps into the case folder.

\- Added an 8th MATLAB script, plot\_roi\_results\_singlecase.m, reads the results CSV and produces the ROI figures.

**Defects found along the way**

These were latent in the inherited code — most only surfaced because the pipeline was being run somewhere new, or run twice.

| What was wrong | Consequence if left |
| :---- | :---- |
| designer\_heal.sh had a hard-coded subject path that silently overrode the one passed in | Would have processed a different subject than requested, without any error |
| rename\_folders.py moved every sub-folder into the DICOM folder | Swallowed the log folder; on any re-run would have swallowed nifti/, derivatives/ and maps/ — real data loss |
| mri\_synthstrip was called with a flag FreeSurfer 7.4.1 rejects | Brain mask never written, so all five eddy corrections failed — this blocked the whole run |
| The RPBM script never put the NIfTI reader folder on the MATLAB path | Step 6 could not load any image |
| The bundled NIfTI toolbox called a helper file that is not in the folder, and had no .gz support | Step 6 failed twice more — first no reader, then mangled .nii.gz filenames |
| resample\_roi.py named its output from the segmentation filename rather than the subject folder | Downstream steps look it up by subject name — a mismatch would break steps 6 and 7 |

**The DTI fitting step had to be rewritten**

This is the one substantive change to the science code, so it is worth explaining properly.

The MATLAB tensor fitting depends on two compiled helper files, mex\_wls and mex\_dti\_eig. Both exist only as macOS binaries, and the C source is nowhere in the repository — so step 4 simply could not run on a Linux machine. I asked around the lab; the consensus was that this class is outdated and there are current alternatives.

## **Why not just use MRtrix**

MRtrix’s dwi2tensor and dipy both accept a gradient table only — a direction plus a scalar b-value. Our pipeline runs a dedicated step to pull the scanner’s effective b-matrix out of the DICOM headers, including the off-diagonal cross terms, and the MATLAB fit used all six elements. Switching to either tool would mean discarding the very thing that earlier step exists to produce.

So I reimplemented the fit directly in Python, matching the original settings exactly: the same neighbourhood smoothing, the same dropping of the leading and trailing b0, the same full b-matrix design, weighted linear least squares with one reweighting step, eigenvalues sorted descending, and the same scaling and clipping of the output maps.

## **How it was checked**

| Test | Agreement |
| :---- | :---- |
| Known tensor recovered through the real b-matrix | 3 x 10^-14 um^2/ms |
| Eigenvalues against a direct decomposition | 3 x 10^-17 |
| Smoothing against a literal transcription of the MATLAB loop | 2 x 10^-16 (1 ulp) |
| Which neighbouring voxels the smoothing selects | identical |

## 

In other words, the two implementations agree to the limits of floating-point arithmetic on synthetic data. They have not been compared on a real case — see below.

**Reproducibility: what could differ from Gabrielle’s results**

Ordered by how much I would want to check them before trusting a number.

**Nothing has been validated against a known-good case   \[CHECK FIRST\]**

The pipeline runs end to end and produces plausible output, but no case has been processed both ways and compared. Gabrielle already processed the full WCMyofascial cohort, so a direct comparison is available and is the obvious next step.

**The brain-masking flag was changed   \[CHECK FIRST\]**

The original called SynthStrip with a \-t 1 flag that FreeSurfer 7.4.1 rejects outright, so no mask was produced and everything downstream failed. In versions that accept it, \-t sets the thread count and has no effect on the mask — but I cannot confirm which version she used. This mask feeds every eddy correction, so it is worth verifying visually.

**Tensor fitting moved from MATLAB to Python   \[VERIFY ON REAL DATA\]**

Verified equivalent to machine precision on synthetic data, but never compared against her actual maps. This would be settled by the same validation run.

**Tool versions on the cluster differ from her Mac   \[UNKNOWN\]**

FSL, MRtrix, the designer container, FreeSurfer and MATLAB are all whatever is installed on server 07\. We do not have a record of the versions she used, so this cannot be quantified — only bounded by the validation comparison.

**Output images are now written by a different library   \[LOW CONCERN\]**

Step 4 writes NIfTI files through Python rather than MATLAB. Header and affine handling should match, but subtle differences are possible.

**NIfTI reader now skips intensity rescaling   \[VERIFIED NO-OP\]**

The repaired reader returns raw stored values rather than applying the header scale factor. I checked every file step 6 opens — all have a scale of 1 and an offset of 0, so the two paths return identical data.

**ROI file naming rule changed   \[VERIFIED NO-OP\]**

The resampled ROI is now named from the subject folder rather than the segmentation filename. Identical whenever the segmentation is named after the subject, which it is in every case we have.

## **Confirmed unchanged**

Every processing parameter was checked against the originals and is byte-for-byte identical:

* Denoising and Gibbs-correction settings

* Topup configuration and all eddy flags

* Smoothing kernel and threshold, b0 trimming, weighted-least-squares iteration count

* Diffusivity scaling and clipping limits

* RPBM training size, diameter factor, polynomial degree, parameter bounds and noise scaling

* Registration transform and interpolation

* ROI outlier trimming percentiles

**Open items**

Roughly in the order I would tackle them.

1. **Run a validation case.** Process one already-published WCMyofascial subject and compare maps and ROI values against Gabrielle’s. This settles most of the reproducibility questions above at once.

2. **A slope inconsistency in the ROI step.** A condition that can never be true means the diffusivity-versus-time slopes are fit on unregistered maps, while the RPBM step and the ROI readout in the same notebook both use registered ones. One-line fix, but it changes the numbers, so I have not applied it unilaterally.

3. **Outlier detection is computed and then discarded.** The original DTI script detects outlier volumes and never passes them to the fit. Pre-existing, and carried into the Python version deliberately so behaviour did not change silently. Worth deciding whether it was meant to be used.

4. **The RPBM fit may be running into its bounds.** Permeability comes out at 0.044-0.049 against an upper bound of 0.05 in all four ROIs, and diameter at 47-55 against a bound of 60\. Worth checking whether the priors need widening.

5. **The bundled NIfTI toolbox is not the stock version.** Two separate pieces were missing from it. Better to drop in the upstream release than keep patching this copy.

6. **C2 cohort segmentations do not exist yet.** Confirmed none are on lab storage. Waiting on the Cornell segmentation model to be published.

7. **Two steps still run as notebooks.** Converting them to plain scripts removes a whole class of failure — one already cost a run.

8. **Server 07’s shared home directory is full.** 40 GB across all users, at 100%. It caused three separate failures. Worked around by redirecting caches, but it needs an admin fix and affects everyone on that machine.

Code and full change history: github.com/omniahass/phaseTransition\_expData\_dt. Every modification to the inherited scripts is marked in-line and described in the commit history.
