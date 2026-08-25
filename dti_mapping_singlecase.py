#!/usr/bin/env python3
#
# DTI mapping - python replacement for dti_mapping_jelle_singlecase.m
#
# The MATLAB version relied on two compiled MEX files (mex_wls, mex_dti_eig) that
# only exist as macOS binaries, so it cannot run on the linux cluster. This script
# reproduces the same fit with numpy:
#
#   - same input   : derivatives/all/<delta>/dwiec.nii  +  bmatrix/<delta>.txt
#   - same output  : maps/<delta>_{RD,AD,L2,L3,FA}.nii.gz
#   - same settings: NLM smoothing (kernel 3, threshold 0.1), first and last volume
#                    dropped, full b matrix design, weighted linear least squares
#                    with one reweighting step, eigenvalues sorted descending,
#                    diffusivities scaled by 1000 and clipped to [0, 3.5]
#
# The b matrix from the scanner is used directly (including the cross terms), which
# is why this is not done with mrtrix dwi2tensor - that only accepts [x y z b].

import os
import numpy as np
import nibabel as nib

diffusion_times = ['0022ms', '0042ms', '0081ms', '0156ms', '0300ms']

# SET PARENT FOLDER HERE
# run_heal_pipeline.sh sets the HEAL_PARENT environment variable, otherwise the
# hard-coded default below is used
parent_folder = os.environ.get(
    'HEAL_PARENT',
    '/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial2075').rstrip('/')


def nlm_smoothing(dwi, kernel=3, threshold=0.1):
    """Average each voxel's signal vector with its most similar neighbours.

    Same as dwi_nlm_smoothing in the MATLAB version: a kernel^3 neighbourhood,
    circular edge wrapping, and the round(threshold*kernel^3) closest patches kept.
    """
    k = (kernel - 1) // 2
    nend = max(1, int(round(threshold * kernel ** 3)))

    # offsets in the same order MATLAB's reshape produces (first axis fastest)
    offsets = [(di, dj, dk)
               for dk in range(-k, k + 1)
               for dj in range(-k, k + 1)
               for di in range(-k, k + 1)]

    dist = np.empty((len(offsets),) + dwi.shape[:3])
    for p, (di, dj, dk) in enumerate(offsets):
        shifted = np.roll(dwi, shift=(-di, -dj, -dk), axis=(0, 1, 2))
        dist[p] = np.mean((shifted - dwi) ** 2, axis=3)

    # stable sort so ties break the same way MATLAB's sort does
    order = np.argsort(dist, axis=0, kind='stable')[:nend]
    del dist

    total = np.zeros_like(dwi)
    count = np.zeros(dwi.shape[:3])
    for p, (di, dj, dk) in enumerate(offsets):
        selected = np.any(order == p, axis=0)
        if not selected.any():
            continue
        shifted = np.roll(dwi, shift=(-di, -dj, -dk), axis=(0, 1, 2))
        total[selected] += shifted[selected]
        count[selected] += 1

    return total / count[..., None]


def fit_tensor(signal, bmatrix, chunk=20000):
    """Weighted linear least squares tensor fit.

    signal   [nvol, nvox]
    bmatrix  [nvol, 6] as [bxx bxy bxz byy byz bzz]

    Model: log(S) = log(S0) - (bxx*Dxx + 2bxy*Dxy + 2bxz*Dxz
                               + byy*Dyy + 2byz*Dyz + bzz*Dzz)
    which is the design matrix the DTI class builds in setGrad.
    """
    nvol, nvox = signal.shape
    design = np.column_stack([np.ones(nvol),
                              -bmatrix[:, 0], -2 * bmatrix[:, 1], -2 * bmatrix[:, 2],
                              -bmatrix[:, 3], -2 * bmatrix[:, 4], -bmatrix[:, 5]])

    # same rescaling and clamping the DTI class applies before fitting
    signal = signal * (1000.0 / signal.max())
    signal[signal <= 0] = np.finfo(float).eps
    logsignal = np.log(signal)

    # unweighted fit first, exactly like dt_lls
    params = np.linalg.pinv(design) @ logsignal

    # one reweighting step with the predicted signal as weights, like dt_wlls
    weights = np.exp(design @ params) ** 2

    for start in range(0, nvox, chunk):
        stop = min(start + chunk, nvox)
        w = weights[:, start:stop]
        y = logsignal[:, start:stop]
        lhs = np.einsum('ni,nv,nj->vij', design, w, design)
        rhs = np.einsum('ni,nv,nv->vi', design, w, y)
        try:
            params[:, start:stop] = np.linalg.solve(lhs, rhs[..., None])[..., 0].T
        except np.linalg.LinAlgError:  # singular voxels, fall back
            params[:, start:stop] = (np.linalg.pinv(lhs) @ rhs[..., None])[..., 0].T

    return params[1:7]  # drop the log(S0) row -> [Dxx Dxy Dxz Dyy Dyz Dzz]


def eigenvalues(tensor):
    """Descending eigenvalues of the diffusion tensor. tensor is [6, nvox]."""
    dxx, dxy, dxz, dyy, dyz, dzz = tensor
    nvox = tensor.shape[1]
    matrices = np.empty((nvox, 3, 3))
    matrices[:, 0, 0] = dxx
    matrices[:, 1, 1] = dyy
    matrices[:, 2, 2] = dzz
    matrices[:, 0, 1] = matrices[:, 1, 0] = dxy
    matrices[:, 0, 2] = matrices[:, 2, 0] = dxz
    matrices[:, 1, 2] = matrices[:, 2, 1] = dyz
    matrices[~np.isfinite(matrices)] = 0
    return np.linalg.eigvalsh(matrices)[:, ::-1].T  # eigvalsh is ascending


def main():
    print('Parent folder:', parent_folder)
    maps_folder = os.path.join(parent_folder, 'maps')
    os.makedirs(maps_folder, exist_ok=True)

    for diffusion_time in diffusion_times:
        nifti_folder = os.path.join(parent_folder, 'derivatives', 'all', diffusion_time)
        if not os.path.isdir(nifti_folder):
            continue
        dwi_file = os.path.join(nifti_folder, 'dwiec.nii')
        bmat_file = os.path.join(parent_folder, 'bmatrix', diffusion_time + '.txt')
        if not os.path.exists(dwi_file) or not os.path.exists(bmat_file):
            print('skipping %s, missing dwiec.nii or bmatrix' % diffusion_time)
            continue

        print(diffusion_time)
        image = nib.load(dwi_file)
        dwi = np.asarray(image.dataobj, dtype=np.float64)
        bmatrix = np.loadtxt(bmat_file)

        if dwi.shape[3] != bmatrix.shape[0]:
            print('  WARNING: %d volumes but %d b matrix rows'
                  % (dwi.shape[3], bmatrix.shape[0]))

        dwi = nlm_smoothing(dwi, 3, 0.1)

        # remove the b0 at the beginning and end
        dwi = dwi[:, :, :, 1:-1]
        bmatrix = bmatrix[1:-1, :]

        shape = dwi.shape[:3]
        signal = dwi.reshape(-1, dwi.shape[3]).T
        tensor = fit_tensor(signal, bmatrix)
        eigval = eigenvalues(tensor)

        ad = eigval[0]
        l2 = eigval[1]
        l3 = eigval[2]
        rd = (l2 + l3) / 2.0
        fa = (np.sqrt(0.5) * np.sqrt((ad - l2) ** 2 + (l2 - l3) ** 2 + (l3 - ad) ** 2)
              / np.sqrt(ad ** 2 + l2 ** 2 + l3 ** 2))

        def save(values, name, scale=True):
            values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
            if scale:
                values = np.clip(1000 * np.abs(values), 0, 3.5)
            out = nib.Nifti1Image(values.reshape(shape), image.affine, image.header)
            out.header.set_data_dtype(np.float64)
            filename = os.path.join(maps_folder, '%s_%s.nii.gz' % (diffusion_time, name))
            nib.save(out, filename)
            print('  wrote', filename)

        save(rd, 'RD')
        save(ad, 'AD')
        save(l2, 'L2')
        save(l3, 'L3')
        save(fa, 'FA', scale=False)


if __name__ == '__main__':
    main()
