#!/usr/bin/env python3
#
# Compare two sets of DTI maps voxel by voxel
#
# Used to check the python tensor fit (dti_mapping_singlecase.py) against the
# MATLAB one (dti_mapping_jelle_singlecase.m), which needs the compiled MEX files.
#
# Usage:
#   python3 compare_dti_maps.py -folder_a <maps_python> -folder_b <maps_matlab>
#
# Only the step 4 maps are compared, not the _withinscanreg ones, since those
# come from the registration step.

import os
import argparse
import numpy as np
import nibabel as nib

DIFFUSION_TIMES = ['0022ms', '0042ms', '0081ms', '0156ms', '0300ms']
MAPS = ['AD', 'RD', 'FA', 'L2', 'L3']


def load(path):
    return np.asarray(nib.load(path).dataobj, dtype=np.float64)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-folder_a', required=True, help='first maps folder')
    parser.add_argument('-folder_b', required=True, help='second maps folder')
    parser.add_argument('-tolerance', type=float, default=1e-4,
                        help='count voxels differing by more than this (default 1e-4)')
    args = parser.parse_args()

    label_a = os.path.basename(args.folder_a.rstrip('/'))
    label_b = os.path.basename(args.folder_b.rstrip('/'))

    print('A: %s' % args.folder_a)
    print('B: %s' % args.folder_b)
    print()
    print('%-6s %-4s %11s %11s %11s %9s %8s' %
          ('time', 'map', 'max|A-B|', 'mean|A-B|', 'p99|A-B|', 'corr', '>tol %'))
    print('-' * 68)

    worst = 0.0
    worst_where = ''
    all_corr = []
    missing = []
    per_map = {m: {'max': 0.0, 'corr': 1.0, 'where': ''} for m in MAPS}

    for dt in DIFFUSION_TIMES:
        for m in MAPS:
            name = '%s_%s.nii.gz' % (dt, m)
            pa = os.path.join(args.folder_a, name)
            pb = os.path.join(args.folder_b, name)
            if not (os.path.exists(pa) and os.path.exists(pb)):
                missing.append(name)
                continue

            a, b = load(pa), load(pb)
            if a.shape != b.shape:
                print('%-6s %-4s  SHAPE MISMATCH %s vs %s' % (dt, m, a.shape, b.shape))
                continue

            # compare where either implementation put signal - background is all
            # zeros in both and would only dilute the statistics
            mask = (a != 0) | (b != 0)
            if not mask.any():
                print('%-6s %-4s  both all zero' % (dt, m))
                continue

            d = np.abs(a[mask] - b[mask])
            mx, mean, p99 = d.max(), d.mean(), np.percentile(d, 99)
            frac = 100.0 * np.count_nonzero(d > args.tolerance) / d.size

            va, vb = a[mask], b[mask]
            if va.std() > 0 and vb.std() > 0:
                corr = np.corrcoef(va, vb)[0, 1]
            else:
                corr = np.nan
            all_corr.append(corr)

            if mx > worst:
                worst, worst_where = mx, '%s %s' % (dt, m)
            if mx > per_map[m]['max']:
                per_map[m]['max'] = mx
                per_map[m]['where'] = dt
            if not np.isnan(corr):
                per_map[m]['corr'] = min(per_map[m]['corr'], corr)

            print('%-6s %-4s %11.3e %11.3e %11.3e %9.6f %8.3f' %
                  (dt, m, mx, mean, p99, corr, frac))
        print()

    print('-' * 68)
    if missing:
        print('missing from one folder     : %s' % ', '.join(missing))

    # AD and RD are clipped to [0, 3.5] on write, so agreement there can hide a
    # disagreement that both implementations clipped away - say so explicitly
    print('note: AD/RD/L2/L3 are clipped to [0, 3.5] before saving, so voxels that')
    print('      saturate in both implementations agree by construction. FA is not')
    print('      clipped, so it is the more informative comparison.')

    # ----------------------------------------------------------------- summary
    # compact enough to screenshot straight onto a slide
    scale = ('below single-precision rounding' if worst < 1e-7 else
             'far below measurement precision' if worst < 1e-3 else
             'VISIBLE - worth investigating')

    print()
    print('=' * 52)
    print('  %s  vs  %s' % (label_a, label_b))
    print('=' * 52)
    print('  %-5s %14s %12s %10s' % ('map', 'largest diff', 'worst at', 'min corr'))
    print('  ' + '-' * 48)
    for m in MAPS:
        info = per_map[m]
        if not info['where']:
            continue
        print('  %-5s %14.2e %12s %10.6f'
              % (m, info['max'], info['where'], info['corr']))
    print('  ' + '-' * 48)
    print('  %-5s %14.2e %12s' % ('ALL', worst, worst_where))
    print()
    print('  Largest disagreement anywhere: %.1e' % worst)
    print('  (%s)' % scale)
    print('=' * 52)


if __name__ == '__main__':
    main()
