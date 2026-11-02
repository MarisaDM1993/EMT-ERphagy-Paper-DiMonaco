#!/usr/bin/env python3
"""
Quantify GR-COL4A1 reporter puncta over a folder of confocal images.

This is the general-purpose entry point to the pipeline in `puncta_core.py`.
It replaces the per-panel processing scripts: the quantification is identical,
and which fields belong to which condition and replicate is supplied as a
manifest rather than hard-coded.

Two modes
---------

1. With `--manifest`, the named CSV decides which files are processed and how
   they are labelled. This is how the published panels are reproduced; the
   manifests for them are in `manifests/`. Expected columns:

       filename   file name within --images
       condition  condition label
       replicate  biological replicate (optional, may be blank)
       batch      imaging batch label  (optional, may be blank)
       include    yes/no -- rows marked no are skipped and reported
       note       why a row is excluded, or any other provenance note

   Only `filename` and `condition` are required. A missing `include` column is
   treated as yes throughout.

2. Without `--manifest`, every readable image in `--images` is processed.
   Condition and replicate are then taken from `--condition-regex` if given,
   and otherwise left blank.

Examples
--------

    # reproduce the published HA-P4HA2 panel
    python quantify_puncta.py --images /path/to/raw \\
        --manifest manifests/hap4ha2_panel.csv --out results/hap4ha2

    # quantify a new folder, condition taken from the filename prefix
    python quantify_puncta.py --images /path/to/new_experiment \\
        --condition-regex '^(?P<condition>[^_]+)_' --out results/new --qc

Outputs `<out>/per_image.csv` and `<out>/per_cell.csv`, and with `--qc` an
overlay PNG per field in `<out>/QC/`.

Detection parameters default to the published values. Changing them changes
the results; see the README for the sensitivity analysis over the two
threshold multipliers.
"""
import argparse
import csv
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puncta_core import (CH_IDX, load_image, subtract_background,  # noqa: E402
                         detect_puncta_split, classify_coloc, segment_nuclei,
                         segment_cells, assign_puncta_to_cells, scale01,
                         green_magenta_merge, build_qc_overlay)

warnings.filterwarnings('ignore')

IMAGE_EXTENSIONS = ('.lsm', '.czi', '.tif', '.tiff')

PER_IMAGE_COLUMNS = [
    'filename', 'condition', 'replicate', 'batch',
    'total_RFP_puncta', 'colocalized_n', 'redonly_n',
    'frac_redonly', 'frac_colocalized',
    'nuclei_n', 'cells_segmented_n', 'cells_real_n', 'cells_sliver_n',
    'puncta_unassigned_to_cell',
    'total_per_nucleus', 'coloc_per_nucleus', 'redonly_per_nucleus',
    'n_blobs_split', 'n_subpuncta_recovered',
]

PER_CELL_COLUMNS = [
    'filename', 'condition', 'replicate', 'batch',
    'cell_id', 'cell_area_px',
    'total_RFP_puncta', 'colocalized_n', 'redonly_n',
    'frac_redonly', 'frac_colocalized',
]


def read_image(path):
    """Load a multi-channel field. `.lsm` and `.czi` go through puncta_core;
    plain TIFF is read with the same reader `.lsm` uses, since both are TIFF."""
    low = path.lower()
    if low.endswith(('.lsm', '.czi')):
        return load_image(path)
    if low.endswith(('.tif', '.tiff')):
        import tifffile
        arr = np.squeeze(tifffile.imread(path))
        if arr.ndim != 3:
            raise ValueError(f'expected a 3-channel image, got shape {arr.shape}')
        if arr.shape[0] > arr.shape[-1]:      # channels-last -> channels-first
            arr = np.moveaxis(arr, -1, 0)
        return {name: arr[idx].astype(np.float32) for name, idx in CH_IDX.items()}
    raise ValueError(f'unsupported file type: {path}')


def read_manifest(path):
    with open(path, newline='', encoding='utf-8-sig') as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f'manifest is empty: {path}')
    missing = {'filename', 'condition'} - set(rows[0])
    if missing:
        sys.exit(f'manifest {path} is missing required column(s): {sorted(missing)}')
    return rows


def discover_images(images_dir, condition_regex):
    pattern = re.compile(condition_regex) if condition_regex else None
    rows = []
    for name in sorted(os.listdir(images_dir)):
        if not name.lower().endswith(IMAGE_EXTENSIONS):
            continue
        condition, replicate = '', ''
        if pattern:
            m = pattern.search(name)
            if m:
                condition = m.groupdict().get('condition', '') or ''
                replicate = m.groupdict().get('replicate', '') or ''
        rows.append({'filename': name, 'condition': condition,
                     'replicate': replicate, 'batch': '', 'include': 'yes', 'note': ''})
    return rows


def quantify_field(path, args):
    ch = read_image(path)
    gfp, rfp, dapi = ch['GFP'], ch['RFP'], ch['DAPI']

    rfp_sub = subtract_background(rfp)
    gfp_sub = subtract_background(gfp)
    puncta, n_blobs_split, n_subpuncta = detect_puncta_split(
        rfp_sub, k=args.k_rfp, min_area=args.min_area, max_area=args.max_area)
    coloc, redonly, _ = classify_coloc(puncta, gfp_sub, k_gfp=args.k_gfp)

    nuc_labels, nuc_markers = segment_nuclei(dapi)
    n_nuclei = int(nuc_labels.max())
    cell_labels = segment_cells(gfp, rfp, nuc_markers)
    n_cells = int(cell_labels.max())

    coloc_ids = assign_puncta_to_cells(coloc, cell_labels)
    redonly_ids = assign_puncta_to_cells(redonly, cell_labels)
    n_unassigned = sum(1 for cid in coloc_ids + redonly_ids if cid == 0)

    tally = {cid: {'coloc': 0, 'redonly': 0} for cid in range(1, n_cells + 1)}
    for cid in coloc_ids:
        if cid > 0:
            tally[cid]['coloc'] += 1
    for cid in redonly_ids:
        if cid > 0:
            tally[cid]['redonly'] += 1

    cells, n_slivers = [], 0
    for cid, counts in tally.items():
        cell_area = int((cell_labels == cid).sum())
        if cell_area < args.min_cell_area:
            n_slivers += 1
            continue
        total_c = counts['coloc'] + counts['redonly']
        cells.append({
            'cell_id': cid, 'cell_area_px': cell_area,
            'total_RFP_puncta': total_c,
            'colocalized_n': counts['coloc'], 'redonly_n': counts['redonly'],
            'frac_redonly': counts['redonly'] / total_c if total_c else np.nan,
            'frac_colocalized': counts['coloc'] / total_c if total_c else np.nan,
        })

    total = len(puncta)
    image = {
        'total_RFP_puncta': total, 'colocalized_n': len(coloc), 'redonly_n': len(redonly),
        'frac_redonly': len(redonly) / total if total else np.nan,
        'frac_colocalized': len(coloc) / total if total else np.nan,
        'nuclei_n': n_nuclei, 'cells_segmented_n': n_cells,
        'cells_real_n': n_cells - n_slivers, 'cells_sliver_n': n_slivers,
        'puncta_unassigned_to_cell': n_unassigned,
        'total_per_nucleus': total / n_nuclei if n_nuclei else np.nan,
        'coloc_per_nucleus': len(coloc) / n_nuclei if n_nuclei else np.nan,
        'redonly_per_nucleus': len(redonly) / n_nuclei if n_nuclei else np.nan,
        'n_blobs_split': n_blobs_split, 'n_subpuncta_recovered': n_subpuncta,
    }
    return image, cells, (gfp, rfp, coloc, redonly, cell_labels)


def main():
    p = argparse.ArgumentParser(
        description='Quantify GR-COL4A1 reporter puncta over a folder of images.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--images', required=True, metavar='DIR',
                   help='folder holding the raw .lsm/.czi/.tif fields')
    p.add_argument('--manifest', metavar='CSV',
                   help='CSV naming which files to process and how to label them')
    p.add_argument('--out', default='results', metavar='DIR',
                   help='output folder (default: results)')
    p.add_argument('--condition-regex', metavar='RE',
                   help="regex with named groups 'condition' and optionally "
                        "'replicate', applied to the filename when no manifest is given")
    p.add_argument('--qc', action='store_true', help='write a QC overlay PNG per field')
    p.add_argument('--k-rfp', type=float, default=20.0, help='mCherry detection multiplier (default: 20)')
    p.add_argument('--k-gfp', type=float, default=8.0, help='eGFP colocalisation multiplier (default: 8)')
    p.add_argument('--min-area', type=int, default=6, help='minimum punctum area, px² (default: 6)')
    p.add_argument('--max-area', type=int, default=250, help='blob-splitting threshold, px² (default: 250)')
    p.add_argument('--min-cell-area', type=int, default=300,
                   help='minimum cell footprint for per-cell output, px² (default: 300)')
    args = p.parse_args()

    if not os.path.isdir(args.images):
        sys.exit(f'--images is not a folder: {args.images}')

    if args.manifest:
        rows = read_manifest(args.manifest)
        source = f'manifest {args.manifest}'
    else:
        rows = discover_images(args.images, args.condition_regex)
        source = f'folder scan of {args.images}'
        if not rows:
            sys.exit(f'no {"/".join(IMAGE_EXTENSIONS)} files found in {args.images}')

    todo = [r for r in rows if (r.get('include') or 'yes').strip().lower() in ('yes', 'y', 'true', '1')]
    skipped = [r for r in rows if r not in todo]

    print(f'Source: {source}')
    if args.manifest:
        print(f'{len(rows)} row(s); {len(todo)} to process, '
              f'{len(skipped)} excluded by the manifest')
    else:
        print(f'{len(todo)} image(s) found')
    for r in skipped:
        print(f"  excluded: {r['filename']}" + (f" -- {r['note']}" if r.get('note') else ''))

    os.makedirs(args.out, exist_ok=True)
    qc_dir = os.path.join(args.out, 'QC')
    if args.qc:
        os.makedirs(qc_dir, exist_ok=True)

    image_rows, cell_rows, missing = [], [], []
    for i, r in enumerate(todo, 1):
        fname = r['filename']
        path = os.path.join(args.images, fname)
        if not os.path.isfile(path):
            missing.append(fname)
            print(f'[{i}/{len(todo)}] MISSING {fname}')
            continue
        label = {'filename': fname, 'condition': r.get('condition', ''),
                 'replicate': r.get('replicate', ''), 'batch': r.get('batch', '')}
        image, cells, qc = quantify_field(path, args)
        image_rows.append({**label, **image})
        cell_rows.extend({**label, **c} for c in cells)
        print(f"[{i}/{len(todo)}] {fname}  condition={label['condition'] or '-'}  "
              f"puncta={image['total_RFP_puncta']}  nuclei={image['nuclei_n']}  "
              f"cells={image['cells_real_n']}")

        if args.qc:
            gfp, rfp, coloc, redonly, cell_labels = qc
            merge = green_magenta_merge(scale01(gfp, np.percentile(gfp, 99.5) or 1),
                                        scale01(rfp, np.percentile(rfp, 99.5) or 1))
            title = f"{fname}  {label['condition']}".strip()
            build_qc_overlay(merge, coloc, redonly, title, cell_labels=cell_labels).save(
                os.path.join(qc_dir, os.path.splitext(fname)[0] + '_QC.png'))

    if not image_rows:
        sys.exit('nothing was processed -- check --images and the manifest filenames')

    img_df = pd.DataFrame(image_rows, columns=PER_IMAGE_COLUMNS)
    cell_df = pd.DataFrame(cell_rows, columns=PER_CELL_COLUMNS)
    sort_cols = [c for c in ('condition', 'replicate', 'batch', 'filename') if c in img_df]
    img_df = img_df.sort_values(sort_cols).reset_index(drop=True)
    if not cell_df.empty:
        cell_df = cell_df.sort_values(sort_cols + ['cell_id']).reset_index(drop=True)

    img_path = os.path.join(args.out, 'per_image.csv')
    cell_path = os.path.join(args.out, 'per_cell.csv')
    img_df.to_csv(img_path, index=False)
    cell_df.to_csv(cell_path, index=False)

    print(f'\n{len(img_df)} field(s), {len(cell_df)} cell(s) written')
    print(f'  {img_path}\n  {cell_path}')
    if args.qc:
        print(f'  {qc_dir}/')
    if missing:
        print(f'\nWARNING: {len(missing)} manifest file(s) not found in {args.images}:')
        for m in missing:
            print(f'  {m}')
        sys.exit(1)


if __name__ == '__main__':
    main()
