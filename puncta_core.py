"""
Core functions for GR-COL4A1 puncta quantification:
- background subtraction (rolling-ball-like via large-sigma Gaussian)
- robust per-image adaptive thresholding (median + k*MAD)
- RFP punctum detection (connected components, area filter)
- GFP-at-punctum colocalization classification
- DAPI nucleus counting
- QC overlay rendering
"""
import re
import numpy as np
import tifffile
import czifile
from scipy import ndimage as ndi
from skimage import measure, morphology, filters
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from PIL import Image, ImageDraw, ImageFont

CH_IDX = {'RFP': 0, 'GFP': 1, 'DAPI': 2}
PIXEL_SIZE_NM = 70.58
SCALE_BAR_UM = 20
SCALE_BAR_PX = int(round(SCALE_BAR_UM * 1000 / PIXEL_SIZE_NM))

# Fonts are used only to annotate QC overlays (scale bar label, header text);
# no quantification depends on them. Previously these two lines called
# ImageFont.truetype() on a hard-coded Debian/Ubuntu path at module scope, so
# `import puncta_core` raised OSError outright on macOS, Windows, and any Linux
# distribution that installs DejaVu elsewhere -- before a single quantification
# function could be reached. Resolved per-platform below, with a guaranteed
# fallback so that importing this module can never fail on account of a missing
# font. If the fallback is used, overlay text renders in PIL's small default
# bitmap font and ignores the requested size; overlays remain readable and all
# numeric output is identical.
_FONT_CANDIDATES = (
    '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',   # Debian / Ubuntu
    '/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf',            # Fedora / RHEL
    '/usr/share/fonts/TTF/DejaVuSans-Bold.ttf',               # Arch
    '/System/Library/Fonts/Supplemental/Arial Bold.ttf',      # macOS
    '/Library/Fonts/Arial Bold.ttf',                          # macOS (older)
    'C:\\Windows\\Fonts\\arialbd.ttf',                        # Windows
    'DejaVuSans-Bold.ttf',                                    # Pillow's own search path
    'Arial Bold.ttf',
)


FONT_SMALL_SIZE = 26
FONT_HDR_SIZE = 34

_FONT_CACHE = {}


def _load_font(size):
    """Return a bold sans-serif PIL font at `size`, falling back to PIL's
    built-in default rather than raising if no TrueType font can be found.
    Resolved on first use and cached, never at import time."""
    if size in _FONT_CACHE:
        return _FONT_CACHE[size]
    font = None
    for path in _FONT_CANDIDATES:
        try:
            font = ImageFont.truetype(path, size)
            break
        except OSError:
            continue
    if font is None:
        try:
            # Pillow >= 10.1 accepts a size, but satisfies it via an embedded
            # TrueType font, so this path still needs FreeType support.
            font = ImageFont.load_default(size)
        except (TypeError, OSError, ImportError):
            # Older Pillow signature, or a Pillow built without FreeType: the
            # sizeless default is a bitmap font with no such dependency.
            font = ImageFont.load_default()
    _FONT_CACHE[size] = font
    return font


def __getattr__(name):
    """Keep `puncta_core.FONT_SMALL` / `.FONT_HDR` working for any caller that
    referenced them when they were eager module-level constants."""
    if name == 'FONT_SMALL':
        return _load_font(FONT_SMALL_SIZE)
    if name == 'FONT_HDR':
        return _load_font(FONT_HDR_SIZE)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def parse_condition(fname):
    if fname.lower().endswith('.lsm'):
        base = re.sub(r'^[0-9a-f]{8}-(\d+_)?GR-COL4A1_', '', fname)
        base = base.replace('.lsm', '').rstrip('_')
        return base, 'LSM'
    else:
        base = re.sub(r'^[0-9a-f]{8}-', '', fname)
        base = base.replace('.czi', '').replace('GR-COL4A1_', '')
        m = re.match(r'^(.*_2hr)\.(\d+)_czi$', base)
        if m:
            base = f"{m.group(1)}_{m.group(2)}"
        base = base.rstrip('_')
        return base, 'CZI'


def load_image(f):
    if f.lower().endswith('.lsm'):
        arr = tifffile.imread(f)
    else:
        with czifile.CziFile(f) as czi:
            arr = np.squeeze(czi.asarray())
    return {name: arr[idx].astype(np.float32) for name, idx in CH_IDX.items()}


def subtract_background(img, sigma=15, denoise_sigma=1.5):
    """Light denoise (small-sigma Gaussian, to suppress single/few-pixel shot
    noise while preserving true multi-pixel puncta) followed by large-scale
    background estimation and subtraction. Returns the UNCLIPPED residual
    (background pixels fluctuate symmetrically around 0) so that robust noise
    statistics computed on it are meaningful. Clip to >=0 only at the point of
    use (display / intensity sampling)."""
    den = ndi.gaussian_filter(img, sigma=denoise_sigma) if denoise_sigma else img
    bg = ndi.gaussian_filter(den, sigma=sigma)
    return den - bg


def noise_sigma(residual):
    """Robust noise SD estimated from the NEGATIVE half of a background-
    subtracted residual only. Background pixels fluctuate symmetrically
    around 0, but real puncta only ever push the residual positive, so using
    only residual<0 values gives a noise estimate uncontaminated by real
    signal (unlike a plain median/MAD over the whole image, which collapses
    to a near-zero estimate when most pixels are true zero)."""
    neg = residual[residual < 0]
    if neg.size == 0:
        return 1e-6
    return np.median(np.abs(neg)) / 0.6745  # median(|half-normal|) -> sigma


def detect_puncta(rfp_sub, k=20.0, min_area=6, max_area=250):
    sigma_n = noise_sigma(rfp_sub)
    thresh = k * sigma_n
    mask = rfp_sub > thresh
    mask = morphology.remove_small_objects(mask, min_size=min_area)
    labeled = measure.label(mask)
    rfp_clipped = np.clip(rfp_sub, 0, None)
    props = measure.regionprops(labeled, intensity_image=rfp_clipped)
    puncta = [p for p in props if min_area <= p.area <= max_area]
    return puncta, thresh, sigma_n


class _SplitPunctum:
    """Lightweight stand-in for a skimage regionprops object, used for the
    sub-regions produced when an oversized merged blob is split (see
    detect_puncta_split). Only exposes the attributes downstream code
    (classify_coloc, assign_puncta_to_cells, build_qc_overlay) actually
    reads: .coords and .centroid."""
    __slots__ = ('coords', 'centroid', 'area')


def detect_puncta_split(rfp_sub, k=20.0, min_area=6, max_area=250,
                         split_min_distance=4, split_dist_smooth=1.0,
                         min_subpunctum_area=2):
    """Same detection as detect_puncta (robust adaptive threshold + connected
    components), but instead of discarding connected components whose area
    exceeds max_area, splits them into individual puncta via local-maxima
    watershed on the (lightly smoothed) RFP intensity restricted to that
    blob's footprint -- the same strategy already used for touching nuclei
    in segment_nuclei.

    Rationale (2026-10-12 QC): in dense conditions (observed in HA-P4HA2 and
    in pilot conditions not reported in the paper) multiple real, separate
    puncta sit close enough together that their above-threshold RFP signal
    forms one connected component exceeding max_area=250px. detect_puncta's
    `min_area <= area <= max_area` filter then drops that WHOLE component --
    not just caps its size -- which silently discards every punctum inside
    it. In a representative HA-P4HA2 field this amounted to 30 oversized
    components (11,731 px of above-threshold RFP signal, ~31% of all signal
    pixels in the field) with zero representation in the quantification; a
    second representative dense field lost ~40% of signal pixels the same
    way. Visual QC overlay confirmed
    these merged blobs sit on genuine, visually distinct punctate structures,
    and splitting them recovers markers that track that texture well (see
    project notes, 2026-10-12 QC-of-method section).

    A blob that cannot be split into >=2 local maxima (e.g. truly uniform,
    or smaller sub-peaks are not resolvable) is kept as a single oversized
    punctum rather than dropped, so nothing above-threshold is ever silently
    discarded by this function.

    Returns (regions, n_blobs_split, n_subpuncta_recovered) where `regions`
    is a list of punctum-like objects (each with .coords and .centroid,
    usable anywhere a detect_puncta() punctum is used)."""
    sigma_n = noise_sigma(rfp_sub)
    thresh = k * sigma_n
    mask = rfp_sub > thresh
    mask = morphology.remove_small_objects(mask, min_size=min_area)
    labeled = measure.label(mask)
    rfp_clipped = np.clip(rfp_sub, 0, None)
    props = measure.regionprops(labeled, intensity_image=rfp_clipped)

    regions = []
    n_blobs_split = 0
    n_subpuncta = 0
    for p in props:
        if p.area < min_area:
            continue
        if p.area <= max_area:
            regions.append(p)
            continue
        minr, minc, maxr, maxc = p.bbox
        sub_mask = p.image
        sub_int = rfp_clipped[minr:maxr, minc:maxc] * sub_mask
        sub_int_sm = ndi.gaussian_filter(sub_int, sigma=split_dist_smooth) if split_dist_smooth else sub_int
        coords = peak_local_max(sub_int_sm, min_distance=split_min_distance, labels=sub_mask)
        if len(coords) < 2:
            regions.append(p)
            continue
        peak_mask = np.zeros(sub_mask.shape, dtype=bool)
        peak_mask[tuple(coords.T)] = True
        markers = measure.label(peak_mask)
        ws = watershed(-sub_int_sm, markers, mask=sub_mask)
        sub_props = measure.regionprops(ws)
        n_blobs_split += 1
        for sp in sub_props:
            if sp.area < min_subpunctum_area:
                continue
            fp = _SplitPunctum()
            fp.coords = sp.coords + np.array([minr, minc])
            fp.centroid = (sp.centroid[0] + minr, sp.centroid[1] + minc)
            fp.area = sp.area
            regions.append(fp)
            n_subpuncta += 1
    return regions, n_blobs_split, n_subpuncta


def classify_coloc(puncta, gfp_sub, k_gfp=8.0):
    sigma_n = noise_sigma(gfp_sub)
    gfp_thresh = k_gfp * sigma_n
    coloc, redonly = [], []
    for p in puncta:
        coords = p.coords
        mean_gfp = gfp_sub[coords[:, 0], coords[:, 1]].mean()
        (coloc if mean_gfp > gfp_thresh else redonly).append(p)
    return coloc, redonly, gfp_thresh


def segment_nuclei(dapi, smooth_sigma=3, min_area=150, min_distance=30, dist_smooth=4):
    """Multi-Otsu threshold (3 classes: background / dim nuclei / bright
    nuclei, keeping the lower background-vs-everything-else cut) + smoothed-
    distance-transform watershed to separate touching nuclei.

    NOTE (2026-10-08): this replaces a plain global-Otsu threshold, which
    failed badly whenever a field contained one or two much brighter nuclei
    alongside many normal-brightness ones (common in this dataset) -- a
    single bimodal Otsu split calibrates to "bright vs everything else"
    rather than "background vs any nucleus", silently dropping most dim-but-
    real nuclei from the mask entirely (confirmed by visual inspection: e.g.
    sgA_TGFb_1 and sgNT_TGFb_2 in the sg-knockout panel went from 6 and 9
    detected nuclei, respectively -- missing the majority of clearly-focused,
    clearly-DAPI+ nuclei in frame -- to 47 and 32 with this fix, matching the
    image by eye). Multi-Otsu's lower cut is far more robust to that skew
    since the 1-2 bright nuclei form their own (upper) class instead of
    dragging the background/foreground split upward. The distance transform
    is also lightly smoothed before peak-finding (dist_smooth) to curb
    watershed over-fragmentation from internal nuclear texture (nucleoli /
    heterochromatin patches), and min_distance raised 20->30 for the same
    reason. Returns (labeled_image, marker_image) so callers can reuse the
    nucleus seeds (marker_image) to seed whole-cell segmentation."""
    sm = ndi.gaussian_filter(dapi, sigma=smooth_sigma)
    try:
        thr = filters.threshold_multiotsu(sm, classes=3)[0]
    except ValueError:
        # falls back if the field is too uniform/sparse for 3 well-separated classes
        thr = filters.threshold_otsu(sm)
    mask = sm > thr
    mask = morphology.remove_small_objects(mask, min_size=min_area)
    mask = morphology.binary_closing(mask, morphology.disk(5))
    if not mask.any():
        empty = np.zeros_like(dapi, dtype=np.int32)
        return empty, empty
    dist = ndi.distance_transform_edt(mask)
    dist_sm = ndi.gaussian_filter(dist, sigma=dist_smooth) if dist_smooth else dist
    coords = peak_local_max(dist_sm, min_distance=min_distance, labels=mask)
    peak_mask = np.zeros(dist.shape, dtype=bool)
    if len(coords):
        peak_mask[tuple(coords.T)] = True
    markers = measure.label(peak_mask)
    labels_ws = watershed(-dist_sm, markers, mask=mask)
    return labels_ws, markers


def count_nuclei(dapi, **kwargs):
    labels_ws, _ = segment_nuclei(dapi, **kwargs)
    return int(labels_ws.max())


def segment_cells(gfp, rfp, nuc_markers, cyto_sigma=10, fg_frac=0.25, min_area=500):
    """Approximate single-cell territories using the diffuse ER (GFP+RFP)
    signal as a cytoplasm-fill channel, seeded from nucleus markers and
    grown outward by watershed (a 'propagate from nuclei' style
    segmentation, standard when a dedicated membrane/cytoplasm dye is not
    available). Returns a labeled image matching nuc_markers' cell IDs
    (0 = not assigned to any cell)."""
    cyto = ndi.gaussian_filter(gfp.astype(np.float32) + rfp.astype(np.float32), sigma=cyto_sigma)
    if not np.any(nuc_markers):
        return np.zeros_like(gfp, dtype=np.int32)
    cthr = filters.threshold_otsu(cyto)
    fg = cyto > (cthr * fg_frac)
    fg = fg | (nuc_markers > 0)
    fg = morphology.binary_closing(fg, morphology.disk(5))
    fg = morphology.remove_small_objects(fg, min_size=min_area)
    cell_labels = watershed(-cyto, nuc_markers, mask=fg)
    return cell_labels


def assign_puncta_to_cells(puncta, cell_labels):
    """Return a list of cell IDs (one per punctum, same order as `puncta`;
    0 means the punctum's centroid did not fall inside any segmented cell
    territory, e.g. an isolated/ambiguous detection near the field edge)."""
    ids = []
    h, w = cell_labels.shape
    for p in puncta:
        y, x = p.centroid
        yi, xi = int(round(y)), int(round(x))
        yi = min(max(yi, 0), h - 1)
        xi = min(max(xi, 0), w - 1)
        ids.append(int(cell_labels[yi, xi]))
    return ids


def scale01(arr, wp):
    return np.clip(arr.astype(np.float32), 0, wp) / wp


def green_magenta_merge(gfp01, rfp01):
    h, w = gfp01.shape
    rgb = np.zeros((h, w, 3), dtype=np.float32)
    rgb[:, :, 1] += gfp01
    rgb[:, :, 0] += rfp01
    rgb[:, :, 2] += rfp01
    return (np.clip(rgb, 0, 1) * 255).astype(np.uint8)


def add_scale_bar(img, px_len=SCALE_BAR_PX, label=f"{SCALE_BAR_UM} µm"):
    img = img.convert('RGB')
    draw = ImageDraw.Draw(img)
    w, h = img.size
    margin = 30
    bar_y = h - margin
    x1 = w - margin - px_len
    x2 = w - margin
    draw.rectangle([x1 - 2, bar_y - 6, x2 + 2, bar_y + 2], fill=(0, 0, 0))
    draw.rectangle([x1, bar_y - 4, x2, bar_y], fill=(255, 255, 255))
    tw = draw.textlength(label, font=_load_font(FONT_SMALL_SIZE))
    tx = x1 + (px_len - tw) / 2
    ty = bar_y - 34
    for dx, dy in [(-1, -1), (1, -1), (-1, 1), (1, 1)]:
        draw.text((tx + dx, ty + dy), label, font=_load_font(FONT_SMALL_SIZE), fill=(0, 0, 0))
    draw.text((tx, ty), label, font=_load_font(FONT_SMALL_SIZE), fill=(255, 255, 255))
    return img


def build_qc_overlay(merge01_rgb, coloc, redonly, title, cell_labels=None):
    """merge01_rgb: HxWx3 uint8 green/magenta merge. Draw markers + legend,
    and (if cell_labels is given) faint white cell-territory boundaries."""
    arr = merge01_rgb.copy()
    if cell_labels is not None:
        from skimage.segmentation import find_boundaries
        bounds = find_boundaries(cell_labels, mode='outer')
        arr[bounds] = [120, 120, 120]
    img = Image.fromarray(arr, 'RGB')
    img = add_scale_bar(img)
    draw = ImageDraw.Draw(img)
    r = 6
    for p in coloc:
        y, x = p.centroid
        draw.ellipse([x - r, y - r, x + r, y + r], outline=(0, 255, 255), width=2)
    for p in redonly:
        y, x = p.centroid
        draw.ellipse([x - r, y - r, x + r, y + r], outline=(255, 255, 0), width=2)
    # header bar
    w, h = img.size
    hdr = Image.new('RGB', (w, 54), (18, 18, 18))
    hdraw = ImageDraw.Draw(hdr)
    hdraw.text((10, 8), title, font=_load_font(FONT_HDR_SIZE), fill=(255, 255, 255))
    legend = f"cyan=coloc (n={len(coloc)})   yellow=red-only (n={len(redonly)})   gray=cell boundary"
    tw = hdraw.textlength(legend, font=_load_font(FONT_SMALL_SIZE))
    hdraw.text((w - tw - 15, 12), legend, font=_load_font(FONT_SMALL_SIZE), fill=(255, 255, 255))
    out = Image.new('RGB', (w, h + 54), (18, 18, 18))
    out.paste(hdr, (0, 0))
    out.paste(img, (0, 54))
    return out
