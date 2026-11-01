// Split-GFP / mCherry-COL4A1 colocalization analysis -- PER CELL (one row per nucleus)
// Sec24D(nG1-10) x FAM134A/B/C(WT/LIRmut) or ARL6IP1(nG11) split-GFP reconstitution
// vs mCherry-COL4A1, DAPI nuclei. Single confocal slices.
//
// This is the per-cell counterpart of split_gfp_colocalization.ijm (the per-field
// version). Both exist side by side -- use whichever matches how you want to present
// results; they measure the same underlying thing at different resolutions.
//
// HOW THIS DIFFERS FROM THE PER-FIELD VERSION:
//   - Nuclei are segmented from DAPI as before, but each nucleus is then expanded into an
//     approximate "cell region" by a FIXED RADIUS (cellExpandRadiusUm in CONFIG, default
//     10 um) using ImageJ's standard run("Enlarge...") ROI command -- one simple, reliable
//     operation per nucleus. (Two more ambitious "uncapped tessellation" approaches --
//     first the core run("Voronoi") command, then a hand-built label-dilation replacement --
//     were tried and abandoned after repeated failures; see VALIDATED note below. This fixed-
//     radius approach is simpler and far more reliable, at the cost of two closely-packed
//     cells' regions being able to overlap -- see the KNOWN CAVEAT below.)
//   - GFP/mCherry are still thresholded ONCE PER FIELD (Otsu on the whole field), not
//     separately per cell. This is deliberate: a single cell's pixel sample is small and
//     noisy, and re-running Otsu on that little a sample risks the same threshold
//     instability that was found to badly distort the per-field puncta counts (see the
//     puncta note below). Computing the threshold on the full field, then measuring each
//     cell's Manders/PCC against that shared, more stable threshold, avoids reintroducing
//     that problem.
//   - Manders' M1/M2 and Pearson's r are reported PER CELL. GFP puncta metrics are
//     DELIBERATELY NOT included in this per-cell version yet -- n_gfp_puncta was found to
//     correlate with the per-field Otsu threshold at Spearman r=-0.85 (p<0.0001) in your
//     first batch, meaning raw counts were dominated by threshold instability rather than
//     real punctum density. That needs a threshold floor (or a brightness-independent
//     detector) before puncta numbers should be trusted at any resolution, field or cell.
//     Manders/PCC did NOT show this problem (both had |r|<0.15, p>0.5 against threshold
//     and against nucleus count), which is why they're the metrics reported here.
//
// KNOWN CAVEAT: for two nuclei closer together than ~2x cellExpandRadiusUm (e.g. a recently
// divided pair, or a crowded field), the expanded regions can overlap -- inspect the QC
// overlay for each field and treat overlapping/crowded cells with appropriate caution rather
// than trusting the automatic expansion blindly. Tune cellExpandRadiusUm in CONFIG if your QC
// overlays show regions that look consistently too small (clipping real signal) or too large
// (regularly overlapping neighbours).
//
// VALIDATED: the Manders/PCC computation reuses formulas already cross-checked against Python
// (scikit-image) and, for PCC, against numpy.corrcoef to 7 decimal places. Two earlier, more
// ambitious versions of this .ijm file tried to build an UNCAPPED nearest-nucleus tessellation
// (every pixel in the field assigned to its nearest nucleus, no radius limit) instead of a
// fixed radius: first via the core run("Voronoi") command, which -- despite producing the
// expected thin-line pixel pattern on this Fiji install (confirmed via a diagnostic print) --
// consistently produced 0 usable regions downstream of it on every real file; then via a
// hand-built label-dilation replacement, which surfaced a real "wrong nucleus wins" bug caught
// before shipping and then ran too slowly to finish a batch. Given repeated failures across
// both approaches, this reverted to the simpler fixed-radius method (your original choice from
// when this analysis was first scoped). A separate, unrelated real bug was also found and
// fixed along the way: the first nucleus in every field was being measured with no ROI applied
// at all (window/selection order was reversed), giving it a bogus ~whole-field area and
// duplicate/missing rows -- fixed by reordering. Please run this on a small test folder first
// and check the Log window, QC overlays, and CSV before batch-processing everything.
//
// HOW TO RUN: same as the per-field version -- Plugins > Macros > Run..., pick this
// file, then pick your input folder and output folder when prompted. Loops over every
// matching file in the input folder automatically.
// =====================================================================================


// ============================== CONFIG -- EDIT THIS ================================

fileExt = ".lsm";

// Channel indices as they appear after Bio-Formats import (1-based).
chMCherry = 1;
chGFP     = 2;
chDAPI    = 3;

// Pixel size fallback, only used if Bio-Formats fails to import calibration.
DEFAULT_PIXEL_SIZE_UM = 0.07058;

// --- nucleus segmentation ---
nucGaussianSigmaPx = 2;
nucMinDiamUm       = 8;
nucMaxDiamUm       = 30;
excludeEdgeNuclei  = true;

// --- cell region: how far to expand each nucleus outward to approximate the cell body ---
// Typical distance from nuclear envelope to plasma membrane for cultured cells is roughly
// 8-15 um; adjust based on what you see in the QC overlays (QC_overlays/ in your output
// folder) -- too small clips real juxtanuclear signal, too large starts overlapping neighbours.
cellExpandRadiusUm = 10;

// --- GFP / mCherry preprocessing ---
medianRadiusPx      = 1;
rollingBallRadiusPx = 25;

// --- output ---
saveOverlayImages = true;

// Filename -> condition label lookup (substring match, case-sensitive; matches your
// actual naming scheme with spaces, e.g. "Split GFP v2 A LIRm 1.lsm").
conditionKeys   = newArray("A WT", "A LIRm", "B WT", "B LIRm", "C WT", "C LIRm", "A6IPI", "ARL6IP1");
conditionLabels = newArray("FAM134A_WT", "FAM134A_LIRmut", "FAM134B_WT", "FAM134B_LIRmut", "FAM134C_WT", "FAM134C_LIRmut", "ARL6IP1_negctrl", "ARL6IP1_negctrl");

// ============================== END CONFIG ==========================================


// ---------------------------- setup ----------------------------
inputDir = getDirectory("Select the folder containing your .lsm files");
if (inputDir == "") exit("No input folder selected -- macro cancelled.");

outputDir = getDirectory("Select (or create, via the dialog's New Folder button) the output folder for results");
if (outputDir == "") exit("No output folder selected -- macro cancelled.");

overlayDir = outputDir + "QC_overlays" + File.separator;
if (saveOverlayImages && !File.exists(overlayDir)) File.makeDirectory(overlayDir);

setBatchMode(true);
setOption("BlackBackground", true);
run("Set Measurements...", "area mean integrated centroid redirect=None decimal=4");

resultsPath = outputDir + "percell_results.csv";
header = "filename,condition,cell_id,n_cells_in_field,nucleus_area_um2,cell_region_area_um2,region_ok,"
       + "gfp_mean_bgsub,mch_mean_bgsub,gfp_threshold,mch_threshold,"
       + "PCC,Manders_M1_GFPoverlapMCH,Manders_M2_MCHoverlapGFP\n";
File.saveString(header, resultsPath);

fileListAll = getFileList(inputDir);
nProcessed = 0;
for (f = 0; f < fileListAll.length; f++) {
    fname = fileListAll[f];
    if (!endsWith(fname, fileExt)) continue;
    print("=== Processing " + fname + " ===");
    processImage(inputDir + fname, fname);
    nProcessed++;
}
print("Done. Processed " + nProcessed + " files. Results: " + resultsPath);
setBatchMode(false);


// ==================================================================================
// MAIN PER-IMAGE FUNCTION
// ==================================================================================
function processImage(path, fname) {

    run("Bio-Formats Importer", "open=[" + path + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
    srcTitle = getTitle();

    getPixelSize(unit, pw, ph);
    if (pw == 1 && ph == 1) {
        print("  WARNING: no calibration found in file, applying DEFAULT_PIXEL_SIZE_UM=" + DEFAULT_PIXEL_SIZE_UM);
        run("Properties...", "unit=um pixel_width=" + DEFAULT_PIXEL_SIZE_UM + " pixel_height=" + DEFAULT_PIXEL_SIZE_UM);
        getPixelSize(unit, pw, ph);
    }
    W = getWidth(); H = getHeight();

    selectWindow(srcTitle);
    run("Duplicate...", "duplicate channels=" + chDAPI);
    rename("DAPI_raw");
    selectWindow(srcTitle);
    run("Duplicate...", "duplicate channels=" + chGFP);
    rename("GFP_raw");
    selectWindow(srcTitle);
    run("Duplicate...", "duplicate channels=" + chMCherry);
    rename("MCH_raw");
    selectWindow(srcTitle);
    close();

    selectWindow("DAPI_raw"); getStatistics(nul, meanD);
    selectWindow("GFP_raw");  getStatistics(nul, meanG);
    selectWindow("MCH_raw");  getStatistics(nul, meanM);
    print("  channel means -- DAPI(ch" + chDAPI + ")=" + d2s(meanD,2) + "  GFP(ch" + chGFP + ")=" + d2s(meanG,2) + "  mCherry(ch" + chMCherry + ")=" + d2s(meanM,2));
    if (meanD < meanG || meanD < meanM) {
        print("  WARNING: DAPI channel mean is not the highest -- double check chDAPI/chGFP/chMCherry are correct for this file!");
    }

    selectWindow("GFP_raw");
    run("Duplicate...", "title=GFP_proc");
    run("Median...", "radius=" + medianRadiusPx);
    run("Subtract Background...", "rolling=" + rollingBallRadiusPx);

    selectWindow("MCH_raw");
    run("Duplicate...", "title=MCH_proc");
    run("Median...", "radius=" + medianRadiusPx);
    run("Subtract Background...", "rolling=" + rollingBallRadiusPx);

    // ============================== NUCLEUS SEGMENTATION ==============================
    selectWindow("DAPI_raw");
    run("Duplicate...", "title=DAPI_mask");
    run("Gaussian Blur...", "sigma=" + nucGaussianSigmaPx);
    setAutoThreshold("Otsu dark");
    run("Convert to Mask");
    run("Fill Holes");
    run("Watershed");

    minAreaUm2 = PI * pow(nucMinDiamUm / 2, 2);
    maxAreaUm2 = PI * pow(nucMaxDiamUm / 2, 2);
    excludeStr = "";
    if (excludeEdgeNuclei) excludeStr = "exclude";

    // DIAGNOSTIC: some files are coming back with 0 nuclei despite visibly having them. Before
    // applying both filters together, check separately how many candidate blobs exist at each
    // stage, so the Log tells us whether this is thresholding/segmentation (finds nothing at
    // all), the size filter (finds blobs, but outside the configured um diameter range -- e.g.
    // a calibration mismatch for that file), or edge exclusion (correctly-sized blobs that all
    // happen to touch the image border in a tightly-framed field).
    selectWindow("DAPI_mask");
    roiManager("reset");
    run("Analyze Particles...", "size=0-Infinity clear add");
    nRaw = roiManager("count");
    selectWindow("DAPI_mask");
    roiManager("reset");
    run("Analyze Particles...", "size=" + minAreaUm2 + "-" + maxAreaUm2 + " clear add");
    nSizeOnly = roiManager("count");
    selectWindow("DAPI_mask");
    roiManager("reset");
    run("Analyze Particles...", "size=0-Infinity " + excludeStr + " clear add");
    nEdgeOnly = roiManager("count");
    print("  nucleus segmentation diagnostic: " + nRaw + " raw blobs found (any size), " + nSizeOnly
        + " pass the " + nucMinDiamUm + "-" + nucMaxDiamUm + " um diameter size filter alone, " + nEdgeOnly
        + " pass edge-exclusion alone  [size filter = " + d2s(minAreaUm2,1) + "-" + d2s(maxAreaUm2,1) + " um^2; pixel size pw=" + pw + " " + unit + "]");

    selectWindow("DAPI_mask");
    roiManager("reset");
    run("Analyze Particles...", "size=" + minAreaUm2 + "-" + maxAreaUm2 + " " + excludeStr + " clear add");
    nNuclei = roiManager("count");
    print("  " + nNuclei + " nuclei retained after size/edge filtering");
    if (nNuclei == 0) {
        print("  WARNING: no nuclei retained -- skipping file");
        cleanupWindows();
        return;
    }

    nucArea = newArray(nNuclei);
    nucCX   = newArray(nNuclei);
    nucCY   = newArray(nNuclei);
    for (i = 0; i < nNuclei; i++) {
        // IMPORTANT: selectWindow() MUST come before roiManager("select", ...). Each image
        // window keeps its own independent selection state in ImageJ; roiManager("select")
        // only applies the ROI to whichever window is active AT THAT MOMENT. The previous
        // version called roiManager("select", i) first and selectWindow("DAPI_raw") second,
        // so the ROI landed on the wrong (still-active) window and DAPI_raw was measured with
        // whatever selection it already had -- none, for the first nucleus in every field,
        // which measured the WHOLE IMAGE instead (explaining the ~18000 um2 "nucleus" and the
        // missing/duplicated cell_id=1 rows in your last CSV). Fixed by switching windows first.
        selectWindow("DAPI_raw");
        roiManager("select", i);
        roiManager("rename", "Nuc_" + (i + 1));
        run("Measure");
        // NOTE: previously used List.setMeasurements()+List.getValue("X"/"Y") here to get the
        // centroid, which is a second, separate measurement call and was the likely source of
        // nuclei failing to match any region (a bad/mismatched centroid would land outside every
        // region, however large the fallback search radius). Read X/Y straight from the Results
        // table instead, exactly like Area is read on the line above -- same row, no ambiguity.
        nucArea[i] = getResult("Area", nResults - 1);
        nucCX[i] = getResult("X", nResults - 1) / pw;
        nucCY[i] = getResult("Y", nResults - 1) / ph;
    }
    for (i = 0; i < nNuclei; i++) {
        print("  nucleus " + (i + 1) + " centroid: (" + d2s(nucCX[i],1) + ", " + d2s(nucCY[i],1) + ") px  [image is " + W + " x " + H + " px]");
    }

    // ============================== CELL REGIONS: FIXED-RADIUS NUCLEUS EXPANSION ==============
    // Switched away from the uncapped tessellation entirely (both run("Voronoi") and the
    // custom label-dilation replacement caused real trouble across several rounds -- an
    // unexplained 0-regions failure, then a genuine "wrong nucleus wins" bug caught before
    // shipping, then a performance problem). This is back to your ORIGINAL choice from when we
    // first scoped this analysis: expand each nucleus's own outline outward by a fixed radius
    // to approximate the cell body, using ImageJ's standard run("Enlarge...") ROI command --
    // one simple, extremely well-tested operation per nucleus, no loops, no thresholding, no
    // custom pixel-growing logic.
    // CAVEAT: two nuclei closer together than ~2x cellExpandRadiusUm will get overlapping
    // regions (both cells' expanded outlines clip into each other) -- inspect the QC overlay
    // and treat overlapping/crowded cells with appropriate caution, same as the existing caveat
    // for touching/dividing nuclei.
    cellExpandPx = round(cellExpandRadiusUm / pw);
    // roiManager currently holds the nNuclei nucleus ROIs at indices 0..nNuclei-1 (added during
    // nucleus segmentation, untouched since). Enlarge a copy of each and append as new entries
    // at nNuclei..2*nNuclei-1, so regionForNucleus[i] = nNuclei + i.
    for (i = 0; i < nNuclei; i++) {
        selectWindow("DAPI_raw");
        roiManager("select", i);
        // explicit " pixel" suffix -- without it, a calibrated image can make Enlarge... interpret
        // the number in calibrated units instead of pixels, which is exactly the kind of silent
        // unit-ambiguity that caused trouble earlier in this file. cellExpandPx is already
        // converted from cellExpandRadiusUm to pixels above, so force pixel interpretation here.
        run("Enlarge...", "enlarge=" + cellExpandPx + " pixel");
        roiManager("add");
    }
    run("Select None");
    nRegions = nNuclei;
    regionOk = 1;   // fixed-radius expansion always succeeds; kept as a column for consistency/QC

    regionForNucleus = newArray(nNuclei);
    for (i = 0; i < nNuclei; i++) regionForNucleus[i] = nNuclei + i;
    // roiManager list is now STABLE for the rest of this image: indices 0..nNuclei-1 are the
    // raw nucleus outlines (unused from here on), nNuclei..2*nNuclei-1 are the expanded cell
    // regions actually used for measurement (regionForNucleus[i] points to the right one).

    if (saveOverlayImages) {
        selectWindow("MCH_raw");
        run("Duplicate...", "title=QC_overlay");
        run("Enhance Contrast", "saturated=0.5");
        run("Red");
        run("RGB Color");
        for (r = 0; r < nRegions; r++) {
            roiManager("select", nNuclei + r);
            run("Add Selection...");
        }
        run("Flatten");
        saveAs("PNG", overlayDir + replace(fname, fileExt, "") + "_QC.png");
        close();
        close("QC_overlay");
    }
    run("Clear Results");

    // ============================== WHOLE-FIELD THRESHOLDS + MASKS (shared by all cells) ==============================
    selectWindow("GFP_proc");
    run("Select None");
    setAutoThreshold("Otsu dark");
    getThreshold(gfpLo, gfpHi);
    resetThreshold();

    selectWindow("MCH_proc");
    run("Select None");
    setAutoThreshold("Otsu dark");
    getThreshold(mchLo, mchHi);
    resetThreshold();

    selectWindow("GFP_proc");
    run("Duplicate...", "title=GFP_pos");
    setThreshold(gfpLo, gfpHi);
    run("Convert to Mask");

    selectWindow("MCH_proc");
    run("Duplicate...", "title=MCH_pos");
    setThreshold(mchLo, mchHi);
    run("Convert to Mask");

    // product images for the fast, ImageCalculator-based Pearson's r (built once, reused per cell)
    imageCalculator("Multiply create 32-bit", "GFP_proc", "MCH_proc");
    rename("XY_prod");
    imageCalculator("Multiply create 32-bit", "GFP_proc", "GFP_proc");
    rename("X2_prod");
    imageCalculator("Multiply create 32-bit", "MCH_proc", "MCH_proc");
    rename("Y2_prod");

    condition = classifyCondition(fname);

    // ============================== PER-CELL MEASUREMENTS ==============================
    for (i = 0; i < nNuclei; i++) {
        regIdx = regionForNucleus[i];
        if (regIdx == -1) {
            // dead code with the fixed-radius method (regionForNucleus[i] is always nNuclei+i,
            // never -1) -- left in defensively in case this section is edited again later.
            print("  WARNING: nucleus " + (i + 1) + " has no assigned cell region, skipping");
            continue;
        }

        selectWindow("GFP_proc");
        roiManager("select", regIdx);
        run("Measure");
        cellAreaUm2    = getResult("Area", nResults - 1);
        gfpMean        = getResult("Mean", nResults - 1);
        gfpIntDenTotal = getResult("IntDen", nResults - 1);

        selectWindow("MCH_proc");
        roiManager("select", regIdx);
        run("Measure");
        mchMean        = getResult("Mean", nResults - 1);
        mchIntDenTotal = getResult("IntDen", nResults - 1);

        // --- Manders M1: IntDen(GFP) on (region AND MCH-positive) / IntDen(GFP) over whole region ---
        selectWindow("MCH_pos");
        setThreshold(255, 255);
        run("Create Selection");
        resetThreshold();
        roiManager("add");
        mchPosIdx = roiManager("count") - 1;
        roiManager("select", newArray(regIdx, mchPosIdx));
        roiManager("AND");
        selectWindow("GFP_proc");
        run("Restore Selection");
        run("Measure");
        gfpIntDenOverlap = getResult("IntDen", nResults - 1);
        run("Select None");
        roiManager("select", mchPosIdx);
        roiManager("delete");
        M1 = NaN;
        if (gfpIntDenTotal > 0) M1 = gfpIntDenOverlap / gfpIntDenTotal;

        // --- Manders M2: IntDen(mCherry) on (region AND GFP-positive) / IntDen(mCherry) over whole region ---
        selectWindow("GFP_pos");
        setThreshold(255, 255);
        run("Create Selection");
        resetThreshold();
        roiManager("add");
        gfpPosIdx = roiManager("count") - 1;
        roiManager("select", newArray(regIdx, gfpPosIdx));
        roiManager("AND");
        selectWindow("MCH_proc");
        run("Restore Selection");
        run("Measure");
        mchIntDenOverlap = getResult("IntDen", nResults - 1);
        run("Select None");
        roiManager("select", gfpPosIdx);
        roiManager("delete");
        M2 = NaN;
        if (mchIntDenTotal > 0) M2 = mchIntDenOverlap / mchIntDenTotal;

        // --- Pearson's r within this cell's region, via the pre-built product images ---
        selectWindow("XY_prod"); roiManager("select", regIdx); getStatistics(nul, meanXY);
        selectWindow("X2_prod"); roiManager("select", regIdx); getStatistics(nul, meanX2);
        selectWindow("Y2_prod"); roiManager("select", regIdx); getStatistics(nul, meanY2);
        run("Select None");
        covXY = meanXY - gfpMean * mchMean;
        varX  = meanX2 - gfpMean * gfpMean;
        varY  = meanY2 - mchMean * mchMean;
        PCC = NaN;
        if (varX > 0 && varY > 0) PCC = covXY / sqrt(varX * varY);

        row = fname + "," + condition + "," + (i + 1) + "," + nNuclei + "," + nucArea[i] + "," + cellAreaUm2 + "," + regionOk + ","
            + gfpMean + "," + mchMean + "," + gfpLo + "," + mchLo + ","
            + PCC + "," + M1 + "," + M2 + "\n";
        File.append(row, resultsPath);
    }

    close("XY_prod"); close("X2_prod"); close("Y2_prod");
    close("GFP_pos"); close("MCH_pos");

    cleanupWindows();
}

function classifyCondition(fname) {
    for (k = 0; k < conditionKeys.length; k++) {
        if (indexOf(fname, conditionKeys[k]) >= 0) return conditionLabels[k];
    }
    print("  WARNING: filename did not match any entry in conditionKeys -- " + fname);
    return "UNCLASSIFIED_CHECK_FILENAME";
}

function cleanupWindows() {
    imgList = getList("image.titles");
    for (i = 0; i < imgList.length; i++) {
        selectWindow(imgList[i]);
        close();
    }
    roiManager("reset");
    run("Clear Results");
}