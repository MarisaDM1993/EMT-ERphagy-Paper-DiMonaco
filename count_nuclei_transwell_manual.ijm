// Count nuclei in a folder of transwell images, with MANUAL threshold per image
//
// For each image the macro pauses, shows the (background-subtracted) image with a
// suggested threshold, and waits for you to adjust it in the Threshold window.
// Everything else (masking, watershed, debris filtering, counting, QC, saving) is automatic.
//
// Blinding: images are shown in random order and titled "Image k of n" so you set
// thresholds without knowing the condition. Filenames are only revealed in the saved CSVs.
//
// Outputs, saved in the image folder:
//   nuclei_counts.csv  - per image: nuclei kept, objects rejected, threshold you used
//   nuclei_objects.csv - per object: area, signal, solidity, kept (1/0)
//   QC/                - each image with kept objects in green, rejected in red

// ---------------- SETTINGS ----------------
var darkBackground = true;  // true = DAPI/Hoechst on dark; false = crystal violet on light
var rollingBall   = 50;     // background subtraction radius (pixels); 0 = off
var blurSigma     = 1.0;    // smoothing before thresholding; 0 = off
var minSize       = 80;     // pixel^2 (um^2 if calibrated)
var maxSize       = 2000;
var minCirc       = 0.50;
var minSolidity   = 0.85;   // set to 0 to switch the solidity filter off
var minSignal     = 40;     // set to 0 to switch the intensity filter off
var excludeEdges  = true;
var useWatershed  = true;
var saveQC        = true;
var blind         = true;   // random order + coded window titles
// ------------------------------------------

// Globals used inside the function (ImageJ functions can only see variables declared with var)
var qcDir = "";
var lastLower = -1;
var lastUpper = -1;
var cImage = 0; var cNuclei = 0; var cRejected = 0; var cLower = 0; var cUpper = 0;
var oImage = 0; var oArea = 0; var oSignal = 0; var oSol = 0; var oKept = 0;

inputDirectory = getDirectory("Choose a Directory of Images");
fileList = getFileList(inputDirectory);
qcDir = inputDirectory + "QC" + File.separator;
if (saveQC) File.makeDirectory(qcDir);

// Collect image files only
images = newArray(0);
for (i = 0; i < fileList.length; i++) {
    lower = toLowerCase(fileList[i]);
    if (endsWith(fileList[i], "/")) continue;
    if (endsWith(lower, ".tif") || endsWith(lower, ".tiff") || endsWith(lower, ".png")
        || endsWith(lower, ".jpg") || endsWith(lower, ".jpeg"))
        images = Array.concat(images, fileList[i]);
}
n = images.length;
if (n == 0) exit("No image files found in:\n" + inputDirectory);

// Shuffle for blinding
if (blind) {
    for (i = n - 1; i > 0; i--) {
        j = floor(random() * (i + 1));
        tmp = images[i]; images[i] = images[j]; images[j] = tmp;
    }
}

cImage = newArray(0); cNuclei = newArray(0); cRejected = newArray(0); cLower = newArray(0); cUpper = newArray(0);
oImage = newArray(0); oArea = newArray(0); oSignal = newArray(0); oSol = newArray(0); oKept = newArray(0);

setOption("BlackBackground", true);
run("Clear Results");
roiManager("reset");

for (k = 0; k < n; k++)
    processImage(inputDirectory, images[k], k + 1, n);

setBatchMode(false);
roiManager("reset");

// Build and save the tables only at the end, so filenames stay hidden while thresholding
Table.create("Counts");
Table.setColumn("Image", cImage);
Table.setColumn("Nuclei", cNuclei);
Table.setColumn("Rejected", cRejected);
Table.setColumn("ThresholdLower", cLower);
Table.setColumn("ThresholdUpper", cUpper);
Table.save(inputDirectory + "nuclei_counts.csv", "Counts");

if (oImage.length > 0) {
    Table.create("Objects");
    Table.setColumn("Image", oImage);
    Table.setColumn("Area", oArea);
    Table.setColumn("Signal", oSignal);
    Table.setColumn("Solidity", oSol);
    Table.setColumn("Kept", oKept);
    Table.save(inputDirectory + "nuclei_objects.csv", "Objects");
}

showMessage("Done", n + " images processed.\nCounts saved to:\n" + inputDirectory + "nuclei_counts.csv");


function processImage(dir, name, k, n) {
    // Open and pre-process hidden, so the filename never appears on screen
    setBatchMode(true);
    open(dir + name);
    title = getTitle();
    if (nSlices > 1) print("WARNING: image " + k + " is a stack/multichannel - only the first slice is used");
    run("8-bit");
    if (rollingBall > 0) {
        if (darkBackground) run("Subtract Background...", "rolling=" + rollingBall);
        else run("Subtract Background...", "rolling=" + rollingBall + " light");
    }
    rename("original");

    maskTitle = "Image " + k + " of " + n;
    run("Duplicate...", "title=[" + maskTitle + "]");
    if (blurSigma > 0) run("Gaussian Blur...", "sigma=" + blurSigma);

    // Starting threshold: the one you used on the previous image, otherwise auto
    if (lastLower >= 0) setThreshold(lastLower, lastUpper);
    else if (darkBackground) setAutoThreshold("Default dark");
    else setAutoThreshold("Default");

    // ---- Manual step ----
    setBatchMode("exit and display");
    selectImage(maskTitle);
    run("Threshold...");
    waitForUser("Set threshold (" + k + " of " + n + ")",
        "Adjust the sliders in the Threshold window so nuclei are red\n" +
        "and background/debris are not. Then click OK here.\n \n" +
        "Do NOT click 'Apply' in the Threshold window.");
    selectImage(maskTitle);
    getThreshold(lower, upper);
    if (lower == -1) {
        print("Image " + k + ": no threshold set, auto threshold used");
        if (darkBackground) setAutoThreshold("Default dark");
        else setAutoThreshold("Default");
        getThreshold(lower, upper);
    }
    lastLower = lower;
    lastUpper = upper;

    // ---- Automatic from here ----
    setBatchMode(true);
    selectImage(maskTitle);
    run("Convert to Mask");
    if (useWatershed) run("Watershed");

    run("Set Measurements...", "area mean shape redirect=original decimal=3");
    run("Clear Results");
    roiManager("reset");
    opts = "size=" + minSize + "-" + maxSize + " circularity=" + minCirc + "-1.00 show=Nothing display clear add";
    if (excludeEdges) opts = opts + " exclude";
    run("Analyze Particles...", opts);

    nObj = nResults;
    kept = 0;
    selectImage("original");
    if (saveQC) run("RGB Color");
    for (j = 0; j < nObj; j++) {
        area = getResult("Area", j);
        mean = getResult("Mean", j);
        sol  = getResult("Solidity", j);
        if (darkBackground) signal = mean;
        else signal = 255 - mean;
        keep = (sol >= minSolidity) && (signal >= minSignal);
        if (keep) kept++;

        oImage = Array.concat(oImage, title);
        oArea = Array.concat(oArea, area);
        oSignal = Array.concat(oSignal, signal);
        oSol = Array.concat(oSol, sol);
        oKept = Array.concat(oKept, keep);

        if (saveQC) {
            selectImage("original");
            roiManager("select", j);
            if (keep) Overlay.addSelection("green");
            else Overlay.addSelection("red");
        }
    }

    if (saveQC) {
        selectImage("original");
        run("Select None");
        run("Flatten");
        saveAs("PNG", qcDir + "QC_" + title);
    }

    cImage = Array.concat(cImage, title);
    cNuclei = Array.concat(cNuclei, kept);
    cRejected = Array.concat(cRejected, nObj - kept);
    cLower = Array.concat(cLower, lower);
    cUpper = Array.concat(cUpper, upper);

    run("Clear Results");
    roiManager("reset");
    close("*");
}
