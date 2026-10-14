# ncorr_2D_matlab

This fork derives from the [upstream Ncorr repository](https://github.com/justinblaber/ncorr_2D_matlab) for:

```
Ncorr: open-source 2D digital image correlation matlab software
J Blaber, B Adair, A Antoniou
Experimental Mechanics 55 (6), 1105-1122
```

Please cite this paper if you use this software in your research. 

## Changes in this fork

- Supports large (including 8K × 8K) images by using overflow-safe indexing in the region-growing DIC calculation.
- Improves seed setup: seed locations remain visible in large-image previews, and seeds can be added or adjusted using typed full-resolution X/Y coordinates. Entered coordinates snap automatically to the nearest valid DIC grid point.

## Installation

This is a source-only distribution. MATLAB, Image Processing Toolbox and a
supported C++ compiler are required; generated MEX binaries are not included.
MATLAB R2026a with Microsoft Visual C++ 2022 has been runtime-verified.

```sh
git clone --branch matlab_2026a https://github.com/FloZr/ncorr_2D_matlab_8k_images.git
```

In MATLAB, open the cloned folder and run:

```matlab
cd('path/to/ncorr_2D_matlab_8k_images');
addpath(pwd,'-begin');
mex -setup C++
build_info = ncorr_build_mex(true,4); % OpenMP, four GUI worker threads
handles_ncorr = ncorr;
```

MATLAB's current folder outranks `addpath`. Do not launch from a folder
containing another Ncorr copy.

Optional synthetic regression checks, after building:

```matlab
test_ncorr_loadimgs();
test_ncorr_rgdic_spacing();
```

These tests create their own temporary inputs; no measurement dataset or
machine-specific paths are required.

## Building and MATLAB R2026a compatibility

`ncorr_build_mex(false,1)` builds without OpenMP. The builder uses rooted
source/output paths and does not change the working directory. All builds pin
the supported separate-complex API, `-R2017b`, and retain double precision and
compiler optimization defaults. OpenMP flags are mapped for GNU and Microsoft
compilers; other compilers can build serial but cannot silently substitute a
serial build for a requested OpenMP build.

The GUI uses `ncorr_buildinfo.mat` plus `ncorr_installinfo.txt` to validate the
MATLAB release, compiler, API, options, source metadata and required artifacts.
Missing or stale build metadata requests installation. Failed
builds invalidate the stamp. Reinstall works from another working directory
and still asks before discarding the active analysis.

R2026a removes `opengl`; the obsolete software-renderer override is removed.
Existing Image Processing Toolbox ROI tools are retained. The RG-DIC
displacement-jump gate now follows the current call's spacing rather than the
first call since loading the MEX. `test_ncorr_rgdic_spacing` exercises both
mixed-spacing orders with a known affine stretch, coverage and accuracy checks.

### Loading numbered image sequences

Multi-image selections must share a prefix and extension and end in a numeric
suffix, for example `DIC_B_02.tif` through `DIC_B_05.tif`. Prefixes may contain
multiple underscores; the last underscore separates the image number. Images
are sorted numerically in both **Load All** and **Load Lazy**, not by selection
order.

R2026a rejects nonscalar operands to `:`. The old filename parser supplied all
underscore positions and misreported the resulting exception as running out of
memory. The parser now uses scalar positions and `fileparts`; unexpected
multi-image loading exceptions show their actual identifier and message.
`test_ncorr_loadimgs` exercises both modes with real temporary TIFF files and
checks image content in numeric order.

After updating an already-running session, dismiss any existing load-error
dialog and run `clear ncorr_util_loadimgs` before retrying. This MATLAB-only
fix does not require rebuilding MEX files or reinstalling Ncorr.


## Public package and local-only files

- Root: reviewed Ncorr MATLAB/C++ sources, headers, generic regression tests,
  documentation and license.
- `local_private/`: local analysis scripts, measurements, results, diagnostics
  and validation captures. This folder is ignored and is not distributed.
- Generated MEX/object files and `ncorr_buildinfo.mat`/`ncorr_installinfo.txt`
  remain beside the sources for execution, but are ignored and rebuilt locally.

`.gitignore` is an explicit public-file allowlist: new files and directories
are ignored unless deliberately added to it. Keep dataset-specific scripts
and absolute local paths out of the public sources and documentation. Do not
force-add private files. Ignoring a file does not remove earlier published
copies or Git history.

For API details see the
[R2026a MEX documentation](https://www.mathworks.com/help/releases/R2026a/matlab/ref/mex.html).


For convienience, a direct link to the manual is available [here](http://www.ncorr.com/download/ncorrmanual_v1_2_2.pdf) and an instructional step-by-step youtube video is linked below:

[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/cXfeiBXjN18/0.jpg)](http://www.youtube.com/watch?v=cXfeiBXjN18)

