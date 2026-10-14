# Test report {#test_report}

\htmlonly
<div class="test-report">
<section class="test-run">
<p class="test-summary"><strong>198 tests collected</strong>, 198 passed ✅, 0 failed ❌, 0 skipped ⏭️ in 0:01:25s on 14/10/2026 at 19:50:23</p>
<details class="test-file" open><summary><strong title="test_fit_AstigFit.cpp">Fit Astig Fit</strong><span>14/14 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Calibrate</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Estimate</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Estimate Batch Within Limits</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.011 s</td></tr>
<tr><td>Calibrate Bad</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Estimate Bad</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Axial Limits / Estimate Within Limits <small>(9 cases)</small></span><span>✅ 9/9</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Below Lower Bound</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>At Lower Bound</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Near Lower Bound</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Negative Between Grid Points</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>At Center</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Positive Between Grid Points</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Near Upper Bound</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>At Upper Bound</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Above Upper Bound</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_Palm.cpp">Palm</strong><span>33/33 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Black Stack</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Selected Planes</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.108 s</td></tr>
<tr><td>Quadrants</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.144 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Stack <small>(10 cases)</small></span><span>✅ 10/10</span><span>1.523 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>No Fit Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.132 s</td></tr>
<tr><td>No Fit With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.057 s</td></tr>
<tr><td>Gaussian XY Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.171 s</td></tr>
<tr><td>Gaussian XY With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.095 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.138 s</td></tr>
<tr><td>Gaussian Sigma With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.139 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.161 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.161 s</td></tr>
<tr><td>Gaussian Theta Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.232 s</td></tr>
<tr><td>Gaussian Theta With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.237 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Image <small>(10 cases)</small></span><span>✅ 10/10</span><span>0.198 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>No Fit Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.012 s</td></tr>
<tr><td>No Fit With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.008 s</td></tr>
<tr><td>Gaussian XY Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.017 s</td></tr>
<tr><td>Gaussian XY With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.014 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.022 s</td></tr>
<tr><td>Gaussian Sigma With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.018 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.024 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.021 s</td></tr>
<tr><td>Gaussian Theta Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.032 s</td></tr>
<tr><td>Gaussian Theta With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.030 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Value <small>(10 cases)</small></span><span>✅ 10/10</span><span>0.173 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Plane 1</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.018 s</td></tr>
<tr><td>Plane 2</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.018 s</td></tr>
<tr><td>Plane 3</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.016 s</td></tr>
<tr><td>Plane 4</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.022 s</td></tr>
<tr><td>Plane 5</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.017 s</td></tr>
<tr><td>Plane 6</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.017 s</td></tr>
<tr><td>Plane 7</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.019 s</td></tr>
<tr><td>Plane 8</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.017 s</td></tr>
<tr><td>Plane 9</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.012 s</td></tr>
<tr><td>Plane 10</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.017 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_Alignment.cpp">SP Alignment</strong><span>14/14 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Blank Call With Reference</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Blank Call With Factor</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Identity And Reuse</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Integer Translation</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Transpose Pixels</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Test Suites / Saves And Reloads Transformed Stack <small>(9 cases)</small></span><span>✅ 9/9</span><span>5.980 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Collapsed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.141 s</td></tr>
<tr><td>Identity</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.446 s</td></tr>
<tr><td>Translate Quarter</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.690 s</td></tr>
<tr><td>Translate Half</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.722 s</td></tr>
<tr><td>Translate Three Quarters</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.853 s</td></tr>
<tr><td>Translate One</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.457 s</td></tr>
<tr><td>Zoom</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.283 s</td></tr>
<tr><td>Transpose</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.253 s</td></tr>
<tr><td>Translate Half Upsampled</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>2.135 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_BSpline.cpp">SP B Spline</strong><span>1/1 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Dirac Impulse</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_FFT.cpp">SP FFT</strong><span>13/13 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Cross Power Spectrum Basic</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Volumes / Matches Direct Discrete Transform <small>(6 cases)</small></span><span>✅ 6/6</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Single Voxel Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Single Voxel 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Odd Rectangle Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Odd Rectangle 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cube Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cube 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Volumes / Round Trip And Parseval <small>(6 cases)</small></span><span>✅ 6/6</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Single Voxel Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Single Voxel 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Odd Rectangle Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Odd Rectangle 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cube Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Cube 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_Interpolation.cpp">SP Interpolation</strong><span>15/15 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Zoom Factor</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Shift Inside Bounds</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Shift Out Of Bounds</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Interpolate Shifted</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
<tr><td>Subpixel Max Interp 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.019 s</td></tr>
<tr><td>Get Shift 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.024 s</td></tr>
<tr><td>Upsample</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Modes / Non Linear Matrix <small>(8 cases)</small></span><span>✅ 8/8</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Trilinear Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Trilinear Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Weight Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Weight Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Spline Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Spline Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_Matrix.cpp">SP Matrix</strong><span>22/22 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Print</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Send Irregular</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Has Na N</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Has Invalid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Is Invalid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Get Valid Index</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Count Valid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Max Ignores Na N</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Sum Valid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Mean Valid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Mean With Mask</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Median</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Quantile</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Rms Difference</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Normalisation With Shift</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Crop Volume</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Extract Sub Volume</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Center Of Mass</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cholesky Success</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cholesky Fail</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Inverse Gauss Jordan</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Inverse LU</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_Spline.cpp">SP Spline</strong><span>10/10 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Set Members</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.029 s</td></tr>
<tr><td>Extract Beads</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.028 s</td></tr>
<tr><td>Align</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>15.995 s</td></tr>
<tr><td>Get Overlap Minimal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Get Overlap</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.163 s</td></tr>
<tr><td>Soft Normalization Minimal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Soft Normalization</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.173 s</td></tr>
<tr><td>Get PSF Minimal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.338 s</td></tr>
<tr><td>Get PSF</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>38.729 s</td></tr>
<tr><td>Calibrate</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>22.050 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_tr_Point.cpp">TR Point</strong><span>1/1 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Base</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_tr_Tracks.cpp">TR Tracks</strong><span>35/35 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Tracking With Gaps</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Tracking Empty Input</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Blinking Empty State</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>MSD With Gaps</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Instant Diffusion</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Analysis Empty State</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Links <small>(20 cases)</small></span><span>✅ 20/20</span><span>0.123 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>No Fit Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>No Fit Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
<tr><td>No Fit With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>No Fit With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Gaussian XY Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Gaussian XY Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Gaussian XY With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
<tr><td>Gaussian XY With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.007 s</td></tr>
<tr><td>Gaussian Sigma With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.007 s</td></tr>
<tr><td>Gaussian Sigma With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.007 s</td></tr>
<tr><td>Gaussian Theta Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.007 s</td></tr>
<tr><td>Gaussian Theta Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.008 s</td></tr>
<tr><td>Gaussian Theta With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.009 s</td></tr>
<tr><td>Gaussian Theta With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Reconnection <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.002 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Stationary</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Diffusive</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Linear</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">MSD <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Planar 2D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Spatial 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Track Analysis / Fit <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Planar 2D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Spatial 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Track Fit / Fit <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Power</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Exponential</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_loc_ATrous.cpp">LOC A Trous</strong><span>12/12 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Shapes / Reconstructs Input <small>(4 cases)</small></span><span>✅ 4/4</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Single Pixel</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Single Row</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Single Column</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Rectangle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Shapes / First Level Matches Direct Convolution <small>(4 cases)</small></span><span>✅ 4/4</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Single Pixel</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Single Row</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Single Column</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Rectangle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Intensities / Details Vanish And Residual Is Preserved <small>(4 cases)</small></span><span>✅ 4/4</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Intensity 0</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Intensity 1</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Intensity 1024</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Intensity 65535</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_PSF.cpp">SP PSF</strong><span>28/28 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Positions / Isotropic Value <small>(7 cases)</small></span><span>✅ 7/7</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Center</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>X Axis</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Y Axis</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Diagonal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>One Degree</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Rotated</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Quarter Turn</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Positions / Circular Anisotropic Matches Isotropic <small>(7 cases)</small></span><span>✅ 7/7</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Center</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>X Axis</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Y Axis</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Diagonal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>One Degree</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Rotated</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Quarter Turn</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Positions / Isotropic Gradient Matches Finite Differences <small>(7 cases)</small></span><span>✅ 7/7</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Center</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>X Axis</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Y Axis</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Diagonal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>One Degree</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Rotated</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Quarter Turn</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Positions / Anisotropic Gradient Matches Finite Differences <small>(7 cases)</small></span><span>✅ 7/7</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Center</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>X Axis</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Y Axis</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Diagonal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>One Degree</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Rotated</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Quarter Turn</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
</section>
</div>
\endhtmlonly
