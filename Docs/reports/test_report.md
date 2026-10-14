# Test report {#test_report}

\htmlonly
<div class="test-report">
<section class="test-run">
<h2>Windows/Release.xml</h2>
<p class="test-summary"><strong>252 tests collected</strong>, 252 passed ✅, 0 failed ❌, 0 skipped ⏭️ in 0:00:22s on 15/10/2026 at 11:46:00</p>
<details class="test-file" open><summary><strong title="test_fit_AstigFit.cpp">Fit Astig Fit</strong><span>14/14 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Calibrate</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Estimate</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Estimate Batch Within Limits</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.011 s</td></tr>
<tr><td>Calibrate Bad</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Estimate Bad</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
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
<details class="test-file" open><summary><strong title="test_fit_Fitting.cpp">Fit Fitting</strong><span>16/16 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Extracts Roi</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Localization Fields <small>(5 cases)</small></span><span>✅ 5/5</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>No Fit</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Gaussian XY</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Gaussian Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Gaussian Sigma XY</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Gaussian Theta</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Recovers Parameters <small>(5 cases)</small></span><span>✅ 5/5</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Fixed Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Free Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Fixed Angle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>One Degree</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Free Angle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Score Matches Squared Correlation <small>(5 cases)</small></span><span>✅ 5/5</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Fixed Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Free Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Fixed Angle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>One Degree</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Free Angle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_loc_Segmentation2D.cpp">LOC Segmentation 2D</strong><span>11/11 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Default State And Reset</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Fields / Region Measurements <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Intensity Ramp</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Constant Plateau</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Thresholds / Strict Threshold <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Below Plateau</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Equal Plateau</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Above Plateau</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Surfaces / Minimum Surface <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Below Minimum</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>At Minimum</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Above Minimum</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Modes / Connected Peaks <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Separated</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Merged</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_tr_TracksFit.cpp">TR Tracks Fit</strong><span>13/13 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Reconfigure And Reuse</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Preserves Missing Lag Positions</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Solves Zero Leading Pivot</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Determination Score</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Recovers Synthetic Curve <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Linear</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Power</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Exponential</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Model Values <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Linear</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Power</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Exponential</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Gradient Matches Finite Differences <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Linear</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Power</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Exponential</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_Palm.cpp">Palm</strong><span>22/22 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Black Stack</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Quadrants</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.070 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Stack <small>(10 cases)</small></span><span>✅ 10/10</span><span>0.459 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>No Fit Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.020 s</td></tr>
<tr><td>No Fit With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.020 s</td></tr>
<tr><td>Gaussian XY Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.033 s</td></tr>
<tr><td>Gaussian XY With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.033 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.046 s</td></tr>
<tr><td>Gaussian Sigma With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.047 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.053 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.054 s</td></tr>
<tr><td>Gaussian Theta Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.075 s</td></tr>
<tr><td>Gaussian Theta With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.078 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Value <small>(10 cases)</small></span><span>✅ 10/10</span><span>0.057 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Plane 1</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Plane 2</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Plane 3</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
<tr><td>Plane 4</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.007 s</td></tr>
<tr><td>Plane 5</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Plane 6</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Plane 7</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Plane 8</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Plane 9</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.004 s</td></tr>
<tr><td>Plane 10</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_TiffManagement.cpp">Tiff Management</strong><span>1/1 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Rejects Wrong Depth And Closes File</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_Alignment.cpp">SP Alignment</strong><span>37/37 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Blank Call With Reference</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Blank Call With Factor</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.004 s</td></tr>
<tr><td>Identity And Reuse</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Integer Translation</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Transpose Pixels</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Fractional Translation</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Coefficients / Coordinate Monomials <small>(20 cases)</small></span><span>✅ 20/20</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Horizontal X Cubed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical X Cubed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal Y Cubed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical Y Cubed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal X Squared Y</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical X Squared Y</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal XY Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical XY Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal X Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical X Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal Y Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical Y Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal XY</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical XY</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal X</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical X</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal Y</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical Y</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal Constant</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical Constant</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Warps / Polynomial Pixel Overlap <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Quadratic</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Test Suites / Saves And Reloads Transformed Stack <small>(9 cases)</small></span><span>✅ 9/9</span><span>9.015 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Collapsed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.043 s</td></tr>
<tr><td>Identity</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.142 s</td></tr>
<tr><td>Translate Quarter</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.221 s</td></tr>
<tr><td>Translate Half</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.978 s</td></tr>
<tr><td>Translate Three Quarters</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>1.362 s</td></tr>
<tr><td>Translate One</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.864 s</td></tr>
<tr><td>Zoom</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.536 s</td></tr>
<tr><td>Transpose</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.493 s</td></tr>
<tr><td>Translate Half Upsampled</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>4.376 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_BSpline.cpp">SP B Spline</strong><span>2/2 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Moves Buffers</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Dirac Impulse</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
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
<tr><td>Odd Rectangle Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Odd Rectangle 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cube Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Cube 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Volumes / Round Trip And Parseval <small>(6 cases)</small></span><span>✅ 6/6</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Single Voxel Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Single Voxel 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Odd Rectangle Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Odd Rectangle 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Cube Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cube 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_Interpolation.cpp">SP Interpolation</strong><span>15/15 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Zoom Factor</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Shift Inside Bounds</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Shift Out Of Bounds</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Interpolate Shifted</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Subpixel Max Interp 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Get Shift 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.007 s</td></tr>
<tr><td>Upsample</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Modes / Non Linear Matrix <small>(8 cases)</small></span><span>✅ 8/8</span><span>0.003 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Trilinear Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Trilinear Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Weight Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Weight Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Spline Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Spline Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
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
<tr><td>Set Members</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.013 s</td></tr>
<tr><td>Extract Beads</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.014 s</td></tr>
<tr><td>Align</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>3.443 s</td></tr>
<tr><td>Get Overlap Minimal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Get Overlap</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.021 s</td></tr>
<tr><td>Soft Normalization Minimal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Soft Normalization</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.030 s</td></tr>
<tr><td>Get PSF Minimal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.051 s</td></tr>
<tr><td>Get PSF</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>4.592 s</td></tr>
<tr><td>Calibrate</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>4.921 s</td></tr>
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
<tr><td>Instant Diffusion</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Analysis Empty State</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Links <small>(20 cases)</small></span><span>✅ 20/20</span><span>0.228 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>No Fit Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.017 s</td></tr>
<tr><td>No Fit Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.011 s</td></tr>
<tr><td>No Fit With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.010 s</td></tr>
<tr><td>No Fit With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.010 s</td></tr>
<tr><td>Gaussian XY Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.010 s</td></tr>
<tr><td>Gaussian XY Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.010 s</td></tr>
<tr><td>Gaussian XY With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.011 s</td></tr>
<tr><td>Gaussian XY With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.008 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.011 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.010 s</td></tr>
<tr><td>Gaussian Sigma With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.012 s</td></tr>
<tr><td>Gaussian Sigma With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.011 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.014 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.012 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.011 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.014 s</td></tr>
<tr><td>Gaussian Theta Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.013 s</td></tr>
<tr><td>Gaussian Theta Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.014 s</td></tr>
<tr><td>Gaussian Theta With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.011 s</td></tr>
<tr><td>Gaussian Theta With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.008 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Reconnection <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.003 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Stationary</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Diffusive</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Linear</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">MSD <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Planar 2D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Spatial 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Track Analysis / Fit <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.002 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Planar 2D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Spatial 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Track Fit / Fit <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.002 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Power</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Exponential</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_loc_ATrous.cpp">LOC A Trous</strong><span>12/12 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Shapes / Reconstructs Input <small>(4 cases)</small></span><span>✅ 4/4</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Single Pixel</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
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
<section class="test-run">
<h2>Windows/Debug.xml</h2>
<p class="test-summary"><strong>252 tests collected</strong>, 252 passed ✅, 0 failed ❌, 0 skipped ⏭️ in 0:12:15s on 15/10/2026 at 11:27:59</p>
<details class="test-file" open><summary><strong title="test_fit_AstigFit.cpp">Fit Astig Fit</strong><span>14/14 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Calibrate</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.023 s</td></tr>
<tr><td>Estimate</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
<tr><td>Estimate Batch Within Limits</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.011 s</td></tr>
<tr><td>Calibrate Bad</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Estimate Bad</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.013 s</td></tr>
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
<details class="test-file" open><summary><strong title="test_fit_Fitting.cpp">Fit Fitting</strong><span>16/16 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Extracts Roi</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Localization Fields <small>(5 cases)</small></span><span>✅ 5/5</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>No Fit</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Gaussian XY</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Gaussian Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Gaussian Sigma XY</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Gaussian Theta</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Recovers Parameters <small>(5 cases)</small></span><span>✅ 5/5</span><span>0.005 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Fixed Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Free Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Fixed Angle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>One Degree</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Free Angle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Score Matches Squared Correlation <small>(5 cases)</small></span><span>✅ 5/5</span><span>0.003 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Fixed Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Free Sigma</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Fixed Angle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>One Degree</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Free Angle</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_loc_Segmentation2D.cpp">LOC Segmentation 2D</strong><span>11/11 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Default State And Reset</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Fields / Region Measurements <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Intensity Ramp</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Constant Plateau</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Thresholds / Strict Threshold <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Below Plateau</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Equal Plateau</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Above Plateau</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Surfaces / Minimum Surface <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Below Minimum</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>At Minimum</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Above Minimum</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Modes / Connected Peaks <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.003 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Separated</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Merged</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_tr_TracksFit.cpp">TR Tracks Fit</strong><span>13/13 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Reconfigure And Reuse</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.038 s</td></tr>
<tr><td>Preserves Missing Lag Positions</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
<tr><td>Solves Zero Leading Pivot</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Determination Score</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Recovers Synthetic Curve <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Linear</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Power</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Exponential</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Model Values <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Linear</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Power</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Exponential</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Models / Gradient Matches Finite Differences <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.000 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Linear</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Power</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Exponential</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_Palm.cpp">Palm</strong><span>22/22 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Black Stack</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.009 s</td></tr>
<tr><td>Quadrants</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>1.284 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Stack <small>(10 cases)</small></span><span>✅ 10/10</span><span>9.881 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>No Fit Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.704 s</td></tr>
<tr><td>No Fit With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.648 s</td></tr>
<tr><td>Gaussian XY Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.876 s</td></tr>
<tr><td>Gaussian XY With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.829 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>1.042 s</td></tr>
<tr><td>Gaussian Sigma With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.918 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>1.117 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>1.071 s</td></tr>
<tr><td>Gaussian Theta Without Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>1.328 s</td></tr>
<tr><td>Gaussian Theta With Watershed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>1.348 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Value <small>(10 cases)</small></span><span>✅ 10/10</span><span>1.109 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Plane 1</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.125 s</td></tr>
<tr><td>Plane 2</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.110 s</td></tr>
<tr><td>Plane 3</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.112 s</td></tr>
<tr><td>Plane 4</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.144 s</td></tr>
<tr><td>Plane 5</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.109 s</td></tr>
<tr><td>Plane 6</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.108 s</td></tr>
<tr><td>Plane 7</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.109 s</td></tr>
<tr><td>Plane 8</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.110 s</td></tr>
<tr><td>Plane 9</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.075 s</td></tr>
<tr><td>Plane 10</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.107 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_TiffManagement.cpp">Tiff Management</strong><span>1/1 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Rejects Wrong Depth And Closes File</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_Alignment.cpp">SP Alignment</strong><span>37/37 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Blank Call With Reference</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Blank Call With Factor</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>1.668 s</td></tr>
<tr><td>Identity And Reuse</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Integer Translation</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Transpose Pixels</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Fractional Translation</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Coefficients / Coordinate Monomials <small>(20 cases)</small></span><span>✅ 20/20</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Horizontal X Cubed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical X Cubed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal Y Cubed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical Y Cubed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal X Squared Y</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical X Squared Y</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Horizontal XY Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical XY Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal X Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical X Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal Y Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical Y Squared</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal XY</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical XY</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal X</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical X</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal Y</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical Y</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Horizontal Constant</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Vertical Constant</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Warps / Polynomial Pixel Overlap <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.002 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Quadratic</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Cubic</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Test Suites / Saves And Reloads Transformed Stack <small>(9 cases)</small></span><span>✅ 9/9</span><span>503.367 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Collapsed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>121.330 s</td></tr>
<tr><td>Identity</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>29.622 s</td></tr>
<tr><td>Translate Quarter</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>61.218 s</td></tr>
<tr><td>Translate Half</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>45.154 s</td></tr>
<tr><td>Translate Three Quarters</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>25.002 s</td></tr>
<tr><td>Translate One</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>21.777 s</td></tr>
<tr><td>Zoom</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>28.736 s</td></tr>
<tr><td>Transpose</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>35.199 s</td></tr>
<tr><td>Translate Half Upsampled</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>135.329 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_BSpline.cpp">SP B Spline</strong><span>2/2 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Moves Buffers</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
<tr><td>Dirac Impulse</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.130 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_FFT.cpp">SP FFT</strong><span>13/13 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Cross Power Spectrum Basic</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Volumes / Matches Direct Discrete Transform <small>(6 cases)</small></span><span>✅ 6/6</span><span>0.003 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Single Voxel Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Single Voxel 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Odd Rectangle Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Odd Rectangle 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cube Separable</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cube 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
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
<tr><td>Zoom Factor</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.035 s</td></tr>
<tr><td>Shift Inside Bounds</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Shift Out Of Bounds</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Interpolate Shifted</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.088 s</td></tr>
<tr><td>Subpixel Max Interp 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.344 s</td></tr>
<tr><td>Get Shift 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.468 s</td></tr>
<tr><td>Upsample</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Modes / Non Linear Matrix <small>(8 cases)</small></span><span>✅ 8/8</span><span>0.007 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Trilinear Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Trilinear Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Weight Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Cubic Weight Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Spline Keep Borders</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.004 s</td></tr>
<tr><td>Spline Crop</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
</tbody></table></details></td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_Matrix.cpp">SP Matrix</strong><span>22/22 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Print</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Send Irregular</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
<tr><td>Has Na N</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Has Invalid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Is Invalid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Get Valid Index</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Count Valid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Max Ignores Na N</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Sum Valid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Mean Valid</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Mean With Mask</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Median</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Quantile</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
<tr><td>Rms Difference</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Normalisation With Shift</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
<tr><td>Crop Volume</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
<tr><td>Extract Sub Volume</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Center Of Mass</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.004 s</td></tr>
<tr><td>Cholesky Success</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
<tr><td>Cholesky Fail</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Inverse Gauss Jordan</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Inverse LU</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.006 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_sp_Spline.cpp">SP Spline</strong><span>10/10 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Set Members</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.094 s</td></tr>
<tr><td>Extract Beads</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.123 s</td></tr>
<tr><td>Align</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>103.243 s</td></tr>
<tr><td>Get Overlap Minimal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.020 s</td></tr>
<tr><td>Get Overlap</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>2.722 s</td></tr>
<tr><td>Soft Normalization Minimal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
<tr><td>Soft Normalization</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.765 s</td></tr>
<tr><td>Get PSF Minimal</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>3.092 s</td></tr>
<tr><td>Get PSF</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>54.519 s</td></tr>
<tr><td>Calibrate</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>51.274 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_tr_Point.cpp">TR Point</strong><span>1/1 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Base</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
</tbody></table></div></details>
<details class="test-file" open><summary><strong title="test_tr_Tracks.cpp">TR Tracks</strong><span>35/35 passed</span></summary>
<div class="test-table-scroll"><table class="test-results">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Tracking With Gaps</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.004 s</td></tr>
<tr><td>Tracking Empty Input</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Blinking Empty State</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>MSD With Gaps</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
<tr><td>Instant Diffusion</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Analysis Empty State</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Links <small>(20 cases)</small></span><span>✅ 20/20</span><span>1.125 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>No Fit Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.055 s</td></tr>
<tr><td>No Fit Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.055 s</td></tr>
<tr><td>No Fit With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.059 s</td></tr>
<tr><td>No Fit With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.060 s</td></tr>
<tr><td>Gaussian XY Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.056 s</td></tr>
<tr><td>Gaussian XY Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.057 s</td></tr>
<tr><td>Gaussian XY With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.059 s</td></tr>
<tr><td>Gaussian XY With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.059 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.057 s</td></tr>
<tr><td>Gaussian Sigma Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.056 s</td></tr>
<tr><td>Gaussian Sigma With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.060 s</td></tr>
<tr><td>Gaussian Sigma With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.061 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.057 s</td></tr>
<tr><td>Gaussian Sigma XY Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.057 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.058 s</td></tr>
<tr><td>Gaussian Sigma XY With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.052 s</td></tr>
<tr><td>Gaussian Theta Without Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.051 s</td></tr>
<tr><td>Gaussian Theta Without Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.050 s</td></tr>
<tr><td>Gaussian Theta With Watershed Ordered</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.054 s</td></tr>
<tr><td>Gaussian Theta With Watershed Reversed</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.052 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Reconnection <small>(3 cases)</small></span><span>✅ 3/3</span><span>0.009 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Stationary</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.005 s</td></tr>
<tr><td>Diffusive</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Linear</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">MSD <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Planar 2D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
<tr><td>Spatial 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Track Analysis / Fit <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.002 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Planar 2D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Spatial 3D</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
</tbody></table></details></td></tr>
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Track Fit / Fit <small>(2 cases)</small></span><span>✅ 2/2</span><span>0.005 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Power</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.002 s</td></tr>
<tr><td>Exponential</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.003 s</td></tr>
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
<tr class="test-group"><td colspan="3"><details><summary><span class="test-group-name">Intensities / Details Vanish And Residual Is Preserved <small>(4 cases)</small></span><span>✅ 4/4</span><span>0.001 s</span></summary>
<table class="test-results test-variants">
<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>
<tbody>
<tr><td>Intensity 0</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Intensity 1</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.000 s</td></tr>
<tr><td>Intensity 1024</td><td><span class="test-status" title="Passed">✅ Passed</span></td><td>0.001 s</td></tr>
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
