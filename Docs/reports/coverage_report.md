# Coverage report {#coverage_report}

\htmlonly
<div class="test-report coverage-report">
<p class="test-summary"><strong>87.62% line coverage</strong> · 2,738/3,125 lines covered · 387 uncovered · 34 files<br>Local Windows / MSVC 19.44 / Debug (/Od, NDEBUG) / OpenCppCoverage 0.9.9 · 15/10/2026 at 11:40:15 +0200</p>
<p>CPU library only (<code>lib/</code>, excluding <code>lib/cuda/</code>). Counts refer to instrumented executable lines, not all source lines. Branch coverage is not measured. No target threshold.</p>
<div class="test-table-scroll"><table class="test-results coverage-results">
<thead><tr><th>File / Uncovered Lines</th><th>Line Coverage</th><th>Covered / Total</th></tr></thead><tbody>
<tr><td><strong>lib/fit/AstigFit.cpp</strong>
<details class="coverage-lines"><summary>4 uncovered lines</summary><pre>
<span><b>  191</b>  		else { lambda *= 10.0; }</span>
<span><b>  393</b>  	return Vec2D{.x = x, .y = fx};</span>
<span><b>  446</b>  			return bestZ;</span>
<span><b>  453</b>  	return bestZ;</span>
</pre></details>
</td><td><span>98.51%</span><progress max="100" value="98.5130" aria-label="lib/fit/AstigFit.cpp line coverage">98.51%</progress></td><td>265 / 269</td></tr>
<tr><td><strong>lib/fit/AstigFit.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/fit/AstigFit.hpp line coverage">100.00%</progress></td><td>1 / 1</td></tr>
<tr><td><strong>lib/fit/Fitting.cpp</strong>
<details class="coverage-lines"><summary>9 uncovered lines</summary><pre>
<span><b>   25</b>  		case fit_mode::spline_3d_scmos: nParams_ = 5; // offset, amplitude, x, y, z</span>
<span><b>   79</b>  		case fit_mode::spline_3d_scmos: point[7] = intensity;</span>
<span><b>  154</b>  			const double RoiC = (fitParams_[0] - 1) / 2.0,</span>
<span><b>  155</b>  			             RoiD = fitParams_[0] / 4.0;</span>
<span><b>  156</b>  			const double Dmin = RoiC - RoiD,</span>
<span><b>  157</b>  			             Dmax = RoiC + RoiD;</span>
<span><b>  159</b>  			p[0] = {-0.01, -1, Dmin, Dmin, 0};</span>
<span><b>  160</b>  			p[1] = {minI, maxI - minI, x, y, z};</span>
<span><b>  161</b>  			p[2] = {0.01, 1, Dmax, Dmax, 1}; // Le max pour Z est borné par la taille des PSF spline (voir paramètres).</span>
</pre></details>
</td><td><span>90.00%</span><progress max="100" value="90.0000" aria-label="lib/fit/Fitting.cpp line coverage">90.00%</progress></td><td>81 / 90</td></tr>
<tr><td><strong>lib/fit/Fitting.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/fit/Fitting.hpp line coverage">100.00%</progress></td><td>12 / 12</td></tr>
<tr><td><strong>lib/fit/GaussianFit.cpp</strong>
<details class="coverage-lines"><summary>8 uncovered lines</summary><pre>
<span><b>   33</b>  		fill(infos, point, infos.x_, infos.y_, infos.z_);</span>
<span><b>   34</b>  		return;</span>
<span><b>   77</b>  			if (lambda &gt;= LAMBDA_MAX) { break; }               // On est déjà au max, on arrête</span>
<span><b>   78</b>  			lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de λ.</span>
<span><b>   79</b>  			continue;                                          // On repart au début de la boucle</span>
<span><b>  141</b>  		case fit_mode::gaussian_2d_fixed: break;</span>
<span><b>  263</b>  		case fit_mode::spline_3d_scmos: return 0;</span>
<span><b>  265</b>  	return 0;</span>
</pre></details>
</td><td><span>95.09%</span><progress max="100" value="95.0920" aria-label="lib/fit/GaussianFit.cpp line coverage">95.09%</progress></td><td>155 / 163</td></tr>
<tr><td><strong>lib/fit/GaussianFit.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/fit/GaussianFit.hpp line coverage">100.00%</progress></td><td>20 / 20</td></tr>
<tr><td><strong>lib/fit/SplineFit.cpp</strong>
<details class="coverage-lines"><summary>200 uncovered lines</summary><pre>
<span><b>   28</b>  {</span>
<span><b>   29</b>  	const std::complex Z(value, 0.0);</span>
<span><b>   30</b>  	return std::real(std::sqrt(Z));</span>
<span><b>   31</b>  }</span>
<span><b>   35</b>  {</span>
<span><b>   36</b>  	size_t idx = 0;</span>
<span><b>   37</b>  	roi_       = static_cast&lt;size_t&gt;(fitParams_[idx++]);</span>
<span><b>   38</b>  	for (size_t i = 0; i &lt; 3; ++i) { splineDim_[i] = static_cast&lt;size_t&gt;(fitParams_[idx++]); }</span>
<span><b>   39</b>  	splineSizeXY_ = splineDim_.GetPlaneSize();</span>
<span><b>   40</b>  	splineSize_   = splineDim_.GetSize();</span>
<span><b>   41</b>  	initZ_        = static_cast&lt;double&gt;(splineDim_.z) / 2.0;</span>
<span><b>   42</b>  	deltaZ_       = fitParams_[idx++];</span>
<span><b>   43</b>  	coeffs_.resize(splineSize_ * SPLINE_PATCH_SIZE);</span>
<span><b>   44</b>  	for (double&amp; coeff : coeffs_) { coeff = fitParams_[idx++]; }</span>
<span><b>   45</b>  	const double CenterCoeff = coeffs_[32 * splineSize_ + splineDim_.z / 2 * splineSizeXY_ + splineDim_.y / 2 * splineDim_.x + splineDim_.x / 2];</span>
<span><b>   46</b>  	intensityRatio_          = std::abs(CenterCoeff) &lt; std::numeric_limits&lt;double&gt;::epsilon() ? 1.0 : 4.0 / CenterCoeff;</span>
<span><b>   47</b>  	splineOffset_            = static_cast&lt;size_t&gt;(std::floor((static_cast&lt;double&gt;(splineDim_.x) + 1.0 - static_cast&lt;double&gt;(roi_)) / 2));</span>
<span><b>   48</b>  }</span>
<span><b>   52</b>  {</span>
<span><b>   54</b>  	const auto Ux = static_cast&lt;size_t&gt;(infos.x_), Uy = static_cast&lt;size_t&gt;(infos.y_);</span>
<span><b>   56</b>  	if (roiHalf_ &gt; Ux || Ux &gt;= dim_.x - roiHalf_ || roiHalf_ &gt; Uy || Uy &gt;= dim_.y - roiHalf_) {</span>
<span><b>   57</b>  		fill(infos, point, infos.x_, infos.y_, infos.z_);</span>
<span><b>   58</b>  		return;</span>
<span><b>   61</b>  	double x = infos.x_, y = infos.y_, z = initZ_;</span>
<span><b>   64</b>  	double lambda = LAMBDA_0, minI = 0, maxI = 0, model = 0, err = std::numeric_limits&lt;double&gt;::max(), mu = 1 + lambda;</span>
<span><b>   65</b>  	std::vector&lt;double&gt; hessian(nParams2_), l(nParams2_), u(nParams2_),</span>
<span><b>   66</b>  	                    dudt(nParams_), jacobian(nParams_), jump(nParams_), update(nParams_),</span>
<span><b>   67</b>  	                    f(SPLINE_PATCH_SIZE);</span>
<span><b>   68</b>  	std::array df{std::vector&lt;double&gt;(SPLINE_PATCH_SIZE), std::vector&lt;double&gt;(SPLINE_PATCH_SIZE), std::vector&lt;double&gt;(SPLINE_PATCH_SIZE)};</span>
<span><b>   69</b>  	std::vector&lt;Pixel&gt; grid(n_);</span>
<span><b>   70</b>  	std::array&lt;std::vector&lt;double&gt;, 3&gt; p;</span>
<span><b>   73</b>  	preCompute(Ux, Uy, grid, minI, maxI);</span>
<span><b>   77</b>  	matrix::CenterOfMass(grid, roi_, x, y);</span>
<span><b>   79</b>  	initP(p, x, y, z, minI, maxI);</span>
<span><b>   80</b>  	p[1][1] *= intensityRatio_;                                                                // Change l&#x27;intensité pour la position actuelle.</span>
<span><b>   81</b>  	p[2][4] = static_cast&lt;double&gt;(splineDim_.z);                                               // Change le maximum pour Z.</span>
<span><b>   82</b>  	for (size_t i = 0; i &lt; nParams_; ++i) { p[1][i] = std::clamp(p[1][i], p[0][i], p[2][i]); } // Clamp en cas de problème.</span>
<span><b>   83</b>  	jump = {std::max(20.0, p[1][1]), std::max(100.0, p[1][1]), 1.0, 1.0, std::max(static_cast&lt;double&gt;(splineDim_.z) / 3.0, 2.0)};</span>
<span><b>   85</b>  	bool stop        = false;</span>
<span><b>   86</b>  	size_t iteration = 0;</span>
<span><b>   87</b>  	while (!stop) { // Main Loop</span>
<span><b>   88</b>  		iteration++;</span>
<span><b>   90</b>  		const double OldErr                 = err;</span>
<span><b>   91</b>  		const std::vector&lt;double&gt; OldP      = p[1];</span>
<span><b>   92</b>  		const std::vector&lt;double&gt; OldUpdate = update;</span>
<span><b>   93</b>  		const double OldLambda              = lambda;</span>
<span><b>   95</b>  		std::ranges::fill(jacobian, 0.0);</span>
<span><b>   96</b>  		std::ranges::fill(hessian, 0.0);</span>
<span><b>   97</b>  		std::ranges::fill(update, 0.0);</span>
<span><b>   98</b>  		std::ranges::fill(l, 0.0);</span>
<span><b>   99</b>  		std::ranges::fill(u, 0.0);</span>
<span><b>  103</b>  		const double Xr = -1 * (p[1][2] - static_cast&lt;double&gt;(roiHalf_) + 0.5),</span>
<span><b>  104</b>  		             Yr = -1 * (p[1][3] - static_cast&lt;double&gt;(roiHalf_) + 0.5);</span>
<span><b>  106</b>  		const int Xi = static_cast&lt;int&gt;(std::floor(Xr)),</span>
<span><b>  107</b>  		          Yi = static_cast&lt;int&gt;(std::floor(Yr)),</span>
<span><b>  108</b>  		          Zi = static_cast&lt;int&gt;(std::floor(p[1][4]));</span>
<span><b>  110</b>  		const double Xc = Xr - Xi, Yc = Yr - Yi, Zc = p[1][4] - Zi;</span>
<span><b>  111</b>  		delta3D(Xc, Yc, Zc, f, df);</span>
<span><b>  115</b>  		Vec3I pos = {.x = 0, .y = 0, .z = static_cast&lt;size_t&gt;(Zi)};</span>
<span><b>  116</b>  		for (size_t i = 0, idx = 0; i &lt; roi_; ++i) {</span>
<span><b>  117</b>  			pos.x = Xi + i + splineOffset_;</span>
<span><b>  118</b>  			for (size_t j = 0; j &lt; roi_; ++j, ++idx) {</span>
<span><b>  119</b>  				pos.y = Yi + j + splineOffset_;</span>
<span><b>  120</b>  				derivativeSpline(pos, f, df, p[1], dudt, model);</span>
<span><b>  121</b>  				updateJacobianHessian(grid[idx].v, model, err, jacobian, hessian, dudt);</span>
<span><b>  122</b>  			}</span>
<span><b>  123</b>  		}</span>
<span><b>  126</b>  		for (size_t i = 0; i &lt; nParams_; ++i) { hessian[i * nParams_ + i] = hessian[i * nParams_ + i] * mu; }</span>
<span><b>  129</b>  		if (matrix::Cholesky(hessian, nParams_, l, u)) { // Si la Factorisation de Cholesky a réussi</span>
<span><b>  131</b>  			luEvaluate(l, u, jacobian, update);</span>
<span><b>  133</b>  			for (size_t i = 0; i &lt; nParams_; ++i) {</span>
<span><b>  134</b>  				if (std::abs(OldUpdate[i]) &gt; EPSILON &amp;&amp; update[i] / OldUpdate[i] &lt; -0.5) { jump[i] *= 0.5; }</span>
<span><b>  135</b>  				const double Scale = update[i] / (1.0 + std::abs(update[i]) / jump[i]); // Pas de division par une grandeur instable</span>
<span><b>  136</b>  				p[1][i] -= Scale;</span>
<span><b>  137</b>  				p[1][i]   = std::clamp(p[1][i], p[0][i], p[2][i]);</span>
<span><b>  138</b>  				update[i] = Scale;</span>
<span><b>  139</b>  			}</span>
<span><b>  140</b>  			if (err &gt; OldErr) { // si erreur plus importante, restauration des valeurs et augmentation du lambda</span>
<span><b>  141</b>  				err    = OldErr;</span>
<span><b>  142</b>  				p[1]   = OldP;</span>
<span><b>  143</b>  				update = OldUpdate;</span>
<span><b>  144</b>  				lambda = OldLambda;</span>
<span><b>  145</b>  				updateMuLambda(mu, lambda, true);</span>
<span><b>  146</b>  			}</span>
<span><b>  147</b>  			else { updateMuLambda(mu, lambda, false); } // sinon réduction</span>
<span><b>  148</b>  		}</span>
<span><b>  150</b>  			err    = OldErr;</span>
<span><b>  151</b>  			p[1]   = OldP;</span>
<span><b>  152</b>  			update = OldUpdate;</span>
<span><b>  153</b>  			lambda = OldLambda;</span>
<span><b>  154</b>  			updateMuLambda(mu, lambda, true);</span>
<span><b>  158</b>  		if (iteration &gt;= N_ITER) { stop = true; } // Si Nombre d&#x27;itérations maximum atteint, on arrête.</span>
<span><b>  159</b>  		if (iteration &gt; 2 &amp;&amp; (err &lt; EPSILON_1 || std::abs(err - OldErr) / err &lt; EPSILON_1)) { stop = true; } // Si Convergence, on arrête.</span>
<span><b>  160</b>  	}</span>
<span><b>  163</b>  	std::vector&lt;double&gt; m(nParams2_, 0.0), mInv(nParams2_, 0.0), crlb(nParams_);</span>
<span><b>  164</b>  	std::ranges::fill(dudt, 0.0);</span>
<span><b>  166</b>  	const double Xr = -1 * (p[1][2] - static_cast&lt;double&gt;(roiHalf_) + 0.5),</span>
<span><b>  167</b>  	             Yr = -1 * (p[1][3] - static_cast&lt;double&gt;(roiHalf_) + 0.5);</span>
<span><b>  169</b>  	const int Xi = static_cast&lt;int&gt;(std::floor(Xr)),</span>
<span><b>  170</b>  	          Yi = static_cast&lt;int&gt;(std::floor(Yr)),</span>
<span><b>  171</b>  	          Zi = static_cast&lt;int&gt;(std::floor(p[1][4]));</span>
<span><b>  173</b>  	const double Xc = Xr - Xi, Yc = Yr - Yi, Zc = p[1][4] - Zi;</span>
<span><b>  174</b>  	delta3D(Xc, Yc, Zc, f, df);</span>
<span><b>  178</b>  	Vec3I pos = {.x = 0, .y = 0, .z = static_cast&lt;size_t&gt;(Zi)};</span>
<span><b>  179</b>  	for (size_t i = 0, idx = 0; i &lt; roi_; ++i) {</span>
<span><b>  180</b>  		pos.x = Xi + i + splineOffset_;</span>
<span><b>  181</b>  		for (size_t j = 0; j &lt; roi_; ++j, ++idx) {</span>
<span><b>  182</b>  			pos.y = Yi + j + splineOffset_;</span>
<span><b>  183</b>  			derivativeSpline(pos, f, df, p[1], dudt, model);</span>
<span><b>  185</b>  			for (size_t i2 = 0; i2 &lt; nParams_; ++i2) {</span>
<span><b>  186</b>  				for (size_t j2 = i2; j2 &lt; nParams_; ++j2) {</span>
<span><b>  187</b>  					m[i2 * nParams_ + j2] += dudt[j2] * dudt[i2] / model;</span>
<span><b>  188</b>  					m[j2 * nParams_ + i2] = m[i2 * nParams_ + j2];</span>
<span><b>  189</b>  				}</span>
<span><b>  190</b>  			}</span>
<span><b>  191</b>  		}</span>
<span><b>  192</b>  	}</span>
<span><b>  196</b>  	for (size_t i = 0; i &lt; nParams_; ++i) { if (std::abs(m[i * nParams_ + i]) &lt; std::numeric_limits&lt;double&gt;::epsilon()) { m[i * nParams_ + i] = EPSILON_2; } }</span>
<span><b>  199</b>  	matrix::InverseLU(m, mInv, nParams_);</span>
<span><b>  201</b>  	for (size_t i = 0; i &lt; nParams_; ++i) { crlb[i] = SqrtReal(mInv[i * nParams_ + i]); }</span>
<span><b>  204</b>  	const double Offset = crlb[0], Intensity             = crlb[1],</span>
<span><b>  205</b>  	             RSS_XY = (crlb[2] + crlb[3]) / 2, RSS_Z = deltaZ_ * crlb[4];</span>
<span><b>  207</b>  	x = p[1][2] - static_cast&lt;double&gt;(roiHalf_) + infos.x_;</span>
<span><b>  208</b>  	y = p[1][3] - static_cast&lt;double&gt;(roiHalf_) + infos.y_;</span>
<span><b>  209</b>  	z = deltaZ_ * (p[1][4] - initZ_); //p[1][4]; //deltaZ_ * (p[1][4] - initZ_);</span>
<span><b>  211</b>  	fill(infos, point, x, y, z, Offset, Intensity, 0, 0, 0, RSS_XY, RSS_Z);</span>
<span><b>  212</b>  }</span>
<span><b>  217</b>  {</span>
<span><b>  218</b>  	pos = pos.Clamp(splineDim_ - 1);</span>
<span><b>  220</b>  	const size_t Pos = pos.GetIndex(splineDim_);</span>
<span><b>  221</b>  	std::ranges::fill(dudt, 0.0); // Initialise toutes les dérivées à zéro</span>
<span><b>  224</b>  	for (size_t i = 0; i &lt; SPLINE_PATCH_SIZE; ++i) {</span>
<span><b>  225</b>  		const size_t Id = i * splineSize_ + Pos;</span>
<span><b>  226</b>  		dudt[1] += f[i] * coeffs_[Id];     // Contribution à l&#x27;intensité.</span>
<span><b>  227</b>  		dudt[2] += df[0][i] * coeffs_[Id]; // Dérivée en x.</span>
<span><b>  228</b>  		dudt[3] += df[1][i] * coeffs_[Id]; // Dérivée en y.</span>
<span><b>  229</b>  		dudt[4] += df[2][i] * coeffs_[Id]; // Dérivée en z.</span>
<span><b>  230</b>  	}</span>
<span><b>  231</b>  	dudt[0] = 1;                   // Dérivée par rapport à l&#x27;offset</span>
<span><b>  232</b>  	dudt[2] *= -1 * p[1];          // Ajustement en x.</span>
<span><b>  233</b>  	dudt[3] *= -1 * p[1];          // Ajustement en y.</span>
<span><b>  234</b>  	dudt[4] *= p[1];               // Ajustement en z.</span>
<span><b>  235</b>  	model = p[0] + p[1] * dudt[1]; // Évaluation du modèle spline : offset + amplitude * contribution</span>
<span><b>  236</b>  }</span>
<span><b>  241</b>  {</span>
<span><b>  243</b>  	if (data &gt; 0) { err = err + 2 * (model - data - data * log(model / data)); }</span>
<span><b>  244</b>  	else { err = err + 2 * model; }</span>
<span><b>  246</b>  	const double T1 = 1 - data / model;</span>
<span><b>  247</b>  	for (size_t i = 0; i &lt; nParams_; ++i) { jacobian[i] += T1 * dudt[i]; } // ∂L/∂θ_i</span>
<span><b>  249</b>  	const double T2 = data / (model * model);</span>
<span><b>  250</b>  	for (size_t i = 0; i &lt; nParams_; ++i) {</span>
<span><b>  251</b>  		for (size_t j = i; j &lt; nParams_; ++j) {</span>
<span><b>  252</b>  			hessian[i * nParams_ + j] += T2 * dudt[i] * dudt[j];</span>
<span><b>  253</b>  			hessian[j * nParams_ + i] = hessian[i * nParams_ + j]; // Hessienne symétrique</span>
<span><b>  254</b>  		}</span>
<span><b>  255</b>  	}</span>
<span><b>  256</b>  }</span>
<span><b>  260</b>  {</span>
<span><b>  261</b>  	if (up) { // Augmentation du lambda : on recule dans la descente</span>
<span><b>  262</b>  		mu     = std::max((1 + lambda * LAMBDA_UP) / (1 + lambda), 1.3);</span>
<span><b>  263</b>  		lambda = LAMBDA_UP * lambda;</span>
<span><b>  264</b>  	}</span>
<span><b>  266</b>  		lambda = LAMBDA_DN * lambda;</span>
<span><b>  267</b>  		mu     = 1 + lambda;</span>
<span><b>  269</b>  }</span>
<span><b>  273</b>  {</span>
<span><b>  275</b>  	std::vector y(nParams_, 0.0);</span>
<span><b>  278</b>  	for (size_t i = 0; i &lt; nParams_; ++i) {</span>
<span><b>  279</b>  		y[i] = b[i];</span>
<span><b>  280</b>  		for (size_t j = 0; j &lt; i; ++j) { y[i] -= l[j * nParams_ + i] * y[j]; }</span>
<span><b>  281</b>  		y[i] /= l[i * nParams_ + i];</span>
<span><b>  282</b>  	}</span>
<span><b>  285</b>  	for (int i = static_cast&lt;int&gt;(nParams_ - 1); i &gt;= 0; --i) {</span>
<span><b>  286</b>  		x[i] = y[i];</span>
<span><b>  287</b>  		for (size_t j = i + 1; j &lt; nParams_; j++) { x[i] -= u[j * nParams_ + i] * x[j]; }</span>
<span><b>  288</b>  		x[i] /= u[i * nParams_ + i];</span>
<span><b>  289</b>  	}</span>
<span><b>  290</b>  }</span>
<span><b>  293</b>  {</span>
<span><b>  294</b>  	std::ranges::fill(f, 0.0);</span>
<span><b>  295</b>  	for (size_t i = 0; i &lt; 3; ++i) { std::ranges::fill(df[i], 0.0); }</span>
<span><b>  296</b>  	double cz = 1.0;</span>
<span><b>  297</b>  	for (size_t i = 0; i &lt; SPLINE_PATCH_SIDE; ++i) {</span>
<span><b>  298</b>  		double cy = 1.0;</span>
<span><b>  299</b>  		for (size_t j = 0; j &lt; SPLINE_PATCH_SIDE; ++j) {</span>
<span><b>  300</b>  			double cx = 1.0;</span>
<span><b>  301</b>  			for (size_t k = 0; k &lt; SPLINE_PATCH_SIDE; ++k) {</span>
<span><b>  302</b>  				const size_t Id   = i * SPLINE_PATCH_PLANE + j * SPLINE_PATCH_SIDE + k;</span>
<span><b>  303</b>  				const double Coef = cz * cy * cx;</span>
<span><b>  304</b>  				f[Id]             = Coef;</span>
<span><b>  305</b>  				if (k + 1 &lt; SPLINE_PATCH_SIDE) { df[0][Id + 1] = static_cast&lt;double&gt;(k + 1) * Coef; }</span>
<span><b>  306</b>  				if (j + 1 &lt; SPLINE_PATCH_SIDE) { df[1][Id + SPLINE_PATCH_SIDE] = static_cast&lt;double&gt;(j + 1) * Coef; }</span>
<span><b>  307</b>  				if (i + 1 &lt; SPLINE_PATCH_SIDE) { df[2][Id + SPLINE_PATCH_PLANE] = static_cast&lt;double&gt;(i + 1) * Coef; }</span>
<span><b>  308</b>  				cx *= x;</span>
<span><b>  309</b>  			}</span>
<span><b>  310</b>  			cy *= y;</span>
<span><b>  311</b>  		}</span>
<span><b>  312</b>  		cz *= z;</span>
<span><b>  313</b>  	}</span>
<span><b>  314</b>  }</span>
</pre></details>
</td><td><span>0.00%</span><progress max="100" value="0.0000" aria-label="lib/fit/SplineFit.cpp line coverage">0.00%</progress></td><td>0 / 200</td></tr>
<tr><td><strong>lib/fit/SplineFit.hpp</strong>
<details class="coverage-lines"><summary>9 uncovered lines</summary><pre>
<span><b>   22</b>  	explicit SplineFit(const fit_mode mode, const uint16_t* image, const Vec2I&amp; dim, const size_t maxPoints,</span>
<span><b>   23</b>  	                   const double* fitParams) : Fitting(mode, image, dim, maxPoints, fitParams) { initialize(); }</span>
<span><b>   34</b>  	Vec3I splineDim_{};          ///&lt; Dimensions de la spline (X, Y, Z)</span>
<span><b>   35</b>  	size_t splineSizeXY_{},      ///&lt; Taille totale d&#x27;un plan XY. Utilisé pour accéder rapidement aux coefficients.</span>
<span><b>   36</b>  	       splineSize_{},        ///&lt; Taille totale de l&#x27;ensemble du volume spline.</span>
<span><b>   37</b>  	       splineOffset_{};      ///&lt; Offset utilisé pour indexer ou accéder efficacement aux sous-blocs ou aux patchs de la spline.</span>
<span><b>   38</b>  	double initZ_{},             ///&lt; Position initiale en Z (origine physique ou de modélisation de la spline en profondeur).</span>
<span><b>   39</b>  	       deltaZ_{},            ///&lt; Pas ou résolution entre deux couches Z successives dans la spline.</span>
<span><b>   40</b>  	       intensityRatio_{};    ///&lt; Rapport d&#x27;intensité entre l&#x27;image et le modèle spline (pour corriger une mise à l&#x27;échelle ou un gain).</span>
</pre></details>
</td><td><span>0.00%</span><progress max="100" value="0.0000" aria-label="lib/fit/SplineFit.hpp line coverage">0.00%</progress></td><td>0 / 9</td></tr>
<tr><td><strong>lib/localization/ATrous.cpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/localization/ATrous.cpp line coverage">100.00%</progress></td><td>34 / 34</td></tr>
<tr><td><strong>lib/localization/ATrous.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/localization/ATrous.hpp line coverage">100.00%</progress></td><td>3 / 3</td></tr>
<tr><td><strong>lib/localization/Segmentation2D.cpp</strong>
<details class="coverage-lines"><summary>9 uncovered lines</summary><pre>
<span><b>   80</b>  		for (size_t i = 0; i &lt; dim_.y; ++i) {</span>
<span><b>   81</b>  			const auto Idx = i * dim_.x;</span>
<span><b>   82</b>  			if (map_[Idx] != 0) { bounds[map_[Idx]] = true; }</span>
<span><b>   83</b>  			if (map_[Idx + dim_.x - 1] != 0) { bounds[map_[Idx + dim_.x - 1]] = true; }</span>
<span><b>   84</b>  		}</span>
<span><b>   85</b>  		for (size_t i = 1; i &lt; dim_.x - 1; ++i) {</span>
<span><b>   86</b>  			if (map_[i] != 0) { bounds[map_[i]] = true; }</span>
<span><b>   87</b>  			if (map_[(dim_.y - 1) * dim_.x + i] != 0) { bounds[map_[(dim_.y - 1) * dim_.x + i]] = true; }</span>
<span><b>   88</b>  		}</span>
</pre></details>
</td><td><span>96.05%</span><progress max="100" value="96.0526" aria-label="lib/localization/Segmentation2D.cpp line coverage">96.05%</progress></td><td>219 / 228</td></tr>
<tr><td><strong>lib/localization/Segmentation2D.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/localization/Segmentation2D.hpp line coverage">100.00%</progress></td><td>2 / 2</td></tr>
<tr><td><strong>lib/localization/SegmentationInfo.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/localization/SegmentationInfo.hpp line coverage">100.00%</progress></td><td>17 / 17</td></tr>
<tr><td><strong>lib/Palm.cpp</strong>
<details class="coverage-lines"><summary>9 uncovered lines</summary><pre>
<span><b>  121</b>  		}</span>
<span><b>  127</b>  		}</span>
<span><b>  133</b>  		}</span>
<span><b>  139</b>  		}</span>
<span><b>  142</b>  			const fit::SplineFit Sfit(fit::fit_mode::spline_3d_emccd, image, Dim, nPointsMaxPerPlane_, fitParams_);</span>
<span><b>  143</b>  			Sfit.Fit(Segmentation, points);</span>
<span><b>  144</b>  			break;</span>
<span><b>  145</b>  		}</span>
<span><b>  151</b>  		}</span>
</pre></details>
</td><td><span>90.32%</span><progress max="100" value="90.3226" aria-label="lib/Palm.cpp line coverage">90.32%</progress></td><td>84 / 93</td></tr>
<tr><td><strong>lib/signal_processing/Alignment.cpp</strong>
<details class="coverage-lines"><summary>4 uncovered lines</summary><pre>
<span><b>   25</b>  	align();</span>
<span><b>   26</b>  }</span>
<span><b>   51</b>  }</span>
<span><b>  166</b>  			default: return false;</span>
</pre></details>
</td><td><span>97.26%</span><progress max="100" value="97.2603" aria-label="lib/signal_processing/Alignment.cpp line coverage">97.26%</progress></td><td>142 / 146</td></tr>
<tr><td><strong>lib/signal_processing/Alignment.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/signal_processing/Alignment.hpp line coverage">100.00%</progress></td><td>13 / 13</td></tr>
<tr><td><strong>lib/signal_processing/BSpline3D.cpp</strong>
<details class="coverage-lines"><summary>97 uncovered lines</summary><pre>
<span><b>   53</b>  			if (x &gt;= -0.5 &amp;&amp; x &lt; 0.5) { return 1.0; }</span>
<span><b>   54</b>  			break;</span>
<span><b>   58</b>  			if (x &gt;= -1 &amp;&amp; x &lt; 0) { return 1.0 + x; }</span>
<span><b>   59</b>  			if (x &gt;= 0 &amp;&amp; x &lt; 1) { return 1.0 - x; }</span>
<span><b>   60</b>  			break;</span>
<span><b>   63</b>  			constexpr auto Edge   = 1.125, // = 9 / 8 = 3^2 / 2^3</span>
<span><b>   64</b>  			               Center = 0.75,  // = 3 / 4 = 3   / 2^2</span>
<span><b>   65</b>  			               A      = 1.5,</span>
<span><b>   66</b>  			               B      = 0.5;</span>
<span><b>   68</b>  			const double X2 = x * x;</span>
<span><b>   69</b>  			if (x &gt;= -A &amp;&amp; x &lt; -B) { return Edge + A * x + B * X2; }</span>
<span><b>   70</b>  			if (x &gt;= -B &amp;&amp; x &lt; B) { return Center - x * x; }</span>
<span><b>   71</b>  			if (x &gt;= B &amp;&amp; x &lt; A) { return Edge - A * x + B * X2; }</span>
<span><b>   72</b>  			break;</span>
<span><b>   87</b>  			break;</span>
<span><b>   90</b>  			const double HalfDegP1 = (static_cast&lt;double&gt;(deg) + 1.0) / 2.0;</span>
<span><b>   91</b>  			return ((x + HalfDegP1) * Eval(x + 0.5, deg - 1) + (HalfDegP1 - x) * Eval(x - 0.5, deg - 1)) / static_cast&lt;double&gt;(deg);</span>
<span><b>   95</b>  	return 0;</span>
<span><b>  177</b>  				if (lambda_[i] &lt;= std::numeric_limits&lt;double&gt;::epsilon()) { c0_[i] = 1.0; }</span>
<span><b>  178</b>  				else { throw std::invalid_argument(&quot;Even degree 0 with nonzero lambda_[i] not supported&quot;); }</span>
<span><b>  179</b>  				break;</span>
<span><b>  182</b>  				if (lambda_[i] &lt;= std::numeric_limits&lt;double&gt;::epsilon()) { c0_[i] = 1.0; }</span>
<span><b>  184</b>  					const double R = 1.0 + 1.0 / (2.0 * lambda_[i]) - std::sqrt(1.0 + 4.0 * lambda_[i]) / (2.0 * lambda_[i]);</span>
<span><b>  185</b>  					filters_[i].emplace_back(R);</span>
<span><b>  186</b>  					c0_[i] = -1.0 / lambda_[i];</span>
<span><b>  188</b>  				break;</span>
<span><b>  191</b>  				if (lambda_[i] &lt;= std::numeric_limits&lt;double&gt;::epsilon()) {</span>
<span><b>  192</b>  					filters_[i].emplace_back(2.0 * std::numbers::sqrt2 - 3.0);</span>
<span><b>  193</b>  					c0_[i] = 8.0;</span>
<span><b>  194</b>  				}</span>
<span><b>  195</b>  				else { throw std::invalid_argument(&quot;Even degree 2 with nonzero lambda_[i] not supported&quot;); }</span>
<span><b>  196</b>  				break;</span>
<span><b>  217</b>  				if (lambda_[i] &lt;= std::numeric_limits&lt;double&gt;::epsilon()) {</span>
<span><b>  218</b>  					filters_[i] = {{-0.36134122590022}, {-0.0137254292973391}};</span>
<span><b>  219</b>  					c0_[i]      = 384.0;</span>
<span><b>  220</b>  				}</span>
<span><b>  221</b>  				else { throw std::invalid_argument(&quot;Even degree 4 with nonzero lambda_[i] not supported&quot;); }</span>
<span><b>  222</b>  				break;</span>
<span><b>  225</b>  				if (lambda_[i] &lt;= std::numeric_limits&lt;double&gt;::epsilon()) {</span>
<span><b>  226</b>  					filters_[i] = {{-0.430575347099973}, {-0.0430962882032647}};</span>
<span><b>  227</b>  					c0_[i]      = 120.0;</span>
<span><b>  228</b>  				}</span>
<span><b>  230</b>  					std::vector&lt;double&gt; p{1, -6, 15, -20, 15, -6, 1};</span>
<span><b>  231</b>  					p[1] -= 1.0 / (120.0 * lambda_[i]);</span>
<span><b>  232</b>  					p[2] -= 13.0 / (60.0 * lambda_[i]);</span>
<span><b>  233</b>  					p[3] -= 11.0 / (20.0 * lambda_[i]);</span>
<span><b>  234</b>  					p[4] -= 13.0 / (60.0 * lambda_[i]);</span>
<span><b>  235</b>  					p[5] -= 1.0 / (120.0 * lambda_[i]);</span>
<span><b>  237</b>  					auto roots = solvePolynomial(p);</span>
<span><b>  238</b>  					std::ranges::copy_if(roots, std::back_inserter(filters_[i]), [](auto z) { return std::abs(z) &lt;= 1.0; });</span>
<span><b>  240</b>  					c0_[i] = -1.0 / lambda_[i];</span>
<span><b>  241</b>  				}</span>
<span><b>  242</b>  				break;</span>
<span><b>  245</b>  				if (lambda_[i] &lt;= std::numeric_limits&lt;double&gt;::epsilon()) {</span>
<span><b>  246</b>  					filters_[i] = {{-0.488294589303046}, {-0.0816792710762375}, {-0.00141415180832582}};</span>
<span><b>  247</b>  					c0_[i]      = 46080.0;</span>
<span><b>  248</b>  				}</span>
<span><b>  249</b>  				else { throw std::invalid_argument(&quot;Even degree 6 with nonzero lambda_[i] not supported&quot;); }</span>
<span><b>  250</b>  				break;</span>
<span><b>  253</b>  				if (lambda_[i] &lt;= std::numeric_limits&lt;double&gt;::epsilon()) {</span>
<span><b>  254</b>  					filters_[i] = {{-0.535280430796439}, {-0.122554615192327}, {-0.00914869480960828}};</span>
<span><b>  255</b>  					c0_[i]      = 5040.0;</span>
<span><b>  256</b>  				}</span>
<span><b>  258</b>  					std::vector&lt;double&gt; p{1, -8, 28, -56, 70, -56, 28, -8, 1};</span>
<span><b>  259</b>  					p[1] += 1.0 / (5040.0 * lambda_[i]);</span>
<span><b>  260</b>  					p[2] += 1.0 / (42.0 * lambda_[i]);</span>
<span><b>  261</b>  					p[3] += 397.0 / (1680.0 * lambda_[i]);</span>
<span><b>  262</b>  					p[4] += 151.0 / (315.0 * lambda_[i]);</span>
<span><b>  263</b>  					p[5] += 397.0 / (1680.0 * lambda_[i]);</span>
<span><b>  264</b>  					p[6] += 1.0 / (42.0 * lambda_[i]);</span>
<span><b>  265</b>  					p[7] += 1.0 / (5040.0 * lambda_[i]);</span>
<span><b>  267</b>  					auto roots = solvePolynomial(p);</span>
<span><b>  268</b>  					std::ranges::copy_if(roots, std::back_inserter(filters_[i]), [](auto z) { return std::abs(z) &lt;= 1.0; });</span>
<span><b>  270</b>  					c0_[i] = 1.0 / lambda_[i];</span>
<span><b>  271</b>  				}</span>
<span><b>  272</b>  				break;</span>
<span><b>  274</b>  			default: throw std::invalid_argument(&quot;Unsupported B-spline degree&quot;);</span>
<span><b>  386</b>  		const double A1 = coeffs[1] / coeffs[0],</span>
<span><b>  387</b>  		             A2 = coeffs[2] / coeffs[0],</span>
<span><b>  388</b>  		             A3 = coeffs[3] / coeffs[0];</span>
<span><b>  390</b>  		const double Q = (3.0 * A2 - A1 * A1) / 9.0,</span>
<span><b>  391</b>  		             R = (9.0 * A1 * A2 - 27.0 * A3 - 2.0 * A1 * A1 * A1) / 54.0;</span>
<span><b>  392</b>  		const double D = Q * Q * Q + R * R;</span>
<span><b>  394</b>  		std::vector&lt;std::complex&lt;double&gt;&gt; roots;</span>
<span><b>  396</b>  		if (D &gt;= 0) {</span>
<span><b>  397</b>  			const double S = std::cbrt(R + std::sqrt(D));</span>
<span><b>  398</b>  			const double T = std::cbrt(R - std::sqrt(D));</span>
<span><b>  399</b>  			roots.emplace_back(-A1 / 3.0 + (S + T));</span>
<span><b>  400</b>  			roots.emplace_back(-A1 / 3.0 - (S + T) / 2.0 + std::complex&lt;double&gt;(0, std::numbers::sqrt3 * (S - T) / 2.0));</span>
<span><b>  401</b>  			roots.emplace_back(-A1 / 3.0 - (S + T) / 2.0 - std::complex&lt;double&gt;(0, std::numbers::sqrt3 * (S - T) / 2.0));</span>
<span><b>  402</b>  		}</span>
<span><b>  404</b>  			const double Theta = std::acos(R / std::sqrt(-Q * Q * Q));</span>
<span><b>  405</b>  			const double X     = 2.0 * std::sqrt(-Q);</span>
<span><b>  406</b>  			roots.emplace_back(X * std::cos(Theta / 3.0) - A1 / 3.0);</span>
<span><b>  407</b>  			roots.emplace_back(X * std::cos((Theta + 2.0 * std::numbers::pi) / 3.0) - A1 / 3.0);</span>
<span><b>  408</b>  			roots.emplace_back(X * std::cos((Theta + 4.0 * std::numbers::pi) / 3.0) - A1 / 3.0);</span>
<span><b>  410</b>  		return roots;</span>
</pre></details>
</td><td><span>69.97%</span><progress max="100" value="69.9690" aria-label="lib/signal_processing/BSpline3D.cpp line coverage">69.97%</progress></td><td>226 / 323</td></tr>
<tr><td><strong>lib/signal_processing/BSpline3D.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/signal_processing/BSpline3D.hpp line coverage">100.00%</progress></td><td>13 / 13</td></tr>
<tr><td><strong>lib/signal_processing/FFT.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/signal_processing/FFT.hpp line coverage">100.00%</progress></td><td>89 / 89</td></tr>
<tr><td><strong>lib/signal_processing/Interpolation.cpp</strong>
<details class="coverage-lines"><summary>4 uncovered lines</summary><pre>
<span><b>  244</b>  		case interpolation_mode::cubic: interp = &amp;cubic3D;</span>
<span><b>  245</b>  			break;</span>
<span><b>  248</b>  		case interpolation_mode::spline: interp = &amp;spline3D;</span>
<span><b>  269</b>  		return {.x = 0.0, .y = 0.0, .z = 0.0};</span>
</pre></details>
</td><td><span>98.08%</span><progress max="100" value="98.0769" aria-label="lib/signal_processing/Interpolation.cpp line coverage">98.08%</progress></td><td>204 / 208</td></tr>
<tr><td><strong>lib/signal_processing/Interpolation.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/signal_processing/Interpolation.hpp line coverage">100.00%</progress></td><td>4 / 4</td></tr>
<tr><td><strong>lib/signal_processing/Matrix.cpp</strong>
<details class="coverage-lines"><summary>6 uncovered lines</summary><pre>
<span><b>  266</b>  		x = static_cast&lt;double&gt;(n) / 2.0;</span>
<span><b>  267</b>  		y = static_cast&lt;double&gt;(n) / 2.0;</span>
<span><b>  347</b>  				for (size_t k = 0; k &lt;= j - 1; ++k) { tmp1 += in[i + k * n] * in[k + j * n]; }</span>
<span><b>  348</b>  				in[i + j * n] = (1 / in[j + j * n]) * (in[i + j * n] - tmp1);</span>
<span><b>  349</b>  				tmp1          = 0;</span>
<span><b>  350</b>  			}</span>
</pre></details>
</td><td><span>97.67%</span><progress max="100" value="97.6744" aria-label="lib/signal_processing/Matrix.cpp line coverage">97.67%</progress></td><td>252 / 258</td></tr>
<tr><td><strong>lib/signal_processing/Matrix.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/signal_processing/Matrix.hpp line coverage">100.00%</progress></td><td>5 / 5</td></tr>
<tr><td><strong>lib/signal_processing/PSF.cpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/signal_processing/PSF.cpp line coverage">100.00%</progress></td><td>38 / 38</td></tr>
<tr><td><strong>lib/signal_processing/Spline.cpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/signal_processing/Spline.cpp line coverage">100.00%</progress></td><td>252 / 252</td></tr>
<tr><td><strong>lib/signal_processing/Spline.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/signal_processing/Spline.hpp line coverage">100.00%</progress></td><td>14 / 14</td></tr>
<tr><td><strong>lib/Structs.hpp</strong>
<details class="coverage-lines"><summary>10 uncovered lines</summary><pre>
<span><b>  102</b>  	{</span>
<span><b>  103</b>  		return {.x = std::clamp&lt;size_t&gt;(x, 0, other.x), .y = std::clamp&lt;size_t&gt;(y, 0, other.y), .z = std::clamp&lt;size_t&gt;(z, 0, other.z)};</span>
<span><b>  104</b>  	}</span>
<span><b>  131</b>  	{</span>
<span><b>  132</b>  		switch (index) {</span>
<span><b>  133</b>  			case 0: return x;</span>
<span><b>  134</b>  			case 1: return y;</span>
<span><b>  135</b>  			case 2: return z;</span>
<span><b>  136</b>  			default: throw std::out_of_range(&quot;Vec3D index must be 0, 1 or 2.&quot;);</span>
<span><b>  138</b>  	}</span>
</pre></details>
</td><td><span>76.74%</span><progress max="100" value="76.7442" aria-label="lib/Structs.hpp line coverage">76.74%</progress></td><td>33 / 43</td></tr>
<tr><td><strong>lib/tracking/ComputesParams.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/tracking/ComputesParams.hpp line coverage">100.00%</progress></td><td>4 / 4</td></tr>
<tr><td><strong>lib/tracking/Point.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/tracking/Point.hpp line coverage">100.00%</progress></td><td>14 / 14</td></tr>
<tr><td><strong>lib/tracking/TrackingPoint.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/tracking/TrackingPoint.hpp line coverage">100.00%</progress></td><td>4 / 4</td></tr>
<tr><td><strong>lib/tracking/Tracks.cpp</strong>
<details class="coverage-lines"><summary>12 uncovered lines</summary><pre>
<span><b>  365</b>  	return 0;</span>
<span><b>  450</b>  			return a.to &lt; b.to;</span>
<span><b>  495</b>  {</span>
<span><b>  496</b>  	os &lt;&lt; &quot;[\n&quot;;</span>
<span><b>  497</b>  	for (const auto&amp; T : tracks.t_) {</span>
<span><b>  498</b>  		os &lt;&lt; &quot;  Track &quot; &lt;&lt; T.id &lt;&lt; &quot;: | &quot;;</span>
<span><b>  499</b>  		for (const auto&amp; P : T.points) { os &lt;&lt; P &lt;&lt; &quot; | &quot;; }</span>
<span><b>  500</b>  		os &lt;&lt; &quot;\n&quot;;</span>
<span><b>  501</b>  	}</span>
<span><b>  502</b>  	os &lt;&lt; &quot;]&quot;;</span>
<span><b>  503</b>  	return os;</span>
<span><b>  504</b>  }</span>
</pre></details>
</td><td><span>96.41%</span><progress max="100" value="96.4072" aria-label="lib/tracking/Tracks.cpp line coverage">96.41%</progress></td><td>322 / 334</td></tr>
<tr><td><strong>lib/tracking/Tracks.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/tracking/Tracks.hpp line coverage">100.00%</progress></td><td>2 / 2</td></tr>
<tr><td><strong>lib/tracking/TracksFit.cpp</strong>
<details class="coverage-lines"><summary>6 uncovered lines</summary><pre>
<span><b>  112</b>  	return 0;</span>
<span><b>  121</b>  		return false;</span>
<span><b>  127</b>  		return false;</span>
<span><b>  132</b>  		case fit_mode::linear: { break; }</span>
<span><b>  137</b>  				return false;</span>
<span><b>  195</b>  			else { p_[2] = 1.0; }</span>
</pre></details>
</td><td><span>97.20%</span><progress max="100" value="97.1963" aria-label="lib/tracking/TracksFit.cpp line coverage">97.20%</progress></td><td>208 / 214</td></tr>
<tr><td><strong>lib/tracking/TracksFit.hpp</strong>
</td><td><span>100.00%</span><progress max="100" value="100.0000" aria-label="lib/tracking/TracksFit.hpp line coverage">100.00%</progress></td><td>6 / 6</td></tr>
</tbody></table></div>
</div>
\endhtmlonly
