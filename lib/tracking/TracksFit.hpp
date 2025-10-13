///-------------------------------------------------------------------------------------------------
///
/// \file TracksFit.hpp
/// \brief Classe d’ajustement non linéaire (Levenberg–Marquardt) pour séries 1D.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <cstdint>
#include <vector>

namespace palm_tracer {
/// <summary>Mode du modèle d'ajustement pour les calculs statistiques des trajectoires.
/// L'ajustement cherche le vecteur de paramètres <c>p</c> qui minimise :
/// \f[ \Phi(p) = \sum_{i=0}^{N-1} \big(y_i - f(x_i; p)\big)^2 \f]
/// (moindres carrés non pondérés). La qualité peut être évaluée via <c>computeR2()</c>.
/// </summary>
/// <remarks>Plus de cas avec le domaine spectral et les histogrammes seront implémantés plus tard :
/// Spectral : Bi-dimensional Polynomial, Translation, Affine, Affine (Translation), Polywarp.
/// Histogramm : Gaussian, Bi-Gaussian (donc Histo_Gaussian, Histo_BiGaussian).
///	Autres méthodes appelées directement : Third Order Polynomial, 3D Astygmatism, NeNa, Gaussian.
/// </remarks>
enum class fit_mode : uint8_t {
	linear,     ///< \f$p = [a, b]\f$ pour optimiser la fonction \f$ f(x) = a \times x + b\f$
	power,      ///< \f$p = [p, a, b]\f$ pour optimiser la fonction \f$ f(x) = a \times x^{p} + b\f$
	exponential ///< \f$p = [x_0, \tau, a, b]\f$ pour optimiser la fonction \f$ f(x) = (b - a) \times \exp\Big(-\tfrac{x - x_0}{\tau}\Big) + a \f$
};

/// <summary>
/// Moteur d’ajustement non linéaire (Levenberg–Marquardt) pour séries 1D (ex. MSD), supportant trois modèles : linéaire, puissance et exponentiel.
/// Le but est de trouver le vecteur de paramètres <c>p</c> qui minimise la somme des carrés des résidus
/// \f[
///   \Phi(p)=\sum_{i=0}^{N-1}\big(y_i - f(x_i; p)\big)^2,
/// \f]
/// puis d’évaluer la qualité via <c>computeR2()</c>.
/// </summary>
class TracksFit {
public:
	/// <summary>Constructeur par défaut de la classe <see cref="TracksFit"/>.</summary>
	/// <param name="mode">Le mode d'ajustement (voir <see cref="fit_mode"/>).</param>
	/// <param name="fitParams">Paramètres de configuration de l'ajustement.
	/// Actuellement, <c>fitParams[0]</c> est interprété comme la longueur de fenêtre <c>fitLength_</c> (nombre d’échantillons utilisés).</param>
	explicit TracksFit(fit_mode mode, const std::vector<double>& fitParams);

	/// <summary>Configure le mode d’ajustement et initialise les structures internes en conséquence.\n
	/// Appelle <c>SetFitLength(...)</c> pour allouer les tableaux dépendant de la fenêtre, puis redimensionne les structures liées au nombre de paramètres (<c>pSize_</c>).</summary>
	/// <param name="mode">Le mode d'ajustement (voir <see cref="fit_mode"/>).</param>
	/// <param name="fitParams">Paramètres de configuration de l'ajustement.
	/// Actuellement, <c>fitParams[0]</c> est interprété comme la longueur de fenêtre <c>fitLength_</c> (nombre d’échantillons utilisés).</param>
	/// <remarks>
	/// - Met à jour <c>mode_</c> puis <c>pSize_</c> via <c>getPSize()</c>.\n
	/// - Appelle <c>SetFitLength(...)</c> : <c>x_</c>, <c>y_</c>, <c>f_</c>, <c>df_</c> sont redimensionnés.\n
	/// - Alloue/resize <c>p_</c>, <c>alpha_</c> (taille \f[pSize_ \times (pSize_+1)\f]), <c>alphaR_</c> (\f[pSize_ \times pSize_\f]), <c>beta_</c>, <c>pTry_</c>.\n
	/// - Ne valide pas le contenu de <c>fitParams</c> (au-delà de l’index 0) : à compléter si d’autres options sont ajoutées.
	/// </remarks>
	void SetMode(fit_mode mode, const std::vector<double>& fitParams);

	/// <summary>Définit la longueur de fenêtre <c>fitLength_</c> et alloue les buffers de données associés.</summary>
	/// <param name="fitLength">Longueur de fenêtre (nombre d’échantillons) utilisée par l'ajustement.</param>
	/// <remarks>
	/// - Redimensionne  <c>x_</c>, <c>y_</c>, <c>f_</c>, à <c>fitLength_</c>.\n
	/// - Redimensionne <c>df_</c> à \f[ \text{fitLength\_} \times pSize_ \f].\n
	/// - Ne modifie pas <c>pSize_</c> (dépend du <see cref="fit_mode"/>). Appeler <c>SetMode(...)</c> pour changer <c>pSize_</c>.\n
	/// - Complexité : \f$\mathcal{O}\!\left(fitLength \times pSize\right)\f$.
	/// </remarks>
	void SetFitLength(size_t fitLength);

	/// <summary>Indique si un scintillement a été détecté.</summary>
	/// <returns><c>true</c> si un scintillement a été détecté; sinon, <c>false</c>.</returns>
	[[nodiscard]] bool HasBlink() const { return fitLength_ != nY_; }

	/// <summary>Récupère le premier paramètre d'un fit de MSD.</summary>
	/// <returns>Le premier paramètre du Fit qui correspond à :
	/// - <c>fit_mode::linear</c> : la pente
	/// - <c>fit_mode::power</c> : l'exposant
	/// - <c>fit_mode::exponential</c> : le décalage horizontal.
	/// </returns>
	[[nodiscard]] double GetP0() const { return p_[0]; }

	/// <summary>Exécute l'ajustement non linéaire (Levenberg–Marquardt) pour le mode courant sur le vecteur en entrée.</summary>
	/// <param name="input">Série d’entrée (ex. MSD), de longueur ≥ <c>fitLength_</c>.</param>
	/// <returns>R² au point final <c>p_</c> (via <c>computeR2(p_)</c>). Renvoie <c>-1</c> si l’entrée est invalide.</returns>
	double Fit(const std::vector<double>& input);

protected:
	fit_mode mode_ = fit_mode::linear;        ///< Mode de fitting utilisé
	size_t pSize_{},                          ///< Nombre total de paramètres du modèle.
	       fitLength_{},                      ///< Taille de la fenêtre d'analyse
	       nY_{};                             ///< Nombre d'éléments valide dans la fenêtre initiale
	std::vector<double> x_,                   ///< Compteur de 1 à fitLength_.
	                    y_,                   ///< Valeurs valides dans la fenêtre initiale.
	                    p_,                   ///< Paramètres du modèle
	                    pTry_,                ///< Paramètres temporaires du modèle lors des essais.
	                    beta_,                ///< Vecteur d’incrément des paramètres Δp (taille pSize_).
	                    f_;                   ///< Résultat pour chaque élément de la fonction de fitting avec les paramètres du modèle.
	std::vector<std::vector<double>> alpha_,  ///< Matrice augmentée de LM.
	                                 alphaR_, ///< Matrice des multiplicateurs d’élimination gaussienne.
	                                 df_;     ///< Résultat pour chaque élément de la dérivée de la fonction de fitting avec les paramètres du modèle.

	static constexpr double LAMBDA_FACTOR_LOW  = 2,          ///< factor multiplicatif de lambda minimal
	                        LAMBDA_FACTOR_HIGH = 10,         ///< factor multiplicatif de lambda maximal
	                        EPSILON_ERR        = 2 * 1e-6,   ///< Seuil de convergence.
	                        EPSILON_MU         = 10,         ///< Seuil pour mu.
	                        EPSILON_RHO        = 0.001;      ///< Seuil pour Rho.
	static constexpr size_t N_ITER             = 500;        ///< Nombre maximal d'itérations autorisées dans l'algorithme d'optimisation.
	static constexpr size_t N_K                = N_ITER / 4; ///< Nombre maximal d'itérations ayant entrainé une mise à jour du modèle.

	/// <summary>Retourne le nombre de paramètres <c>pSize_</c> attendu pour le <see cref="fit_mode"/> courant.</summary>
	/// <returns>
	/// \f[
	/// \begin{cases}
	/// 2 & \text{pour } \c fit\_mode::linear \quad (p=[a,b])\\[2pt]
	/// 3 & \text{pour } \c fit\_mode::power \quad (p=[p,a,b])\\[2pt]
	/// 4 & \text{pour } \c fit\_mode::exponential \quad (p=[x_0,a,\tau,b])
	/// \end{cases}
	/// \f]
	/// </returns>
	/// <remarks>
	/// - En cas de  <c>mode_</c> non reconnu, retourne <c>0</c> (à éviter en pratique).
	/// </remarks>
	[[nodiscard]] size_t getPSize() const;

	/// <summary>Vérifie la validité minimale de l’entrée <c>input</c> vis-à-vis du <c>fit_mode</c> et de <c>fitLength_</c>.\n
	/// Écrit un diagnostic sur <c>std::cerr</c> en cas d’échec.</summary>
	/// <param name="input">Série d’observations (ex. MSD). Seuls les <c>fitLength_</c> premiers éléments sont requis.</param>
	/// <returns><c>true</c> si l’entrée est compatible ; <c>false</c> sinon.</returns>
	/// <remarks>
	/// Conditions vérifiées :\n
	/// - <c>fitLength_</c> \f$\ge\f$ 2.\n
	/// - <c>input.size()</c> \f$\ge\f$ <c>fitLength_</c>.\n
	/// - <b>Power</b> : <c>input[0] < 0 </c> et <c>input[1] < 0 </c> (pour \f$\log(\mathrm{MSD}(1)/\mathrm{MSD}(0))\f$).\n
	/// - <b>Linear</b>, <b>Exponential</b> : pas d’exigence spécifique supplémentaire ici.\n
	/// <b>Note</b> : cette fonction ne vérifie pas la finitude (NaN/Inf) ni d’autres contraintes métier ; ajouter si nécessaire.
	/// </remarks>
	[[nodiscard]] bool checkInputValidity(const std::vector<double>& input) const;

	/// <summary>
	/// Remplit <c>X</c> et <c>Y</c> à partir de <c>input</c> (ex. MSD) en conservant les valeurs strictement supérieures à -1.\n
	/// Et en « tassant » <c>Y</c> (et implicitement <c>X</c>) pour supprimer les trous (blinking).</summary>
	/// <param name="input">Série d’entrée de longueur supérieure ou égale à <c>fitLength_</c>. Les valeurs suéprieure à -1 sont considérées valides.
	/// </param>
	/// <remarks>
	/// - Pour <c>i=0 .. fitLength_ - 1</c> : \f$X[i] \leftarrow i+1\f$ ; si \f$input[i] > -1\f$ alors \f$Y[write] \leftarrow input[i]\f$, sinon l’élément est sauté.\n
	/// - Met à jour <c>nY_ = fitLength_ - (nombre de sauts)/c>.\n
	/// - <b>Important</b> : cette implémentation ne redimensionne pas les buffers ( <c>x_</c>, <c>y_</c>, <c>f_</c>, <c>df_</c>).\n
	///   Les éléments au-delà de <c>nY_</c> peuvent rester à 0 ; certaines routines VB parcourent <c>Y.Length</c> (équivalent <c>fitLength_</c> ici).\n
	///   // TODO : clarifier la convention finale « <c>fitLength_</c> vs <c>nY_</c> » pour les boucles d’évaluation (base/deriv/R²/alpha).\n
	/// - Complexité : \f$\mathcal{O}\!\left(fitLength\right)\f$.
	/// </remarks>
	void fillXY(const std::vector<double>& input);

	/// <summary>Initialise le vecteur de paramètres <c>>p_</c> en fonction du <c>fit_mode</c> courant (heuristiques d’amorçage).</summary>
	/// <param name="input">Série d’observations de longueur attendue <c>fitLength_</c>.\n Utilisée pour estimer les valeurs initiales.</param>
	/// <remarks>
	/// - <b>Linear</b> : \f$p=[a, b]\f$\n
	/// - <b>Power</b> : \f$p=[p, a, b]\f$\n
	/// - <b>Exponential</b> : \f$p = [x_0, \tau, a, b]\f$
	/// - Aucune vérification forte n’est faite sur les zéros/divisions : on réplique l’esprit du VB mais ça sent pas bon.
	/// </remarks>
	void initP(const std::vector<double>& input);

	/// <summary>Évalue la fonction de base du modèle courant et remplit <c>f_</c> pour chaque abscisse <c>x_[i]</c>.</summary>
	/// <param name="p">Paramètres du modèle :
	///   - <c>fit_mode::linear</c> : \f$coeffs = [a, b]\f$ et \f$f(x)=a\,x+b\f$
	///   - <c>fit_mode::power</c> : \f$coeffs = [p, A, c]\f$ et \f$f(x)=A\,x^{p}+c\f$
	///   - <c>fit_mode::exponential</c> : \f$coeffs = [x_0, c_1, \tau, c_3]\f$ et \f$f(x)=(c_3-c_1)\,\exp\!\Big(-\frac{x-x_0}{\tau}\Big)+c_1\f$
	/// </param>
	/// <remarks>
	/// - La taille de sortie est \f$f_.size() == x_.size()\f$.
	/// - Si <c>coeffs</c> ne contient pas assez d’éléments pour le mode choisi, la fonction n’écrit rien et retourne immédiatement.
	/// - <b>Power</b> : pour \f$x<0\f$ et un exposant non entier, <c>std::pow</c> peut retourner <c>NaN</c> (comportement identique à C++ standard).
	/// - <b>Exponential</b> : si \f$tau == 0\f$, l’expression est singulière, aucun traitement spécifique n’est appliqué ici.
	/// </remarks>
	void base(const std::vector<double>& p);

	/// <summary>Calcule la matrice jacobienne <c>df_(n × m)</c> des dérivées partielles \f$\partial f / \partial c_j\f$.
	/// pour chaque point <c>x_[i]</c> et chaque paramètre <c>p_[j]</c>, selon le <c>fit_mode</c> courant.</summary>
	/// <remarks>
	/// - Dimensions : <c>df_.size() == x_.size()</c> et <c>df_[i].size() == pSize_</c>.
	/// - <b>Linear</b> : \f$f(x)=a\,x+b,\quad p=[a,b]\f$.\n
	///   \f[\partial f/\partial a = x,\quad \partial f/\partial b = 1\f]\n
	/// - <b>Power</b> : \f$f(x)=A\,x^{p}+c,\quad p=[p,A,c]\f$\n
	///   \f[\partial f/\partial p = A\,\ln(x)\,x^{p},\quad \partial f/\partial A = x^{p},\quad \partial f/\partial c = 1\f]\n
	///   Pour <c>x≤0</c> et <c>p</c> non entier, <c>std::pow</c> peut produire <c>NaN</c>.\n
	/// - <b>Exponential</b> : \f[f(x)=(c_3-c_1)\exp(-\tfrac{x-c_0}{c_2})+c_1,\quad p=[c_0,c_1,c_2,c_3]\f]\n
	///   \f[\begin{aligned}
	///   \partial f/\partial c_0 &= (c_3-c_1)\,\frac{1}{c_2}\,E,\\
	///   \partial f/\partial c_1 &= 1 - E,\\
	///   \partial f/\partial c_2 &= (c_3-c_1)\,\frac{x-c_0}{c_2^2}\,E,\\
	///   \partial f/\partial c_3 &= E,\quad E=\exp\!\big(-\tfrac{x-c_0}{c_2}\big)
	///   \end{aligned}\f]\n
	/// - Aucun <c>try/except</c> : en cas de <c>c2==0</c>, le modèle est singulier.
	/// </remarks>
	void deriv();

	/// <summary>Construit la matrice augmentée <c>alpha_</c> pour un pas de Levenberg–Marquardt, à partir de la prédiction <c>f_</c> et du Jacobien <c>df_</c>.</summary>
	/// <param name="lambda">Paramètre de Levenberg–Marquardt (\f$\lambda \ge 0\f$). La diagonale est multipliée par \f$(1+\lambda)\f$.</param>
	/// <remarks>On forme la matrice de taille <c>pSize_ × (pSize_ + 1)</c> :
	/// \f[
	/// \alpha_{jk} \mathrel{+}= \sum_i \frac{\partial f_i}{\partial p_j}\,\frac{\partial f_i}{\partial p_k}
	/// \quad\text{et}\quad
	/// \alpha_{j,\;pSize_} \mathrel{+}= \sum_i (y_i - f_i)\,\frac{\partial f_i}{\partial p_j}.
	/// \f]
	/// La partie inférieure est symétrisée : \f$\alpha_{jk} \leftarrow \alpha_{kj}\f$ pour \f$k&lt;j\f$.\n
	/// Enfin, LM : \f$\alpha_{jj}\leftarrow (1+\lambda)\,\alpha_{jj}\f$.\n
	/// <b>Bornes de boucle :</b> par défaut on itère sur <c>fitLength_</c> (VB emploie <c>Y.Length</c>).
	/// Si l’on souhaite n’utiliser que les points « valides » compactés, remplacer par <c>nY_</c>.\n
	/// TODO: valider la convention finale.
	/// </remarks>
	void computeAlpha(double lambda);

	/// <summary>Résout le système triangulaire après élimination gaussienne avec pivot partiel échelonné
	/// sur la matrice augmentée <c>alpha_</c> (taille <c>pSize_ × (pSize_ + 1))</c> et remplit <c>beta_</c>.</summary>
	/// <returns><c>true</c> si la factorisation + remontée réussissent, <c>false</c> sinon (dans ce cas <c>beta_</c> est nul).</returns>
	/// <remarks> Singularités traitées : pivot quasi nul, diagonale non positive après LM.\n
	/// En cas d’échec, <c>beta_</c> ← 0 (pas nul) ; l’appelant peut augmenter <c>lambda</c> et retenter.</remarks>
	bool computeBeta();

	/// <summary>Effectue un pivoteo parcial escalado (pivot partiel échelonné) sur la matrice augmentée <c>alpha_</c> pour le système linéaire courant, en prenant <c>index</c> comme ligne de départ.</summary>
	/// <param name="index">Ligne de départ. Doit vérifier <c>index < pSize_</c>.</param>
	/// <remarks>La matrice <c>alpha_</c> est supposée augmentée : <c>alpha_.cols() == pSize_ + 1</c>.\n
	/// L’algorithme :\n
	/// 1) Pour chaque ligne <c>i ∈ [index, pSize_-1]</c>, on calcule <c>m[i-index] = max_j |alpha_(i,j)|</c> pour <c>j ∈ [index, pSize_-1]</c>.\n
	/// 2) Si <c>m[...] ≠ 0</c>, on normalise <c>m[i-index] ← |alpha_(i,index)| / m[i-index]</c>.\n
	/// 3) On amène en ligne <c>index</c> la ligne ayant la plus grande valeur <c>m[...]</c> via un échange de lignes,
	///    en échangeant toutes les colonnes <c>j ∈ [index, pSize_]</c> (incluant le membre de droite).</remarks>
	void pivot(size_t index);

	/// <summary>Calcule le coefficient de détermination R² entre les observations <c>y_</c> et la prédiction <c>f_(p)</c> calculée via <c>base(p)</c>.\n
	/// R² est défini par : \f[R^2 = 1 - \frac{\sum_i (y_i - f_i)^2}{\sum_i (y_i - \bar{y})^2}\f]  avec \f$\bar{y}\f$ la moyenne des \f$y_i\f$.</summary>
	/// <param name="p">Paramètres du modèle utilisés par <c>base(...)</c> pour produire la prédiction <c>f_</c>. </param>
	/// <returns>La valeur de R².</returns>
	/// <remarks>
	/// - <b>Interprétation</b> : <c>1.0</c> = ajustement parfait, peut être <b>négatif</b> si le modèle est moins bon que la moyenne.\n
	/// - Ce n’est pas un « chi-deux réduit ». Pour un \f$\chi^2\f$ ou une log-vraisemblance (p.ex. fit gaussien), utiliser une autre métrique.
	/// - Si <c>y_</c> est vide, la fonction renvoie <c>NaN</c>.
	/// - Si la variance de <c>y_</c> est nulle (tous les <c>y_i</c> égaux), on renvoie <c>0.0</c> pour éviter une division par zéro.
	/// - <c>y_</c> et <c>f_</c> doivent avoir la même taille après l’appel à <c>base(p)</c>.
	/// </remarks>
	[[nodiscard]] double computeR2(const std::vector<double>& p);
};
} // namespace palm_tracer
