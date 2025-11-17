///-------------------------------------------------------------------------------------------------
///
/// \file TracksFit.hpp
/// \brief Classe d'ajustement non linéaire (Levenberg–Marquardt) pour séries 1D.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 13/10/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <cstdint>
#include <vector>

#include "ComputesParams.hpp"

namespace palm_tracer::tracking {
/// <summary>
/// Moteur d'ajustement non linéaire (Levenberg–Marquardt) pour séries 1D (ex. MSD), supportant trois modèles : linéaire, puissance et exponentiel.
/// Le but est de trouver le vecteur de paramètres <c>p</c> qui minimise la somme des carrés des résidus
/// \f[
///   \Phi(p)=\sum_{i=0}^{N-1}\big(y_i - f(x_i, p)\big)^2
/// \f]
/// puis d'évaluer la qualité via <c>computeR2()</c>.
/// </summary>
class TracksFit {
public:
	/// <summary>Constructeur par défaut de la classe <see cref="TracksFit"/>.</summary>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	explicit TracksFit(const ComputesParams& p);

	/// <summary>Configure le mode d'ajustement et initialise les structures internes en conséquence.\n
	/// Appelle <c><see cref="SetFitLength"/>(...)</c> pour allouer les tableaux dépendant de la fenêtre, puis redimensionne les structures liées au nombre de paramètres (<see cref="pSize_"/>).</summary>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	/// <remarks>
	/// - Met à jour <see cref="mode_"/> puis <see cref="pSize_"/> via <c>getPSize()</c>.\n
	/// - Appelle <c><see cref="SetFitLength"/>(...)</c> : <see cref="x_"/>, <see cref="y_"/>, <see cref="f_"/>, <see cref="df_"/> sont redimensionnés.\n
	/// - Alloue/resize <see cref="p_"/>, <see cref="alpha_"/> (taille \f$p_{Size} \times (p_{Size}+1)\f$), <see cref="alphaR_"/> (\f$p_{Size} \times p_{Size}\f$), <see cref="beta_"/>, <see cref="pTry_"/>.\n
	/// </remarks>
	void SetMode(const ComputesParams& p);

	/// <summary>Définit la longueur de fenêtre <see cref="fitLength_"/> et alloue les buffers de données associés.</summary>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	/// <remarks>
	/// - Redimensionne  <see cref="x_"/>, <see cref="y_"/>, <see cref="f_"/>, à <see cref="fitLength_"/>.\n
	/// - Redimensionne <see cref="df_"/> à \f$ \text{fit}_{Length} \times p_{Size}\f$.\n
	/// - Ne modifie pas <see cref="pSize_"/> (dépend du <see cref="fit_mode"/>). Appeler <c><see cref="SetMode"/>(...)</c> pour changer <see cref="pSize_"/>.\n
	/// - Complexité : \f$\mathcal{O}\!\left(\text{fit}_{Length} \times p_{Size}\right)\f$.
	/// </remarks>
	void SetFitLength(const ComputesParams& p);

	/// <summary>Indique si un scintillement a été détecté.</summary>
	/// <returns><c>True</c> si un scintillement a été détecté; sinon, <c>false</c>.</returns>
	[[nodiscard]] bool HasBlink() const { return hasBlink_; }

	/// <summary>Récupère les paramètres d'un fit.</summary>
	/// <returns>Les paramètres ajustés.</returns>
	[[nodiscard]] std::vector<double> GetP() const { return p_; }

	/// <summary>Exécute l'ajustement non linéaire (Levenberg–Marquardt) pour le mode courant sur le vecteur en entrée.</summary>
	/// <param name="input">Série d'entrées (ex. MSD), de longueur supérieure ou égale à <see cref="fitLength_"/>.</param>
	/// <returns>R² au point final <see cref="p_"/> (via <c><see cref="computeR2"/>(p_)</c>). Renvoie <c>-1</c> si l'entrée est invalide.</returns>
	double Fit(const std::vector<double>& input);

protected:
	fit_mode mode_ = fit_mode::linear; ///< Mode de fitting utilisé.
	size_t pSize_         = 0,                       ///< Nombre total de paramètres du modèle.
	       fitLength_     = 0;                       ///< Taille de la fenêtre d'analyse.
	bool hasBlink_        = false;                   ///< Indique si un scintillement a été détecté.
	std::vector<double> x_,                          ///< Compteur de 1 à <see cref="fitLength_"/>.
	                    y_,                          ///< Valeurs valides dans la fenêtre initiale.
	                    p_,                          ///< Paramètres du modèle.
	                    pTry_,                       ///< Paramètres temporaires du modèle lors des essais.
	                    beta_,                       ///< Vecteur d'incrément des paramètres \f$\Delta p\f$ (taille <see cref="pSize_"/>).
	                    f_;                          ///< Résultat pour chaque élément de la fonction de fitting avec les paramètres du modèle.
	std::vector<std::vector<double>> alpha_,         ///< Matrice augmentée de LM.
	                                 alphaR_,        ///< Matrice des multiplicateurs d'élimination gaussienne.
	                                 df_;            ///< Résultat pour chaque élément de la dérivée de la fonction de fitting avec les paramètres du modèle.

	static constexpr double LAMBDA_FACTOR_LOW  = 2,          ///< factor multiplicatif de lambda minimal.
	                        LAMBDA_FACTOR_HIGH = 10,         ///< factor multiplicatif de lambda maximal.
	                        EPSILON_ERR        = 2 * 1e-6,   ///< Seuil de convergence.
	                        EPSILON_MU         = 10,         ///< Seuil pour \f$\mu\f$.
	                        EPSILON_RHO        = 0.001,      ///< Seuil pour \f$\rho\f$.
	                        EPSILON_STEP       = 1e-6,       ///< Seuil pour un pas trop faible.
	                        EPSILON_ERR_STOP   = 1e-6;       ///< Seuil pour une nouvelle erreur trop faible.
	static constexpr size_t N_ITER             = 500;        ///< Nombre maximal d'itérations autorisées dans l'algorithme d'optimisation.
	static constexpr size_t N_K                = N_ITER / 4; ///< Nombre maximal d'itérations ayant entrainé une mise à jour du modèle.

	/// <summary>Retourne le nombre de paramètres <see cref="pSize_"/> attendu pour le <see cref="fit_mode"/> courant.</summary>
	/// <returns>
	/// \f[
	/// \begin{cases}
	/// 2 & \text{pour fit_mode::linear} \quad & (p=[a, b])\\
	/// 3 & \text{pour fit_mode::power} \quad & (p=[\alpha, a, b])\\
	/// 4 & \text{pour fit_mode::exponential} \quad & (p=[x_0, a, \tau, b])
	/// \end{cases}
	/// \f]
	/// </returns>
	/// <remarks>En cas de  <see cref="mode_"/> non reconnu, retourne <c>0</c> (à éviter en pratique).</remarks>
	[[nodiscard]] size_t getPSize() const;

	/// <summary>Vérifie la validité minimale de l'entrée <paramref name="input"/> vis-à-vis du <c>fit_mode</c> et de <see cref="fitLength_"/>.\n
	/// L'écriture d'un diagnostic sur <c>std::cerr</c> en cas d'échec est commenté, car trop bruyante sur la sortie.</summary>
	/// <param name="input">Série d'observations (ex. MSD). Seuls les <see cref="fitLength_"/> premiers éléments sont requis.</param>
	/// <returns><c>true</c> si l'entrée est compatible ; <c>false</c> sinon.</returns>
	/// <remarks>
	/// Conditions vérifiées :\n
	/// - <see cref="fitLength_"/> \f$\ge\f$ 2.\n
	/// - <c>input.size()</c> \f$\ge\f$ <see cref="fitLength_"/>.\n
	/// - <b>Power</b> : \f$\text{input}_0 \ge 0\f$ et \f$\text{input}_1 \ge 0\f$ (pour \f$\log\big(\frac{\mathrm{MSD}(1)}{\mathrm{MSD}(0)}\big)\f$).\n
	/// - <b>Linear</b>, <b>Exponential</b> : pas d'exigence spécifique supplémentaire ici.\n
	/// <b>Note</b> : cette fonction ne vérifie pas la finitude (<c>NaN/Inf</c>) ni d'autres contraintes.
	/// </remarks>
	[[nodiscard]] bool checkInputValidity(const std::vector<double>& input) const;

	/// <summary>Remplit <see cref="x_"/> et <see cref="y_"/> à partir de l'<paramref name="input"/>  (ex. MSD) en conservant les valeurs strictement supérieures à <c>-1</c>.\n
	/// Et en « tassant » <see cref="y_"/> (et implicitement <see cref="x_"/>) pour supprimer les trous (blinking).</summary>
	/// <param name="input">Série d'entrées de longueur supérieure ou égale à <see cref="fitLength_"/>. Les valeurs supérieures à <c>-1</c> sont considérées valides.</param>
	/// <remarks>
	/// - Pour <c>i=0 .. fitLength_ - 1</c> : \f$X[i] \leftarrow i+1\f$ ; si \f$input[i] > -1\f$ alors \f$Y[write] \leftarrow input[i]\f$, sinon l'élément est sauté.\n
	/// - Met à jour <c>nY_ = fitLength_ - (nombre de sauts)</c>.\n
	/// - Complexité : \f$\mathcal{O}(\text{fit}_{Length})\f$.
	/// </remarks>
	void fillXY(const std::vector<double>& input);

	/// <summary>Initialise le vecteur de paramètres <see cref="p_"/> en fonction du <c>fit_mode</c> courant (heuristiques d'amorçage).</summary>
	/// <remarks>
	/// - <b>Linear</b> : \f$p=[a, b]\f$\n
	/// - <b>Power</b> : \f$p=[\alpha, a, b]\f$\n
	/// - <b>Exponential</b> : \f$p = [x_0, a, \tau, b]\f$
	/// </remarks>
	void initP();

	/// <summary>Évalue la fonction de base du modèle courant et remplit <see cref="f_"/> pour chaque abscisse \f$x_i\f$.</summary>
	/// <param name="p">Paramètres du modèle :
	///   - <c>fit_mode::linear</c> : \f$\mathcal{C} = [a, b]\f$ et \f$f(x)=a \times x + b\f$
	///   - <c>fit_mode::power</c> : \f$\mathcal{C} = [\alpha, a, b]\f$ et \f$f(x)=a \times x^{\alpha} + b\f$
	///   - <c>fit_mode::exponential</c> : \f$\mathcal{C} = [x_0, a, \tau, b]\f$ et \f$f(x) = (b - a) \times \exp\!\Big(-\frac{x - x_0}{\tau}\Big) + a\f$
	/// </param>
	/// <remarks>
	/// - La taille de sortie est <c>f_.size() == x_.size()</c>.
	/// - Si <c>coeffs</c> ne contient pas assez d'éléments pour le mode choisi, la fonction n'écrit rien et retourne immédiatement.
	/// - <b>Power</b> : pour \f$x<0\f$ et un exposant non entier, <c>std::pow</c> peut retourner <c>NaN</c> (comportement identique à C++ standard).
	/// - <b>Exponential</b> : si \f$\tau == 0\f$, l'expression est singulière, aucun traitement spécifique n'est appliqué ici.
	/// </remarks>
	void base(const std::vector<double>& p);

	/// <summary>Calcule la matrice jacobienne <c>df_(n × m)</c> des dérivées partielles \f$\frac{\partial f}{\partial c_j}\f$.
	/// pour chaque point \f$x_i\f$ et chaque paramètre \f$p_j\f$, selon le <c>fit_mode</c> courant.</summary>
	/// <remarks>
	/// - Dimensions : <c>df_.size() == x_.size()</c> et <c>df_[i].size() == pSize_</c>.
	/// - <b>Linear</b> : \f$f(x)=a\times x + b,\quad p=[a, b]\f$.\n
	///   \f[\frac{\partial f}{\partial a} = x,\quad\frac{\partial f}{\partial b} = 1\f]\n
	/// - <b>Power</b> : \f$f(x)=a \times x^{\alpha}+b,\quad p=[\alpha, a, b]\f$\n
	///   \f[\frac{\partial f}{\partial p} = a \times \ln(x) \times x^{\alpha},\quad \frac{\partial f}{\partial a} = x^{\alpha},\quad \frac{\partial f}{\partial b} = 1\f]\n
	///   Pour \f$x \le 0\f$ et \f$p\f$ non entier, <c>std::pow</c> peut produire <c>NaN</c>.\n
	/// - <b>Exponential</b> : \f[f(x)=(b - a)\exp(-\tfrac{x-x_0}{\tau}) + a,\quad p=[x_0, a, \tau, b]\f]\n
	///   \f[\begin{aligned}
	///   \frac{\partial f}{\partial x_0} &= (b-a) \times \frac{1}{\tau} \times E,\quad & \frac{\partial f}{\partial a} &= 1 - E, \\
	///   \frac{\partial f}{\partial \tau} &= (b-a) \times \frac{x-x_0}{\tau^2} \times E,\quad & \frac{\partial f}{\partial b} &= E, \\
	///   E&=\exp\!\big(-\tfrac{x-x_0}{\tau}\big)
	///   \end{aligned}\f]\n
	/// - Aucun <c>try/except</c> : en cas de \f$\tau==0\f$, le modèle est singulier.
	/// </remarks>
	void deriv();

	/// <summary>Construit la matrice augmentée <see cref="alpha_"/> pour un pas de Levenberg–Marquardt, à partir de la prédiction <see cref="f_"/> et du Jacobien <see cref="df_"/>.</summary>
	/// <param name="lambda">Paramètre de Levenberg–Marquardt (\f$\lambda \ge 0\f$). La diagonale est multipliée par \f$(1+\lambda)\f$.</param>
	/// <remarks>On forme la matrice de taille <c><see cref="pSize_"/> × (<see cref="pSize_"/> + 1)</c> :
	/// \f[
	/// \alpha_{jk} \mathrel{+}= \sum_i \frac{\partial f_i}{\partial p_j} \times \frac{\partial f_i}{\partial p_k}
	/// \quad\text{et}\quad
	/// \alpha_{j,pSize} \mathrel{+}= \sum_i (y_i - f_i) \times \frac{\partial f_i}{\partial p_j}.
	/// \f]
	/// La partie inférieure est symétrisée : \f$\alpha_{j,k} \leftarrow \alpha_{k,j}\f$ pour \f$k&lt;j\f$.\n
	/// Enfin, Levenberg–Marquardt : \f$\alpha_{j,j}\leftarrow (1+\lambda) \times \alpha_{j,j}\f$.\n
	/// </remarks>
	void computeAlpha(double lambda);

	/// <summary>Résout le système triangulaire après élimination gaussienne avec pivot partiel échelonné
	/// sur la matrice augmentée <see cref="alpha_"/> (taille <c><see cref="pSize_"/> × (<see cref="pSize_"/> + 1)</c> et remplit <see cref="beta_"/>.</summary>
	/// <returns><c>true</c> si la factorisation et la remontée réussissent, <c>false</c> sinon (dans ce cas <see cref="beta_"/> est nul).</returns>
	/// <remarks> Singularités traitées : pivot quasi nul, diagonale non positive après Levenberg–Marquardt.\n
	/// En cas d'échec, <see cref="beta_"/> ← 0 (pas nul) ; l'appelant peut augmenter <c>lambda</c> et retenter.</remarks>
	bool computeBeta();

	/// <summary>Effectue un pivot partiel échelonné sur la matrice augmentée <see cref="alpha_"/> pour le système linéaire courant, en prenant <paramref name="index"/> comme ligne de départ.</summary>
	/// <param name="index">Ligne de départ. Doit vérifier <c><paramref name="index"/> < <see cref="pSize_"/></c>.</param>
	/// <remarks>La matrice <see cref="alpha_"/> est supposée augmentée : <c>alpha_.cols() == pSize_ + 1</c>.\n
	/// L'algorithme :\n
	/// 1) Pour chaque ligne \f$i \in [\text{index}, p_{Size}-1]\f$, on calcule  \f$m_{i-\text{index}} = \max_j \lvert\alpha_{i, j}\rvert\f$ pour \f$j \in [\text{index}, p_{Size}-1]\f$.\n
	/// 2) Si \f$m[...] \neq 0\f$, on normalise \f$m_{i-\text{index}} = \frac{\lvert\alpha_{i, \text{index}}\rvert}{m_{i-\text{index}}}\f$.\n
	/// 3) On amène en ligne <paramref name="index"/> la ligne ayant la plus grande valeur \f$m[...]\f$ via un échange de lignes,
	///    en échangeant toutes les colonnes \f$j \in [\text{index}, p_{Size}]\f$ (incluant le membre de droite).</remarks>
	void pivot(size_t index);

	/// <summary>Calcule le coefficient de détermination R² entre les observations <see cref="y_"/> et la prédiction <c>f_(p)</c> calculée via <c><see cref="base"/>(<paramref name="p"/>)</c>.\n
	/// R² est défini par : \f[R^2 = 1 - \frac{\sum_i (y_i - f_i)^2}{\sum_i (y_i - \bar{y})^2}\f]  avec \f$\bar{y}\f$ la moyenne des \f$y_i\f$.</summary>
	/// <param name="p">Paramètres du modèle utilisés par <c><see cref="base"/>(...)</c> pour produire la prédiction <see cref="f_"/>. </param>
	/// <returns>La valeur de R².</returns>
	/// <remarks>
	/// - <b>Interprétation</b> : \f$1.0\f$ = ajustement parfait, peut être <b>négatif</b> si le modèle est moins bon que la moyenne.\n
	/// - Ce n'est pas un « chi-deux réduit ». Pour un \f$\chi^2\f$ ou une log-vraisemblance (p.ex. fit gaussien), utiliser une autre métrique.
	/// - Si <see cref="y_"/> est vide, la fonction renvoie <c>NaN</c>.
	/// - Si la variance de <see cref="y_"/> est nulle (tous les \f$y_i\f$ égaux), on renvoie \f$0.0\f$ pour éviter une division par zéro.
	/// - <see cref="y_"/> et <see cref="f_"/> doivent avoir la même taille après l'appel à <c><see cref="base"/>(p)</c>.
	/// </remarks>
	[[nodiscard]] double computeR2(const std::vector<double>& p);
};
} // namespace palm_tracer::tracking
