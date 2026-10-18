# Changelog

Toutes les modifications notables de ce projet seront documentées dans ce fichier.

Le format est basé sur [Keep a Changelog](https://keepachangelog.com/fr/1.0.0/)
et ce projet adhère à [Semantic Versioning](https://semver.org/lang/fr/).

## [Unreleased]

<!-- run command git log v1.3.0..HEAD --pretty=format:"- %s"-->

## [1.3.0] - 2025-10-29

### Ajouté

- Ajout de plusieurs helper d'affichage de matrices pour le debug.
- Ajout d'une classe dédiée à l'analyse des trajectoires.
- Ajout d'une fonction de reconnexion de trajectoires suite au scintillements. ([#27](https://github.com/tmonseigne/PALMTracer-cpp/pull/27))
- Ajout des calculs de différentes métriques sur les trajectoires. ([#29](https://github.com/tmonseigne/PALMTracer-cpp/pull/29))
- Ajout de nombreux tests
- Exécution des tests sur Windows, Mac et Unix

### Modifié

- Unification des noms des structures (vecteur de dimension et de coordonnées) et de leur utilisation.
- Amélioration du calcul du Gaussian Fit. ([#24](https://github.com/tmonseigne/PALMTracer-cpp/pull/24))
- Suppression d'éléments inutiles dans le code.

## [1.2.1] - 2025-07-28

### Ajouté

- Ajout d'une première calibration pour l'ajustement par Splines à améliorer avec des tests poussés. ([#18](https://github.com/tmonseigne/PALMTracer-cpp/pull/18))
- Ajout de nombreuses fonctions de traitement d'images : Interpolation, FFT, opérations sur matrices
- Ajout du calcul automatique du seuil en interne.
- Ajout de tests unitaires simples.

### Modifié

- Première correction pour l'ajustement par Splines, descente de gradient à améliorer. ([#21](https://github.com/tmonseigne/PALMTracer-cpp/pull/21))
- Corrections mineures sur les Clang code Issues

## [1.2.0] - 2025-06-18

### Ajouté

- Ajout d'un premier ajustement par Splines à améliorer avec des tests poussés. ([#15](https://github.com/tmonseigne/PALMTracer-cpp/pull/15))

### Modifié

- Correction d'une erreur critique sur la segmentation pouvant faire crasher le traitement dans certains cas.
- Corrections mineures sur les Clang code Issues

## [1.1.0] - 2025-06-05

### Ajouté

- Ajout d’un système d’**alignement d’images** incluant suréchantillonnage (`Upsampling`) et transformation avec pré-calcul de la grille. ([#14](https://github.com/tmonseigne/PALMTracer-cpp/pull/14))
- Application de traitements sur l'ensemble de la surface d'analyse.
- Ajout de **tests unitaires** basiques et avancés incluant des translations pour l'alignement.
- Nouvelle **classe de fitting** pour ajustement de données segmentées en vue de l'ajout des Splines. ([#12](https://github.com/tmonseigne/PALMTracer-cpp/pull/12))

### Modifié

- Appels au processus modifiés pour utiliser une **pile d’images** au lieu d’une seule. ([#9](https://github.com/tmonseigne/PALMTracer-cpp/pull/9))
- Passage des paramètres de configuration à l’appel du processus via `settings`. ([#10](https://github.com/tmonseigne/PALMTracer-cpp/pull/10))
- Réorganisation des constantes dans le traitement.
- Optimisations diverses du code, notamment via **OpenMP** pour le parallélisme.

## [1.0.0] - 2025-05-19

### Ajouté

- Première version stable de PALM.
- Intégration continue avec GitHub Actions.
- Documentation générée avec Doxygen, disponible dans le dossier `Docs/`.

### Modifié

- Ajout de la configuration `pre-commit` pour assurer la qualité du code.

### Supprimé

- N/A
