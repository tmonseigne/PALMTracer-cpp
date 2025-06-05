# Changelog

Toutes les modifications notables de ce projet seront documentées dans ce fichier.

Le format est basé sur [Keep a Changelog](https://keepachangelog.com/fr/1.0.0/)
et ce projet adhère à [Semantic Versioning](https://semver.org/lang/fr/).

## [Unreleased]

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
