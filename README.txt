# Projet Analyse Numérique II - Équations de prédation de Lotka-Volterra

## Auteurs
- Rayane Troudi
- Jassiem Zouga
- Florent Gerbaud

## Date
- 18 avril 2023

## Description
Ce projet se concentre sur la résolution et l'analyse des équations de Lotka-Volterra, qui modélisent les interactions entre prédateurs et proies dans un écosystème. Nous utilisons différentes méthodes numériques, notamment les schémas d'Euler implicite et de Runge-Kutta d'ordre 4, pour simuler ces interactions et visualiser les comportements des populations.

Les équations de Lotka-Volterra traitées ici incluent deux configurations :
1. **Système d'ordre 2** : Modélise l'interaction entre deux espèces (une proie et un prédateur).
2. **Système d'ordre 3** : Introduit une troisième espèce, ajoutant de la complexité au modèle avec une interaction prédateur-prédateur.

## Structure du Projet
### 1. Schémas d'Euler Implicite
- **Lotka-Volterra d'ordre 2** : Présentation du schéma d'Euler pour un système avec deux variables.
- **Lotka-Volterra d'ordre 3** : Extension du modèle d'ordre 2 avec une troisième espèce.

### 2. Schéma de Runge-Kutta d'Ordre 4
- Amélioration de la précision de la simulation par rapport au schéma d'Euler avec un calcul plus complexe, mais plus précis.

### 3. Équilibres
- Calcul des points d'équilibre pour les systèmes d'ordre 2 et 3.
- Analyse de la stabilité des équilibres, avec la détermination des points cols, centres et les conditions de stabilité.

### 4. Interprétations et Visualisations
- Étude des paramètres du système d'équations de Lotka-Volterra pour comprendre leur impact sur les dynamiques des populations.
- Visualisation des résultats avec des graphiques des populations de proies et prédateurs en fonction du temps et de leurs interactions.

### 5. Améliorations du Modèle
- Introduction de nouveaux comportements dans le système d'équations, avec la chasse de plusieurs proies par un même prédateur. Ceci permet de modéliser des interactions plus complexes.

### 6. Limites du Modèle
- Les modèles de Lotka-Volterra présentés ont des limites, notamment dans les cas où la croissance des populations est exponentielle, ce qui n'est pas réaliste dans des environnements où les ressources sont limitées.

## Installation
1. Clonez ce dépôt :
   ```bash
   git clone <url-du-repo>

## utilisations

lodka-volterra-3sys.py -> Euler Implicite LV3
lodka-volterra.py -> Euler Implicite LV2
LV3-RK4.py -> Runge-Kutta 4 LV2 et LV3


ATTENTION : changer la variable "choixSysEspece" à la valeur 2 si LV2. 3 
si LV3.

