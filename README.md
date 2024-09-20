# Numerical Analysis II Project - Lotka-Volterra Predation Equations


## Table of Contents
- [Authors](#authors)
- [Date](#date)
- [Description](#description)
- [Project Structure](#project-structure)
  - [1. Implicit Euler Schemes](#1-implicit-euler-schemes)
  - [2. Fourth-Order Runge-Kutta Scheme](#2-fourth-order-runge-kutta-scheme)
  - [3. Equilibria](#3-equilibria)
  - [4. Interpretations and Visualizations](#4-interpretations-and-visualizations)
  - [5. Model Improvements](#5-model-improvements)
  - [6. Model Limitations](#6-model-limitations)
- [utilisations](#utilisations)
- 

## Authors
- Rayane Troudi
- Jassiem Zouga
- Florent Gerbaud

## Date
- April 18, 2023

## Description
This project focuses on solving and analyzing the Lotka-Volterra equations, which model the interactions between predators and prey in an ecosystem. We use various numerical methods, including the implicit Euler scheme and the fourth-order Runge-Kutta method, to simulate these interactions and visualize the behaviors of the populations.

The Lotka-Volterra equations addressed here include two configurations:
1. **Second-Order System**: Models the interaction between two species (one prey and one predator).
2. **Third-Order System**: Introduces a third species, adding complexity to the model with predator-predator interaction.

## Project Structure
### 1. Implicit Euler Schemes
- **Second-Order Lotka-Volterra**: Presentation of the Euler scheme for a system with two variables.
- **Third-Order Lotka-Volterra**: Extension of the second-order model with a third species.

### 2. Fourth-Order Runge-Kutta Scheme
- Improvement of the simulation accuracy compared to the Euler scheme with a more complex but more precise calculation.

### 3. Equilibria
- Calculation of equilibrium points for the second-order and third-order systems.
- Analysis of the stability of equilibria, determining saddle points, centers, and stability conditions.

### 4. Interpretations and Visualizations
- Study of the parameters of the Lotka-Volterra system of equations to understand their impact on population dynamics.
- Visualization of results with graphs of prey and predator populations over time and their interactions.

### 5. Model Improvements
- Introduction of new behaviors in the system of equations, with the hunting of multiple prey by the same predator. This allows for modeling more complex interactions.

### 6. Model Limitations
- The presented Lotka-Volterra models have limitations, particularly in cases where population growth is exponential, which is not realistic in environments with limited resources.

## utilisations

1. change the parameters of the simulation in the code
   - you have some exemples of parameters that is working in the " Valeurs de Test et Sorties " section
2. run the code

   - lodka-volterra-3sys.py -> Euler Implicite LV3
   - lodka-volterra.py -> Euler Implicite LV2
   - LV3-RK4.py -> Runge-Kutta 4 LV2 et LV3


WARNING: change the “choixSysEspece” variable to value 2 if LV2. 3 
if LV3.
