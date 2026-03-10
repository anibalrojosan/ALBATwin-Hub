# Product Requirements Document (PRD) - BioProcess-Twin Hub

## Table of Contents
1. [Executive Summary](#executive-summary)
2. [Objectives](#objectives)
3. [User Stories](#user-stories)
4. [Functional Requirements](#functional-requirements)
5. [Non-Functional Requirements](#non-functional-requirements)
6. [Tech Stack](#tech-stack)

---

## Executive Summary
**BioProcess-Twin Hub** is a high-fidelity Digital Twin project designed to simulate and optimize the ALBA (Algae-Bacteria) wastewater treatment model. It integrates a rigorous mechanistic model (ODEs/DAEs) with a Physics-Informed Neural Network (PINN). The system aims to simulate the complex interactions between microalgae and bacteria in an open raceway pond, providing a tool for researchers and operators to predict system behavior under varying environmental conditions and optimize operational parameters.

## Objectives

### Scientific Objectives
*   **Validate the ALBA Model:** Accurately reproduce the biological and chemical dynamics described in Casagli et al. (2021).
*   **Hybrid Modeling:** Demonstrate the effectiveness of PINNs in bioprocess engineering by enforcing physical laws (mass balance) within the neural network loss function.
*   **Understanding Interactions:** Quantify the symbiotic relationship between algae ($O_2$ producers) and bacteria ($CO_2$ producers) under dynamic conditions.

### Business/Operational Objectives
*   **Process Optimization:** Reduce energy consumption (aeration) while maintaining treatment efficiency.
*   **Scenario Analysis:** Enable "What-if" simulations for extreme weather events or operational failures.
*   **Synthetic Data Generation:** Create a robust dataset ("Ground Truth") to train AI models where real-world data is scarce or noisy.

## User Stories

### Researcher
*   "As a researcher, I want to modify the incident light intensity ($I_0$) and temperature profiles to study their impact on algal growth rates."
*   "As a researcher, I want to access the internal state variables (e.g., specific growth rates $\mu$, inhibition factors) to debug the biological kinetics."
*   "As a researcher, I want to compare the mechanistic model outputs against the PINN predictions to validate the hybrid approach."
*   "As a researcher, I want to assess how the geometry of the reactor affects the performance of the system."

### Operator
*   "As an operator, I want to simulate the effect of turning off the paddle wheel at night to save energy without crashing the oxygen levels."
*   "As an operator, I want to visualize the pH and Dissolved Oxygen ($DO$) dynamics over a 24-hour cycle to prevent inhibition."
*   "As an operator, I want to input real weather forecasts (CSV) to predict the reactor's performance for the next week."

## Functional Requirements

### Core Simulation Engine
*   **ODE Solver:** Must solve 19 coupled differential equations representing biological processes.
*   **Algebraic Solver:** Must solve the pH system (chemical equilibria) at each time step.
*   **Dynamic Inputs:** Must accept time-series data for light, temperature, evaporation, and rainfall.
*   **Mass Balance:** Must strictly conserve mass for Carbon, Nitrogen, Phosphorus, Hydrogen, and Oxygen.

### Hybrid Module (PINN)
*   **Physics Loss:** The neural network training must include a loss term representing the residual of the ALBA differential equations.
*   **Data Ingestion:** Must be able to train on both synthetic data (from the engine) and potentially real experimental data.

### Interface
*   **Configuration:** Simulation parameters (geometry, initial conditions, kinetic constants) must be configurable via configuration files (YAML/JSON) or UI.
*   **Visualization:** A dashboard must display time-series plots for biomass, substrates, and gases.

## Non-Functional Requirements
*   **Performance:** A 1-year simulation should complete in under 5 minutes on a standard workstation.
*   **Accuracy:** Mass balance error must be less than 0.1% throughout the simulation.
*   **Maintainability:** Code must follow modular architecture principles, separating physics, biology, chemistry, and numerical methods.
*   **Reproducibility:** The environment must be strictly managed using `uv` to ensure consistent dependency versions.

## Tech Stack
*   **Language:** Python 3.10+
*   **Package Management:** `uv`
*   **Scientific Computing:** `numpy`, `scipy` (integration), `pandas` (data handling)
*   **Machine Learning:** `pytorch` or `DeepXDE`
*   **Visualization:** `streamlit`, `matplotlib`/`plotly`
*   **Containerization:** `Docker`

---
*Last updated: March 10, 2026 by Anibal Rojo*
