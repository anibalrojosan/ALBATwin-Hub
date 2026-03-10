# Project Roadmap - BioProcess-Twin Hub

## Table of Contents
1. [Phase 0: Foundation & Documentation](#phase-0-foundation--documentation)
2. [Phase 1: Mechanistic Engine (ALBATwin)](#phase-1-mechanistic-engine-albatwin)
3. [Phase 2: Synthetic Data Generation](#phase-2-synthetic-data-generation)
4. [Phase 3: Hybrid Modeling (PINN)](#phase-3-hybrid-modeling-pinn)
5. [Phase 4: User-Ready Product](#phase-4-user-ready-product)

---

## Phase 0: Foundation & Documentation
**Goal:** Establish the theoretical and architectural groundwork.
*   [x] **Requirement Analysis:** Define PRD and User Stories.
*   [x] **Architecture Design:** Modular system design (Architecture.md).
*   [x] **Mathematical Specification:** Transcribe ALBA model equations (Math_Model.md).
*   [ ] **Environment Setup:** Initialize `uv` project, `pyproject.toml`, and git repository.

## Phase 1: Mechanistic Engine (ALBATwin)
**Goal:** Build a robust Python simulator that reproduces the results of Casagli et al. (2021).
*   [ ] **Core Implementation:**
    *   Implement `ReactorConfig` and `StateVector`.
    *   Implement `BioKinetics` (Stoichiometry + Rates).
    *   Implement `HydroChemistry` (pH solver).
*   [ ] **Solver Integration:** Connect components with `scipy.integrate`.
*   [ ] **Testing (TDD):**
    *   Unit tests for individual kinetic functions.
    *   Mass balance verification tests (Critical).
*   [ ] **Validation:** Compare simulation outputs with paper benchmarks (Figure 2 & 3 of Casagli 2021).

## Phase 2: Synthetic Data Generation
**Goal:** Use the validated engine to create a "Ground Truth" dataset for AI training.
*   [ ] **Scenario Definition:** Define ranges for Light, Temperature, and Influent Loads.
*   [ ] **Batch Simulation:** Run parallel simulations covering the operational space.
*   [ ] **Data Processing:** Clean, normalize, and structure the output into training tensors.
*   [ ] **Noise Injection:** Add Gaussian noise to simulate sensor inaccuracy (optional, for robustness testing).

## Phase 3: Hybrid Modeling (PINN)
**Goal:** Develop a Physics-Informed Neural Network that learns from data while respecting physical laws.
*   [ ] **Network Architecture:** Design the neural network (Input: Env/State -> Output: dX/dt or Next State).
*   [ ] **Loss Function Engineering:**
    *   $L_{total} = L_{data} + \lambda \cdot L_{physics}$
    *   Implement $L_{physics}$ using the residual of the Phase 1 differential equations.
*   [ ] **Training:** Train the model using the synthetic dataset.
*   [ ] **Evaluation:** Compare PINN performance vs. Pure Data-Driven models vs. Mechanistic Model.

## Phase 4: User-Ready Product
**Goal:** Deploy the tool for end-users (researchers/operators).
*   [ ] **Dashboard:** Build a Streamlit app for visualization and interaction.
*   [ ] **Dockerization:** Create a `Dockerfile` for easy deployment.
*   [ ] **Documentation:** User manual and API reference.
*   [ ] **Publication:** Prepare technical report or paper on the Hybrid Twin implementation.

---
*Last updated: March 10, 2026 by Anibal Rojo*
