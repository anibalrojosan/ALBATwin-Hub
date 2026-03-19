# DEVLOG

This document is a log of the development process of the project. It is used to track the progress of the project and to document the decisions made during the development process and the learned concepts.

## Index

- [2026-03-19 - Sprint 2: Petersen Matrix and Stoichiometry (Mass Balance Verification)](#2026-03-19---sprint-2-petersen-matrix-and-stoichiometry-mass-balance-verification)
- [2026-03-17 - Sprint 1: Implementation of State Vector](#2026-03-17---sprint-1-implementation-of-state-vector)
- [2026-03-16 - Sprint 1: Implementation of Reactor Configuration](#2026-03-16---sprint-1-implementation-of-reactor-configuration)
- [2026-03-13 - Sprint 0: CI/CD and Infrastructure Setup](#2026-03-13---sprint-0-cicd-and-infrastructure-setup)
- [2026-03-12 - Phase 1: Technical Specification & Architecture Definition](#2026-03-12---phase-1-technical-specification--architecture-definition)
- [2026-03-12 - Phase 1: ALBA Model Analysis & Data Digitization](#2026-03-12---phase-1-alba-model-analysis--data-digitization)
- [2026-03-10 - Phase 0: Project Initialization and Foundation](#2026-03-10---phase-0-project-initialization-and-foundation)

---

## [2026-03-19] - Sprint 2: Petersen Matrix and Stoichiometry (Mass Balance Verification)

### Context & Goals
The goal was to implement the **Stoichiometry** module for the ALBA model (`phase1-02`), providing the Petersen matrix (19×17) and composition matrix (6×17) required for mass balance validation. This work builds on the `StateVector` structure and establishes the biological stoichiometric backbone for the simulation engine.

### Technical Implementation
- **Stoichiometry Module (`src/bioprocess_twin/models/stoichiometry.py`):**
  - Implemented `get_composition_matrix()` returning a 6×17 array (rows: COD, O, C, N, P, H).
  - Implemented `get_petersen_matrix()` with 19 biological process rows.
  - Corrected composition matrix values to match supplementary material: `I_H_S_IC=0`, `I_H_S_ND=0.14`, `I_H_S_NH=0.22`, and `S_N2` N=1.
- **Documentation (`docs/STOICHIOMETRY.md`):**
  - Fixed typos in the H row of the composition matrix to align with SI.3.1.
- **Unit Tests (`tests/unit/test_stoichiometry.py`):**
  - Added dimension tests for both matrices.
  - Implemented mass balance conservation test for all 6 elements (COD, O, C, N, P, H) with per-species error reporting.
  - Spot-checked Petersen matrix structure (phototrophic growth, decay, urea hydrolysis).
- **Mass Balance Status:**
  - C, N, P balances close within numerical precision.
  - COD, O, and H show residuals (max ~5.56 for O in AOB process); `atol=6.0` used temporarily pending root cause analysis.

### 💡 Deep Dive: Mass Balance Verification via Matrix Multiplication
A fundamental physics law is the mass conservation law, that states that the mass of a system is conserved (There's no creation or destruction of mass). This law can be summarized as follows: 

The stoichiometric mass balance for each process $i$ and element $k$ is 
$$\sum_j S_{ij} \cdot I_{kj} \approx 0$$

Where $S_{ij}$ is the stoichiometric coefficient of the $i$-th process and the $j$-th element, and $I_{kj}$ is the composition of the $k$-th element in the $j$-th process.

In other words, the sum of the products of the stoichiometric coefficients and the composition matrix should be zero. This expression can be verified by computing `balance = S @ I.T`, which is a matrix multiplication between the Petersen matrix and the composition matrix, where a 19×6 matrix is obtained. Each entry of this matrix should be zero if mass is conserved. This algebraic check catches transcription errors and composition mismatches before they propagate into simulations.

### Next Steps
- **Investigate why mass balances do not close** for COD, O, and H across all processes; root cause may relate to ALBA composition conventions or matrix construction.
- Implement the BioKinetics module with the 19 process rate equations.
- Integrate the Petersen matrix into the ODE right-hand side for dynamic simulation.

---

## [2026-03-17] - Sprint 1: Implementation of State Vector

### Context & Goals
The goal was to implement the `StateVector` component, which serves as the central data structure for the ALBA model. This component must bridge the gap between human-readable biological variables and the numerical arrays required by the ODE solver, while enforcing physical constraints (non-negative concentrations). This completes the core data structures for `phase1-01`.

### Technical Implementation
- **State Vector Schema (`src/bioprocess_twin/core/state.py`):**
  - Implemented `StateVector` using `pydantic` to manage exactly 17 state variables.
  - Enforced non-negative constraints on all variables using `Field(..., ge=0.0)`.
  - Applied `frozen=True` to ensure state immutability at any given time step.
- **Numerical Integration Support:**
  - Added `to_array()` method to export states as NumPy arrays in the strict order required by the stoichiometric matrix.
  - Added `from_array()` class method to reconstruct the `StateVector` from solver outputs.
- **Testing & Quality:**
  - Developed `tests/unit/test_state.py` covering instantiation, array mapping, physical validation, and immutability.
  - Fixed Pydantic V2 deprecation warnings regarding `model_fields` access.
  - Verified 100% test pass rate and code formatting with `ruff`.

### 💡 Deep Dive: Name-to-Index Mapping in Bioprocess Models
In large-scale models like ALBA (17 variables), managing array indices manually is a major source of bugs (e.g., swapping $S_{NH}$ with $S_{NO3}$). By using a Pydantic-based `StateVector`, the "Source of Truth" can be centralized for variable ordering. The `to_array` and `from_array` methods act as a translation layer, allowing the simulation logic to use descriptive names (e.g., `state.X_ALG`) while the numerical engine receives the raw vectors it expects.

### Next Steps
- Begin **Sprint 2** (`phase1-02`): Implementation of the `Stoichiometry` module.
- Define the dynamic construction of the Petersen Matrix using the parameters from `MATH_MODEL.md`.
- Implement mass balance verification tests for the stoichiometric coefficients.

---

## [2026-03-16] - Sprint 1: Implementation of Reactor Configuration

### Context & Goals
The objective was to implement the physical and operational configuration layer for the ALBA model. This layer serves as the foundation for the Digital Twin, ensuring that all physical boundaries (geometry, volume, depth) are correctly validated and immutable during simulation. This work addresses the first part of issue `phase1-01`.

### Technical Implementation
- **Configuration Schema (`src/bioprocess_twin/core/reactor.py`):**
  - Implemented `ReactorConfig` using `pydantic` for strict data validation.
  - Enforced immutability using `ConfigDict(frozen=True)` to prevent accidental changes to the reactor's physical design.
  - Created nested models for `GeometryConfig`, `OperationalConfig`, and `LocationConfig`.
- **Process Logic & Engineering Methods:**
  - Added `calculate_hrt` to compute Hydraulic Retention Time dynamically.
  - Implemented `theoretical_volume` and `volume_discrepancy` properties to handle the difference between design specifications and ideal geometry ($Area \times Depth$).
- **Data Integration:**
  - Created `config/reactor_setup.yaml` with the exact specifications of the Narbonne HRABP reactor (56 m², 17 m³, 0.3 m depth).
  - Added a robust `from_yaml` loader with root-key validation.
- **Testing (TDD):**
  - Developed `tests/unit/test_reactor.py` covering 7 test cases: loading, validation, missing fields, HRT calculation, file-not-found handling, immutability, and volume logic.

### 💡 Deep Dive: Immutability with Pydantic Frozen Models
In a complex bioprocess simulation, maintaining the integrity of physical constants is critical. By using `frozen=True` in our Pydantic models, we ensure that once the `ReactorConfig` is loaded, it cannot be modified by any other part of the system. This prevents a whole class of "silent bugs" where a developer might accidentally overwrite a geometric parameter (like surface area) mid-simulation, which would invalidate all subsequent mass balance and light penetration calculations.

### Next Steps
- Complete **Sprint 1** (`phase1-01`): Implement the `StateVector` class.
- Define the 17 state variables in `src/bioprocess_twin/core/state.py`.
- Develop unit tests to verify name-to-index mapping and NumPy array conversion for the ODE solver.

---

## [2026-03-13] - Sprint 0: CI/CD and Infrastructure Setup

### Context & Goals
To ensure the reliability and maintainability of the BioProcess-Twin Hub, I established a robust development environment before implementing complex biological logic. The goal was to configure a "Modular Monolith" architecture where backend and frontend coexist in a single installable package, backed by automated quality gates (linting and testing).

### Technical Implementation
- **Package Structure:** Created `src/bioprocess_twin` as the main package, containing `core`, `models`, and `ui` submodules to support the unified architecture.
- **Build System:** Configured `pyproject.toml` using `hatchling` as the build backend to support editable installs (`pip install -e .`).
- **Quality Tools:**
  - **Testing:** Configured `pytest` with `pytest-cov` to enforce code coverage. Added `tests/conftest.py` and a sanity check suite.
  - **Linting:** Configured `ruff` to enforce PEP 8 standards and import sorting.
- **CI/CD Pipeline:** Implemented `.github/workflows/ci.yml` using `uv` for ultra-fast dependency resolution and testing on every push.

### 💡 Deep Dive: Modular Monolith with Hatchling
A "Modular Monolith" approach was chosen over separating backend and frontend repositories. By using `hatchling` and defining `src/bioprocess_twin` as a single installable package, the Streamlit UI (`src/bioprocess_twin/ui/app.py`) can import core logic (`src/bioprocess_twin/core/`) directly as a library. This eliminates the need for complex REST APIs for internal communication and simplifies the deployment pipeline to a single artifact.

### Next Steps
- Begin **Sprint 1** (`phase1-01`): Implement `StateVector` and `ReactorConfig`.
- Translate the mathematical constants from `MATH_MODEL.md` into `src/config/constants.yaml`.
- Develop the initial unit tests for the reactor geometry and the initial state vector.

---

## [2026-03-12] - Phase 1: Technical Specification & Architecture Definition

### Context & Goals
The project transitioned from initial research to concrete technical specification. The primary goal was to establish a Single Source of Truth (SSOT) for the ALBA mathematical model to prevent ambiguity during implementation. Additionally, I needed to formalize architectural decisions regarding the simulation engine and data handling to ensure the system can handle the stiffness of bioprocess kinetics.

### Technical Implementation
- **Mathematical Specification (`docs/MATH_MODEL.md`):**
  - Consolidated all kinetic and stoichiometric parameters from supplementary materials into a unified configuration table.
  - Defined the state vector explicitly (17 variables).
  - Wrote out the 19 specific process rate equations ($\rho$) including Liebig's Law and inhibition terms.
- **Architecture Update (`docs/ARCHITECTURE.md`):**
  - Added the `Stoichiometry` module to decouple matrix construction from kinetic rate calculations.
  - Refined the data flow between Solver, BioKinetics, and HydroChemistry.
- **Architectural Decision Records (`docs/adrs/`):**
  - Created 6 ADRs to document critical engineering choices:
    1.  **Engine:** Python + SciPy + Numba (JIT) for performance.
    2.  **Stiffness:** Implicit Solvers (LSODA/Radau) for stability.
    3.  **DAE Strategy:** Nested Algebraic Solver for pH.
    4.  **Config:** YAML + Pydantic for validation.
    5.  **Data:** Apache Parquet for efficient I/O.
    6.  **Units:** Strict SI Base ($g, m^3, d$) standard.

### 💡 Deep Dive: Handling Stiff Bioprocess Systems
Bioprocess models are notoriously "stiff" because they couple phenomena occurring at vastly different time scales: chemical equilibria (milliseconds), gas transfer (minutes), and biomass growth (days). Using standard explicit solvers (like RK45) forces the time step to be infinitesimally small to maintain stability, making long-term simulations impractical. By adopting implicit solvers (LSODA) and decoupling the algebraic pH equations (Nested Solver approach), we ensure the Digital Twin can simulate weeks of operation in seconds.

### Next Steps
- Implement the configuration loader using Pydantic (`src/config/`).
- Create the `Stoichiometry` module to generate the Petersen Matrix dynamically.
- Implement the `BioKinetics` module with Numba-optimized rate equations.

---

## [2026-03-12] - Phase 1: ALBA Model Analysis & Data Digitization

### Context & Goals
To ensure the fidelity of the **ALBATwinHub** simulation engine, it is essential to ground the mechanistic model in verified data. I performed a comprehensive review of the original ALBA model (Casagli et al., 2021), focusing on the mathematical backbone of the system. The goal was to digitize these complex kinetic expressions and stoichiometric data from static PDF formats into structured, machine-readable Markdown, laying the foundation for the first iteration of the simulation engine.

### Technical Implementation
- **Data Digitization:** Extracted and structured critical model parameters from the supplementary material, specifically:
    - **Growth Kinetics:** Digitized the mathematical expressions governing algal growth rates (SI.1.1).
    - **Stoichiometry:** Transcribed the stoichiometric matrix and coefficient values (SI.3) into Markdown tables.
- **Model Verification:** Reviewed the system's stoichiometry to ensure mass balance consistency before implementation.
- **Documentation:** Created structured Markdown files to serve as the "source of truth" for the Python implementation, facilitating easier review and version control compared to the original PDF documents.

### 💡 Deep Dive: Structured Data for Mechanistic Models
Transcribing mathematical models from academic papers to code is a critical step often prone to error. By first converting the raw PDF data into structured Markdown tables (rather than jumping straight to code), we create an intermediate verification layer. This allows for easier peer review of the constants and equations against the original paper, ensuring that the `ALBATwin` engine is built on accurate, traceable definitions of the biological system.

### Next Steps
- Update the `MATH_MODEL.md` file with the verified kinetic equations and stoichiometric matrix.
- Review how to translatate this mathematical expressions into Python code.

---

## 2026-03-10 - Phase 0: Project Initialization and Foundation

### Context & Goals
The project was kicked off with a focus on establishing a professional software engineering foundation before writing functional code. The primary goal was to translate high-level scientific research (Casagli et al., 2021) into a structured technical plan and initialize a reproducible, modern Python development environment.

### Technical Implementation
- **Documentation Suite:** Created core technical documents in the `docs/` directory:
    - `PRD.md`: Defined product vision, user stories, and functional requirements.
    - `ARCHITECTURE.md`: Designed a modular system separating physics, biology, and chemistry using Mermaid.js diagrams.
    - `MATH_MODEL.md`: Transcribed 19 biological processes and pH equilibria into LaTeX format.
    - `ROADMAP.md`: Established a 5-phase execution plan.
- **Environment Management:** Initialized the project using `uv`. Configured `pyproject.toml` with a comprehensive scientific stack: `numpy`, `scipy`, `pandas`, `pyyaml` for the mechanistic engine; `torch` for the PINN module; and `streamlit` for the visualization layer.
- **Git Configuration:** Implemented a robust `.gitignore` to handle large binary PDFs and generated datasets. Performed index cleanup (`git rm --cached`) to ensure reference papers are not tracked by version control while remaining available locally.

### 💡 Deep Dive: Environment Reproducibility with `uv`
I chose `uv` as the primary package manager to solve the common "research code" problem of broken dependencies. By using `uv.lock`, I ensure that the exact versions of `scipy` (for ODE solving) and `torch` (for neural networks) are locked. This is critical for the Hybrid Twin approach, as small changes in underlying library versions can lead to numerical instability in high-dimensional biological simulations.

### Next Steps
The phase 0 is complete. The project is now ready to start the Phase 1: Mechanistic Engine (ALBATwin).
