# DEVLOG

This document is a log of the development process of the project. It is used to track the progress of the project and to document the decisions made during the development process and the learned concepts.

## Index

- [2026-03-12 - Phase 1: ALBA Model Analysis & Data Digitization](#2026-03-12---phase-1-alba-model-analysis--data-digitization)
- [2026-03-10 - Phase 0: Project Initialization and Foundation](#2026-03-10---phase-0-project-initialization-and-foundation)

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
