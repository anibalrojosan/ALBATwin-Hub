# DEVLOG

This document is a log of the development process of the project. It is used to track the progress of the project and to document the decisions made during the development process and the learned concepts.

## Index

- [2026-03-10 - Phase 0: Project Initialization and Foundation](#2026-03-10---phase-0-project-initialization-and-foundation)

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