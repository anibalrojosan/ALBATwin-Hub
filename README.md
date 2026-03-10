# ALBATwin-Hub 🌿🦠

> This project is a work in progress. Currently in the **Phase 1: Mechanistic Engine (ALBATwin)**.

**ALBATwin-Hub** is a high-fidelity Digital Twin for microalgae-bacteria consortia in High Rate Algal Ponds (HRAP). It implements the **ALBA model** (Casagli et al., 2021) using a hybrid approach that combines rigorous mechanistic equations with Physics-Informed Neural Networks (PINNs).

## Index

- [✨ Key Features](#key-features)
- [📂 Project Structure](#project-structure)
- [🚀 Roadmap](#roadmap)
- [📚 References](#references)

## Key Features
- **Mechanistic Engine:** 19 coupled ODEs/DAEs simulating biological and chemical dynamics (pH, Gas Transfer, Light).
- **Hybrid Modeling:** Integration of physical laws into Deep Learning via PINNs for robust predictions under noise.
- **Modern Stack:** Managed with `uv` for lightning-fast dependency handling and reproducibility.
- **Scenario Analysis:** Interactive simulation of "What-if" scenarios (weather impact, aeration control).

For a more detailed description of the system, please refer to the [ARCHITECTURE](docs/ARCHITECTURE.md) and [MATH_MODEL](docs/MATH_MODEL.md) pages.

## Project Structure
```
ALBATwin-Hub/
├── docs/
│   ├── development/
│   │   └── DEVLOG.md       # 📝 Chronological log of technical decisions and daily progress
│   ├── ARCHITECTURE.md     # 🏗️ System design, component interaction, and data flow diagrams
│   ├── MATH_MODEL.md       # ➗ Detailed ODE/DAE equations and biological kinetics (ALBA model)
│   ├── PRD.md              # 📋 Core objectives, user stories, and functional requirements
│   ├── REFERENCES.md       # 📚 Scientific bibliography and data sources
│   └── ROADMAP.md          # 🗺️ High-level project phases and long-term goals
├── .gitignore               
├── pyproject.toml          # ⚙️ Project metadata and dependency management
├── uv.lock                 # 🔒 Locked dependencies
└── README.md               # 📖 Project overview
```

## Roadmap
To see the development phases, please refer to the [ROADMAP](docs/ROADMAP.md).

## References
The main reference paper for the MVP of the project is:
- Casagli et al. (2021). **ALBA: A comprehensive growth model to optimize algae-bacteria wastewater treatment in raceway ponds**. *Water Research*, 190, 116734. [https://doi.org/10.1016/j.watres.2020.116734](https://doi.org/10.1016/j.watres.2020.116734) 

For the complete list of references, please refer to the [REFERENCES](docs/REFERENCES.md) page.

---
> Made with ❤️ by [Anibal Rojo](https://github.com/anibalrojosan)