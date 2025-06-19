# iDC-cortex-model
This is a morphological approximation model that investigates how ionic direct current (iDC) stimulation modulates the membrane potential in rat cortical neurons. 
The code generates:
1. **Extracellular potential $V_e$** produced by a disk electrode on the pia mater and the corresponding **mirror‐estimate $\Delta V_m$** for multiple neuron types **(Fig. 2A–B)**.
2. **Population‐level heatmaps** via weighted‐sum interpolation of soma potentials, and **layer‐wise averages** as a function of lateral distance **(Fig. 2C–E)**.

---

## Prerequisites

* **MATLAB** R2020a or later (no external toolboxes required).
* Download or clone this repository locally.
* Ensure the working directory contains the **`iDC_paper_model.m`** script.

---

## Script Overview

All code is contained in a single MATLAB script (iDC_paper_model.m) with two independent sections. You can run each section in isolation by selecting its code block and pressing **Run Section** (in the MATLAB Editor) or by copying it into the Command Window.

### Section 1: Ve & ΔVm for Neuron Types (Fig. 2A–B)

This section computes and plots:

* **Extracellular potential** (Ve) on a 2D cortical cross‐section (±1.1 mm laterally, 0–2.2 mm depth) using a point‐source summation from a 250 µm‐diameter disk electrode.
* **Mirror‐estimate** $\Delta V_m$ for 12 neuron morphologies (8 excitatory, 4 inhibitory) represented as vertical rods.  Each rod’s length covers its dendritic span and a black marker indicates the soma.

**To run**:

1. Open the script in MATLAB.
2. Highlight the code under **Section 1** (lines 1–108).
3. Press **Run Section**.
4. Examine figures titled “V\_e (V)” and each neuron’s $\Delta V_m$ subplot.

### Section 2: Weighted 2D Interpolation & Layer‐wise Averages (Fig. 2C–E)

This section simulates a population of rods randomly scattered according to known layer‐specific densities, computes each soma’s $\Delta V_m$, and generates:

* A 2D cross-section with randomly scattered neuron rods.
* A heatmap of **weighted‐sum soma potential change** of all scattered neurons on a finer grid (200×200) using a Gaussian kernel.
* **Layer‐wise average** potential changes at four lateral distances (0.2, 0.55, 0.9, 1.25 mm).

**To run**:

1. Highlight the code under **Section 2** (lines 109–end).
2. Press **Run Section**.
3. Examine the figures:
   * “Combined Randomly Scattered Neurons”
   * “2D Weighted Sum of Soma Potentials”
   * “Depth‑Layer Averages of Weighted Soma Potentials”

---

## Adjusting Parameters

* **iDC electrode current** (`Icenter_amp`): default anodic +20 µA.
* **iDC electrode radius** (`electrode_radius`): default 0.125mm. 
* **iDC electrode position** (`electrode_center`): default [0,0,0] mm.
* **Grid resolution** (`dx`, `dy`): finer resolution increases accuracy at cost of speed. 
* **Neuron types**: edit `neuronTypes` to reflect alternative neuron types/morphologies.
* **Neuron densities**: edit `densities_percent` to reflect alternative relative neuron densities.
* **Gaussian kernel** (`sigma`): controls spatial smoothing of weighted-sum interpolation.


---

## License

This project is released under the MIT License—see [LICENSE](LICENSE) for details.  Feel free to adapt and extend for your own research.


