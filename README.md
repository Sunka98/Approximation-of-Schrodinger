
# Approximation of the Schrödinger Equation of 2D Materials
This repository provides the foundation for solving the Schrödinger equation of 2D materials. 
The approach is based on Bloch theory, with the potential term modeled using a generalized Kronig-Penney model.


The Hamiltonian of the Schrodinger Equation Hψ = Eψ is H = -∂²/∂x² - ∂²/∂z² + V(x, z),

where the potential based on Kronig Penny Model is given by:

V(x, z) = Σᵢ Σₗₙ V₀ⁱ δₑ(x + bᵢ - nᵢ) exp(-(z + Llᵢ)²ηᵢ).
The Gaussian function is defined as δₑ(x) = 1 / √(2π) e⁻(x² / 2ε²).

The matrix equation is given by:

ℋⱼ'ₘ',ⱼₘ cⱼₘ = Eₙₖ cⱼ'ₘ'

where ℋⱼ'ₘ',ⱼₘ represents the matrix elements of the Hamiltonian. Solving this matrix equation provides the eigenvalues Eₙₖ (the energy) and the eigenvectors cⱼₘ.

The matrix elements of the Hamiltonian are expressed as:

ℋⱼ'ₘ',ⱼₘ = [(k - 2πm')² + (4π²j'²) / L²] δⱼⱼ' δₘₘ'
          + (√(π / ηᵢ))(V₀ⁱ εᵢ / L) e^(-π²(|j - j'|)² / (ηᵢL²)) 
            e^(-2π²(m - m')² εᵢ²) e^(2πi(m - m')bᵢ)

By taking the truncations as:

|m| ≤ Nₓ, |j| ≤ N𝓏, Nₖ, k ∈ [0, 2π]

| ![Band Structure for [0,0.4,0.6]](Band_structure[0,0.4,0.6].png) | ![Band Structure for [0,0.33,0.66]](Bandstructure_[0,0.33,0.66].png) |
|:-------------------------------------------------------:|:--------------------------------------------------------:|
| Band Structure  [0,0.4,0.6]                                    | Band Structure  [0,0.33,0.66]                                   |






# Observations:
If the centers of the Gaussian functions are taken as [0, 0.4, 0.6], there are three isolated bands below the Fermi level, each separated by a band gap. However, we obtain one isolated and two entangled bands if the centers are taken at one-third ratios, i.e., [0, 0.33, 0.66]. Since all these bands are truncated within a certain value, we observe band folding in the given range.

# Acknowledgements
This research was done under Prof. Daniel Massatt- Louisiana State University.
# References
Electronic Structure Theory of Weakly Interacting Bilayers
Authors: Shiang Fang, Efthimios Kaxiras
Journal: Phys. Rev. B
Volume: 93, Issue 23, Pages: 235153 (2016).
DOI:10.1103/PhysRevB.93.235153

Exponential Decay Properties of Wannier Functions and Related Quantities
Authors: Lixin He, David Vanderbilt
Journal: Phys. Rev. Lett.
Volume: 86, Issue 23, Pages: 5341–5344 (2001).
DOI:10.1103/PhysRevLett.86.5341

Electronic Density of States for Incommensurate Layers
Authors: Daniel Massatt, Mitchell Luskin, Christoph Ortner
Journal: Multiscale Modeling & Simulation
Volume: 15, Issue 1, Pages: 476–499 (2017).
DOI:10.1137/16M1088363

Convergence of the Planewave Approximations for Quantum Incommensurate Systems
Authors: Ting Wang et al.
Year: 2022
arXiv:2204.00994

Honeycomb Lattice Potentials and Dirac Points
Authors: Charles L. Fefferman, Michael I. Weinstein
Journal: Journal of the American Mathematical Society
Volume: 25, Issue 4, Pages: 1169–1220 (2012).
DOI:10.1090/S0894-0347-2012-00745-0
