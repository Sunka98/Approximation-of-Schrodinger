
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
1. Kohn, W. "Analytic Properties of Bloch Waves and Wannier Functions." Physical Review, vol. 115, no. 4, 1959, pp. 809–821. DOI: 10.1103/PhysRev.115.809.

2. Damle, Anil, Antoine Levitt, and Lin Lin. "Variational Formulation for Wannier Functions with Entangled Band Structure." 2018. arXiv, arXiv:1801.08572.

3. Marzari, Nicola, et al. "Maximally Localized Wannier Functions: Theory and Applications." Reviews of Modern Physics, vol. 84, no. 4, 2012, pp. 1419–1475. DOI: 10.1103/revmodphys.84.1419.

4. Kuchment, Peter. "Tight Frames of Exponentially Decaying Wannier Functions." Journal of Physics A: Mathematical and Theoretical, vol. 42, no. 2, 2008, p. 025203. DOI: 10.1088/1751-8113/42/2/025203.

5. Fang, Shiang, and Efthimios Kaxiras. "Electronic Structure Theory of Weakly Interacting Bilayers." Physical Review B, vol. 93, no. 23, 2016, p. 235153. DOI: 10.1103/PhysRevB.93.235153.

6. He, Lixin, and David Vanderbilt. "Exponential Decay Properties of Wannier Functions and Related Quantities." Physical Review Letters, vol. 86, no. 23, 2001, pp. 5341–5344. DOI: 10.1103/PhysRevLett.86.5341.

7. Massatt, Daniel, Mitchell Luskin, and Christoph Ortner. "Electronic Density of States for Incommensurate Layers." Multiscale Modeling & Simulation, vol. 15, no. 1, 2017, pp. 476–499. DOI: 10.1137/16M1088363.

8. Wang, Ting, et al. "Convergence of the Planewave Approximations for Quantum Incommensurate Systems." 2022. arXiv, arXiv:2204.00994.

9. Fefferman, Charles L., and Michael I. Weinstein. "Honeycomb Lattice Potentials and Dirac Points." Journal of the American Mathematical Society, vol. 25, no. 4, 2012, pp. 1169–1220. DOI: 10.1090/S0894-0347-2012-00745-0.



