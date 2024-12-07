% Parameters for the Lattice
params = struct();
params.N_x = 20; % Range of m values
params.N_z = 14; % Range of j values
params.Nk = 301; % Number of k points in the Brillouin zone
params.k_values = linspace(0, 2 * pi, params.Nk); % Brillouin zone discretization
params.num_bands = 15; % Number of bands
params.L = 5; % Lattice parameter

% Parameters for the Potentials
params.V0 = -350.0; params.eta = 2; params.epsilon = 0.02; % Potential 1
params.V0_prime = -380.0; params.epsilon_prime = 0.03; params.b = 0.4; params.eta_prime = 2; % Potential 2
params.V0_pri = -400.0; params.c = 0.6; params.epsilon_pri = 0.03; params.eta_pri = 2; % Potential 3

% Constants for Laplacian and Potentials
params.c1 = 4 * pi^2 / params.L^2;
params = calculate_potential_constants(params);

% Meshgrid for Lattice Indices
[J, M, Index] = create_lattice_indices(params.N_x, params.N_z);

% Preallocate for Eigenvalues and Eigenvectors
eigenvalue_matrix = zeros(params.Nk, params.num_bands);
eigenvector_matrix = zeros(params.Nk, numel(Index), params.num_bands);

% Calculate Eigenvalues and Eigenvectors
[eigenvalue_matrix, eigenvector_matrix] = compute_eigenvalues(params, J, M, Index);

% Plot Eigenvector for Band
plot_eigenvector(eigenvector_matrix, params, 1, -100:0.01:100);

% Define and Plot Combined Potential
plot_combined_potential(params);

% Plot Band Structure
plot_band_structure(params.k_values, eigenvalue_matrix, params.num_bands);

%% Functions

function params = calculate_potential_constants(params)
    % Calculate constants for potentials
    params.c2 = sqrt(pi / params.eta) * params.V0 * params.epsilon / params.L;
    params.c3 = pi^2 / params.eta / params.L^2;
    params.c4 = 2 * pi^2 * params.epsilon^2;
    params.c5 = sqrt(pi / params.eta_prime) * params.V0_prime * params.epsilon_prime / params.L;
    params.c6 = pi^2 / params.eta_prime / params.L^2;
    params.c7 = 2 * pi^2 * params.epsilon_prime^2;
    params.c8 = sqrt(pi / params.eta_pri) * params.V0_pri * params.epsilon_pri / params.L;
    params.c9 = pi^2 / params.eta_pri / params.L^2;
    params.c10 = 2 * pi^2 * params.epsilon_pri^2;
end

function [J, M, Index] = create_lattice_indices(N_x, N_z)
    % Create lattice indices for m and j
    [J, M] = meshgrid(-N_z:N_z, -N_x:N_x);
    Index = (J + N_z) * (2 * N_x + 1) + (M + N_x) + 1;
    J = J(:); M = M(:); Index = Index(:);
end

function [eigenvalue_matrix, eigenvector_matrix] = compute_eigenvalues(params, J, M, Index)
    % Compute eigenvalues and eigenvectors for the Hamiltonian
    
    N_index = numel(Index);
    laplacian_term_fn = @(j_prime, m_prime, k) ...
        (k - 2 * pi * m_prime).^2 + params.c1 * j_prime.^2;
    potential_term_fn = @(j, j_prime, m, m_prime) ...
        params.c2 * exp(-params.c3 * (j - j_prime).^2) .* exp(-params.c4 * (m - m_prime).^2) + ...
        params.c5 * exp(-params.c6 * (j - j_prime).^2) .* exp(-params.c7 * (m - m_prime).^2) .* exp(1i * 2 * pi * (m - m_prime) * params.b) + ...
        params.c8 * exp(-params.c9 * (j - j_prime).^2) .* exp(-params.c10 * (m - m_prime).^2) .* exp(1i * 2 * pi * (m - m_prime) * params.c);
    
    parfor i = 1:params.Nk
        fprintf('Processing k-point %d / %d (%.2f%% complete)\n', i, params.Nk, (i / params.Nk) * 100);
        H = zeros(N_index, N_index);
        K = diag(laplacian_term_fn(J, M, params.k_values(i)));
        for n = 1:N_index
            H(Index(n), Index) = potential_term_fn(J(n), J, M(n), M);
        end
        H = H + K;
        [eigvecs, eigvals] = eig(H);
        [eigvals, I] = sort(real(diag(eigvals)));
        eigvecs = eigvecs(:, I);
        eigenvalue_matrix(i, :) = eigvals(1:params.num_bands);
        eigenvector_matrix(i, :, :) = eigvecs(:, 1:params.num_bands);
    end
end

function plot_eigenvector(eigenvector_matrix, params, band, XI)
    % Plot the eigenvector for a selected band
    psi = zeros(1, numel(XI));
    for l = 1:numel(XI)
        xi = XI(l);
        k = mod(xi, 2 * pi);
        k_int = min(floor(params.Nk * k / (2 * pi) + 1), params.Nk);
        m = round((xi - k) / (2 * pi));
        psi(l) = eigenvector_matrix(k_int, m + params.N_x + 1, band);
    end
    figure;
    plot(XI, real(psi), 'b');
    xlabel('\xi');
    ylabel('\psi(\xi)');
    title(['Eigenvector for Band ', num2str(band)]);
end

function plot_combined_potential(params)
    % Define and plot the combined potential
    x = linspace(0, 1, 100);
    z = linspace(-params.L / 2, params.L / 2, 100);
    [X, Z] = meshgrid(x, z);
    
    V_one = @(x, z) params.V0/sqrt(2*pi) * ...
        (exp(-(x - (-1)).^2 / (2*params.epsilon^2)) + ...
         exp(-(x - 0).^2 / (2*params.epsilon^2)) + ...
         exp(-(x - 1).^2 / (2*params.epsilon^2))) .* ...
        exp(-(z + params.L * 0).^2 * params.eta);
    
    V_two = @(x, z) params.V0_prime/sqrt(2*pi) * ...
        (exp(-(x + params.b - (-1)).^2 / (2*params.epsilon_prime^2)) + ...
         exp(-(x + params.b - 0).^2 / (2*params.epsilon_prime^2)) + ...
         exp(-(x + params.b - 1).^2 / (2*params.epsilon_prime^2))) .* ...
        exp(-(z + params.L * 0).^2 * params.eta_prime);
    
    V_three = @(x, z) params.V0_pri/sqrt(2*pi) * ...
        (exp(-(x + params.c - (-1)).^2 / (2*params.epsilon_pri^2)) + ...
         exp(-(x + params.c - 0).^2 / (2*params.epsilon_pri^2)) + ...
         exp(-(x + params.c - 1).^2 / (2*params.epsilon_pri^2))) .* ...
        exp(-(z + params.L * 0).^2 * params.eta_pri);
    
    V_tilde = @(x, z) V_one(x, z) + V_two(x, z) + V_three(x, z);
    V_grid_multi = V_tilde(X, Z);
    
    figure;
    surf(X, Z, V_grid_multi);
    xlabel('x');
    ylabel('z');
    zlabel('Potential V(x,z)');
    title('Combined Potential V(x,z)');
    colorbar;
end

function plot_band_structure(k_values, eigenvalue_matrix, num_bands)
    % Plot the band structure for selected bands
    figure;
    hold on;
    for n = 1:num_bands
        plot(k_values, eigenvalue_matrix(:, n), 'LineWidth', 1.5);
    end
    hold off;
    xlabel('k');
    ylabel('Eigenvalue');
    title('Band Structure');
    grid on;
    legend(arrayfun(@(n) ['n = ', num2str(n)], 1:num_bands, 'UniformOutput', false));
end
