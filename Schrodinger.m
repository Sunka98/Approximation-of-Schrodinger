% Define the parameters
N_x = 20; % Range of m values from -N_x to N_x
N_z = 14; % Range of j values from -N_z to N_z
Nk = 301; % Number of k points in the Brillouin zone
k_values = linspace(0, 2 * pi, Nk); % Brillouin zone discretization
%eps = pi/4;
num_bands = 15; % Number of bands
L = 5; % Set the value of L, Lattice 
%parameters of potential1
eta = 2; % Set the value of eta, -rate at which it decays 
V0 = -350.0; % Set the value of V0, depth of the potential
epsilon = 0.02 ;% Set the value of epsilon, width 
%parameters of potential2
V0_prime = -380.0; % Set the value of V0'
epsilon_prime = 0.03; % Set the value of epsilon'
b = 0.4; %shift invariance in potential
eta_prime = 2;
%parameters of potential3
V0_pri=-400;
c=0.6;
epsilon_pri=0.03;
eta_pri= 2;



% Define constants for potential / laplacian operators
c1 = 4*pi^2/L^2;
%potential_one
c2 =  sqrt(pi / eta) * V0 * epsilon / L;
c3 = (pi^2) / eta / L^2;
c4 = 2 * pi^2 * epsilon^2;
%potential_two
c5 =  sqrt(pi / eta_prime) * V0 * epsilon_prime / L;
c6 = (pi^2) / eta_prime / L^2;
c7 = 2 * pi^2 * epsilon_prime^2;
%potential three
c8 =  sqrt(pi / eta_pri) * V0 * epsilon_pri / L;
c9 = (pi^2) / eta_pri / L^2;
c10 = 2 * pi^2 * epsilon_pri^2;

% Define the Laplacian term function
laplacian_term_fn = @(j_prime, m_prime, k) ((k - 2 * pi * m_prime).^2 + c1 * j_prime.^2);

% Define the potential term function
potential_term_fn = @(j, j_prime, m, m_prime) ...
    c2 * exp(-c3 * (j - j_prime).^2) .* exp(-c4 * (m - m_prime).^2) + ...
    c5 * exp(-c6 * (j - j_prime).^2) .* exp(-c7 * (m - m_prime).^2) .* exp(1i*2 * pi * (m - m_prime) * b)+...
    c8 * exp(-c9 * (j - j_prime).^2) .* exp(-c10 * (m - m_prime).^2) .* exp(1i*2 * pi * (m - m_prime) * c);

[J,M] = meshgrid(-N_z:N_z,-N_x:N_x);
Index = (J+N_z)*(2*N_x+1) + (M+N_x)+1;
J = J(:);
M = M(:);
Index = Index(:);
N_index = size(Index,1);
eigenvalue_matrix = zeros(Nk, num_bands);
eigenvector_matrix = zeros(Nk, N_index, num_bands);

% Calculate eigenvalues and eigenvectors
parfor i = 1:Nk
    H = zeros((2 * N_x + 1) * (2 * N_z + 1), (2 * N_x + 1) * (2 * N_z + 1));
    fprintf('%d / %d\n', i, Nk)
    K = diag(laplacian_term_fn(J, M, k_values(i)));
    for n = 1:N_index
        H(Index(n), Index) = potential_term_fn(J(n), J, M(n), M);
    end
    H = H + K;
    [eigenvectors, eigenvalues] = eig(H); 
    [eigenvalues, I] = sort(real(diag(eigenvalues)));
    eigenvectors = eigenvectors(:, I);
    eigenvalue_matrix(i, :) = eigenvalues(1:num_bands);
    eigenvector_matrix(i, :, :) = eigenvectors(:, 1:num_bands);
end

band = 1;
j=1;
XI = -100:.01:100;
psi = zeros(1,size(XI,2));
for l = 1:size(XI,2)
xi = XI(l);
k = mod(xi,2*pi);
k_int = min(floor(Nk*k/(2*pi)+1),Nk);
m = round((xi - k)/(2*pi));
psi(l) = eigenvector_matrix(k_int,(j+N_z)*(2*N_x+1) + (m+N_x)+1,band);
end
plot(XI,real(psi),'b')
hold on
%plot(XI, imag(psi),'r')
hold off
%
% Define potential functions V_one and V_two and V_three

V_one = @(x, z) V0/sqrt(2*pi) * (exp(-(x - (-1)).^2 / (2*epsilon^2)) + ...
                                   exp(-(x - 0).^2 / (2*epsilon^2)) + ...
                                   exp(-(x - 1).^2 / (2*epsilon^2)) ).* exp(-(z + L*0).^2 * eta);
V_two = @(x, z) V0_prime/sqrt(2*pi) * (exp(-(x+b - (-1)).^2 / (2*epsilon_prime^2)) + ...
                                   exp(-(x+b - 0).^2 / (2*epsilon_prime^2)) + ...
                                   exp(-(x+b - 1).^2 / (2*epsilon_prime^2))) .* exp(-(z + L*0).^2 * eta_prime);
V_three = @(x, z) V0_pri/sqrt(2*pi) * (exp(-(x+c - (-1)).^2 / (2*epsilon_pri^2)) + ...
                                   exp(-(x+c - 0).^2 / (2*epsilon_pri^2)) + ...
                                   exp(-(x+c - 1).^2 / (2*epsilon_pri^2))) .* exp(-(z + L*0).^2 * eta_pri);

% Define combined potential V_tilde
V_tilde = @(x, z) V_one(x, z) + V_two(x, z)+ V_three(x, z);


% Define grid for plotting
x = linspace(0, 1, 100); 
z = linspace(-L/2, L/2, 100); 
[X, Z] = meshgrid(x, z);
V_grid_multi = V_tilde(X, Z);

% Plot the combined potential
figure;
surf(X, Z, V_grid_multi);
xlabel('x');
ylabel('z');
zlabel('Potential V(x,z)');
title('Combined Potential V(x,z)');
colorbar; % Show colorbar for the potential values
% Save the plot at 300 DPI
%filename = 'Potential_4.png'; % Specify the filename and format
%print(filename, '-dpng', '-r300'); % Save as PNG with 300 DPI

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
legend('n = 1', 'n = 2', 'n = 3', 'n = 4', 'n = 5', 'n = 6', 'n = 7');
% Save the plot at 300 DPI
%filename = 'band_structure.png'; % Specify the filename and format
%print(filename, '-dpng', '-r300'); % Save as PNG with 300 DPI

% Initialize the matrix to store u_j(k)
%%
u_matrix = zeros(size(eigenvector_matrix));

% Loop over each k-point

for j = 1:num_bands  
    [U,U_I] = sort(abs(eigenvector_matrix(1,:,j)));
    index = U_I(size(U,2));
    for i = 1:Nk
    % Loop over each band
        % Get the eigenvector for the current k-point and band
        v_jk = eigenvector_matrix(i, :, j);
        

        % Compute the real and imaginary parts of v_j(k)
        alpha = real(v_jk(index));
        beta = imag(v_jk(index));

        
        % Compute e^{-i\theta}
        theta = atan2(beta, alpha); % Compute the angle
        exp_neg_i_theta = (alpha - 1i * beta) / sqrt(alpha^2 + beta^2);
        
        % Compute u_j(k)
        u_jk = exp_neg_i_theta * v_jk;
        
        % Store the result in the u_matrix
        u_matrix(i, :, j) = u_jk;
    end
end

% Define x values for psi_j(x)

x_values = linspace(-20, 20, 500); 
psi_jx = zeros(1, length(x_values), band);
band = 1; 

% Loop over x values to compute psi_j(x)
for xi = 1:length(x_values)
    x = x_values(xi);
    integral_value = 0;
    
    % Loop over k points
    for k_idx = 1:Nk
        k = k_values(k_idx);
        
        % Sum over lattice points m
        for m = -N_x:N_x
            % Find the corresponding index in J and M
            idx = find(J == 0 & M == m);
            
            % Extract  eigenvector component from u_matrix
            c_jm = u_matrix(k_idx, idx, band);
            
            % Compute the integral value summation term
            integral_value = integral_value + c_jm * exp(1i * (k - 2 * pi * m) * x);
        end
    end
    
    % Normalize the integral value
    psi_jx(xi) = integral_value * (2 * pi / Nk); % Normalizing factor
end

% % Plot the abs part of psi_j(x)
% figure;
% plot(x_values, abs(psi_jx), 'b');
% xlabel('x');
% ylabel('abs(\psi_j(x))');
% title('Absolute value of \psi_j(x)');
% grid on;


% Define the fixed j and m values
j = 1;
m = 1;
band=1;

% Find the index corresponding to j = 1 and m = 1 in the J and M arrays
idx_jm = find(J == j & M == m);

% Initialize the array to store c_{jm} values for each k
c_jm_values = zeros(1, Nk);

% Loop over k points to extract c_{jm}
for k_idx = 1:Nk
    % Extract the corrected eigenvector component from u_matrix
    c_jm_values(k_idx) = u_matrix(k_idx, idx_jm, band);
end

% % Plot the real part, imaginary part, and absolute value of c_{jm} against k-values
% figure;
% 
% plot(k_values, abs(c_jm_values), 'g', 'DisplayName', 'Absolute value');
% xlabel('k');
% ylabel('c_{jm}');
% title(['c_{jm} for j = ' num2str(j) ' and m = ' num2str(m)]);
% legend('show');
% grid on;
% hold off;


% 
% 
% 
% 
% 
% Define z values and compute w_n(x, z)
x_values = linspace(-40, 40, 500);
z_values = linspace(-5, 5, 500); % input as needed
w_nxz = zeros(length(x_values), length(z_values), band); % Storing w_n(x, z) values

% Loop over x and z values to compute w_n(x, z)
for xi = 1:length(x_values)
    x = x_values(xi);

    for zi = 1:length(z_values)
        z = z_values(zi);

        % Sum over j to compute w_n(x, z)
        w_nxz(xi, zi) = sum(exp(-2 * pi * 1i * J * z / L) .* psi_jx(:, xi)) / sqrt(L);
    end
end


% Define finer grids for interpolation
x_values_fine = linspace(-40,40, 1000);
z_values_fine = linspace(-L/2, L/2, 1000);
[X_fine, Z_fine] = meshgrid(x_values_fine, z_values_fine);

% Pad w_nxz with zeros to minimize boundary effects
padding = 50; % Amount of zero padding
w_nxz_padded = padarray(abs(w_nxz), [padding, padding], 0);

% Adjust the x and z values to match the padded array dimensions
x_values_padded = linspace(min(x_values), max(x_values), length(x_values) + 2*padding);
z_values_padded = linspace(min(z_values), max(z_values), length(z_values) + 2*padding);

% Define finer grids for interpolation
x_values_fine = linspace(min(x_values_padded), max(x_values_padded), 1000);
z_values_fine = linspace(min(z_values_padded), max(z_values_padded), 1000);
[X_fine, Z_fine] = meshgrid(x_values_fine, z_values_fine);

% Interpolate w_n(x, z) using spline interpolation
w_nxz_interp = interp2(z_values_padded, x_values_padded, w_nxz_padded, Z_fine, X_fine, 'spline');

% Plot the interpolated absolute value of w_n(x, z)
figure;
imagesc(z_values_fine, x_values_fine, w_nxz_interp);
colorbar;
xlabel('z');
ylabel('x');
title('Surface plot of Absolute value of w_n(x, z)');
grid on;
% Optional save the plot at 300 DPI
filename = 'ABS_wan.png'; 
print(filename, '-dpng', '-r300'); 


% Find the index for z = 0
[~, z_index] = min(abs(z_values_fine)); % Find the index closest to z = 0

% Extract the slice at z = 0
w_nx_z0 = w_nxz_interp(:, z_index);

% Plot log|w_n(x, z = 0)|
figure;
plot(x_values_fine, log(abs(w_nx_z0)), 'LineWidth', 1.5);
xlabel('x');
ylabel('log|w_n(x, z = 0)|');
title('log|w_n(x, z = 0)|');
grid on;

% % Optional save the plot at 300 DPI
% filename = 'log_w_n_x_z0.png'; 
% print(filename, '-dpng', '-r300'); 