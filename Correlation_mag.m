tic
n_l = L/2; % Center site
times = 0:0.1:30; % Time evolution steps

magOps = cell(1, L);
for m = 1:L
    magOps{m} = M' * zTerm{m} * M;
end
for j = 1:L
    magCorrOps{j} = magOps{n_l} * magOps{j};
end

correlator = zeros(L, length(times));
expPhase = exp(-1i * E * times);


for tInd = 1:numel(times)
    psi_t = expPhase(:, tInd) .* psi0Rotated;
    
    % Expectation values
    n_t = zeros(1, L);
    for m = 1:L
        n_t(m) = real(psi_t' * (magOps{m} * psi_t));
    end
    
    % Connected correlator
    for j = 1:L
        corr = real(psi_t' * (magCorrOps{j} * psi_t));
        correlator(j, tInd) = corr - n_t(n_l) * n_t(j);
    end
end

% Transpose for correct orientation in imagesc
correlator = correlator';  
toc
sharedMap = parula(256);
% Create the density plot
figure;
imagesc(1:L, times, correlator);
colormap(sharedMap);
colorbar;
xlabel('Lattice Site');
ylabel('Time');
title('Correlation function C_j^c (t)');
set(gca, 'YDir', 'normal'); % Ensures time runs from bottom to top
% colormap(jet); % Adjust colormap for better visualization
clim([-1, 1]);