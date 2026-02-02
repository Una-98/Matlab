
%build three site interaction local Hamiltonian
[S,I] = getLocalSpace('Spin',1/2);
Aprev = 1;
N=3;
delta=1;

for itN = (1:N)
    Anow = getIdentity(Aprev,2,I,2,[1 3 2]);

    % Sz*(S*S) interaction
    if itN ==2 
        Hpm = updateLeft(Sprev,3,Anow, ...
            permute(conj(S(:,:,1)),[2 1 3]),3,Anow);
        Hmp = updateLeft(Sprev,3,Anow, ...
            permute(conj(S(:,:,3)),[2 1 3]),3,Anow);
        Hzz = updateLeft(Sprev,3,Anow, ...
            permute(conj(S(:,:,2)),[2 1 3]),3,Anow);
    end

    if itN ==3 
        Hpm = updateLeft(Hpm,3,Anow, ...
            permute(conj(S(:,:,3)),[2 1 3]),3,Anow);
        Hmp = updateLeft(Hmp,3,Anow, ...
            permute(conj(S(:,:,1)),[2 1 3]),3,Anow);
        Hzz = updateLeft(Hzz,3,Anow, ...
            permute(conj(S(:,:,2)),[2 1 3]),3,Anow);
    end
    
    Sprev = updateLeft([],[],Anow,S(:,:,2),3,Anow);
    Aprev = Anow;
end

H=Hpm+Hmp+delta*Hzz;

% iTEBD parameters
% Nkeep = [];
Nkeep = [];
Sz = S(:,:,2);

% tau_ini = 1; % initial imaginary time step size
% tau_fin = 1e-6; % final imaginary time step size
% Nstep = 2e3; % number of imaginary time steps
% taus = tau_ini*((tau_fin/tau_ini).^linspace(0,1,Nstep));

dt = 0.1;
tmax = 30;
taus = dt * ones(1, round(tmax/dt));
t = cumsum(taus);
Nstep=numel(taus);

mu=1;
L=8;

M = cell(1,L);
M1=0.5*I+mu*Sz;
M2=0.5*I-mu*Sz;

for itN = (1:L)
    if itN <= (L/2)
    M{itN} = permute(M1, [3 4 1 2]);
    else
    M{itN} = permute(M2, [3 4 1 2]);
    end
end

[M, Ovals, dw_sum] = TEBD2_rho_par(M,H,Sz,Nkeep,taus,dt);

% Create filename that embeds all simulation parameters
finalSaveName = sprintf('results_L%d_Nkeep%d_delta%.4g_dt%.4g_tmax%.4g_mu%.4g.mat', ...
    L, Nkeep, delta, dt, tmax, mu);

fprintf('Saving final results to: %s\n', finalSaveName);

try
    save(finalSaveName, 'M', 'Ovals', '-v7.3');
    fprintf('Final results saved successfully.\n');
catch ME
    warning('Failed to save final results: %s', ME.message);
end

imagesc([1 L],[t(1) t(end)],real(Ovals));
colorbar;
clim([-0.5 0.5]); 
set(gca,'FontSize',13,'LineWidth',1);
xlabel('Site index');
ylabel('Time');title(sprintf('S_z: L=%d, dt=%.3g, Nkeep=%d, \\mu=%.3g, \\Delta=%.3g, t_{max}=%.2f', ...
    L, dt, Nkeep, mu, delta, tmax));