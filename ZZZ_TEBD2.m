%build three site interaction local Hamiltonian
[S,I] = getLocalSpace('Spin',1/2);
Aprev = 1;
N=3;
delta=0.1;
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


L=16;

M = cell(1,L);
% S=cell(1,L+1);
% for it=1:L+1
%     S{it}=1;
% end

for itN = (1:L)
    if itN <= (L/2)
    M{itN} = permute([1 0],[1 3 2]);
    else
    M{itN} = permute([0 1],[1 3 2]);
    end
end
[M, Ovals, OvalsCon] = TEBD2(M,H,Sz,Nkeep,taus);

imagesc([1 L],[t(1) t(end)],real(Ovals));
colorbar;
clim([-0.5 0.5]); 
set(gca,'FontSize',13,'LineWidth',1);
xlabel('Site index');
ylabel('Time');
title(sprintf('S_z: L=%d, dt=%.3g, Nkeep=%d, \\Delta=%.3g, t_{max}=%.2f', ...
    L, dt, Nkeep, delta, tmax));

% OOvals = Ovals .* Ovals(:,6); % <S_z(l)*S_z(n)>
% 
% sharedMap = parula(256);
% Ovals_2p = OvalsCon - OOvals;
% % disp(max(imag(Ovals_2p(:)))); % should be double precision noise
% Ovals_2p = real(Ovals_2p); % remove noise
% figure;
% imagesc([1 L],[t(1) t(end)],real(Ovals_2p));
% 
% set(gca,'YDir','normal')
% colormap(sharedMap);
% colorbar;
% clim([-0.25 0.25]); 
% set(gca,'FontSize',13,'LineWidth',1);
% xlabel('Site index');
% ylabel('Time');
% title('Two point correlation function Sz');