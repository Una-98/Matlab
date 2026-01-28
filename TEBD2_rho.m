function [M, Ovals] = TEBD2_rho(M,H,O,Nkeep,taus)

% < Description >
%
% [M, Ovals] = TEBD(M,H,O,Nkeep,taus)
%
% The TEBD (time-evolving block decimation) method to find the
% ground state of an finite one-dimensional system, by applying imaginary
% time evolutions.
%
% < Input >
% M : [1 x N cell array] Matrix product state (MPS) as the initial guess of
%       the variational search. Each M{n} is a rank-3 tensor acting on site
%       n. Their legs are ordered as left-right-bottom(=physical). The
%       length M, i.e., numel(Hs), defines the system length N.
% H : [matrix] Three-site interaction Hamiltonian. 
%
% Nkeep : [integer] Maximum bond dimension.
% taus : [numeric] Vector of imaginary time step sizes. Each "outer"
%       iteration consists of two imaginary-time evolutions, the first for
%       odd bonds and the second for even bonds. Both time evolutions
%       within the m-th outer iteration take the times size taus(m).
%
% < Output >
% M : [cell] The final MPS after real-time evolution.
% Eiter : [(numel(taus) x 2 x 2 matrix] Eiter(m,n,k) is the measured energy
%       for an odd (k = 1) or even (k = 2) bond after odd (n = 1) or
%       even (n = 2) bonds are updated, at the m-th "outer" iteration.

% tobj = tic2;
% [~,I] = getLocalSpace('Spin',1/2);
N = numel(M);
Nstep = numel(taus);
ldim = 2; % local space dimension
Skeep = 1e-8;
tTEBD = 0.0;        % total time for time evolution
tObs  = 0.0;        % total time for observable calculation

if ~ismatrix(H)
    error('ERR: H should be rank-2.');
elseif any(ldim^3 ~= [size(H,1) size(H,2)])
    error('ERR: Both the legs of H should have the same size.');
end

tStart = tic;
lastPrint = tic;
printEvery = 1; 

% show information
disptime(['TEBD ground state search : Nkeep = ',sprintf('%i',Nkeep), ...
    ', # of imag. time steps = ',sprintf('%.4g',Nstep)]);

% energy expectation value at each step
% Eiter = zeros(Nstep,2,2);
Ovals = zeros(Nstep,N);
% diagonalize the Hamiltonian to exponentiate
[VH,DH] = eig((H+H')/2);
DH = diag(DH);

for it1 = (1:Nstep)
    % exponentiate the matrix representation of Hamiltonian
    expH = VH*diag(exp(-1i*taus(it1)/2*DH))*VH';
    % reshape matrix -> rank-6 tensor
    % expHconj=expH';
    expH = reshape(expH,ldim*ones(1,6));
    % expHconj=reshape(expHconj,ldim*ones(1,6));
    expHconj=conj(expH);
    % t0 = tic2;   % start TEBD timing
     t0_dt = datetime('now');  
   
    for it2 = (1:3:(N-2))
        M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep);
    end
    for it2 = (2:3:(N-2))
        M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep);
    end
    expH3 = VH*diag(exp(-1i*taus(it1)*DH))*VH';
    % expH3conj=expH3';

    % reshape matrix -> rank-6 tensor
    expH3 = reshape(expH3,ldim*ones(1,6));
    % expH3conj=reshape(expH3conj,ldim*ones(1,6));
    expH3conj=conj(expH3);

    for it2 = (3:3:(N-2))
        M = TEBD_3siteUpdate(M,expH3,it2,Nkeep,Skeep);
    end
    for it2 = (2:3:(N-2))
        M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep);
    end
    for it2 = (1:3:(N-2))
        M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep);
    end

    for it2 = (1:3:(N-2))
        M = TEBD_3siteUpdateConj(M,expHconj,it2,Nkeep,Skeep);
    end
    for it2 = (2:3:(N-2))
        M = TEBD_3siteUpdateConj(M,expHconj,it2,Nkeep,Skeep);
    end

    % expH3conj=conj(expH3);

    for it2 = (3:3:(N-2))
        M = TEBD_3siteUpdateConj(M,expH3conj,it2,Nkeep,Skeep);
    end
    for it2 = (2:3:(N-2))
        M = TEBD_3siteUpdateConj(M,expHconj,it2,Nkeep,Skeep);
    end
    for it2 = (1:3:(N-2))
        M = TEBD_3siteUpdateConj(M,expHconj,it2,Nkeep,Skeep);
    end
     dtTEBD = seconds(datetime('now') - t0_dt);
    % tTEBD = tTEBD + double(toc2(t0));
      tTEBD = tTEBD + double(dtTEBD);
     t1 = tic;   % start observable timing
    Tr = 1;

    for itN = 1:N
        Tr = contract(Tr,2,2,M{itN},4,1);
        Tr = contract(Tr,4,[4 3],getIdentity(Tr,4),2,[1 2]);
    end

    A0=1;
    for itN = (1:N)
        A0 = contract(A0,2,2,M{itN},4,1);
        O1=contract(A0,4,[4 3],O,2,[1 2]);
        for it2=itN+1:N
            O1=contract(O1,2,2,M{it2},4,1);
            O1=contract(O1,4,[4 3],getIdentity(O1,4),2,[1 2]);
        end
        Ovals(it1,itN)=trace(O1)/Tr;
        A0 = contract(A0,4,[4 3],getIdentity(A0,4),2,[1 2]);
    end
    % M{N}=M{N}/Tr;
    % for itN = (1:N)
    %     Tr = contract(Tr,2,2,M{itN},4,1);
    %     O1=contract(Tr,4,[4 3],O,2,[1 2]);
    %     for it2=itN+1:N
    %         O1=contract(O1,2,2,M{it2},4,1);
    %         O1=contract(O1,4,[4 3],getIdentity(O1,4),2,[1 2]);
    %     end
    %     Ovals(it1,itN)=trace(O1);
    %     Tr = contract(Tr,4,[4 3],getIdentity(Tr,4),2,[1 2]);
    % end
        tObs = tObs + toc(t1);
    if mod(it1, printEvery) == 0 || it1 == 1
        elapsed = toc(tStart);
        avgPerStep = elapsed / it1;
        remaining = avgPerStep * (Nstep - it1);

        fprintf(['it1 = %4d / %4d | ' ...
            'elapsed = %6.1f s | ' ...
            'ETA = %6.1f s\n'], ...
            it1, Nstep, elapsed, remaining);
    end

end
% check performance
fprintf('\n================ Timing summary ================\n');
fprintf('Total TEBD evolution time   : %8.2f s\n', tTEBD);
fprintf('Total observable time       : %8.2f s\n', tObs);
fprintf('================================================\n');
% toc2(tobj,'-v');
chkmem;
end

function M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep)

T = contract(M{it2},4,2,M{it2+1},4,1);% contract the first and second site of MPS 
% contract the third site of MPS to make a rank 5 tensor 
T = contract(T,6,4,M{it2+2},4,1,[1 6 2 4 7 3 5 8]);
% where legs 2, 3, 4 are to be contracted with the local time
% evolution operator

% Contract a three-site gate exp(-taus(it1)*H) with the rank-5 ket tensor

eHT = contract(T,8,[3 4 5],expH,6,[4 5 6]);
% eHT = contract(T,8,[6 7 8],expH,6,[1 2 3]);
        
% SVD; truncate singular values smaller than Skeep (= 1e-8 by
% default)
[M{it2},S,V] = svdTr(eHT,8,[1 3 6],Nkeep,Skeep);
S=diag(S);
M{it2} = permute(M{it2},[1 4 3 2]);
MM=contract(S,2,2,V,6,1);
[M{it2+1},S,V] = svdTr(MM,6,[1 3 5],Nkeep,Skeep);
S=diag(S);
M{it2+1} = permute(M{it2+1},[1 4 3 2]);
M{it2+2}=contract(S,2,2,V,4,1,[1 2 4 3]); 
end

function M = TEBD_3siteUpdateConj(M,expH,it2,Nkeep,Skeep)

T = contract(M{it2},4,2,M{it2+1},4,1);% contract the first and second site of MPS 
% contract the third site of MPS to make a rank 5 tensor 
T = contract(T,6,4,M{it2+2},4,1,[1 6 2 4 7 3 5 8]);
% where legs 2, 3, 4 are to be contracted with the local time
% evolution operator

% Contract a three-site gate exp(-taus(it1)*H) with the rank-5 ket tensor
% expH = permute(expH,[4 5 6 1 2 3]);
eHT = contract(T,8,[6 7 8],expH,6,[1 2 3]);
% eHT = contract(T,8,[3 4 5],expH,6,[4 5 6]);
        
% SVD; truncate singular values smaller than Skeep (= 1e-8 by
% default)
[M{it2},S,V] = svdTr(eHT,8,[1 3 6],Nkeep,Skeep);
S=diag(S);
M{it2} = permute(M{it2},[1 4 2 3]);
MM=contract(S,2,2,V,6,1);
[M{it2+1},S,V] = svdTr(MM,6,[1 3 5],Nkeep,Skeep);
S=diag(S);
M{it2+1} = permute(M{it2+1},[1 4 2 3]);
M{it2+2}=contract(S,2,2,V,4,1); 
end
