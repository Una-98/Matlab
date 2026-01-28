function [M, Ovals,OvalsCon] = TEBD2(M,H,O,Nkeep,taus)

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

tobj = tic2;
[~,I] = getLocalSpace('Spin',1/2);
N = numel(M);
Nstep = numel(taus);
ldim = 2; % local space dimension
Skeep = 1e-8;
OvalsCon = zeros(Nstep,N);

if ~ismatrix(H)
    error('ERR: H should be rank-2.');
elseif any(ldim^3 ~= [size(H,1) size(H,2)])
    error('ERR: Both the legs of H should have the same size.');
end

tStart = tic;
lastPrint = tic;
printEvery = 1; 
tTEBD = 0.0;        % total time for time evolution
tObs  = 0.0; 
% show information
disptime(['TEBD ground state search : Nkeep = ',sprintf('%i',Nkeep), ...
    ', # of imag. time steps = ',sprintf('%.4g',Nstep)]);

Ovals = zeros(Nstep,N);
% diagonalize the Hamiltonian to exponentiate
[VH,DH] = eig((H+H')/2);
DH = diag(DH);
for it1 = (1:Nstep)
    % exponentiate the matrix representation of Hamiltonian
    expH = VH*diag(exp(-1i*taus(it1)/2*DH))*VH';
    % reshape matrix -> rank-6 tensor
    expH = reshape(expH,ldim*ones(1,6));
     t0_dt = datetime('now');
        for it2 = (1:3:(N-2))
            M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep);    
        end
        for it2 = (2:3:(N-2))
            M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep); 
        end
    expH3 = VH*diag(exp(-1i*taus(it1)*DH))*VH';
    % reshape matrix -> rank-6 tensor
    expH3 = reshape(expH3,ldim*ones(1,6));
        for it2 = (3:3:(N-2))
            M = TEBD_3siteUpdate(M,expH3,it2,Nkeep,Skeep);
        end
        for it2 = (2:3:(N-2))
            M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep); 
        end
        for it2 = (1:3:(N-2))
            M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep);    
        end
      dtTEBD = seconds(datetime('now') - t0_dt);
      tTEBD = tTEBD + double(dtTEBD);
    [M,~,~] = canonForm (M,numel(M),[],0);

    % [Mp,~,~] = canonForm (M,0,[],0);

    %permute legs to use update left
    % for itN = N:-1:1
    % Mp{itN} = permute(M{itN},[2 1 3]);
    % end
     t1 = tic; 
    Mr = [];

    for itN = (N:-1:1)
        Ovals(it1,itN) = trace(updateRight(Mr,2,M{itN},O,2,M{itN})); % the left part is the identity
        Mr = updateRight(Mr,2,M{itN},[],[],M{itN});
    end
      tObs = tObs + toc(t1);
    % for itN = (N:-1:1)
    %     B = permute(M{itN},[2 1 3]);
    %     Ovals(it1,itN) = trace(updateLeft(Mr,2,B,O,2,B)); % the left part is the identity
    %     Mr = updateLeft(Mr,2,B,[],[],B);
    % end

    % CC = [];
    % C=Mp;
    % ids=N/2;
    % 
    % 
    % for itN = (1:N)
    % 
    %     if itN<ids
    %         CCr=updateLeft(CC,2,C{itN},O,2,C{itN});
    % 
    %         for itC = (itN+1):1:(ids-1)
    %             CCr=updateLeft(CCr,2,C{itC},[],[],C{itC});
    %         end
    %         CCr = updateLeft(CCr,2,C{ids},O,2,C{ids});
    % 
    %         OvalsCon(it1,itN)=trace(CCr);
    %         CC=updateLeft(CC,2,C{itN},[],[],C{itN});
    %     end
    % 
    %     CCl = [];
    % 
    %     if itN>ids
    %         for itC=1:(ids-1)
    %             CCl=updateLeft(CCl,2,C{itC},[],[],C{itC});
    %         end
    %         CCl = updateLeft(CCl,2,C{ids},O,2,C{ids});
    % 
    %         if itN == (ids+1)
    %             CCl=updateLeft(CCl,2,C{itN},O,2,C{itN});
    % 
    %         else
    %             for itC = (ids+1):1:itN-1
    %                 CCl=updateLeft(CCl,2,C{itC},[],[],C{itC});
    %             end
    %             CCl=updateLeft(CCl,2,C{itN},O,2,C{itN});
    % 
    %         end
    %         OvalsCon(it1,itN)=trace(CCl);
    % 
    %     elseif itN==ids
    %         for itC=1:(ids-1)
    %             CCl=updateLeft(CCl,2,C{itC},[],[],C{itC});
    %         end
    % 
    %         OvalsCon(it1,itN)=1/4;
    %     end

    % end
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
toc2(tobj,'-v');
chkmem;
end

function M = TEBD_3siteUpdate(M,expH,it2,Nkeep,Skeep)

T = contract(M{it2},3,2,M{it2+1},3,1,[1 3 2 4]);% contract the first and second site of MPS 
% contract the third site of MPS to make a rank 5 tensor 
T = contract(T,4,2,M{it2+2},3,1,[1 4 2 3 5]);
% where legs 2, 3, 4 are to be contracted with the local time
% evolution operator
        
% Contract a three-site gate exp(-taus(it1)*H) with the rank-5 ket tensor
%eHT = contract(expH,6,[4 5 6],T,5,[2 3 5]);
eHT = contract(T,5,[3 4 5],expH,6,[4 5 6]);
        
% SVD; truncate singular values smaller than Skeep (= 1e-8 by
% default)
[M{it2},S,V] = svdTr(eHT,5,[1 3],Nkeep,Skeep);
S=diag(S);
M{it2} = permute(M{it2},[1 3 2]);
MM=contract(S,2,2,V,4,1);
[M{it2+1},S,V] = svdTr(MM,4,[1 3],Nkeep,Skeep);
S=diag(S);
M{it2+1} = permute(M{it2+1},[1 3 2]);
M{it2+2}=contract(S,2,2,V,3,1); 
end
