function [H0, states,occupation] = ZZZ(L,N,matterConfig,PBC)

% myTimes = [0 1e-2:1e-2:1];
% myTimes(end) = [];
% grid = 10;
% lastPower = 12;

% for cnt = 0:(lastPower-1)
%     ti = 10.0^cnt;
%     tf = ti*10.0;
%     myTimes = [myTimes ti:ti/grid:tf];
% 
%     if cnt < (lastPower-1)
%         myTimes(end) = [];
%     end
% end
myTimes=[1/20:1/20:30];

states=BasisCreatorZZZ(L, N);
HilbertSpaceSize = size(states,1);
FlipFlop = cell(1,L);
EnergyShift = cell(1,L);
Spin = cell(1,L);
zTerm = cell(1,L);


I = sparse(1:HilbertSpaceSize,1:HilbertSpaceSize,...
    ones(1,HilbertSpaceSize),HilbertSpaceSize,HilbertSpaceSize);
H0 = sparse(HilbertSpaceSize,HilbertSpaceSize);
Hm = sparse(HilbertSpaceSize,HilbertSpaceSize);

for m = 1:L
    FlipFlopList = [];
    SpinVals = zeros(1,HilbertSpaceSize);
    zTermVals = zeros(1,HilbertSpaceSize);
    EnergyShiftVals = zeros(1,HilbertSpaceSize);
    %FlipFlop{m}=sparse(HilbertSpaceSize,HilbertSpaceSize);

    for cnt = 1:HilbertSpaceSize
        psi = states(cnt,:);

        SpinVals(cnt) = psi(m);
        zTermVals(cnt) = psi(m);
        EnergyShiftVals(cnt) = psi(m);
        phi = psi;

        if zTermVals(cnt) == 0
            zTermVals(cnt)= -1;
        end
            if psi(mod(m,L)+1)==psi(mod(m+1,L)+1)
               EnergyShiftVals(cnt) = 1;
            else
               EnergyShiftVals(cnt) = -1;

            
                    % phi(mod(m,L)+1)=m+1, phi(mod(m+1,L)+1)=m+2
                if phi(mod(m,L)+1)+phi(mod(m+1,L)+1) == 1
                   phi(mod(m,L)+1) = psi(mod(m+1,L)+1);
                   phi(mod(m+1,L)+1) = psi(mod(m,L)+1);

                    if max(phi) <= 1 && min(phi) >= 0
                       tnc = find(ismember(states,phi,'rows'));
                       %FlipFlop{m}(cnt,tnc)=1;
                       %FlipFlop{m}(tnc, cnt) = 1;
                       FlipFlopList = [FlipFlopList [cnt;tnc]];
                    else
                       keyboard
                    end
                end
            end
        %end
    end
    
    zTerm{m} = sparse(1:HilbertSpaceSize,1:HilbertSpaceSize,...
        zTermVals,HilbertSpaceSize,HilbertSpaceSize);

    EnergyShift{m} = sparse(1:HilbertSpaceSize,1:HilbertSpaceSize,...
        EnergyShiftVals,HilbertSpaceSize,HilbertSpaceSize);

    FlipFlop{m} = sparse(FlipFlopList(1,:),...
        FlipFlopList(2,:),ones(1,numel(FlipFlopList(1,:))),...
        HilbertSpaceSize,HilbertSpaceSize);

end

for m = 1:L
    if PBC || m < L-1
       H0 = H0+1/2*zTerm{m}*(1/2*FlipFlop{m}+1/4*EnergyShift{m});
    end

    %H0 = 1/2*zTerm{m}*H0;
end

coeffs = ones(1,L);

switch matterConfig
    case 'DWr'
        coeffs(1:(L/2)) = -coeffs(1:(L/2));
        
    case 'DWl'
        coeffs((L/2+1):end) = -coeffs((L/2+1):end);
        
end

for m = 1:L
    Hm = Hm+coeffs(m)*zTerm{m};
end

[Mm, Dm] = eig(full(Hm));
[~, idx] = max(diag(Dm));
psi0 = Mm(:, idx);

[M,E] = eig(full(H0),'vector');

psi0Rotated = M'*psi0;
magOpRotated = M'*Hm*M/L;

mag = NaN*ones(size(myTimes));
% occupationOpFirstSite =M'*zTerm{1}*M;
% occupationFirstSite = zeros(size(myTimes));
% occupationOpSecondSite =M'*zTerm{2}*M;
% occupationSecondSite = zeros(size(myTimes));
% occupationOpThirdSite =M'*zTerm{3}*M;
% occupationThirdSite = zeros(size(myTimes));
% occupationOpFourthSite =M'*zTerm{4}*M;
% occupationFourthSite = zeros(size(myTimes));
% occupationOpFifthSite =M'*zTerm{5}*M;
% occupationFifthSite = zeros(size(myTimes));
% occupationOpSixthSite =M'*zTerm{6}*M;
% occupationSixthSite = zeros(size(myTimes));
% occupationOpSeventhSite =M'*zTerm{7}*M;
% occupationSeventhSite = zeros(size(myTimes));
% occupationOpEighthSite =M'*zTerm{8}*M;
% occupationEighthSite = zeros(size(myTimes));
occupation = zeros(numel(myTimes), L);

for tInd = 1:numel(myTimes)
    %tic
    t = myTimes(tInd);
    psi = diag(exp(-1i*E*t))*psi0Rotated;
    for itN = 1:L
        occupation(tInd,itN)=real(psi'*(M'*zTerm{itN}*M)*psi);
     %occupationFirstSite(tInd) = real(psi'*(occupationOpFirstSite*psi));
   %occupationSecondSite(tInd) = real(psi'*(occupationOpSecondSite*psi));
   %occupationThirdSite(tInd) = real(psi'*(occupationOpThirdSite*psi));
    %occupationFourthSite(tInd) = real(psi'*(occupationOpFourthSite*psi));
    end
end
keyboard
imagesc([1 L],myTimes,real(occupation));
colorbar;
clim([-1 1]); 
set(gca,'FontSize',13,'LineWidth',1);
xlabel('Site index');
ylabel('Time');
title('Magnetization');

end