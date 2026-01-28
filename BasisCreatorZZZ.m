function basis = BasisCreatorZZZ(L, N)

% Calculate all possible matter configurations
combinations = nchoosek(1:L, N);
basis = zeros(size(combinations, 1), L);

for i = 1:size(combinations, 1)
    basis(i, combinations(i, :)) = 1;
end