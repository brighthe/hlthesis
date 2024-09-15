function [Ce, J1, vol] = FEA(t, t1, t2, p, Ve, BDY, E, nu)
    % 单元总数
    NT = length(t);
    % 单元刚度矩阵
    KK = zeros(6, 6*NT);
    % 计算实心单元的单元刚度矩阵
    for i = length(t1)+1:NT
        KK(:, 6*i-5:6*i) = GetKe(p(t(i,:), 1), p(t(i,:), 2), E, nu);
    end
    elemDof = zeros(NT, 6);
    elemDof(:, [1 3 5]) = 2*t-1;
    elemDof(:, [2 4 6]) = 2*t;
    iK = reshape(kron(elemDof, ones(6,1))', 36*NT, 1);
    jK = reshape(kron(elemDof, ones(1,6))', 36*NT, 1);
    sK = reshape(KK, 36*NT, 1);
    NK = sparse(iK, jK, sK, 2*length(p), 2*length(p));
    NK = (NK+NK')/2;
    fixedNodes = find(p(:,1) == BDY(1,1));
    forceNodes = find(p(:,1) == BDY(2,1) & p(:,2)==0);
    fixedDof = [2*fixedNodes-1; 2*fixedNodes];
    SolidDof = [2*unique(sort(t2(:)))-1; 2*unique(sort(t2(:)))];
    freeDofs = setdiff(SolidDof, fixedDof);
    U = zeros(2*length(p), 1);
    F = sparse(2*forceNodes, 1, -100, 2*length(p), 1);
    U(freeDofs,:) = NK(freeDofs, freeDofs) \ F(freeDofs, 1);
    for i = length(t1)+1:NT
        Ce(i) = 0.5 .* sum((U(elemDof(i,:))' * KK(:,6*i-5:6*i)) .* U(elemDof(i,:))', 2);
    end
    J1 = sum(Ce);
    vol = 1 - sum(Ve(1:length(t1))) / sum(Ve);
end