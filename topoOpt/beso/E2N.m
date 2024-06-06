function [dN] = E2N(t, p, x, Ve)
    dN = zeros(length(p), 1);
    for i = 1:length(p)
        [row, ~] = find(t==i);
        dN(i) = dot(Ve(row), x(row)) / sum(Ve(row));
    end
end