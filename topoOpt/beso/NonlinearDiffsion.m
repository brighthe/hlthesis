function [nonlindiff] = NonlinearDiffsion(p, t, t1, Ve)
    f = E2N(t, p, sparse(length(t1)+1:lenght(t), 1, 1, length(t), 1), Ve);
    for A0 = 1:length(p)
        Lf = [];
        [row, ~] = find(t==A0);
        ts = t(row,:);
        PT = setdiff(ts, A0);
        for i = 1:length(PT)
            angl = [];
            [secp, ~] = find(ts==PT(i));
            for k = 1:length(secp)
                A1 = PT(i)
            end
        end
end