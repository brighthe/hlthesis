function [nonlindiff] = NonlinearDiffsion(p, t, t1, Ve)
    f = E2N(t, p, sparse(length(t1)+1:length(t), 1, 1, length(t), 1), Ve);
    for A0 = 1:length(p)
        Lf = [];
        [row, ~] = find(t==A0);
        ts = t(row,:);
        PT = setdiff(ts, A0);
        for i = 1:length(PT)
            angl = [];
            [secp, ~] = find(ts==PT(i));
            for k = 1:length(secp)
                A1 = PT(i);
                A2 = setdiff(ts(secp(k),:), [A0, A1]);
                angl(k) = atan2( 2*Ve( secp(k) ), dot( p(A0,:)-p(A2,:), p(A1,:)-p(A2,:) ) );
            end
            sigma = 0.8;
            g = 1/(1+((f(A1)-f(A0))^2/sigma^2));
            if length(angl) == 2
                Lf(i) = (cot(angl(1)) + cot(angl(2))) * (f(A1)-f(A0))*g;
            else
                Lf(i) = (2*cot(angl(1))) * (f(A1)-f(A0))*g;
            end
        end
        nonlindiff(A0) = 1/(2*sum(Ve(row))) * sum(Lf);
    end
end