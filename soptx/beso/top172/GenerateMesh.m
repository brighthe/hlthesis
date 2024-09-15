function [p, t, t1, t2, Ve] = GenerateMesh(xn, yn, dN, maxedge, minedge, BDY, maxiter)
    x = xn;
    x(2:2:end,:) = x(2:2:end,:) + 0.005;
    pi = [x(:), yn(:)];
    C = contour(xn, yn, dN, [0 0]);
    C1 = [];
    dump = find((C(1,:)==0) & (C(2,:)>=2));
    for i = 1:length(dump) - 1
        C1 = [C1 (RedisPoints (C(:,dump(i)+1:dump(i+1)-1), 0.004, 0.01))'];
    end
    last = C(:,dump(end)+1:end);
    C = [C1 (RedisPoints (last, 0.004, 0.01))'];
    d = zeros(size(pi, 1), 1);
    for i = 1:size(pi, 1)
        d(i) = sqrt(min((pi(i, 1) - C(1,:)).^2 + (pi(i, 2) - C(2,:)).^2));
    end
    r0 = 1./min(max(minedge, d), maxedge).^2;
    pfix = [C'; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2);
        BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1), 0];
    p = [pfix; pi(r0./max(r0)>0.5,:)];
    p1 = 0;
    warning off;
    for i = 1:maxiter
        if max(sum((p-1).^2, 2)) > 1e-6
            t = delaunay(p);
            edges = unique(sort([t(:,[1,2]); t(:,[1,3]); t(:,[2,3])], 2), 'rows');
            p1 = p;
        end
        midpoint = (p(edges(:,1),:) + p(edges(:,2),:)) / 2;
        d = zeros(size(midpoint,1), 1);
        for j = 1:size(midpoint, 1)
            d(j) = sqrt(min((midpoint(j,1)-C(1,:)).^2 + (midpoint(j,2)-C(2,:)).^2));
        end
        L = sqrt(sum((p(edges(:,1),:) - p(edges(:,2),:)).^2, 2));
        L1 = min(max(minedge,d), maxedge);
        L0 = 1.2*L1*sqrt(sum(L.^2) / sum(L1.^2));
        Fe = max(L0 - L, 0) ./ L*[1,1] .* (p(edges(:,1),:) - p(edges(:,2),:));
        Fp = full(sparse(edges(:, [1,1,2,2]), ones(size(d))*[1,2,1,2], [Fe,-Fe], size(p,1), 2));
        Fp(1:size(pfix, 1),:) = 0;
        p = p + 0.2*Fp;
        p(:, 1) = min(BDY(2, 1), max(BDY(1, 1), p(:, 1)));
        p(:, 2) = min(BDY(2, 2), max(BDY(1, 2), p(:, 2)));
    end
    [p, t] = UniqueNode(p, t);
    pmid = (p(t(:,1),:) + p(t(:,2),:) + p(t(:,3),:)) / 3;
    dE = interp2(xn, yn, dN, pmid(:,1), pmid(:,2), 'cubic');
    t1 = t(dE < 0,:);
    t2 = t(dE >= 0,:);
    t = [t1; t2];
    for kk = 1:length(t)
        Ve(kk) = 0.5.*det([ones(3, 1) p(t(kk,:),:)]);
    end
end