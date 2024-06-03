function [P, T] = UniqueNode(p, t)
    for kk = 1:length(t)
        Ve(kk) = 0.5.*det([ones(3,1) p(t(kk,:),:)]);
    end
    t((Ve==0), :) = [];
    P = unique(p(unique(sort(t(:))),:), 'rows');
    for i = 1:length(t)
        for j = 1:3
            T(i,j) = find(P(:,1)==p(t(i,j),1)&P(:,2)==p(t(i,j),2));
        end
end