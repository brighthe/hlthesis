%% Contour points disatnce adjust
function [C1] = RedisPoints(C, d1, d2)
    C = C';
    C1 = C;
    CL = sqrt(sum(diff(C1, 1, 1).^2, 2));
    for i = 1:(size(C, 1) - 1)
        if CL(i) < d1
            C1(i, :) = [0; 0];
            C1((i+1), :) = 0.5 * (C(i, :) + C((i+1), :));
        end
    end
    C1(all(C1==0, 2), :) = [];
    CL2 = sqrt(sum(diff(C1, 1, 1).^2, 2));
    Cmid = [];
    for i = 1:(size(C1,1)-1)
        if CL2(i) > d2 
            Cmid = [Cmid; 0.5*(C1(i,:) + C1((i+1),:))];
        end
    end
    if isempty(Cmid) == 0
        C1 = union(C1, Cmid, 'rows');
    end
end