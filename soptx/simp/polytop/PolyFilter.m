function [P] = PolyFilter(fem, R)
    if R < 0, P = speye(fem.NElem); return; end % P is set to identity when R < 0
    ElemCtrd = zeros(fem.NElem, 2);
    for e1 = 1:fem.NElem % Compute the centroids of all the elements
        vx = fem.Node(fem.Element{e1}, 1); 
        vy = fem.Node(fem.Element{e1}, 2);
        temp = vx .* vy([2:end 1]) - vy .* vx([2:end 1]);
        A = 0.5 * sum(temp);
        ElemCtrd(e1, 1) = 1 / (6 * A) * sum((vx + vx([2:end 1])) .* temp);
        ElemCtrd(e1, 2) = 1 / (6 * A) * sum((vy + vy([2:end 1])) .* temp);
    end
    [d] = DistPntSets(ElemCtrd, ElemCtrd, R); % Obtain distance values & indices
    P = sparse(d(:, 1), d(:, 2), 1 - d(:, 3) / R); % Assemble the filtering matrix
    P = spdiags(1 ./ sum(P, 2), 0, fem.NElem, fem.NElem) * P;
end

%------------Compute Distance Between Two Point Sets----------------%
function [d] = DistPntSets(PS1, PS2, R)
    d = cell(size(PS1, 1), 1);
    for e1 = 1:size(PS1, 1) % Compute the distance information
        dist = sqrt((PS1(e1, 1) - PS2(:, 1)).^2 + (PS1(e1, 2) - PS2(:, 2)).^2);
        [I, J] = find(dist < R); % Find the indices for distances less than R
        d{e1} = [I, J + (e1 - 1), dist(I)];
    end
    d = cell2mat(d); % Matrix of indices and distance value
end
%----------------------------------------------------------------------%
