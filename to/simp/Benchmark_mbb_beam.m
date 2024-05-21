nelx = 60;
nely = 20;
volfrac = 0.5;
penal = 3.0;
rmin = 1.5;

x(1:nely, 1:nelx) = volfrac;
loop = 0;
change = 1.;
% Start Iteration
while change > 0.01
    loop = loop + 1;
    xold = x;
    % FE-Analysis
    [U] = FE(nelx, nely, x, penal);
    
    % K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
    % F = sparse(2*(nely+1)*(nelx+1), 1); U = zeros(2*(nely+1)*(nelx+1), 1);
    % for elx = 1:nelx
    %     for ely = 1:nely
    %         n1 = (nely+1)*(elx-1)+ely;
    %         n2 = (nely+1)* elx   +ely;
    %         edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    %         K(edof, edof) = K(edof, edof) + x(ely, elx)^penal*KE;
    %     end
    % end
    % % Define Loads And Suppotrs (Half-MBB-Beam)
    % F(2, 1) = -1;
    % fixeddofs = union([1:2:2*(nely+1)], [2*(nelx+1)*(nely+1)]);
    % alldofs = [1:2*(nely+1)*(nelx+1)];
    % freedofs = setdiff(alldofs, fixeddofs);
    % % SOLVING
    % U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
    % U(fixeddofs,:)= 0;

    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [KE] = lk;
    c = 0.;
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1)+ely;
            n2 = (nely+1)* elx   +ely;
            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
            c = c + x(ely, elx)^penal*Ue'*KE*Ue;
            dc(ely, elx) = -penal*x(ely, elx)^(penal-1)*Ue'*KE*Ue;
        end
    end
    % Filtering Of Sensitivity
    [dc] = check(nelx, nely, rmin, x, dc);
    % Design Update By The Optimality Criteria Method
    [x] = OC(nelx, nely, x, volfrac, dc);
    % Print Results
    change = max(max(abs(x-xold)));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
        ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
    % Plot Densities
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
end

    