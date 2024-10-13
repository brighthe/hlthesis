%%-------------- element matrix function --------------%%
function KE = stiff_ele(E, nu, dx, dy)
    % Gauss Points (4 in total)
    egv = [-1/sqrt(3), 1/sqrt(3)]; %2 gauss point x
    ngv = egv;                    %2 gauss point y
    wg = [1,1];                   %weight
    
    % Constitutive matrix
    matC = consti(E, nu);
    % initialize stiffness matrix
    KE = zeros(8, 8);
    % Gauss quadrature integration
    for eit = 1:2
        eg = egv(eit);
        for nit = 1:2
            ng = ngv(nit);
            % Derivative of the shape functions
            [Be, Jdet] = Der_shape_fun(eg, ng, dx, dy);
            %Derivative of the shape function for displacement in x and y
            BKe=[Be(1,1) 0 Be(2,1) 0 Be(3,1) 0 Be(4,1) 0;...
                   0 Be(1,1) 0 Be(2,1) 0 Be(3,1) 0 Be(4,1);
                    Be(1,2) 0 Be(2,2) 0 Be(3,2) 0 Be(4,2) 0;...
                    0 Be(1,2) 0 Be(2,2) 0 Be(3,2) 0 Be(4,2)];
            %Stifness for element in domain 1
            KE=KE+wg(eit)*wg(nit)*BKe'*matC*BKe*Jdet;
        end
    end
    end