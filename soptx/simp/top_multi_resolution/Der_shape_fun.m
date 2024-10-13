%%--- Derivative of Shape function for [dN1/dx dN2/dx ...;dN1/dy dN2/dy ... ] and determinant of Jacobian ---%%
function [Be, Jdet] = Der_shape_fun(eg, ng, dx, dy)
    % coordinates of the 4 nodes (specific case for rectangle)
    Coorde=[0 dx dx 0;
            0 0 dy dy]; 
    
    % Derivative of the shape function
    DNe=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
        -(1-eg) -(1+eg) 1+eg 1-eg]';
    
    % Jacobian of the element e at location e,n
    J = Coorde*DNe;
    Jdet = det(J);
    
    % Derivative of the shape function
    Be = DNe/J;
end