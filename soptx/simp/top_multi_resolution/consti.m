%%----------------- Constitutive relationship in matrix form ---------------%%
function matC = consti(E, nu)
    % Lame constants
    lam = E*nu/(1+nu)/(1-2*nu);   
    miu = E/(2*(1+nu));
    lam = 2*miu*lam/(lam+2*miu);  %plane stress
    % Constitutive matrix
    matC = [lam+2*miu 0 0 lam;0 miu miu 0;0 miu miu 0;lam 0 0 lam+2*miu];
    end