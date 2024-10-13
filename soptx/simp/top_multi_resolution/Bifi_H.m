nelx = 100;
nely = 40;
volfrac = 0.5;
penal=3;
rmin=2;
inf=0;
theta=0;
mag = 1;
Ey = 1;
ft = 1;

%% MATERIAL PROPERTIES
E0 = 1*Ey;

Emin = 1e-9;
nu = 0.3;
%%----------------------- PREPARE FINITE ELEMENT ANALYSIS ----------------------%%
dx = 1;   %size of the element in x
dy = 1;   %size of the element in y

%---------------- stiffness matrix of one element (rectangular) ----------------%
KE = stiff_ele(E0, nu, dx, dy);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

fixeddofs = [1 2];
for i=1:nelx
    fixeddofs = [fixeddofs 1+i*(nely+1)*2 2+i*(nely+1)*2];
end

%fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%%----------------------- INITIALIZE ITERATION -------------------%%

x = repmat(volfrac, nely, nelx);
x(1:nely/2, nelx/2:nelx) = 0.001;
xPhys = x;
loop = 0;
change = 1;

%%--------------------------- START ITERATION ---------------------------%%
while ((change > 0.01) && (loop<277))
    loop = loop + 1;
    %------------------- FE-ANALYSIS -------------------%
    load data_100.mat;
    load des_quad.mat;
    wquad = w';  
    m = size(wquad, 2);
    xfilt = (H*reshape(xPhys, nelx*nely, 1)) ./ Hs;
    
    %------------------- Bifi -------------------%
    
    mesh_l = 10;
    zrb = 100 / mesh_l;
    
    for i = 1:mesh_l
        for j = 1:mesh_l
            xL(i,j) = sum(sum(xPhys(1+(i-1)*zrb:i*zrb, 1+(j-1)*zrb:j*zrb))) / (zrb^2);
        end
    end
    
    [UL] = randomize_L(mesh_l, mesh_l, 3, 3, mesh_l/2+1, xL);
end