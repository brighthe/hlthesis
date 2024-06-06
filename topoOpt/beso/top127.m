nx = 80;
ny = 50;
ER = 0.05;
Vmax = 0.5;
maxedge = 0.2;
minedge = 0.025;
E = 1e5;
nu = 0.3;
BDY = [-.5*(nx) -.5*(ny); .5*(nx), .5*(ny)] / 100;
[xn, yn] = meshgrid(BDY(1,1):(BDY(2,1)-BDY(1,1)) / (nx):BDY(2,1),BDY(1,2):(BDY(2,2)-BDY(1,2)) / (ny):BDY(2,2));
dN = sin(xn/BDY(2, 1) * 6*pi) .* cos(yn/BDY(2, 1) * 6*pi) + 0.5;

for iterNum = 1
    % 生成网格
    [p, t, t1, t2, Ve] = GenerateMesh(xn, yn, dN, maxedge, minedge, BDY, 80);
    % 非线性扩散项
    [nonlidiff] = NonlinearDiffsion(p, t, t1, Ve);
    % 有限元分析和灵敏度计算
    [Ce, J, vol] = FEA(t, t1, t2, p, Ve, BDY, E, nu);
    V = max(Vmax, vol-ER);
    Cn = E2N(t, p, Ce, Ve) + tau*nonlidiff';
    if iterNum > 1
        Cnlast = griddata(pold(:,1), pold(:,2), Cnold, p(:,1), p(:,2), 'cubic');
    end
    Cnold = Cn;
    pold = p;
    if iterNum > 1
        Cn = 0.5 * (Cn +Cnlast);
    end
end

% 绘制所有单元
figure;
triplot(t, p(:,1), p(:,2));
title('所有生成的网格');
xlabel('x');
ylabel('y');

% 绘制实心单元
figure;
triplot(t1, p(:,1), p(:,2));
title('空心单元网格');
xlabel('x');
ylabel('y');

% 绘制空心单元
figure;
triplot(t2, p(:,1), p(:,2));
title('实心单元网格');
xlabel('x');
ylabel('y');