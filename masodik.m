clc; clear;
[eps_r, eps_0, M, B, NoC, Resolution, phi_0, K, R, c_R, c_B, deltaX, deltaY, r_0, N1, N2, d, h, V] = defineConstants();
cPMat = chargePositionMatrix();

[G1, b1] = Gamma1(cPMat);
[G2, b2] = Gamma2(cPMat);
[G34, b34] = Gamma34(cPMat);
[GR, bR] = GammaR(cPMat);

S = ones(1,NoC);
A = [S;G1;G2;GR;G34];
b = [0;b1;b2;bR;b34];

cVec = (A\b)';

disp("kesz");

chargeSum = sum(cVec);

% Számolgatás

x_vec = linspace(-deltaX/2, deltaX/2, Resolution);
y_vec = linspace(-deltaY/2, deltaY/2, Resolution);
[xx,yy] = meshgrid(x_vec,y_vec); 

Ex = zeros(size(xx));
Ey = zeros(size(xx));

% TESZTX = zeros(size(xx));
% TESZTY = zeros(size(yy));

for vi_x = 1:length(x_vec) % vi = vector index
    for vi_y = 1:length(y_vec)
        mi_x = vi_y; % Matrix index X iranyban
        mi_y = vi_x; % Matrix index Y iranyban

        [e_x, e_y, e_z] = tererosseg(x_vec(vi_x), y_vec(vi_y), cVec, cPMat);
        Ex(mi_x,mi_y) = e_x;
        Ey(mi_x,mi_y) = e_y;
        Ez(mi_x,mi_y) = e_z;

% % %       Csak kód teszteléshez a helyes indexelés megtalálásához.
% %         TESZTX(vi_y,vi_x) = x_vec(vi_x); % Olyan legyen, mint xx 
% %         TESZTY(vi_y,vi_x) = y_vec(vi_y); % Olyan legyen, mint yy
    end
end

E_r = sqrt(Ex.^2+Ey.^2+Ez.^2);

% Potenciál kiszámítása

phi = zeros(size(xx));

for vi_x = 1:length(x_vec) % vi = vector index
    for vi_y = 1:length(y_vec)
        mi_x = vi_y; % Matrix index X iranyban
        mi_y = vi_x; % Matrix index Y iranyban

        phi(mi_x, mi_y) = potencial(x_vec(vi_x), y_vec(vi_y), cVec, cPMat);
    end
end

% A hosszegysegre eso energia
W_p = eps_0*eps_r*deltaX*deltaY/(Resolution^2)*sum(E_r.^2, "all")*10^-18 ; % J/m
% A hosszegysegre eso kapacitas
C_p = 2*N1*N2*W_p/(V^2); % F/m





%% INNEN ÁBRÁZOLÁS

figure('name', 'Töltéseloszlás');
% plot3(cPMat(1,:), cPMat(2,:), cVec, 'o', 'Color', [1 0 0]);
scatter3(cPMat(1,:), cPMat(2,:), cVec, 'filled'); % scatter plot
hold on;
stem3(cPMat(1,:), cPMat(2,:), cVec, 'k-'); % stem plot with black dashed lines
patch(1.1*c_B*[-deltaX/2, +deltaX/2, +deltaX/2, -deltaX/2], ...
    1.1*c_B*[-deltaY/2, -deltaY/2, +deltaY/2, +deltaY/2], [0, 0, 0, 0], 'b', 'FaceAlpha', 0.3);
rectangle('Position', [-deltaX/2, -deltaY/2, deltaX, deltaY], 'EdgeColor', 'k', 'LineWidth', 0.6);
rectangle('Position', [-R, -R, 2*R, 2*R], 'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 0.6);
hold off;
zlabel('Töltés (pC/mm)');
xlabel('x (mm)');ylabel('y (mm)');
title('Töltések eloszlása az XY síkon');
% cVecMax = max(abs(cVec(1:M)));
% zlim([-cVecMax, +cVecMax]);

% Térerősség vektorok ábrázolása
figure('name', 'Tererosseg vektorok');
qui = quiver(x_vec, y_vec, Ex, Ey);
% set(qui, 'AutoScaleFactor', 10);

% Térerősség abszolútérték ábrázolása
figure('name', 'E_abs'); 
  h = surf(xx,yy,E_r );
  set(h, 'edgeColor', 'none', 'faceAlpha', 0.5, 'faceLighting', 'flat');
  xlabel('x'); ylabel('y'); 
  zlabel('E');

% Potenciál ábrázolása
figure_potential = figure('name', 'Potencial'); 
  surf_potential = surf(xx,yy,phi);
  colormap(figure_potential, hot);
  set(surf_potential, 'edgeColor', 'none', 'faceAlpha', 0.7, 'faceLighting', 'flat');
  xlabel('x'); ylabel('y'); 
  zlabel('Phi');

% Gamma1 Potenciál  
figure('name', 'Gamma1 Potencial'); 
plot(yy(:,1), phi(:,1))
title('Gamma1 Potenciál');
xlabel('A tér y koordinátája');
ylabel('Potenciál(V)');
% axis([-0.5*10^-3 0.5*10^-3 -1 1]);

% Gamma2 Potenciál  
figure('name', 'Gamma2 Potencial')
; plot(yy(:,Resolution), phi(:,Resolution))
title('Gamma2 Potenciál');
xlabel('A tér y koordinátája');
ylabel('Potenciál(V)');
% axis([-0.5*10^-3 0.5*10^-3 -1 1]);