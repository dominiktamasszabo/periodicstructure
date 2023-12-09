clc; clear;

tic;

% Változó ciklus
global outsideR;
variables = 5;
variableVector = linspace(0.03, 0.45, variables);

eps_R_vec = zeros(1, variables);
rindex = 1;
for outsideR = variableVector



[eps_r, eps_0, M, B, NoC, Resolution, phi_0, K, R, c_R, c_B, deltaX, deltaY, r_0, N1, N2, d, h, V] = defineConstants();
cPMat = chargePositionMatrix();

[G1, b1] = Gamma1(cPMat);
[G2, b2] = Gamma2(cPMat);
[G34, b34] = Gamma34(cPMat);
[GR, bR] = GammaR(cPMat);


% S = ones(1,NoC);
S = [ones(1,M), zeros(1, NoC-M)];
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
W_p = 1/2 * eps_0*1e-12*eps_r*deltaX*deltaY/(Resolution^2)*sum(E_r.^2, "all") ; % J/m
% A hosszegysegre eso kapacitas
C_p = 2*N1*N2*W_p/(V^2)*1e12; % pF/m

C_ref = eps_0*1e-12*eps_r* h / d * 1e12; % pF/m

eps_R_calculated = C_p/C_ref;


disp(['Progress: ', num2str(rindex), ' / ', num2str(variables)]);
disp(['A fémrúd sugara: ', num2str(R), ' mm']);
disp(['A periodikus struktúrával hangolt kondenzátor hosszegységre eső kapacitása: ', num2str(C_p), ' pF/m']);
disp(['A periodikus struktúra nélküli kondenzátor hosszegységre eső kapacitása: ', num2str(C_ref), ' pF/m']);
disp(['Ezek aránya: ', num2str(eps_R_calculated)]);


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
div = 1;
terer_x = x_vec(1:div:end);
terer_y = y_vec(1:div:end);
terer_Ex = Ex(1:div:end,1:div:end);
terer_Ey = Ey(1:div:end,1:div:end);
qui = quiver(terer_x, terer_y, terer_Ex, terer_Ey);
title('Térerősség vektorok (V/m)', 'FontSize', 12); % Add a title
xlabel('x (mm)');ylabel('y (mm)');
% set(qui, 'AutoScaleFactor', 10);

% Térerősség abszolútérték ábrázolása
figure('name', 'E_abs'); 
  hsurf = surf(xx,yy,E_r );
  set(hsurf, 'edgeColor', 'none', 'faceAlpha', 1, 'faceLighting', 'flat');
  xlabel('x (mm)'); ylabel('y (mm)'); 
  zlabel('E (V/m)');
  title('Térerősség nagysága ', 'FontSize', 12); % Add a title

% Potenciál ábrázolása
figure_potential = figure('name', 'Potencial'); 
  surf_potential = surf(xx,yy,phi);
  colormap(figure_potential, hot);
  set(surf_potential, 'edgeColor', 'none', 'faceAlpha', 0.7, 'faceLighting', 'flat');
  xlabel('x (mm)'); ylabel('y (mm)'); 
  %% 
  zlabel('Potenciál (V)');
  title('Potenciál (V)', 'FontSize', 12);

% Gamma1 Potenciál  
figure('name', 'Gamma1 Potencial'); 
plot(yy(:,1), phi(:,1))
title('Gamma1 Potenciál');
xlabel('A tér y koordinátája');
ylabel('Potenciál(V)');
% axis([-0.5*10^-3 0.5*10^-3 -1 1]);

% Gamma2 Potenciál  
figure('name', 'Gamma2 Potencial');
plot(yy(:,Resolution), phi(:,Resolution))
title('Gamma2 Potenciál');
xlabel('A tér y koordinátája');
ylabel('Potenciál(V)');
% axis([-0.5*10^-3 0.5*10^-3 -1 1]);




eps_R_vec(rindex) = eps_R_calculated;
rindex = rindex + 1;


% Sugár ciklus vége
end

figure('name', 'Relatív kapacitás a fémrúd sugarának függvényében');
hold on;

% plot(radiusVector, eps_R_vec);
% xlabel('Fémrúd sugara (mm)');
% ylabel('Relatív hosszegységre eső kapacitás');


plot(variableVector, eps_R_vec, 'LineWidth', 2, 'Color', [0.2, 0.4, 0.8]); % Adjust line style and color

xlabel('Fémrúd sugara (mm)', 'FontSize', 12); % Increase font size and make it bold
ylabel('Relatív hosszegységre eső kapacitás', 'FontSize', 12);

title('Relatív kapacitás', 'FontSize', 12); % Add a title

grid on; % Show grid lines

toc;



