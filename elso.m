clc; clear;
[eps_r, eps_0, M, phi_0, K, R, c, deltaX, deltaY, r_0, N1, N2, d, h, Resolution, V, chargeWeight] = defineConstants();

q_vec = 1e-9*zeros(1,Resolution);



% [Ex, Ey, Ez] = tererosseg(0.4,0.4, q_vec);

% q_vec = 10^-9*[1,0,0,0,0,1,0,0];


[G1, b1] = Gamma1(q_vec);
[G2, b2] = Gamma2(q_vec);
[G34, b34] = Gamma34(q_vec);
[GR, bR] = GammaR(q_vec);

% Együttható normalizálás
mu1 = mean(abs(G1'))';
mu2 = mean(abs(G2'))';
mu34 = mean(abs(G34'))';
muR = mean(abs(GR'))';

% Próbálkozás: leosztok a sorok abszolutértékének átlagaival 
G1 = G1./mu1;
G2 = G2./mu2;
G34 = G34./mu34;
GR = GR./muR;

b1 = b1 ./ mu1;
b2 = b2 ./ mu2;
b34 = b34 ./ mu34;
bR = bR ./ muR;


A = [G1;G2;G34;GR];
b = [b1;b2;b34;bR];

% 
%  A = [G1;G2];
%  b = [b1;b2];
% 
% A = GR;
% b = bR;

q_vec = chargeWeight*( A\b)';



%% INNEN ÁBRÁZOLÁS

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

        [e_x, e_y, e_z] = tererosseg(x_vec(vi_x), y_vec(vi_y), q_vec);
        Ex(mi_x,mi_y) = e_x;
        Ey(mi_x,mi_y) = e_y;
        Ez(mi_x,mi_y) = e_z;

% % %       Csak kód teszteléshez a helyes indexelés megtalálásához.
% %         TESZTX(vi_y,vi_x) = x_vec(vi_x); % Olyan legyen, mint xx 
% %         TESZTY(vi_y,vi_x) = y_vec(vi_y); % Olyan legyen, mint yy
    end
end

% Potenciál kiszámítása

phi = zeros(size(xx));

for vi_x = 1:length(x_vec) % vi = vector index
    for vi_y = 1:length(y_vec)
        mi_x = vi_y; % Matrix index X iranyban
        mi_y = vi_x; % Matrix index Y iranyban

        phi(mi_x, mi_y) = potencial(x_vec(vi_x), y_vec(vi_y), q_vec);
    end
end


% Térerősség vektorok ábrázolása
figure;
qui = quiver(x_vec, y_vec, Ex, Ey);
% set(qui, 'AutoScaleFactor', 10); 
E_r = sqrt(Ex.^2+Ey.^2+Ez.^2);

% Térerősség abszolútérték ábrázolása
figure; 
  h = surf(xx,yy,E_r );
  set(h, 'edgeColor', 'none', 'faceAlpha', 0.5, 'faceLighting', 'flat');
  xlabel('x'); ylabel('y'); 
  zlabel('E');


% Potenciál ábrázolása
figure_potential = figure; 
  surf_potential = surf(xx,yy,phi);
  colormap(figure_potential, hot);
  set(surf_potential, 'edgeColor', 'none', 'faceAlpha', 0.5, 'faceLighting', 'flat');
  xlabel('x'); ylabel('y'); 
  zlabel('Phi');

  
% Gamma1 Potenciál  
figure; plot(yy(:,1), phi(:,1))
title('Gamma1 Potenciál');
xlabel('A tér y koordinátája');
ylabel('Potenciál(V)');
% axis([-0.5*10^-3 0.5*10^-3 -1 1]);

% Gamma1 Potenciál  
figure; plot(yy(:,Resolution), phi(:,Resolution))
title('Gamma2 Potenciál');
xlabel('A tér y koordinátája');
ylabel('Potenciál(V)');
% axis([-0.5*10^-3 0.5*10^-3 -1 1]);
