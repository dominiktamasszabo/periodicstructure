function [GammaR_Matrix, GammaR_Et] = GammaR(q_vec)

[eps_r, eps_0, M, phi_0, K, R, c, deltaX, deltaY, r_0, N1, N2, d, h, Resolution, V, chargeWeight] = defineConstants();

% Peremfeltétel Gamma R felületen (bal oldal): A potenciál -1/(2*N1)*V

% Vektor a megfelelő szögekkel

angles = linspace(0, 2*pi, Resolution+1);
angles = angles(1, 1:Resolution); % Azokat a szögeket tartalmazza, amelyeknek megfelelő pontokban nézzük a tangenciális komponenst

Gamma1_x = R*cos(angles);
Gamma1_y = R*sin(angles);

GammaR_Matrix = zeros(Resolution, length(q_vec));

% Most végig kell menni minden egyes töltésen (1 Coulomb) és beírni, hogy
% hogyan járul hozzá a potenciálhoz

for posi = 1:Resolution
    for vi = 1:length(q_vec)
        q_testing = chargeWeight*[zeros(1, vi-1), 1, zeros(1, length(q_vec)-vi)];
        

        [e_x, e_y, e_z] = tererosseg(Gamma1_x(posi), Gamma1_y(posi), q_testing);
        e_t = e_x*sin(angles(posi))-e_y*cos(angles(posi));
        e_r = sqrt(e_x^2+e_y^2);
%         disp([e_x, e_y, e_t]);
        GammaR_Matrix(posi, vi) = e_t/e_r; % Ez a tangencialis komponens
    end
end
% A várt potenciál vektor
GammaR_Et = zeros(Resolution, 1);
% Meg is van a Mátrix. Most ebből alapesetben egy túlhatározott
% egyenletrendszerünk van, Ami Ax = b alakú, ezt megoldva a négyzetes
% valamivel


end