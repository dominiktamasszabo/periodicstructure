function [Gamma1_Matrix, Gamma1_potentials] = Gamma1(cPMat)

[eps_r, eps_0, M, B, NoC, Resolution, phi_0, K, R, c_R, c_B, deltaX, deltaY, r_0, N1, N2, d, h, V] = defineConstants();

% Peremfeltétel Gamma 1 felületen (bal oldal): A potenciál -1/(2*N1)*V

Gamma1_x = repmat(-deltaX/2,1, Resolution); % Olyan vektor, ami a perem x koordinátáit tartalmazza
Gamma1_y = linspaceNoCorner(-deltaY/2, +deltaY/2, Resolution); % Olyan vektor, ami a perem y koordinátáit tartalmazza

% Most meg kell nézni, hogy ezekben a pontokban hogyan járulnak hozzá az
% egyes töltések a potenciálokhoz.

% Egy olyan mátrixot szeretnék létrehozni, aminek annyi oszlopa van, mint
% töltés és annyi sora mint diszkrét pont Gamma 1 felületen

Gamma1_Matrix = zeros(Resolution, NoC);

% Most végig kell menni minden egyes töltésen (1 Coulomb) és beírni, hogy
% hogyan járul hozzá a potenciálhoz

for posi = 1:Resolution
    for vi = 1:NoC
        cVec_testing = [zeros(1, vi-1), 1, zeros(1, NoC-vi)];
        Gamma1_Matrix(posi, vi) = potencial(Gamma1_x(posi), Gamma1_y(posi), cVec_testing, cPMat);
    end
end

% A várt potenciál vektor
Gamma1_POT = -1/(2*N1)*V;
Gamma1_potentials = repmat(Gamma1_POT,Resolution,1);
% Meg is van a Mátrix. Most ebből alapesetben egy túlhatározott
% egyenletrendszerünk van, Ami Ax = b alakú, ezt megoldva a négyzetes
% valamivel

% q_fromGamma1 = inv(Gamma1_Matrix' * Gamma1_Matrix) * Gamma1_Matrix' *
% Gamma1_potentials; Ezt még ki se kéne számolni


end