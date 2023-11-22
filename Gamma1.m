function [Gamma1_Matrix, Gamma1_potentials] = Gamma1(q_vec)

[eps_r, eps_0, M, phi_0, K, R, c, deltaX, deltaY, r_0, N1, N2, d, h, Resolution, V, chargeWeight] = defineConstants();

% Peremfeltétel Gamma 1 felületen (bal oldal): A potenciál -1/(2*N1)*V

Gamma1_x = repmat(-deltaX/2,1, Resolution); % Olyan vektor, ami a perem x koordinátáit tartalmazza
Gamma1_y = linspace((deltaY/2*(-1+1/Resolution)),(deltaY/2*(1-1/Resolution)),Resolution); % Olyan vektor, ami a perem y koordinátáit tartalmazza

% Most meg kell nézni, hogy ezekben a pontokban hogyan járulnak hozzá az
% egyes töltések a potenciálokhoz.

% Egy olyan mátrixot szeretnék létrehozni, aminek annyi oszlopa van, mint
% töltés és annyi sora mint diszkrét pont Gamma 1 felületen

Gamma1_Matrix = zeros(Resolution, length(q_vec));

% Most végig kell menni minden egyes töltésen (1 Coulomb) és beírni, hogy
% hogyan járul hozzá a potenciálhoz

for posi = 1:Resolution
    for vi = 1:length(q_vec)
        q_testing = chargeWeight*[zeros(1, vi-1), 1, zeros(1, length(q_vec)-vi)];
    
        Gamma1_Matrix(posi, vi) = potencial(Gamma1_x(posi), Gamma1_y(posi), q_testing);
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