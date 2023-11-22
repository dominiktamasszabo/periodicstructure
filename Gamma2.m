function [Gamma2_Matrix, Gamma2_potentials] = Gamma1(q_vec)

[eps_r, eps_0, M, phi_0, K, R, c, deltaX, deltaY, r_0, N1, N2, d, h, Resolution, V, chargeWeight] = defineConstants();

% Peremfeltétel Gamma 2 felületen (jobb oldal): A potenciál 1/(2*N1)*V

Gamma2_x = repmat(+deltaX/2,1, Resolution); % Olyan vektor, ami a perem x koordinátáit tartalmazza
Gamma2_y = linspace((deltaY/2*(-1+1/Resolution)),(deltaY/2*(1-1/Resolution)),Resolution); % Olyan vektor, ami a perem y koordinátáit tartalmazza

% Most meg kell nézni, hogy ezekben a pontokban hogyan járulnak hozzá az
% egyes töltések a potenciálokhoz.

% Egy olyan mátrixot szeretnék létrehozni, aminek annyi oszlopa van, mint
% töltés és annyi sora mint diszkrét pont Gamma 1 felületen

Gamma2_Matrix = zeros(Resolution, length(q_vec));

% Most végig kell menni minden egyes töltésen (1 Coulomb) és beírni, hogy
% hogyan járul hozzá a potenciálhoz

for posi = 1:Resolution
    for vi = 1:length(q_vec)
        q_testing = chargeWeight*[zeros(1, vi-1), 1, zeros(1, length(q_vec)-vi)];
    
        Gamma2_Matrix(posi, vi) = potencial(Gamma2_x(posi), Gamma2_y(posi), q_testing);
    end
end

% A várt potenciál vektor
Gamma2_POT = +1/(2*N1)*V;
Gamma2_potentials = repmat(Gamma2_POT,Resolution,1);
% Meg is van a Mátrix. Most ebből alapesetben egy túlhatározott
% egyenletrendszerünk van, Ami Ax = b alakú, ezt megoldva a négyzetes
% valamivel


end