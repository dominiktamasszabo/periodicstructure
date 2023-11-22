function [Gamma34_Matrix, Gamma34_Ey] = Gamma34(q_vec)

[eps_r, eps_0, M, phi_0, K, R, c, deltaX, deltaY, r_0, N1, N2, d, h, Resolution, V, chargeWeight] = defineConstants();

% Peremfeltétel Gamma 3,4 felületen (fenti,lenti): A térerősség normális
% komponense 0

osztas = linspace((deltaX/2*(-1+1/Resolution)),(deltaX/2*(1-1/Resolution)),Resolution);
Gamma34_x =  [osztas, osztas]; % Olyan vektor, ami a peremek x koordinátáit tartalmazza
Gamma34_y =  [repmat(-deltaY/2,1, Resolution), repmat(+deltaY/2,1, Resolution)];% Olyan vektor, ami a peremek y koordinátáit tartalmazza

% Most meg kell nézni, hogy ezekben a pontokban hogyan járulnak hozzá az
% egyes töltések a potenciálokhoz.

% Egy olyan mátrixot szeretnék létrehozni, aminek annyi oszlopa van, mint
% töltés és annyi sora mint diszkrét pont Gamma 1 felületen

Gamma34_Matrix = zeros(2*Resolution, length(q_vec));

% Most végig kell menni minden egyes töltésen (1 Coulomb) és beírni, hogy
% hogyan járul hozzá a potenciálhoz

for posi = 1:2*Resolution
    for vi = 1:length(q_vec)
        q_testing = chargeWeight*[zeros(1, vi-1), 1, zeros(1, length(q_vec)-vi)];
        [e_x, e_y, e_z] = tererosseg(Gamma34_x(posi), Gamma34_y(posi), q_testing);
        Gamma34_Matrix(posi, vi) = e_y; % Ez a normalis komponens
    end
end


% A várt y irányú térerősség komponens
Gamma34_Ey = zeros(2*Resolution,1);

end