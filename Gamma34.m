function [Gamma34_Matrix, Gamma34_Ey] = Gamma34(cPMat)

[eps_r, eps_0, M, B, NoC, Resolution, phi_0, K, R, c_R, c_B, deltaX, deltaY, r_0, N1, N2, d, h, V] = defineConstants();

% Peremfeltétel Gamma 3,4 felületen (fenti,lenti): A térerősség normális
% komponense 0

osztas = linspaceNoCorner(-deltaX/2, +deltaX/2, Resolution);
Gamma34_x =  [osztas, osztas]; % Olyan vektor, ami a peremek x koordinátáit tartalmazza
Gamma34_y =  [repmat(-deltaY/2,1, Resolution), repmat(+deltaY/2,1, Resolution)];% Olyan vektor, ami a peremek y koordinátáit tartalmazza

% figure; plot(Gamma34_x,Gamma34_y,'LineWidth',2,'Color',[.6 0 0])

% Most meg kell nézni, hogy ezekben a pontokban hogyan járulnak hozzá az
% egyes töltések a potenciálokhoz.

% Egy olyan mátrixot szeretnék létrehozni, aminek annyi oszlopa van, mint
% töltés és annyi sora mint diszkrét pont Gamma 1 felületen

Gamma34_Matrix = zeros(2*Resolution, NoC);

% Most végig kell menni minden egyes töltésen (1 Coulomb) és beírni, hogy
% hogyan járul hozzá a potenciálhoz

for posi = 1:2*Resolution
    for vi = 1:NoC
        cVec_testing = [zeros(1, vi-1), 1, zeros(1, NoC - vi)];
        [e_x, e_y, ~] = tererosseg(Gamma34_x(posi), Gamma34_y(posi), cVec_testing, cPMat);
        e_r = sqrt(e_x^2+e_y^2);
        Gamma34_Matrix(posi, vi) = e_y/e_r; % Ez a normalis komponens
    end
end

% A várt y irányú térerősség komponens
Gamma34_Ey = zeros(2*Resolution,1);

end