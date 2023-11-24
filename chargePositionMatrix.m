function [cPMat] = chargePositionMatrix()

[eps_r, eps_0, M, B, NoC, Resolution, phi_0, K, R, c_R, c_B, deltaX, deltaY, r_0, N1, N2, d, h, V] = defineConstants();

% Annyi oszlop kell, amennyi töltés lesz (Összesen 4*B+M)
% 2 sor kell az x és y koordinátáknak.

cPMat = zeros(2,4*B+M);

for i = 1:M
    angle = i*2*pi/M + phi_0;
    x = c_R*R*cos(angle);
    y = c_R*R*sin(angle);
    cPMat(:,i) = [x;y];
end

uX = deltaX*c_B/2;
uY = deltaY*c_B/2;

% Gamma1 mentén
cPMat(1,M+1:M+B) = -uX;
cPMat(2,M+1:M+B) = linspaceNoCorner(-uY, +uY, B);

% Gamma2 mentén
cPMat(1,M+B+1:M+2*B) = +uX;
cPMat(2,M+B+1:M+2*B) = linspaceNoCorner(-uY, +uY, B);

% Gamma3 mentén
cPMat(1,M+2*B+1:M+3*B) = linspaceNoCorner(-uX, +uX, B);
cPMat(2,M+2*B+1:M+3*B) = -uY;
% Gamma4 mentén
cPMat(1,M+3*B+1:M+4*B) = linspaceNoCorner(-uX, +uX, B);
cPMat(2,M+3*B+1:M+4*B) = +uY;


% figure('name', 'cPMat');
% title('cPMat');
% xv = cPMat(1,:);
% yv = cPMat(2,:);
% plot(xv', yv', 'o','Color',[0.5 0 0.8]);
% % plot([-deltaX/2, +deltaX/2, +deltaX/2, -deltaX/2, -deltaX/2], ...
% %     [-deltaY/2, -deltaY/2, +deltaY/2, +deltaY/2, -deltaY/2], '-','Color',[0 0 0]);
% axis(1.1*[-c_B*deltaX/2 +c_B*deltaX/2 -c_B*deltaY/2 +c_B*deltaY/2 ])

end
