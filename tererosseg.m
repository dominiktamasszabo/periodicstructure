function [Ex, Ey, Ez] = tererosseg(x, y, cVec, cPMat)

cPMatsize = size(cPMat);
assert(length(cVec) == cPMatsize(1,2), "Nem egyezik a cVec es a cPMat hossza");

[eps_r, eps_0, M, B, NoC, Resolution, phi_0, K, R, c_R, c_B, deltaX, deltaY, r_0, N1, N2, d, h, V] = defineConstants();

Ex = 0;
Ey = 0;
Ez = 0;


if sqrt(x^2 + y^2) < R*0.99 | (x>deltaX/2) | (x<-deltaX/2) | (y>deltaY/2) | (y<-deltaY/2) % ITT AZÉRT VAN 0.95 MERT A GR IS EZT HASZNÁLJA
    % Minden 0
else
    i = 1;
    for q_i = cVec
        x_q = cPMat(1,i);
        y_q = cPMat(2,i);
        
        dx = x-x_q;
        dy = y-y_q;

        r = sqrt(dx^2+dy^2);

        dE = q_i / (eps_r*eps_0*2*pi*r);

        dEx = dE * dx/r; % x komponens
        dEy = dE * dy/r; % y komponens
        % z komponens nincs

        Ex = Ex + dEx;
        Ey = Ey + dEy;
        i = i+1;
    end
end

end
