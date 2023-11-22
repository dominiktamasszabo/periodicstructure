function [Ex, Ey, Ez] = tererosseg(x, y, q)

[eps_r, eps_0, M, phi_0, K, R, c, deltaX, deltaY, r_0, N1, N2, d, h, Resolution, V, chargeWeight] = defineConstants();

Ex = 0;
Ey = 0;
Ez = 0;


if sqrt(x^2 + y^2) < R*0.95 | (x>deltaX/2) | (x<-deltaX/2) | (y>deltaY/2) | (y<-deltaY/2) % ITT AZÉRT VAN 0.95 MERT A GR IS EZT HASZNÁLJA
    % Minden 0
else
    i = 0;
    for q_i = q
        M=length(q); % Így nem használjuk fel az M betűs konstanst!

        angle = i*2*pi/M + phi_0;
        x_q = c*R*cos(angle);
        y_q = c*R*sin(angle);
        
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
