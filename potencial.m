function phi = potencial(x, y, q)

[eps_r, eps_0, M, phi_0, K, R, c, deltaX, deltaY, r_0, N1, N2, d, h, Resolution, V, chargeWeight] = defineConstants();

phi = 0;


% sqrt(x^2 + y^2) <= R | %%%% Ez még benne lehetne, de sokkal menőbb az
% ábra ha nincs.
if sqrt(x^2 + y^2) <= R | (x>deltaX/2) | (x<-deltaX/2) | (y>deltaY/2) | (y<-deltaY/2)
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

        dPhi = q_i / (eps_r*eps_0*2*pi) * log(r_0/r);

        phi = phi + dPhi;

        i = i+1;
    end
end

end
