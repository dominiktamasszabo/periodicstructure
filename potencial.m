function phi = potencial(x, y, cVec, cPMat)

cPMatsize = size(cPMat);
assert(length(cVec) == cPMatsize(1,2), "Nem egyezik a cVec es a cPMat hossza");

[eps_r, eps_0, M, B, NoC, Resolution, phi_0, K, R, c_R, c_B, deltaX, deltaY, r_0, N1, N2, d, h, V] = defineConstants();


phi = 0;

if sqrt(x^2 + y^2) <= R | (x>deltaX/2) | (x<-deltaX/2) | (y>deltaY/2) | (y<-deltaY/2)
    % Minden 0
else
    i = 1;
    for q_i = cVec
        x_q = cPMat(1,i);
        y_q = cPMat(2,i);
        
        dx = x-x_q;
        dy = y-y_q;

        r = sqrt(dx^2+dy^2);

        dPhi = q_i / (eps_r*eps_0*2*pi) * log(r_0/r);

        phi = phi + dPhi;

        i = i+1;
    end
end

end
