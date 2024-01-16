clc;
clear;

eps0 = 8.854;

% Olyan vektor sorozata (matrixok), amely tartalmazza a bázis-szakaszok kezdetét, végét.
[basisSectionBegin, basisSectionEnd] = generateBasisSections();
[potTestSectionBegin, potTestSectionEnd,  forceTestSectionBegin, forceTestSectionEnd, b] = generateTestSections();
% 
% figure;
% hold on;
% plot(basisSectionBegin(1,:),basisSectionBegin(2,:) , 'Color','r');
% plot(basisSectionEnd(1,:),basisSectionEnd(2,:) , 'Color','b');
% plot(testSectionBegin(1,:),testSectionBegin(2,:) , 'Color','g');
% plot(testSectionEnd(1,:),testSectionEnd(2,:) , 'Color','y');
% 

potWeights = zeros(length(potTestSectionBegin), length(basisSectionBegin));
tic;
for i = 1:length(basisSectionBegin)
    for j = 1:length(potTestSectionBegin)
        if(i==j)
%             continue
        end

        disp([i,j]);

        bbegin = basisSectionBegin(:, i);
        bend = basisSectionEnd(:, i);
        bc = (bbegin + bend) / 2;
        bv = (bend - bbegin) / 2;

        tbegin = potTestSectionBegin(:, j);
        tend = potTestSectionEnd(:, j);
        tc = (tbegin + tend) / 2;
        tv = (tend - tbegin) / 2;

        fun = @(bparam,tparam) 1 .* log(sqrt((bc(1)+bparam*bv(1)).^2+(bc(2)+bparam*bv(2)).^2) ./ sqrt(((bc(1)+bparam*bv(1))-(tc(1)+tparam*tv(1))).^2+(((bc(2)+bparam*bv(2))-(tc(2)+tparam*tv(2))).^2))); % triangularPulse
        potWeights(j,i) = integral2(fun,-1,+1,-1,+1, 'RelTol',1e-3)/(2*pi*eps0); %/dp!?

%         bbeginx = basisSectionBegin(1,i);
%         bbeginy = basisSectionBegin(2,i);
%         bendx = basisSectionEnd(1,i);
%         bendy = basisSectionEnd(2,i);
% 
%         cx = (bbeginx+bendx)/2;
%         cy = (bbeginy+bendy)/2;
    end
end
toc;

forceWeights = zeros(length(forceTestSectionBegin), length(basisSectionBegin));
tic;
for i = 1:length(basisSectionBegin)
    for j = 1:length(forceTestSectionBegin)
        if(i==j)
%             continue
        end

        disp([i,j]);

        bbegin = basisSectionBegin(:, i);
        bend = basisSectionEnd(:, i);
        bc = (bbegin + bend) / 2;
        bv = (bend - bbegin) / 2;

        tbegin = forceTestSectionBegin(:, j);
        tend = forceTestSectionEnd(:, j);
        tc = (tbegin + tend) / 2;
        tv = (tend - tbegin) / 2;

        fun = @(bparam,tparam) 1 .* ((bc(2)+bparam*bv(2))-(tc(2)+tparam*tv(2))) ./ (((bc(1)+bparam*bv(1))-(tc(1)+tparam*tv(1))).^2+(((bc(2)+bparam*bv(2))-(tc(2)+tparam*tv(2))).^2)); % triangularPulse
        forceWeights(j,i) = integral2(fun,-1,+1,-1,+1, 'RelTol',1e-3)/(2*pi*eps0); %/dp!?

%         bbeginx = basisSectionBegin(1,i);
%         bbeginy = basisSectionBegin(2,i);
%         bendx = basisSectionEnd(1,i);
%         bendy = basisSectionEnd(2,i);
% 
%         cx = (bbeginx+bendx)/2;
%         cy = (bbeginy+bendy)/2;
    end
end
toc;

A = [potWeights;forceWeights];

c = A \ b;

tic;
num = 200;
vec = linspace(-0.5,0.5, num);
[xx,yy] = meshgrid(vec);
potential = zeros(num, num);

for ix = 1:num
    disp(['Calculating potential at ix = ', num2str(ix)]);
    for iy = 1:num
        x = xx(ix, iy);
        y = yy(ix, iy);
        for i = 1:length(basisSectionBegin)
            bbegin = basisSectionBegin(:, i);
            bend = basisSectionEnd(:, i);
            bc = (bbegin + bend) / 2;
            bv = (bend - bbegin) / 2;
    
            fun = @(bparam) 1 .* log(sqrt((bc(1)+bparam*bv(1)).^2+(bc(2)+bparam*bv(2)).^2) ./ sqrt(((bc(1)+bparam*bv(1))-x).^2+(((bc(2)+bparam*bv(2))-y).^2))); %triangularPulse(bparam)
            add = c(i)*integral(fun, -1, +1, 'RelTol',1e-2)/(2*pi*eps0);%/dp!?
            potential(ix, iy) = potential(ix, iy)+add;
        end
    end
end


force = zeros(num, num);

for ix = 1:num
    disp(['Calculating field strength at ix = ', num2str(ix)]);
    for iy = 1:num
        x = xx(ix, iy);
        y = yy(ix, iy);
        for i = 1:length(basisSectionBegin)
            bbegin = basisSectionBegin(:, i);
            bend = basisSectionEnd(:, i);
            bc = (bbegin + bend) / 2;
            bv = (bend - bbegin) / 2;
    
            fun = @(bparam) 1 .* (1 ./ sqrt(((bc(1)+bparam*bv(1))-x).^2+(((bc(2)+bparam*bv(2))-y).^2))); %triangularPulse(bparam)
            add = c(i)*integral(fun, -1, +1, 'RelTol',1e-2)/(2*pi*eps0);%/dp!?
            force(ix, iy) = force(ix, iy)+add;
        end
    end
end

d = 0.3;

for ix = 1:num
    for iy = 1:num

        x = xx(ix, iy);
        y = yy(ix, iy);
  
        if(abs(x) < d/2 && abs(y) < d/2)
            potential(ix, iy) = 0;
            force(ix, iy) = 0;
        end
    end
end

toc;

figure;
surf( xx, yy,potential);

figure;
surf( xx, yy,force);

N1=1;
N2=N1;

% A hosszegysegre eso energia
W_p = 1/2 * eps0*1e-12*1^2/(num^2)*sum(force.^2, "all") ; % J/m
% A hosszegysegre eso kapacitas
C_p = 2*N1*N2*W_p/(1^2)*1e12; % pF/m

C_ref = eps0*1e-12* 1 / 1 * 1e12; % pF/m

eps_R_calculated = C_p/C_ref;

disp(['A periodikus struktúrával hangolt kondenzátor hosszegységre eső kapacitása: ', num2str(C_p), ' pF/m']);
disp(['A periodikus struktúra nélküli kondenzátor hosszegységre eső kapacitása: ', num2str(C_ref), ' pF/m']);
disp(['Ezek aránya: ', num2str(eps_R_calculated)]);