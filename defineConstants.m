function [eps_r, eps_0, M, B, NoC, Resolution, phi_0, K, R, c_R, c_B, deltaX, deltaY, r_0, N1, N2, d, h, V] = defineConstants()

% Konstansok
eps_r = 1; % Relatív permittivitás [1]
eps_0 = 8.86; % Vákuum permittivitás [F/mm*10^12] !!! EZ NEM BIZTOS HOGY JÓ !!!!
phi_0 = 0*2*pi/16*0.5; % Első vonaltöltés szöge [rad]
K = 40; % Y-ra párhuzamos oldalak (Gamma1, Gamma2) diszkrét pontjainak száma //RES VAN HELYETTE HASZNÁLVA NEM IS KELLENE IDEÍRNI
R = 0.3; % A fémhenger sugara [mm]
c_R = 1/3; % A vonaltöltés-kör sugarának aránya a fémhenger sugarához [1]
c_B = 4/3; % A külső töltések négyzetének oldalának aránya a cella oldalához
deltaX = 1; % Cella X irányú mérete [mm]
deltaY = 1; % Cella Y irányú mérete [mm]
V = 1; % Kondenzátorra kapcsolt feszültség [V]

Resolution = 100; % Hány pontra osztjuk a deltaX hosszt. [1]
N1 = 10; % Cellák száma X irányban [1]
N2 = 10; % Cellák száma Y irányban [1]

M = 16; % Fémrúdon belüli vonaltöltések száma /Metal/
B = 200; % Oldalankénti külső töltések száma /Boundary/

% chargeWeight = 1; % Ennyivel szorozzuk a töltéseket a mátrixszámolásnál

% "Származtatott" konstansok

NoC = 4*B+M; % Number of Charges
r_0 = c_R*R; % A referenciapont távolsága a vonaltöltés körtől [mm]
d = deltaX*N1; % A kapacitás fegyverzeteinek távolsága [mm]
h = deltaY*N2; % A kapacitás "magassága" [mm]

% Tehát más mértékegységek:
% Térerősség:
end