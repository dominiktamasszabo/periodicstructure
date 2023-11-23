function [eps_r, eps_0, M, phi_0, K, R, c, deltaX, deltaY, r_0, N1, N2, d, h, Resolution, V, chargeWeight] = defineConstants()

% Konstansok
eps_r = 1; % Relatív permittivitás [1]
eps_0 = 8.86; % Vákuum permittivitás [F/mm*10^12] !!! EZ NEM BIZTOS HOGY JÓ !!!!
M = 4; % Vonaltöltések száma / cella MOST NEM KELL
phi_0 = 2*pi/17; % Első vonaltöltés szöge [rad]
K = 40; % Y-ra párhuzamos oldalak (Gamma1, Gamma2) diszkrét pontjainak száma //RES VAN HELYETTE HASZNÁLVA 
R = 0.1; % A fémhenger sugara [mm]
c = 1/3; % A vonaltöltés-kör sugarának aránya a fémhenger sugarához [1]
deltaX = 1; % Cella X irányú mérete [mm]
deltaY = 1; % Cella Y irányú mérete [mm]
V = 1; % Kondenzátorra kapcsolt feszültség [V]

Resolution = 15; % Hány pontra osztjuk a deltaX hosszt. [1]
N1 = 1; % Cellák száma X irányban [1]
N2 = 1; % Cellák száma Y irányban [1]

chargeWeight = 1; % Ennyivel szorozzuk a töltéseket a mátrixszámolásnál

% "Származtatott" konstansok

r_0 = c*R; % A referenciapont távolsága a vonaltöltés körtől [mm]
d = deltaX*N1; % A kapacitás fegyverzeteinek távolsága [mm]
h = deltaY*N2; % A kapacitás "magassága" [mm]

% Tehát más mértékegységek:
% Térerősség:
end