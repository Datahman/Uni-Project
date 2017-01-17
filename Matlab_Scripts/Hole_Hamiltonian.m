function [hole_hamilton] = Hole_Hamiltonian(kf, thetak,b_h)
Norbitals=2;



me =  0.511 * 10^12; % mass of the electron * c^2, in micro eV
hbarc = 1.24 * 10^6 / (2 * pi); %micro eV * micro m
mex = 0.21 * me; 
% a = 0.1815 (micro-ev (mu m )^2).
% b_h = 0.3.
%Hamiltonian setup

phasef = cos(thetak) + (1i *  sin(thetak));
phasefc = cos(thetak) - (1i *  sin(thetak));

ham = zeros (Norbitals,Norbitals);
devx = zeros(Norbitals, Norbitals); % derivative w.r.t x
devy = zeros(Norbitals, Norbitals); % derivative w.r.t y
vec = zeros(2,2);
vec(2,1) = 1;
vec(1,2) = 1;
for n=1:Norbitals
hhh1(n,n)= (hbarc * kf)^2/ (2* mex) ;
devx(n,n) = 2* ((hbarc  * kf)^2/ (2* mex))  * cos(thetak) ;
devy(n,n) = 2* ((hbarc  * kf)^2/ (2* mex))  * sin(thetak);
end
hhh1(1,1) = hhh1(1,1);
hhh1(1,2) = (b_h * kf * phasefc) ;

devx(1,2) = b_h ;
devy(1,2) =  -1i * b_h;

hhh1(2,1) = (b_h * kf * phasef) ;

devx(2,1) = b_h ;
devy(2,1) =  1i * b_h;

hhh1(2,2) = hhh1(1,1);

    hole_hamilton = zeros(Norbitals,Norbitals,3);
    hole_hamilton(:,:,1) = vec' * hhh1 * vec ;
    hole_hamilton(:,:,2) = vec' * devx * vec;
    hole_hamilton(:,:,3) =  vec' * devy * vec;
end
