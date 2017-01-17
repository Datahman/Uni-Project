%  The vvv1x, vvv1y are really the Fermi velocity multiplied by hbar. 

% BMagnetic is the perpendicular magnetic field in Tesla
% PRB 88 195309 and Supplementary Material to PRL 110 246403 both set 
% deltab = 0.5  micro eV. deltab breaks rotational symmetry of
% the Fermi surface.
function [excitonHamExt] = excitonHamiltonian(kf, thetak, BMagnetic, deltab) 


   Norbitals = 4;

   % these are parameters of the Hamiltonian
mub = 57.88;  % micro eV / Tesla, the Bohr magneton
me =  0.511 * 10^12; % mass of the electron * c^2, in micro eV
hbarc = 1.24 * 10^6 / (2 * pi); %micro eV * micro m

% value from Supplementary Material to PRL 110 246403, Fig. S2
mex = 0.21 * me;  % exciton mass from PRB 88 125307, in rough agreement with Szymanska and Littlewood

% values from Supplementary Material to PRL 110 246403
% k_{ex} = 15.4 / micro m
% E(k_{ex}) = (hbarc k_{ex})^2 / (2 me) =  9.04 micro eV
% T = 0.1 K
gemub = -0.01 * mub; % the electronic g factor, no units, gemub * 7Tesla = -4.05
ghmub = 0.0085 * mub; % the hole g factor, no units, ghmub * 7Tesla = 3.45
% values from PRB 88 195309 and from Supplementary Material to PRL 110 246403
alphae = 0; % electron Rashba strength, micro eV * micro m
betae = 2.7; %electron Dresselhaus strength, micro eV * micro m, betae * kf = 42 with Butov's numbers
alphah = 0;  % hole Rashba strength, micro eV * micro m
betah = 0.92; % hole Dresselhaus strength, micro eV * micro m, betah * kf = 14 with Butov's numbers
gammah = 0.0;  % hole cubic strength
EbmEd = 5; %  micro eV
deltad = -13; %  micro eV

%todo
% betae = 0;
% betah = 0.3;
%  EbmEd =  0;  % Sx etc. proportional to 1/this between 1 and 0.1, 
 %when smaller it's almost independent
% deltad = 0;
% deltab = 0;  % Sx etc. proportional to this.
% gemub = 0;
% ghmub= 0;
% this is the phase which is associated with the momentum
phasef = cos(thetak) + (1i *  sin(thetak));
phasefc = cos(thetak) - (1i *  sin(thetak));

% set up the Hamiltonian
hhh1 = zeros(Norbitals, Norbitals);
vvv1x = zeros(Norbitals, Norbitals);
vvv1y = zeros(Norbitals, Norbitals);
for n=1:Norbitals
hhh1(n,n)= ((hbarc * kf)^2 / (2 * mex));
vvv1x(n,n) = hbarc * (hbarc * kf) * cos(thetak) / mex;
vvv1y(n,n) = hbarc * (hbarc * kf) * sin(thetak) / mex;
end
% 1 state = +e, -h, -1, bright
% 2 state = -e, +h, +1, bright
% 3 state = -e, -h, -2, dark
% 4 state = +e, +h, +2, dark

% 4->1, 1->2, 2->3, 3->4
% 1 state = +e, +h, +2, dark
% 2 state = +e, -h, -1, bright
% 3 state = -e, +h, +1, bright
% 4 state = -e, -h, -2, dark
% - 1-6, - 3-8, - 6-1, - 8-3  (2,4,5,7 are inactive)

% 1->2, 2->3, 3->6, 4->7, 5->4, 6->1, 7->8, 8->5, 
% 9->10, 10->11, 11->14, 12->15, 13->12, 14->9,  15->16, 16-> 13
% x + (y-1)*4
% 1 =6 = (2,2) = (+-,+-)
% 2 = 1 = (1,1) = (++,++)
% 3= 2 = (2,1) =  (+-,++)
% 4= 5 = (1,2)= (++,+-)
% 
% 5=8 = (4,2) = (--,+-)
% 6 =3 = (3,1) =(-+,++)
% 7=4= (4,1) = (--,++)
% 8=7 = (3,2) = (-+,+-)
% 
% 9=14= (2,4) = (+-,--)
% 10=9= (1,3) = (++,-+)
% 11=10=(2,3) = (+-,-+)
% 12=13=(1,4) = (++,--)
% 
% 13=16=(4,4) = (--,--)
% 14=11=(3,3) = (-+,-+)
% 15=12=(4,3) = (--,-+)
% 16=15=(3,4) = (-+,--)

ttt = zeros(4,4);
ttt(4,1) = 1;
ttt(1,2) = 1;
ttt(2,3) = 1;
ttt(3,4) = 1;
hhh1(1,1) = hhh1(1,1) + (gemub * BMagnetic  / 2);
hhh1(4,4) = hhh1(4,4) + (gemub * BMagnetic  / 2);
hhh1(2,2) = hhh1(2,2) - (gemub * BMagnetic  / 2);
hhh1(3,3) = hhh1(3,3) - (gemub * BMagnetic  / 2);
hhh1(1,1) = hhh1(1,1) - ghmub * BMagnetic /2;
hhh1(2,2) = hhh1(2,2) + ghmub * BMagnetic /2;
hhh1(3,3) = hhh1(3,3) - ghmub * BMagnetic /2;
hhh1(4,4) = hhh1(4,4) + ghmub * BMagnetic /2;
hhh1(1,1) = hhh1(1,1) + EbmEd/2;
hhh1(2,2) = hhh1(2,2) + EbmEd/2;
hhh1(3,3) = hhh1(3,3) - EbmEd/2;
hhh1(4,4) = hhh1(4,4) - EbmEd/2;
hhh1(1,2) = -deltab;
hhh1(2,1) = -deltab;
hhh1(3,4) = -deltad;
hhh1(4,3) = -deltad;

hhh1(1,3) = (alphae * phasef * kf) + (betae * phasefc * kf);
vvv1x(1,3) = alphae + betae;
vvv1y(1,3) = 1i * (alphae - betae);

hhh1(2,4) = (alphae * phasefc * kf) + (betae * phasef * kf);
vvv1x(2,4) = alphae + betae;
vvv1y(2,4) = 1i * (-alphae + betae);

hhh1(3,1) = (alphae * phasefc * kf) + (betae * phasef * kf);
vvv1x(3,1) = alphae + betae;
vvv1y(3,1) = 1i * (-alphae + betae);

hhh1(4,2) = (alphae * phasef * kf) + (betae * phasefc * kf);
vvv1x(4,2) = alphae + betae;
vvv1y(4,2) = 1i * (alphae - betae);

hhh1(1,4) = (alphah * phasef * kf) + (betah * phasefc * kf) + (gammah * (phasef * kf)^3);
vvv1x(1,4) = alphah + betah + (3 * gammah * (phasef * kf)^2);
vvv1y(1,4) = 1i * (alphah - betah + (3 * gammah * (phasef * kf)^2));

hhh1(2,3) = (alphah * phasefc * kf) + (betah * phasef * kf) + (gammah * (phasefc * kf)^3);
vvv1x(2,3) = alphah + betah + (3 * gammah * (phasef * kf)^2);
vvv1y(2,3) = 1i * (-alphah + betah - (3 * gammah * (phasef * kf)^2));

hhh1(3,2) = (alphah * phasef * kf) + (betah * phasefc * kf) + (gammah * (phasef * kf)^3);
vvv1x(3,2) = alphah + betah + (3 * gammah * (phasef * kf)^2);
vvv1y(3,2) = 1i * (alphah - betah + (3 * gammah * (phasef * kf)^2));

hhh1(4,1) = (alphah * phasefc * kf) + (betah * phasef * kf) + (gammah * (phasefc * kf)^3);
vvv1x(4,1) = alphah + betah + (3 * gammah * (phasef * kf)^2);
vvv1y(4,1) = 1i * (-alphah + betah - (3 * gammah * (phasef * kf)^2));

     excitonHamExt = zeros(Norbitals, Norbitals, 3);
     excitonHamExt(:,:,1) = ttt' * hhh1 * ttt;
     excitonHamExt(:,:,2) = ttt' * vvv1x * ttt;
     excitonHamExt(:,:,3) = ttt' * vvv1y * ttt;
    
   end