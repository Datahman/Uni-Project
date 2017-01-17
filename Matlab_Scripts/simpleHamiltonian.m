function [excitonHamExt] = SimpleHamiltonian(kf, thetak, BMagnetic, deltab)
Norbitals = 2;
hhh1 = zeros(Norbitals, Norbitals);
vvv1x = zeros(Norbitals, Norbitals);
vvv1y = zeros(Norbitals, Norbitals);
kx = kf * cos(thetak);
ky = kf * sin(thetak);
hhh1(1,1) = kx * kx + (ky + BMagnetic)* (ky + BMagnetic);
hhh1(2,2) = kx * kx + (ky - BMagnetic)* (ky - BMagnetic) + 2;
vvv1x(1,1) = 2 * kx;
vvv1x (2,2) = 2 * kx;
vvv1y (1,1) = 2 * (ky + BMagnetic);
vvv1y (2,2) = 2 * (ky - BMagnetic);

excitonHamExt = zeros(Norbitals, Norbitals, 3);
     excitonHamExt(:,:,1) = hhh1;
     excitonHamExt(:,:,2) = vvv1x;
     excitonHamExt(:,:,3) = vvv1y;
end

