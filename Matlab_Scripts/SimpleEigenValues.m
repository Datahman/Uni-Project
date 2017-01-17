function [excitonEig] = SimpleEigenValues(kf, thetak, BMagnetic, deltab, WhichEig, Ef)

 Norbitals = 2;
 
 hhh1 = zeros(Norbitals, Norbitals,3);

        hhh1 = simpleHamiltonian(kf, thetak, BMagnetic, deltab);
        % diagonalize the Hamiltonian
        [vv,dd] = eig(hhh1(:,:,1));
        % save the eigenvalues in eigarray
        eigarray = sort(diag(dd));

     excitonEig = eigarray(WhichEig) - Ef;

 end