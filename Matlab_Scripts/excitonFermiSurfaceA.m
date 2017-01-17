% This function returns the count of Fermi surfaces in FermiSurfaceNum
% kfArray's first element is kf, its second element is an index into the
% array of k's, and its third element is the Fermi velocity at this kf and
% thetak. The fourth element of kfArray is an index into the array
% eigenvalues/eigenvectors.
% todo: test this function to see how it does when the spin-orbit splitting
% is very small.
function[FermiSurfaceNum,kfArray] = excitonFermiSurfaceA(Norbitals, thetak,BMagnetic,deltab,Ef,kmax,NumkGridPoints,MaxNumFS)
func = @(kf, thetak, BMagnetic, deltab, WhichEig, Ef) excitonEigenvalues(kf, thetak, BMagnetic, deltab, WhichEig, Ef);
%func = @(kf, thetak, BMagnetic, deltab, WhichEig, Ef) SimpleEigenValues(kf, thetak, BMagnetic, deltab, WhichEig, Ef);
%  Norbitals = 2;

% KValues is a 1-D grid of Fermi momenta, from k = 0 to k = kmax
 KValues = zeros(1,NumkGridPoints);
% EigArray, with 1 in the middle index, contains the eigenvalues of the Hamiltonian at each value of k
% EigArray, with 2 in the middle index, contains the numbering of the
% eigenvalues as produced by the matlab diagonalization routine
 EigArray = zeros(Norbitals, 2, NumkGridPoints);

 % hhh1 will contain the Hamiltonian and the two components of the velocity
 % operator
 hhh1 = zeros(Norbitals, Norbitals,3);
 
 % calculate the momenta and eigenvalues at each point on the grid
for kk=1:NumkGridPoints
    KValues(kk) = kmax * (1- kk)/(1-NumkGridPoints)  ; 
    % get the Hamiltonian and velocity operators
    hhh1 = excitonHamiltonian(KValues(kk), thetak, BMagnetic, deltab);
    %hhh1 = simpleHamiltonian(KValues(kk), thetak, BMagnetic, deltab)

  % get the eigenvalues
[vv,dd] = eig(hhh1(:,:,1))
% sort the eigenvalues
 SortThis = zeros(Norbitals, 2);
 SortThis(:,1) = diag(dd);
 SortThis(:,2) = squeeze([1:4]); % todo check this % was 4 for exciton.
EigArray(:,:,kk) = sort(SortThis);
end

% For each Fermi surface, this will store kf, an index into the Kvalues
% array, and the Fermi velocity, and which eigenvalue/eigenvector pair.
% MaxNumFS is a limit on how many Fermi surfaces it is allowed to find.
kfArray = zeros(MaxNumFS,6);

% initialize the search for Fermi surfaces - have not found any yet
FermiSurfaceNum = 0;


% loop over which eigenvalue we are looking at.  We have sorted them
% already.
for EigenLoop = 1:Norbitals 
    % loop starting at k=0 and going out to k=kmax, looking for intervals
    % where the eigenvalue becomes bigger or smaller than Ef
    for WhichValueK = 1:NumkGridPoints - 1
        
        % value of the eigenvalue at the beginning of the interval
        % also the second element index of the eigenvalue as produced as
        % produced by the diag function.
         x1 = squeeze(EigArray(EigenLoop,:,WhichValueK));
         
         % value of the eigenvalue at the end of the interval
        % also the second element index of the eigenvalue as produced as
        % produced by the diag function.
        x2 = squeeze(EigArray(EigenLoop,:, WhichValueK + 1));   

        % test for crossing Ef
        if (x1(1) < Ef && x2(1) > Ef) || (x1(1) > Ef && x2(1) < Ef)
            % We found a Fermi surface so increment our count.
            FermiSurfaceNum  = FermiSurfaceNum + 1; 
            
            % if the count is too high then crash with an error.
           if FermiSurfaceNum > MaxNumFS
              msgbox 'Too many Fermi Surfaces';
           end
           
           % interval1 is the interval of k inside which the eigenvalue
           % must cross Ef
           interval1 = [KValues(WhichValueK),KValues(WhichValueK+1)];
           
           % This uses matlab's fzero to find the exact value of kf to high
           % precision.
             func = @(kf, thetak, BMagnetic, deltab, WhichEig, Ef) excitonEigenvalues(kf, thetak, BMagnetic, deltab, WhichEig, Ef);
             %func = @(kf, thetak, BMagnetic, deltab, WhichEig, Ef) SimpleEigenValues(kf, thetak, BMagnetic, deltab, WhichEig, Ef);
            func1 = @(kf) func(kf, thetak, BMagnetic, deltab, EigenLoop, Ef);
            kzero = fzero(func1,interval1); 
            
            % kf is stored in the first element of kfArray.
            kfArray(FermiSurfaceNum,1) = kzero;
            
            % an index into the k array is the second element of kfArray
            kfArray(FermiSurfaceNum,2) = WhichValueK;
 
            % obselete: store the Fermi velocity as the third element of kfArray.
             kfArray(FermiSurfaceNum,3) = -1; % thisvf < ---- derivative dE/dk 

            % stores the index of the eigenvalue after sorting the
            % eigenvalues
            kfArray(FermiSurfaceNum,4) = EigenLoop;

        % Stores the index of the eigenvalue as produced as
        % produced by the diag function.
        %  If this is different at the two ends of the interval, we put -1.
        %  maybetodo:  We use the index at the ends of the interval, and
        %  don't check whether they match the index at the actual solution
        %  kzero.
            if x1(2) == x2(2)
                kfArray(FermiSurfaceNum,5) = x1(2);
            else
               kfArray(FermiSurfaceNum,5) = -1;
            end
            
            
      
         % store the degeneracy
         kfArray(FermiSurfaceNum, 6) = 1;
         
        end % end of block where the eigenvalue has gone through Ef
    end % end of loop from k=0 to k=kmax
                 
end % end of loop over eigenvalues

% loop starting at k=0 and going out to k=kmax, looking for degenerate
% Fermi surfaces.  If there is one, increase the duplicate of the first
% entry about that FS, and remove the duplicate from consideration by
% marking it with kfArray(iFS,6) = -1;
for WhichFS = 1:FermiSurfaceNum
    if kfArray(WhichFS,6) < 0
        continue;
    else
        for iFS = WhichFS+1:FermiSurfaceNum
            if abs(kfArray(WhichFS,1) - kfArray(iFS,1)) < 1e-12
                kfArray(iFS,6) = -1;
                kfArray(WhichFS,6) =  kfArray(WhichFS,6) + 1;
            end
        end
    end
end

% make a new kfArray which omits the duplicate degenerate states
kfArray1 = zeros(MaxNumFS,6);
FermiSurfaceNum1 = 0;
for WhichFS = 1:FermiSurfaceNum
    if kfArray(WhichFS,6) > 0
        FermiSurfaceNum1 = FermiSurfaceNum1 + 1;
        kfArray1(FermiSurfaceNum1,:) = kfArray(WhichFS,:);        
    end
end

FermiSurfaceNum = FermiSurfaceNum1;
kfArray = kfArray1;
%size(kfArray)% [7,6] for excitonFermiSurfaceA(4, 100,0,0.5,20,20,100,7)
end  % end of the function
