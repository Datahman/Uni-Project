% this function just does a loop over energy Ef, calling SimpleHIntegrate
% each time.  It puts the results into arrays and plots them.

function[xyz] = dosplot(Norbitals,BMagnetic,deltab,Efmin,Efmax,NoEfs,kmax,NumKGridPoints,MaxNumkf,NumValueTheta,hbarOverTau)
LinearLength=1;
kspacegrid = zeros(LinearLength, LinearLength,3);
DinverseArray = zeros(Norbitals * Norbitals, Norbitals * Norbitals,LinearLength, LinearLength, NoEfs); % Diffuson operator 16x16

EfArray = zeros(NoEfs,1); % Array to store Efs.
DOSAtEfArray = zeros(NoEfs,1); % Array to store DOS at each Ef.
kfatzero=zeros(MaxNumkf,NoEfs);
KfAtThetaPi2=zeros(MaxNumkf,NoEfs);
RhoN = zeros(Norbitals, Norbitals, NoEfs);

RhoNEigs = zeros(Norbitals,  NoEfs);
DinverseEigs = zeros(Norbitals * Norbitals, NoEfs);

for Efindex = 1: NoEfs
    Ef = Efmin + ((Efmax - Efmin) / (NoEfs - 1)) * (Efindex - 1); 
    % Ef = Efmin + ((Efmax - Efmin) / (NoEfs -1)) * (r -1);
 EfArray(Efindex) = Ef;
 
 % this function calculates the scalar and matrix DOS, and Dinverse, after integration
% over the Fermi surface.
% It also calculates the Fermi momenta at theta = 0 and theta = pi/2. 
% It does this at particular values of Ef, tau, and the Hamiltonian parameters.  
 [DOSAtEfArray(Efindex),kfatzero(:,Efindex),KfAtThetaPi2(:,Efindex),RhoN(:,:,Efindex), DinverseArray(:,:,:,:,Efindex)]=  SimpleHintegrate(Norbitals,BMagnetic,deltab,Ef,kmax,NumKGridPoints,MaxNumkf,NumValueTheta,hbarOverTau,kspacegrid, LinearLength);

    kfatzero(:,Efindex)= sort(kfatzero(:,Efindex))
    KfAtThetaPi2(:,Efindex)= sort(KfAtThetaPi2(:,Efindex))
   [vv1, dd1 ] = eig(RhoN(:,:,Efindex));
   RhoNEigs(:,Efindex) = sort(real(diag(dd1)));
   [vv1, dd1 ] = eig(DinverseArray(:,:,1,1,Efindex));
   DinverseEigs(:,Efindex) = sort(real(diag(dd1)));
end

subplot(2,1,1);
%set('defaulttextinterpreter', 'latex');
%plot(EfArray,DOSAtEfArray);
%hold on;
hold off;
for i = 1:Norbitals
plot(EfArray, RhoNEigs(i,:));
hold on;
end
plot(EfArray, DOSAtEfArray/4, 'k');
hold on;

xlabel('Energy \mu eV ');
ylabel('DOS');
hold on;
axis([Efmin, Efmax,0,50]);


subplot(2,1,2);
hold off;
for i = 1:Norbitals*Norbitals
plot(EfArray, DinverseEigs(i,:));
hold on;
end
title('$\tau = 10$','interpreter', 'latex');
xlabel('Energy $\mu$ eV ', 'interpreter', 'latex');
ylabel('Decay rate $D^{-1}$ ', 'interpreter', 'latex');
hold on;
axis([Efmin, Efmax,-0.2, 1.5]);
%Parametrs used.. dosplot(4,0,0.5,10,20,8,20,100,8,100,1) constant DOS..
%expected?
    end

