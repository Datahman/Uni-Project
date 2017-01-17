function[xyz] = dosplot(Norbitals,BMagnetic,deltab,Efmin,Efmax,NoEfs,kmax,NumGridPoint,MaxNumEf,NumValueTheta,tau)
DinverseArray = zeros(Norbitals * Norbitals, Norbitals * Norbitals,NoEfs); % Diffuson operator 16x16

EfArray = zeros(NoEfs,1); % Array to store Efs.
DOSAtEfArray = zeros(NoEfs,1); % Array to store DOS at each Ef.
kfatzero=zeros(MaxNumEf,NoEfs);
kfatpi2=zeros(MaxNumEf,NoEfs);
KfAtThetaPi2=zeros(MaxNumEf,NoEfs);
RhoN = zeros(Norbitals, Norbitals, NoEfs);

RhoNEigs = zeros(Norbitals,  NoEfs);
DinverseEigs = zeros(Norbitals * Norbitals, NoEfs);

for Efindex = 1: NoEfs
    Efindex
    Ef = Efmin + ((Efmax - Efmin) / (NoEfs - 1)) * (Efindex - 1); 
    % Ef = Efmin + ((Efmax - Efmin) / (NoEfs -1)) * (r -1);
 EfArray(Efindex) = Ef;
 
 [DOSAtEfArray(Efindex),kfatzero(:,Efindex),KfAtThetaPi2(:,Efindex),RhoN(:,:,Efindex), DinverseArray(:,:,Efindex)]=  SimpleHintegrate(Norbitals,BMagnetic,deltab,Ef,kmax,NumGridPoint,MaxNumEf,NumValueTheta,tau);
%DinverseArray
 kfatpi2(:,Efindex) = sort(KfAtThetaPi2(:,Efindex));
     % kfatpi2
 %kfatpi(:,r)=sort(kfatpi2(:,r));
   kfatzero(:,Efindex)= sort(kfatzero(:,Efindex));
   [vv1, dd1 ] = eig(RhoN(:,:,Efindex));
   RhoNEigs(:,Efindex) = sort(real(diag(dd1)));
   [vv1, dd1 ] = eig(DinverseArray(:,:,Efindex));
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
plot(EfArray, DOSAtEfArray);
hold on;
plot(EfArray, RhoNEigs(1,:)+RhoNEigs(2,:)+RhoNEigs(3,:)+RhoNEigs(4,:), 'k');
hold on;

xlabel('Energy \mu eV ');
ylabel('DOS');
hold on;
axis([Efmin, Efmax,0,150]);


subplot(2,1,2);
hold off;
for i = 1:Norbitals*Norbitals
plot(EfArray, DinverseEigs(i,:));
hold on;
end
xlabel('Energy \mu eV ');
ylabel('Dinverse');
hold on;
axis([Efmin, Efmax,-0.2,1.2]);


% hold off;
% for i = 1:4
% plot(EfArray,real(kfatpi2(i,:))); %Blue
% hold on;
% end
% for i = 1:4
% plot(EfArray,real(kfatzero(i,:)), 'r'); %Blue
% hold on;
% end
% 
% xlabel('Energy \mu eV ');
% ylabel(' k at $\theta = \frac{\pi}{2}$ ', 'Interpreter', 'Latex');
% axis([Efmin,Efmax,-5,40]);
% 
    end

