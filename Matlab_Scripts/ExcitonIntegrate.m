function ExcitonIntegrate  (Norbitals,BMagnetic,deltab,Ef,kmax,NumGridPoint,MaxNumEf,NumValueTheta)
FSNumArray = zeros(NumValueTheta -1 , 1);
thetaarray = zeros(NumValueTheta -1, 1);
theta = zeros(NumValueTheta -1,1);
r_1 = zeros(NumValueTheta -1,1 );
r_2 = zeros(NumValueTheta -1,1 );
r_3 = zeros(NumValueTheta -1,1 );
r_4 = zeros(NumValueTheta -1,1 );
k_Array = zeros(MaxNumEf, 3, NumValueTheta -1);
for n = 1: NumValueTheta -1 
    dthetak = (2* pi / (NumValueTheta - 1)) * (n-1);
    thetaarray(n) = dthetak;
    [FermiSurfaceNum,kfArray] = excitonFermiSurfaceA(Norbitals, dthetak,BMagnetic,deltab,Ef,kmax,NumGridPoint,MaxNumEf);
    size(kfArray)
    k_Array(:,:,n)  = kfArray;
    theta(n) = thetaarray(n);
    r_1(n) = k_Array(1,1,n);
    r_2(n) = k_Array(2,2,n);
    r_3(n) = k_Array(3,3,n);
    r_4(n) = k_Array(4,4,n);
end

squeeze(size(thetaarray(:,1)))
squeeze(size(k_Array(1,1,:)))
hold off;
polar(theta, r_1, '--r');
hold on;
polar(theta, r_2, '--b');
hold on;
polar(theta, r_3, '--b');
end
