clear;
clc;
close all;

load('kanatani_data.mat');

p = 8;
w = 3;

%-------------------------------------%
% Variance component estimation       %
%-------------------------------------%

sp = w;

Cor(:,:,1) = xyz_2010 -  ones(p,1)*mean(xyz_2010);
Cor(:,:,2) = xyz_2011 -  ones(p,1)*mean(xyz_2011);
Cor(:,:,3) = xyz_2012 -  ones(p,1)*mean(xyz_2012);


Q(:,:,1) = blkdiag(Q0_2010,zeros(3*p),zeros(3*p));
Q(:,:,2) = blkdiag(zeros(3*p),Q0_2011,zeros(3*p));
Q(:,:,3) = blkdiag(zeros(3*p),zeros(3*p),Q0_2012);

iter = 0;
epsilon = 1e-8;
MAX_ITER = 100;

sig = ones(sp,1);

A = [];
y = [];
for i = 1:w-1
    Ai = [kron(eye(3),Cor(:,:,i)) kron(eye(3),ones(p,1))];
    A = blkdiag(A,Ai);
    y = [y; reshape(Cor(:,:,i+1),3*p,1)];
end

SIG = [];
while (1)
    sigtest = sig;
    
    Qxyz = zeros(3*p*w);
    for i = 1:sp
        Qxyz = Qxyz + sig(i) * Q(:,:,i);
    end
    
    [x, EA, ey, ~] = MultiTrans(Cor, Qxyz, epsilon);
    
    RD_a = reshape(x,3,(w-1)*4);


    Ah = A - EA;
    
    en = ey - EA*x;
    
    %en = y - EA*x;
    
    Bq = [];
    for i = 1:w-1
        Ri = RD_a(:, 4*(i-1)+1 : 4*(i-1)+3);
        Bi = [zeros(3,3*(i-1)) -Ri' eye(3) zeros(3,3*(w-1-i))];
        Bq = [Bq;Bi];
    end
    
    Lb = kron(Bq,eye(p));
   
    QB = zeros(3*(w-1)*p);
    
    for i = 1:sp
        Qy(:,:,i) = Lb*Q(:,:,i)*Lb';
        QB = QB + sig(i)*Qy(:,:,i);
    end
    
    iQB = inv(QB);
    Pao = eye(3*p*(w-1)) - Ah*inv(Ah'*iQB*Ah)*Ah'*iQB;
    
    for k=1:sp
        %lp(k,1) = 0.5*en'*iQB*Qy(:,:,k)*iQB*en;
        lp(k,1) = 0.5*en'*iQB*Pao*Qy(:,:,k)*iQB*Pao*en;
        for j=1:sp
            N(k,j) = 0.5*trace(iQB*Pao*Qy(:,:,k)*iQB*Pao*Qy(:,:,j));
        end
    end
    
    sig = inv(N)*lp;
    %[sig,~] = nnls_v(N,lp);
    %[sig,~] = nnls_vc2(N,lp,Q);
    
    %if (max(abs(sigtest-sig))<1e-5) || iter>= MAX_ITER
    if (norm(sigtest-sig)<1e-8) || iter>= MAX_ITER
       break; 
    end
    
    iter = iter + 1;
    
    SIG(iter,:) = sig';
    
    disp(max(abs(sigtest-sig)));
end

% plot VC

SIG = [1 1 1;SIG];
n = size(SIG,1);

pSig = plot(1:n,sqrt(SIG));
xlim([1 n]);

ax = gca;
ax.XTick = unique( round(ax.XTick) );

set(pSig(1), 'color', [0.203921568627451,0.423529411764706,0.682352941176471],...
    'Linewidth',2,'Linestyle','-','Marker','^','Markersize',5);
set(pSig(2), 'color', [0.917647058823529,0.0274509803921569,0.498039215686275],...
    'Linewidth',2,'Linestyle','--','Marker','^','Markersize',5);
set(pSig(3), 'color', [0.396078431372549,0.400000000000000,0.396078431372549],...
    'Linewidth',2,'Linestyle','-.','Marker','^','Markersize',5);
