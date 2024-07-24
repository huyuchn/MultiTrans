
clear;
clc;

load('kanatani_data.mat');

epsilon = 1e-8;

n = 8;

Cor(:,:,1) = xyz_2010 -  ones(n,1)*mean(xyz_2010);
Cor(:,:,2) = xyz_2011 -  ones(n,1)*mean(xyz_2011);
Cor(:,:,3) = xyz_2012 -  ones(n,1)*mean(xyz_2012);

%-------------------------------------%
% Combined Adjustment                 %
%-------------------------------------%

Qc = blkdiag(Q0_2010,Q0_2011,Q0_2012);

[x_c, ~, ~, ecor_c] = MultiTrans(Cor,Qc,epsilon);

RD_c = reshape(x_c,3,8);

Rc1 = RD_c(:, 1:3);
Dc1 = RD_c(:, 4);
Rc2 = RD_c(:, 5:7);
Dc2 = RD_c(:, 8);

%-------------------------------------%
% Separate Adjustment                 %
%-------------------------------------%

[x_s1, ~, ~, ecor_s1] = MultiTrans(Cor(:,:,1:2),...
    blkdiag(Q0_2010, Q0_2011),epsilon);
[x_s2, ~, ~, ecor_s2] = MultiTrans(Cor(:,:,2:3),...
    blkdiag(Q0_2011, Q0_2012),epsilon);

RD_s1 = reshape(x_s1,3,4);
RD_s2 = reshape(x_s2,3,4);

Rs1 = RD_s1(:, 1:3);
Ds1 = RD_s1(:, 4);
Rs2 = RD_s2(:, 1:3);
Ds2 = RD_s2(:, 4);

