function [x,EA,ey,ecor] = MultiTrans(Cor,Q,epsilon)

p = size(Cor, 1);
w = size(Cor, 3);

A = [];
y = [];
for i = 1:w-1
    Ai = [kron(eye(3),Cor(:,:,i)) kron(eye(3),ones(p,1))];
    A = blkdiag(A,Ai);
    y = [y; reshape(Cor(:,:,i+1),3*p,1)];
end


x = inv(A'*A)*A'*y;

RD_a = reshape(x,3,(w-1)*4);

iter = 0;
MAX_ITER = 500;


EA = zeros(3*p*(w-1),12*(w-1));

while (1)
    Bq = [];
    for i = 1:w-1
        Ri = RD_a(:, 4*(i-1)+1 : 4*(i-1)+3);
        Bi = [zeros(3,3*(i-1)) -Ri' eye(3) zeros(3,3*(w-1-i))];
        Bq = [Bq;Bi];
    end
    
    B = kron(Bq,eye(p));
    QB = B*Q*B';
    
    iQB = inv(QB);
    
    ecor = Q*B'*iQB*(y - A*x);
    Ecor = reshape(ecor, p, 3*w);
    
    EA = [];
    ey = [];
    for i = 1:w-1
        EAi = [kron(eye(3), Ecor(:,3*(i-1)+1:3*i)) zeros(3*p,3)];
        EA = blkdiag(EA,EAi);
        ey = [ey; reshape(Ecor(:,3*i+1:3*i+3),3*p,1)];
    end
    
    Ah = A - EA;
    
    N = inv(Ah'*iQB*Ah);
    xnew = N*Ah'*iQB*(y - EA*x);
    
    
    delta = xnew-x;
    
    iter = iter + 1;
    
    x = xnew;
    RD_a = reshape(x,3,(w-1)*4);
    
    if (norm(delta)<epsilon)
        break;
    end
    
    if (iter>MAX_ITER)
        disp('-----------------');
        disp('Multitrans fails!');
        disp(norm(delta));
        disp('-----------------');
        break;
    end
    
end


end