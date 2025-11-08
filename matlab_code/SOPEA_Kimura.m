function F = SOPEA_Kimura(A_0,A_1,A_2,B_0,C_0,C_1,C_2,L_1,L_2,Q,Z)
% Solves for gain matrix F that gives eigenvalues L_1,L_2 for the system
% A_2 ddx + A_1 dx + A_0 x = -B_0 F (C_2 ddx + C_1 dx + C_0 x)
% using parameter sets Q and Z to pick from the available eigenspaces.

if length([L_1(:);L_2(:)]) ~= length(unique([L_1(:);L_2(:)]))
    warning('Warning. The desired eigenvalue sets L_1 and L_2 contain duplicate eigenvalues. This algorithm is written for distinct eigenvalues.')
end

C = [C_0, C_1, C_2];

temp = size(B_0);
n = temp(1);
m = temp(2);
temp = size(C);
p = temp(1);

if rank(C) < p
    warning('Warning. Sensors are linearly dependend. This is probably not going to work.')
end
if rank(B_0) < m
    warning('Warning. Actuators are linearly dependend. This is probably not going to work.')
end

r = p-1;
if length(L_1) ~= r
    warning(['Warning. Wrong length of first eigenvalue set. L_1 needs r=p-1=',num2str(r),' entries.'])
end
rho = 2*n-r;
if length(L_2) ~= rho
    warning(['Warning. Wrong length of second eigenvalue set. L_2 needs rho=2n-r=2n-p+1=',num2str(2*n-r),' entries.'])
end

%% build complex conjugate parameter vectors
temp = size(Q);
if temp(2) < length(L_1)
    idx = abs(imag(L_1))>eps;
    didSomething = 0;
    for i = 1:length(L_1)
        if didSomething
            didSomething = 0;
            continue
        end
        if idx(i)
            Q = [Q(:,1:i),conj(Q(:,i)),Q(:,i+1:end)];
            didSomething = 1;
        end
    end
end

temp = size(Q);
if temp(2) ~= r
    warning(['Warning. Wrong number of parameters q_i. There should be ',num2str(r),' q_i, which should be the number of columns of Q. Beware that for complex conjugate eigenvalues there should be only one q_i, since the second will be derived from the first.'])
else
    if r ~=0 & temp(1) ~= m
    warning(['Warning. Wrong length of parameters q_i. Each q_i needs ',num2str(m),' entries, which should be the number of rows of Q.'])
    end
end

temp = size(Z);
if temp(2) < length(L_2)
    idx = abs(imag(L_2))>eps;
    didSomething = 0;
    for i = 1:length(L_2)
        if didSomething
            didSomething = 0;
            continue
        end
        if idx(i)
            Z = [Z(:,1:i),conj(Z(:,i)),Z(:,i+1:end)];
            didSomething = 1;
        end
    end
end

temp = size(Z);
if temp(2) ~= rho
    warning(['Warning. Wrong number of parameters z_i. There should be ',num2str(rho),' z_i, which should be the number of columns of Z. Beware that for complex conjugate eigenvalues there should be only one z_i, since the second will be derived from the first.'])
else
    if rho ~= 0 & temp(1) ~= p-r
    warning(['Warning. Wrong length of parameters z_i. Each z_i needs ',num2str(p-r),' entries, which should be the number of rows of Z.'])
    end
end

%% algorithm
% build F1
if r == 0
    F1 = 1;
    U1Apostrophe = 1;
    %C1 = C;
else
    RightNull = null([L_1(1)^2*A_2+L_1(1)*A_1+A_0, B_0]);
    Ni = RightNull(1:n,1:m);
    Mi = RightNull(n+1:n+m,1:m);
    qi = Q(:,1);
    Vr = (L_1(1)^2*C_2+L_1(1)*C_1+C_0)*Ni*qi;
    Qr = Mi*qi;
    for i = 2:r
        RightNull = null([L_1(i)^2*A_2+L_1(i)*A_1+A_0, B_0]);
        qi = Q(:,i);
        Ni = RightNull(1:n,1:m);
        Mi = RightNull(n+1:n+m,1:m);
        Vr = [Vr, (L_1(i)^2*C_2+L_1(i)*C_1+C_0)*Ni*qi];
        Qr = [Qr, Mi*qi];
    end
    F1 = Qr*pinv(Vr);
    U1Apostrophe = null((Vr)')';
    %C1 = U1Apostrophe*C;
end

% build F2
if rho == 0
    F2 = 0;
else
    A1_0 = A_0+B_0*F1*C_0;
    A1_1 = A_1+B_0*F1*C_1;
    A1_2 = A_2+B_0*F1*C_2;
    B1_0 = B_0;
    C1_0 = U1Apostrophe*C_0;
    C1_1 = U1Apostrophe*C_1;
    C1_2 = U1Apostrophe*C_2;
    C1 = [C1_0, C1_1, C1_2];
    temp = size(C1);
    p1 = temp(1);
    LeftNull = null([L_2(1)^2*A1_2+L_2(1)*A1_1+A1_0; L_2(1)^2*C1_2+L_2(1)*C1_1+C1_0]')';
    Dj = LeftNull(1:p1,1:n);
    %Dj2 = LeftNull(1:p1,n+1:2*n);
    Ej = LeftNull(1:p1,n+1:n+p1);
    zj = Z(:,1);
    Wrho = Dj'*zj;
    Zrho = Ej'*zj;
    for j = 2:rho
        LeftNull = null([L_2(j)^2*A1_2+L_2(j)*A1_1+A1_0; L_2(j)^2*C1_2+L_2(j)*C1_1+C1_0]')';
        zj = Z(:,j);
        Dj = LeftNull(1:p1,1:n);
        %Dj2 = LeftNull(1:p1,n+1:2*n);
        Ej = LeftNull(1:p1,n+1:n+p1);
        Wrho = [Wrho, Dj'*zj];
        Zrho = [Zrho, Ej'*zj];
    end
    F2 = pinv(Wrho'*B1_0) *Zrho';% - minus?
end

% calculate F
F = F1+F2*U1Apostrophe;
%F = Feedback2Feedthrough(F,-C_2/A_2*B_0);
end
