function F = ParametricEigenvalueAssignment(A_0,A_1,A_2,B_0,C_0,C_1,C_2,L_1,L_2,Q,Z)
% Solves for gain matrix F that gives eigenvalues L_1,L_2 for the system
% A_2 ddx + A_1 dx + A_0 x = -B_0 F (C_2 ddx + C_1 dx + C_0 x)
% using parameter sets Q and Z to pick from the available eigenspaces.

if length([L_1(:);L_2(:)]) ~= length(unique([L_1(:);L_2(:)]))
    warning('Warning. The desired eigenvalue sets L_1 and L_2 contain duplicate eigenvalues. This algorithm is written for distinct eigenvalues.')
end

CT_0 = C_0-C_2/A_2*A_0;
CT_1 = C_1-C_2/A_2*A_1;

temp = size(B_0);
n = temp(1);
m = temp(2);
temp = size([CT_0, CT_1]);
p = temp(1);
if rank([CT_0, CT_1]) < p
    warning('Warning. Sensors are linearly dependend. This is probably not going to work.')
end
r = p-1;
if length(L_1) ~= r
    warning(['Warning. Wrong length of first eigenvalue set. L_1 needs',num2str(r),'entries.'])
end
if length(L_2) ~= 2*n-r
    warning(['Warning. Wrong length of second eigenvalue set. L_2 needs',num2str(2*n-r),'entries.'])
end
rho = 2*n-r;

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
if temp(1) ~= m
    warning(['Warning. Wrong length of parameters q_i. Each q_i needs',num2str(m),'entries, which should be the number of rows of Q.'])
end
if temp(2) ~= r
    warning(['Warning. Wrong number of parameters q_i. There should be',num2str(r),'q_i, which should be the number of columns of Q. Beware that for complex conjugate eigenvalues there should be only one q_i, since the second will be derived from the first.'])
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
if temp(1) ~= p-r
    warning(['Warning. Wrong length of parameters z_i. Each z_i needs',num2str(p-r),'entries, which should be the number of rows of Z.'])
end
if temp(2) ~= rho
    warning(['Warning. Wrong number of parameters z_i. There should be',num2str(rho),'z_i, which should be the number of columns of Z. Beware that for complex conjugate eigenvalues there should be only one z_i, since the second will be derived from the first.'])
end

%% algorithm
% build F0
if r == 0
    F0 = 0;
    U1Apostrophe = 1;
    C1 = [CT_0, CT_1];
else
    RightNull = null([L_1(1)^2*A_2+L_1(1)*A_1+A_0, B_0]);
    Ni = RightNull(1:n,1:m);
    Mi = RightNull(n+1:n+m,1:m);
    qi = Q(:,1);
    Vr = [Ni*qi;L_1(1)*Ni*qi];
    Qr = Mi*qi;
    for i = 2:r
        RightNull = null([L_1(i)^2*A_2+L_1(i)*A_1+A_0, B_0]);
        qi = Q(:,i);
        Ni = RightNull(1:n,1:m);
        Mi = RightNull(n+1:n+m,1:m);
        Vr = [Vr, [Ni*qi;L_1(i)*Ni*qi]];
        Qr = [Qr, Mi*qi];
    end
    F0 = Qr*pinv([CT_0, CT_1]*Vr);
    U1Apostrophe = null(([CT_0, CT_1]*Vr)')';
    C1 = U1Apostrophe*[CT_0, CT_1];
end

% build F1
if rho == 0
    F1 = 0;
else
    A1_0 = A_0+B_0*F0*CT_0;
    A1_1 = A_1+B_0*F0*CT_1;
    A1_2 = A_2;
    B1_0 = B_0;
    C1_0 = U1Apostrophe*CT_0;
    C1_1 = U1Apostrophe*CT_1;
    C1 = [C1_0, C1_1];
    temp = size(C1);
    p1 = temp(1);
    LeftNull = null([-L_2(1)*eye(n), eye(n); -A1_0, -L_2(1)*A1_2-A1_1;C1]')';
    Dj1 = LeftNull(1:p1,1:n);
    Dj2 = LeftNull(1:p1,n+1:2*n);
    Ej = LeftNull(1:p1,2*n+1:2*n+p1);
    zj = Z(:,1);
    Wrho = Dj2'*zj;
    Zrho = Ej'*zj;
    for j = 2:rho
        LeftNull = null([-L_2(j)*eye(n), eye(n); -A1_0, -L_2(j)*A1_2-A1_1;C1]')';
        zj = Z(:,j);
        Dj1 = LeftNull(1:p1,1:n);
        Dj2 = LeftNull(1:p1,n+1:2*n);
        Ej = LeftNull(1:p1,2*n+1:2*n+p1);
        Wrho = [Wrho, Dj2'*zj];
        Zrho = [Zrho, Ej'*zj];
    end
    F1 = -pinv(Wrho'*B1_0) *Zrho';
end

% calculate F
F = F0+F1*U1Apostrophe;
F = Feedback2Feedthrough(F,-C_2/A_2*B_0);
end
