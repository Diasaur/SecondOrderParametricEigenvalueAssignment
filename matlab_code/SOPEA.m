function F = SOPEA(A_0,A_1,A_2,B_0,C_0,C_1,C_2,L_1,L_2,L_3,Q,Z)
% Used for optimization
if all(size(L_3) == 0)
    F = SOPEA_Kimura(A_0,A_1,A_2,B_0,C_0,C_1,C_2,L_1,L_2,Q,Z);
else
    F = SOPEA_KimuraMinus1(A_0,A_1,A_2,B_0,C_0,C_1,C_2,L_1,L_2,L_3,Q,Z);
end
end
