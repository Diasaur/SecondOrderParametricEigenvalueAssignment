function imaginary_part = SOPEA_Imag(A_0,A_1,A_2,B_0,C_0,C_1,C_2,L_1,L_2,L_3,Q,Z)
% Used for optimization
F = SOPEA(A_0,A_1,A_2,B_0,C_0,C_1,C_2,L_1,L_2,L_3,Q,Z);
imaginary_part = max(max(imag(F)));
end
