function gain = ParametricEigenvalueAssignmentGain(A_0,A_1,A_2,B_0,C_0,C_1,C_2,L_1,L_2,Q,Z)
% Used for optimization

F = ParametricEigenvalueAssignment(A_0,A_1,A_2,B_0,C_0,C_1,C_2,L_1,L_2,Q,Z);
gain = norm([B_0*F*C_0, B_0*F*C_1, B_0*F*C_2],'fro');
end
