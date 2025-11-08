function F = Feedback2Feedthrough(FTilde,D)
% Recalculate feedback gain for the case with a feedthrough matrx

temp = size(FTilde);
m = temp(1);
F = (eye(m)+FTilde*D)\FTilde;

end
