function [c, ceq, GradC, GradCeq] = geodesic_constr(vars,Aeq,beq)

c = [];
ceq = Aeq*vars - beq;
GradC = [];
GradCeq = Aeq';

end