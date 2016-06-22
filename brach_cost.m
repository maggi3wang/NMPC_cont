function [J, gradJ] = brach_cost(xt,N)

J = xt(end);
gradJ = [zeros(2*(N+1),1);1];
end