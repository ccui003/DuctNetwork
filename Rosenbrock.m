function [ y, dydx ] = Rosenbrock( x )
%ROSENBROCK Summary of this function goes here
%   Detailed explanation goes here
y = [1-x(1);10*(x(2)-x(1)^2)];
dydx = [-1,0;10,-20*x(1)];
end

