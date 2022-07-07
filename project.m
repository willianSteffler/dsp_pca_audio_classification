function [ux] = project(U,X,i,f)
%PROJECT Summary of this function goes here
%   Detailed explanation goes here
Vi = U';
Xi = X';
Vi = Vi(:,i:f);
Xi = Xi(i:f,:);
ux = Vi*Xi;
ux = ux';
end

