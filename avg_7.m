function [ avg_E ] = avg_7( E )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    
    avg_E = zeros(1,numel(E));
    offset = find(E,1);
    
    X = E(find(E));
    
    for k = 1:numel(X)
        coeff=[];
        for l=-min(4,k-1):min(4,numel(X)-k)
            coeff=[coeff X(k+l)];
        end
        X_out(k)=mean(coeff);
    end
    
    avg_E(find(E)) = X_out;
    
end

