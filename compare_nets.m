function [R] = compare_nets(new,orig)

N = sqrt(size(new,2));

I = find(triu(ones(N),1));

R = zeros(size(new,1),1);

for s = 1:size(new,1)
    R(s) = corr(new(s,I)',orig(s,I)');
end

R = 0.5*log((1+R)./(1-R)); R = mean(R);
