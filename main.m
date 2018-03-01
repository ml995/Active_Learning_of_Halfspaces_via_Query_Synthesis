clear;
noise = 0.1;
n = 50;
orackle = rand(1,n);
orackle = orackle / norm(orackle);
epsilon = 1e-3;

es = zeros(n,n);

for i=1:n
    es(i,i) = 1;
end

t = n;
nQueries = 0;
while(t > 1)
    tt = floor(t/2);
    
    for j=1:tt
        [a1,a2,tot_it] = DC2(es(2*j-1,:), es(2*j,:), orackle, epsilon / n, noise, noise);
        nQueries = nQueries + tot_it;
        es(j,:) = a1 * es(2*j-1,:) + a2* es(2*j,:);
    end
    
    if (mod(t,2) == 1)
        es(tt+1,:) = es(t,:);
        t = tt + 1;
    else
        t = tt;
    end
end
disp('Estimated vector: ')
disp(num2str(squeeze(es(1, :))));
disp('True vector: ')
disp(num2str(squeeze(orackle)));
fprintf('Distance: %d\n', norm(es(1, :) - orackle));
fprintf('Total number of queries: %d\n', nQueries);