function [alpha1 , alpha2]= DFA(RR,ID)

x_hat = mean(RR);
y = cumsum(RR-x_hat);
y_len = length(y);
n_short_term = 4:16;
n_long_term  = 16:64;
count = 1;

for n = n_short_term
    rem = mod(y_len,n);
    y_trunc = y(1:end-rem);
    y_reshape = reshape(y_trunc,n,[])';
    win_index = reshape(1:length(y_trunc),n,[])'; 
    y_n = zeros(size(win_index,1),n);
    
    for i=1:size(win_index,1)
        b = y_reshape(i,:)';
        A = [ones(n,1) , win_index(i,:)'];
        x = A\b;
        y_n(i,:) = A*x;

%         p = polyfit(win_index(i,:)',b,1);
%         y_n(i,:)=polyval(p,win_index(i,:)')';
        
    end
    y_n = reshape(y_n',1,[]);
    F_short(count) = sqrt(1/length(y_trunc)*sum((y_trunc-y_n).^2));
    count = count+1;
end

count = 1;
for n = n_long_term
    rem = mod(y_len,n);
    y_trunc = y(1:end-rem);
    y_reshape = reshape(y_trunc,n,[])';
    win_index = reshape(1:length(y_trunc),n,[])'; 
    y_n = zeros(size(win_index,1),n);
    
    for i=1:size(win_index,1)
        b = y_reshape(i,:)';
        A = [ones(n,1) , win_index(i,:)'];
        x = A\b;
        y_n(i,:) = A*x;

%         p = polyfit(win_index(i,:)',b,1);
%         y_n(i,:)=polyval(p,win_index(i,:)')';

    end
    y_n = reshape(y_n',1,[]);
    F_long(count) = sqrt(1/length(y_trunc)*sum((y_trunc-y_n).^2));
    count = count+1;
end

p = polyfit(log10(n_short_term),log10(F_short),1);
alpha1 = p(1);
p = polyfit(log10(n_long_term),log10(F_long),1);
alpha2 = p(1);

figure,
subplot(1,2,1)
plot(log10(n_short_term),log10(F_short))
xlabel('log(n) (window size, short-term)')
ylabel('log(F(n))')
title(sprintf('alpha1 = %.2f , ID = %s',alpha1,ID))

subplot(1,2,2)
plot(log10(n_long_term),log10(F_long))
xlabel('log(n) (window size, long-term)')
ylabel('log(F(n))')
title(sprintf('alpha2 = %.2f , ID = %s',alpha2,ID))

end