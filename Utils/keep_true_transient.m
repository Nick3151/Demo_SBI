function [Y, counts] = keep_true_transient(X,l,L)
% Keep true transients that lasts for at least l frames but less than L
% frames
% compute the number of transients for each neuron
[N,T] = size(X);
X_pad = zeros(N,T+2);
X_pad(:,2:T+1) = X;
counts = zeros(N,1);
dX = diff(X_pad>0,1,2);
for i=1:N
    f_start = find(dX(i,:)==1);
    f_end = find(dX(i,:)==-1);
    if isempty(f_start) && isempty(f_end)
        counts(i) = 0;
    else
        for j=1:length(f_start)
            if (f_end(j)-f_start(j)<l) || f_end(j)-f_start(j)>=L
                X(i,f_start(j):f_end(j)-1) = 0;
            else
                counts(i) = counts(i)+1;
            end
        end
    end
end

Y = X;
end