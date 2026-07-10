function [ids_start, ids_end, X_filt] = locate_true_transients(X,l,L)
% Locate true transients that lasts for at least l frames but less than L
% frames
% return the start and end locations of transients
% X_filt: filtered data after removing all false transients
[N,T] = size(X);
X_pad = zeros(N,T+2);
X_pad(:,2:T+1) = X;
ids_start = cell(N,1);
ids_end = cell(N,1);
dX = diff(X_pad>0,1,2);

for n=1:N
    id_start = find(dX(n,:)==1);
    id_end = find(dX(n,:)==-1);
    if ~isempty(id_start) && ~isempty(id_end)
        for i=1:length(id_start)
            if (id_end(i)-id_start(i)>=l) && id_end(i)-id_start(i)<=L
                ids_start{n} = [ids_start{n}, id_start(i)];
                ids_end{n} = [ids_end{n}, id_end(i)-1];
            else
                X(n,id_start(i):id_end(i)-1) = 0;
            end
        end
    end
end
X_filt = X;
end