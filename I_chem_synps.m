function [i_synps,s] = I_chem_synps(j,tStep,g,S,N,alpha, beta,dt,E_syn,v_a1,E_syn_inh)

n = size(g,1);
s = zeros(1,n);

i_synps = g(j,:) .* S(j,:,tStep) .* (v_a1 - E_syn(j,:));
i_synps = sum(i_synps,2);

% g_syn = g(j,:) .* S(j,:,tStep);
% g_syn = sum(g_syn,2);
% i_synps = g_syn * ( v_a1 - E_syn(j) );


for i = 1:n
    
    if E_syn(j,i) == E_syn_inh
        alphaa = alpha(2);
        betaa = beta(2);
    else
        alphaa = alpha(1);
        betaa = beta(2);
    end
    
    s(i) = S(j,i,tStep) + dt * ( alphaa * N(tStep,i) *(1-S(j,i,tStep)) - ...
        betaa * S(j,i,tStep));
end

end

