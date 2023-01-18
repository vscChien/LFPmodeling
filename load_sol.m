% load solutions p [4x28] for the 4 recording sites 
%       p(s,1:8)   % connection strengths 
%       p(s,9:10)  % ratio of thalamic inputs
%       p(s,11:12) % short-term plasticity 
%       p(s,13)    % synaptic kernel time constants
%       p(s,14)    % sigmoid function slopes
%       p(s,15:19) % thalamic input decay levels
%       p(s,20:24) % lateral inhibition E2→SOM1,2
%       p(s,25:28) % thalamic input strengths

function p=load_sol()
    if ispc % windows
        path='.\solution\';
    else    % linux
        path='./solution/';
    end
    p=zeros(4,28);
    for s=1:4
      tmp=load([path, sprintf('solution_site%d.mat',s)]);
      p(s,:)=tmp.p;
    end
end
