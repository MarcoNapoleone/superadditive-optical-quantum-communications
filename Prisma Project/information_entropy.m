

%%%
%
% Entropy function of a given probability distribution
%
%%%



function [entropy] = h(x)
    entropy = -x*log2(x);
end


function [entropy] = H(m)
    if sum(m) ~= 1
        error('The sum of the probabilities must be 1');
    end

    entropy = 0;
    for i = 1:length(m)
        if m(i) ~= 0
            entropy = entropy + h(m(i));
        end
    end
    
end


function [entropy] = H2(p)
    entropy = H([p, 1-p])
end


% example of usage
m = [0.001 0.999];


disp("Entropy of a uniform distribution: ");
disp(H(m));



