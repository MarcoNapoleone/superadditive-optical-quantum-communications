function H = entropy(x)
    H=-x.*log2(x)-(1-x).*log2(1-x); 
end