function vec = linspaceNoCorner(a1, a2, N)
    vec = linspace((a1*(1-1/N)),(a2*(1-1/N)),N); 
end