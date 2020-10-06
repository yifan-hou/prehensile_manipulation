function feasible = velocityFeasibilityCheck(J, G, b_G)
JG = [J; G];
b_JG = [zeros(size(J,1), 1); b_G];
if rank([JG b_JG]) > rank(JG)
    % Goal is infeasible
    feasible = false;
else
    feasible = true;
end



end