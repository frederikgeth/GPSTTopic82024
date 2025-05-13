function alpha_fairness(alpha)
    if alpha == 1
        return f = vals -> sum(log.(vals)) 
    else
        return f = vals -> (1/1-alpha)*sum(vals.^(a-alpha)) 
    end
end

W2 = alpha_fairness(2)
W1 = alpha_fairness(1)