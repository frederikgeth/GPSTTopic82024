function alpha_fairness(alpha)
    @assert alpha>=0
    if alpha == 1
        return f = vals -> sum(log.(vals)) 
    else
        return f = vals -> (1/1-alpha)*sum(vals.^(1-alpha)) 
    end
end

W2 = alpha_fairness(2)
W1 = alpha_fairness(1)

n_samples = 10
fair = ones(n_samples)
unfair = collect((1:n_samples)./n_samples)
W1f = W1(fair)
W1u = W1(unfair)

W2f = W2(fair)
W2u = W2(unfair)


#W1
using JuMP
m = JuMP.Model(Ipopt.Optimizer)
# register(m, :W1, 1, W1; autodiff = true)
@variable(m, 0.01 <=x[1:10]<=1.0)
@variable(m, logx[1:10])
@NLconstraint(m, [i=1:10], logx[i] == log(x[i]))
JuMP.set_upper_bound(x[2],.01)
@objective(m, Max, sum(logx) )
optimize!(m)

@show value.(x)
@show value.(logx)
termination_status(m)
scatter!([value.(x)],[value.(y11)],markersize=5 )


#W2
using JuMP
m = JuMP.Model(Ipopt.Optimizer)
alpha = 2
@variable(m, 0.01 <=x[1:10]<=1.0)
@variable(m, xut[1:10])
@NLconstraint(m, [i=1:10], xut[i] == x[i]^(1-alpha))
JuMP.set_upper_bound(x[2],.01)
@objective(m, Max, sum(xut)/(1-alpha) )
optimize!(m)

@show value.(x)
@show value.(logx)
termination_status(m)
scatter!([value.(x)],[value.(y11)],markersize=5 )