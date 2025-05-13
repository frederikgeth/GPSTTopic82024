
"""
    objective_mc_max_pg_equal(pm::AbstractUnbalancedPowerModel)

    equal assignment
"""
function objective_mc_max_pg_equal(pm::_PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    refgens = Dict()

    for (n, nw_ref) in _PMD.nws(pm)
        pg_cost = _PMD.var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in _PMD.ids(pm, n, :gen)], base_name="$(n)_pg_cost",
        )
        report && _IM.sol_component_value(pm, pmd_it_sym, n, :gen, :pg_cost, _PMD.ids(pm, n, :gen), pg_cost)
        for (i, gen) in nw_ref[:gen]
            pg = _PMD.var(pm, n, :pg, i)
            phases = length(gen["connections"][1:end-1])
            JuMP.@NLconstraint(pm.model, pg_cost[i] == gen["cost"]*sum(pg[c] for c in 1:length(phases)))
        end
        gensids = 2:length(nw_ref[:gen]) #1 is reference bus, ignore here
        if !(isempty(gensids))
            refgens[n] = gensids[1]
            for i in gensids[2:end]
                phasesi = length(nw_ref[:gen][i]["connections"][1:end-1])
                phasesref = length(nw_ref[:gen][gensids[1]]["connections"][1:end-1])
                JuMP.@NLconstraint(pm.model, sum(_PMD.var(pm, n, :pg, i)[c] for c in 1:phasesi) == sum(_PMD.var(pm, n, :pg, gensids[1])[c] for c in 1:phasesref))
            end
        end
    end
    return JuMP.@NLobjective(pm.model, Max,
            sum(sum(_PMD.var(pm, n, :pg, refgen)[c] for c in 1:length(_PMD.nws(pm)[n][:gen][refgen]["connections"][1:end-1])) for (n, refgen) in refgens)
    )
end
"""
    objective_mc_max_pg_competitive(pm::AbstractUnbalancedPowerModel)

"""
function objective_mc_max_pg_competitive(pm::_PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    objective_variable_pg_competitive(pm; report=report)
    obj = JuMP.@objective(pm.model, Max,
        sum(
            sum( _PMD.var(pm, n, :pg_cost, i) for (i,gen) in nw_ref[:gen] if i!=1)
        for (n, nw_ref) in _PMD.nws(pm))
    )
    return obj 
end


"""
    objective_variable_pg_competitive(pm::AbstractUnbalancedIVRModel)

adds pg_cost variables and constraints for the IVR formulation
"""
function objective_variable_pg_competitive(pm::_PMD.AbstractUnbalancedIVRModel; report::Bool=report)
    for (n, nw_ref) in _PMD.nws(pm)
        #to avoid function calls inside of @NLconstraint
        pg_cost = _PMD.var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in _PMD.ids(pm, n, :gen)], base_name="$(n)_pg_cost",
        )
        report && _IM.sol_component_value(pm, pmd_it_sym, n, :gen, :pg_cost, _PMD.ids(pm, n, :gen), pg_cost)

        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            pg = _PMD.var(pm, n, :pg, i)
            phases = length(gen["connections"][1:end-1])
            JuMP.@NLconstraint(pm.model, pg_cost[i] == gen["cost"]*sum(pg[c] for c in 1:length(phases)))
        end
    end
end

"""
    objective_mc_max_log_fairness(pm::AbstractUnbalancedPowerModel)

"""
function objective_mc_max_log_fairness(pm::_PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    objective_variable_log_fairness(pm; report=report)
    return JuMP.@objective(pm.model, Max,
        sum(
            sum( _PMD.var(pm, n, :pg_cost_log, i) for (i,gen) in nw_ref[:gen] if i!=1)
        for (n, nw_ref) in _PMD.nws(pm))
    )
end

"""
    objective_variable_log_fairness(pm::AbstractUnbalancedIVRModel)

adds pg_cost variables and constraints for the IVR formulation
"""
function objective_variable_log_fairness(pm::_PMD.AbstractUnbalancedIVRModel; report::Bool=report)
    for (n, nw_ref) in _PMD.nws(pm)
            #to avoid function calls inside of @NLconstraint
        pg_cost = _PMD.var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
                    [i in _PMD.ids(pm, n, :gen)], base_name="$(n)_pg_cost", start = 0.0000001
        ) 
        # if LB is zero AD may evaluate log(0)
        report && _IM.sol_component_value(pm, pmd_it_sym, n, :gen, :pg_cost, _PMD.ids(pm, n, :gen), pg_cost)
        
        pg_cost_log = _PMD.var(pm, n)[:pg_cost_log] = JuMP.@variable(pm.model,
            [i in _PMD.ids(pm, n, :gen)], base_name="$(n)_pg_cost_log", 
        )
        report && _IM.sol_component_value(pm, pmd_it_sym, n, :gen, :pg_cost_log, _PMD.ids(pm, n, :gen), pg_cost_log)

        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            pg = _PMD.var(pm, n, :pg, i)
            if i==1
                phases = length(gen["connections"][1:end-1])
                JuMP.@NLconstraint(pm.model, pg_cost[i] == gen["cost"]*sum(pg[c] for c in 1:length(phases)))
                JuMP.@constraint(pm.model, pg_cost_log[i] == 0)
            else
                phases = length(gen["connections"][1:end-1])
                @assert gen["cost"]>0
                JuMP.@NLconstraint(pm.model, pg_cost[i] == gen["cost"]*sum(pg[c] for c in 1:length(phases)))
                JuMP.@NLconstraint(pm.model, pg_cost_log[i] == log(pg_cost[i]))
            end
        end
    end
end


"""
    objective_mc_fair_pg_mse(pm::AbstractUnbalancedPowerModel)

Fuel cost minimisation and fairness objective, using squared difference from mean
"""
function objective_mc_fair_pg_mse(pm::_PMD.AbstractUnbalancedPowerModel; report::Bool=true, weightcompetitive=.01)
    objective_variable_pg_fair_mse(pm; report=report)
 
    return JuMP.@objective(pm.model, Min,
            sum(
                sum(
                    (_PMD.var(pm, n, :pg_cost, i) - _PMD.var(pm, n, :pg_avg))^2 for (i,gen) in nw_ref[:gen]
                ) for (n, nw_ref) in _PMD.nws(pm)
            ) 
            - weightcompetitive*sum(
                sum( _PMD.var(pm, n, :pg_cost, i) for (i,gen) in nw_ref[:gen])
            for (n, nw_ref) in _PMD.nws(pm))
        )
end

"""
    objective_mc_fair_pg_abs(pm:AbstractUnbalancedPowerModel)

Fuel cost minimisation and fairness objective, using the absolute difference from mean
"""
function objective_mc_fair_pg_abs(pm::_PMD.AbstractUnbalancedPowerModel; report::Bool=true, weightcompetitive=.01)
    objective_variable_pg_fair_abs(pm; report=report)
 
    return JuMP.@objective(pm.model, Min,
            sum(
                sum(
                    (_PMD.var(pm, n, :pg_abs, i)) for (i,gen) in nw_ref[:gen]
                ) for (n, nw_ref) in _PMD.nws(pm)
            ) 
            -weightcompetitive*sum(
                sum( 
                    _PMD.var(pm, n, :pg_cost, i) for (i,gen) in nw_ref[:gen]
                ) for (n, nw_ref) in _PMD.nws(pm))
        )
end

"""
    objective_variable_pg_fair_mse(pm::AbstractUnbalancedIVRModel)

Adds pg_cost and pg_avg variables and constraints for the IVR formulation
"""
function objective_variable_pg_fair_mse(pm::_PMD.AbstractUnbalancedIVRModel; report::Bool=report)
    for (n, nw_ref) in _PMD.nws(pm)
        #to avoid function calls inside of @NLconstraint
        pg_cost = _PMD.var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in _PMD.ids(pm, n, :gen)], base_name="$(n)_pg_cost",
        )
        report && _IM.sol_component_value(pm, pmd_it_sym, n, :gen, :pg_cost, _PMD.ids(pm, n, :gen), pg_cost)

        pg_avg = _PMD.var(pm, n)[:pg_avg] = JuMP.@variable(pm.model, avg, base_name="pg_avg")
        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            phases = length(gen["connections"][1:end-1])
            pg = _PMD.var(pm, n, :pg, i)
            JuMP.@NLconstraint(pm.model, pg_cost[i] == gen["cost"]*sum(pg[c] for c in 1:phases))
        end

        JuMP.@constraint(pm.model, pg_avg == sum(pg_cost)/length(pg_cost))
    end
end

"""
    objective_variable_pg_fair_mse(pm::AbstractUnbalancedIVRModel)

Adds pg_cost, pg_avg and pg_abs variables and constraints for the IVR formulation
"""
function objective_variable_pg_fair_abs(pm::_PMD.AbstractUnbalancedIVRModel; report::Bool=report)
    for (n, nw_ref) in _PMD.nws(pm)
        #to avoid function calls inside of @NLconstraint
        pg_cost = _PMD.var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in _PMD.ids(pm, n, :gen)], base_name="$(n)_pg_cost",
        )
        report && _IM.sol_component_value(pm, pmd_it_sym, n, :gen, :pg_cost, _PMD.ids(pm, n, :gen), pg_cost)

        pg_avg = _PMD.var(pm, n)[:pg_avg] = JuMP.@variable(pm.model, avg, base_name="pg_avg")

        pg_abs = _PMD.var(pm, n)[:pg_abs] = JuMP.@variable(pm.model,
            [i in _PMD.ids(pm, n, :gen)], base_name="$(n)_pg_abs",
        )

        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            phases = length(gen["connections"][1:end-1])
            pg = _PMD.var(pm, n, :pg, i)
            JuMP.@NLconstraint(pm.model, pg_cost[i] == gen["cost"]*sum(pg[c] for c in 1:phases))
        end

        JuMP.@constraint(pm.model, pg_avg == sum(pg_cost)/length(pg_cost))
        
        for (i, _) in nw_ref[:gen]
            JuMP.@NLconstraint(pm.model, pg_abs[i] >= pg_cost[i] - pg_avg)
            JuMP.@NLconstraint(pm.model, pg_abs[i] >= -(pg_cost[i] - pg_avg))
        end
    end
end

