"""
	function solve_mc_vvvw_opf(
		data::Union{Dict{String,<:Any},String},
		model_type::Type,
		solver;
		kwargs...
	)

Solve DOE quantification problem
"""
function solve_mc_vvvw_opf(data::Union{Dict{String,<:Any},String},  solver; kwargs...)
    return _PMD.solve_mc_model(data, _PMD.IVRENPowerModel, solver, build_mc_vvvw_opf; kwargs...)
end

"""
	function build_mc_vvvw_opf(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for DOE quantification in current-voltage variable space with explicit neutrals
"""
function build_mc_vvvw_opf(pm::_PMD.AbstractExplicitNeutralIVRModel)
    # Register volt-var/watt functions
    vv_curve_pu = voltvar_handle()
    vw_curve_pu = voltwatt_handle()
    JuMP.register(pm.model, :vv_curve_pu, 1, vv_curve_pu; autodiff = true)
    JuMP.register(pm.model, :vw_curve_pu, 1, vw_curve_pu; autodiff = true)

    # Variables
    _PMD.variable_mc_bus_voltage(pm)
    _PMD.variable_mc_branch_current(pm)
    _PMD.variable_mc_load_current(pm)
    _PMD.variable_mc_load_power(pm)
    _PMD.variable_mc_generator_current(pm)
    _PMD.variable_mc_generator_power(pm)
    variable_mc_generator_voltage_magnitude(pm)
    _PMD.variable_mc_transformer_current(pm)
    _PMD.variable_mc_transformer_power(pm)
    _PMD.variable_mc_switch_current(pm)

    # Constraints
    for i in _PMD.ids(pm, :bus)

        if i in _PMD.ids(pm, :ref_buses)
            _PMD.constraint_mc_voltage_reference(pm, i)
        end

        _PMD.constraint_mc_voltage_absolute(pm, i)
        _PMD.constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in _PMD.ids(pm, :gen)
        _PMD.constraint_mc_generator_power(pm, id)
        _PMD.constraint_mc_generator_current(pm, id)
        constraint_mc_gen_vpn(pm, id)
        constraint_mc_gen_voltvar(pm, id)
        constraint_mc_gen_voltwatt(pm, id)
    end

    for id in _PMD.ids(pm, :load)
        _PMD.constraint_mc_load_power(pm, id)
        _PMD.constraint_mc_load_current(pm, id)
    end

    for i in _PMD.ids(pm, :transformer)
        _PMD.constraint_mc_transformer_voltage(pm, i)
        _PMD.constraint_mc_transformer_current(pm, i)

        _PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in _PMD.ids(pm, :branch)
        _PMD.constraint_mc_current_from(pm, i)
        _PMD.constraint_mc_current_to(pm, i)
        _PMD.constraint_mc_bus_voltage_drop(pm, i)

        _PMD.constraint_mc_branch_current_limit(pm, i)
        _PMD.constraint_mc_thermal_limit_from(pm, i)
        _PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in _PMD.ids(pm, :switch)
        _PMD.constraint_mc_switch_current(pm, i)
        _PMD.constraint_mc_switch_state(pm, i)

        _PMD.constraint_mc_switch_current_limit(pm, i)
        _PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in _PMD.ids(pm, :bus)
        _PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    _PMD.objective_mc_min_fuel_cost(pm)
    # objective_mc_max_pg_competitive(pm)
end

function rectifier(x,y,a;type="smooth", ϵ=0.0001)
    if type=="nonsmooth"
        return f = vals -> a*max(0, vals - x) + y
    elseif type=="smooth"
        return f = vals -> a*ϵ*StatsFuns.log1pexp((vals-x)/ϵ) + y
    end
end

```
Inputs in PU voltage
outputs in pu power relative to apparent power rating
```
function voltvar_handle(;ϵ=0.0001, type="smooth")
    V_vv = [195; 207; 220; 240; 258; 276]./230
    Q_vv = [44; 44;   0;   0;  -60;  -60]./100

    r1pu = rectifier(207/230,44/100,-44/13*(230/100);type=type, ϵ=ϵ)
    r2pu = rectifier(220/230,0,+44/13*(230/100);type=type, ϵ=ϵ)
    r3pu = rectifier(240/230,0,-60/18*(230/100);type=type, ϵ=ϵ)
    r4pu = rectifier(258/230,0,+60/18*(230/100);type=type, ϵ=ϵ)
    vv_curve_pu(x) = r1pu(x) + r2pu(x) + r3pu(x) + r4pu(x)
end

```
Inputs in PU voltage
outputs in pu power relative to apparent power rating
```
function voltwatt_handle(;ϵ=0.0001, type="smooth")
    V_vw = [195; 253; 260; 276]./230
    P_vw = [100; 100; 20; 20]./100

    relupu = rectifier(253/230,1,-80/7*(230/100);type=type, ϵ=ϵ)
    relu2pu =rectifier(260/230,0,+80/7*(230/100);type=type, ϵ=ϵ)
    vw_curvepu(x) = relupu(x) + relu2pu(x)
end

function constraint_mc_gen_vpn(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    generator = _PMD.ref(pm, nw, :gen, id)
    configuration = generator["configuration"]
    vpn = _PMD.var(pm, nw, :vg_pn, id)
    # bus = _PMD.ref(pm, nw, :bus, generator["gen_bus"])["bus_i"]
    vr = _PMD.var(pm, nw, :vr, generator["gen_bus"])
    vi = _PMD.var(pm, nw, :vi, generator["gen_bus"])
    phases = generator["connections"][1]
    for (idx, p) in enumerate(phases)
        # @show idx, p, id
        JuMP.@constraint(pm.model,  (vr[p]-vr[4])^2 + (vi[p]-vi[4])^2 == vpn[idx]^2)
    end
end


function constraint_mc_gen_voltvar(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    generator = _PMD.ref(pm, nw, :gen, id)
    configuration = generator["configuration"]
    N = length(generator["connections"])
    smax = 5.0*ones(N-1)

    qg = _PMD.var(pm, nw, :qg, id)
    vpn = _PMD.var(pm, nw, :vg_pn, id)
    phases = generator["connections"][1]
    @show "vv  gen $id"
    if id == 1
        # do nothing for source bus
    else
        if configuration==_PMD.WYE
            for (idx, p) in enumerate(phases)
                JuMP.@NLconstraint(pm.model, qg[idx] == vv_curve_pu(vpn[idx])*smax[idx])
                @show "vv added for gen $id"
            end
        else #Delta
            error("delta connections not supported")
        end
    end
end

function constraint_mc_gen_voltwatt(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    generator = _PMD.ref(pm, nw, :gen, id)
    configuration = generator["configuration"]
    N = length(generator["connections"])
    smax = 5.0*ones(N-1)
    pg  = _PMD.var(pm, nw, :pg,    id)
    vpn = _PMD.var(pm, nw, :vg_pn, id)
    phases = generator["connections"][1:end-1]

    if !contains(generator["source_id"], "voltage_source.source")
        if configuration==_PMD.WYE
            for (idx, p) in enumerate(phases)
                JuMP.@NLconstraint(pm.model, pg[idx] <= vw_curve_pu(vpn[idx])*smax[idx])
            end

        else #Delta
            error("delta connections not supported")
        end
    end
end

function variable_mc_generator_voltage_magnitude(pm::_PMD.AbstractUnbalancedIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_generator_voltage_magnitude_phase_neutral(pm; nw=nw, bounded=bounded, report=report)
end


function variable_mc_generator_voltage_magnitude_phase_neutral(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _PMD._infer_int_dim_unit(gen, false) for (i,gen) in _PMD.ref(pm, nw, :gen))
    @show int_dim
    vg_pn = _PMD.var(pm, nw)[:vg_pn] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_vg_pn_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :gen, i), "vg_pn_start", c, 1.0)
        ) for i in _PMD.ids(pm, nw, :gen)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :vg_pn, _PMD.ids(pm, nw, :gen), vg_pn)
end