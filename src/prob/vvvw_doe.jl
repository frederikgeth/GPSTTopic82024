"""
	function solve_mc_vvvw_opf(
		data::Union{Dict{String,<:Any},String},
		model_type::Type,
		solver;
		kwargs...
	)

Solve DOE quantification problem with volt-var/watt response and objective for competitive outcome
"""
function solve_mc_vvvw_doe_competitive(data::Union{Dict{String,<:Any},String},  solver; kwargs...)
    return _PMD.solve_mc_model(data, _PMD.IVRENPowerModel, solver, build_mc_vvvw_doe_competitive; kwargs...)
end

"""
	function build_mc_vvvw_doe_competitive(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for DOE quantification in current-voltage variable space with explicit neutrals and  with volt-var/watt response
    and objective for comeptitive outcome
"""
function build_mc_vvvw_doe_competitive(pm::_PMD.AbstractExplicitNeutralIVRModel)
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
    objective_mc_max_pg_competitive(pm)
end


"""
	function solve_mc_vvvw_doe_mse(
		data::Union{Dict{String,<:Any},String},
		model_type::Type,
		solver;
		kwargs...
	)

Solve DOE quantification problem with volt-var/watt response and objective for mean square differences 
"""
function solve_mc_vvvw_doe_mse(data::Union{Dict{String,<:Any},String},  solver; kwargs...)
    return _PMD.solve_mc_model(data, _PMD.IVRENPowerModel, solver, build_mc_vvvw_doe_fair_mse; kwargs...)
end

"""
	function build_mc_vvvw_doe_fair_mse(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for DOE quantification in current-voltage variable space with explicit neutrals and  with volt-var/watt response
    and objective for mean square differences 
"""
function build_mc_vvvw_doe_fair_mse(pm::_PMD.AbstractExplicitNeutralIVRModel)
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
    objective_mc_fair_pg_mse(pm)
end


"""
	function solve_mc_vvvw_doe_abs(
		data::Union{Dict{String,<:Any},String},
		model_type::Type,
		solver;
		kwargs...
	)

Solve DOE quantification problem with volt-var/watt response and objective for mean absolute value differences 
"""
function solve_mc_vvvw_doe_abs(data::Union{Dict{String,<:Any},String},  solver; kwargs...)
    return _PMD.solve_mc_model(data, _PMD.IVRENPowerModel, solver, build_mc_vvvw_doe_fair_abs; kwargs...)
end

"""
	function build_mc_vvvw_doe_fair_abs(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for DOE quantification in current-voltage variable space with explicit neutrals and  with volt-var/watt response
    and objective for mean absolute value differences 
"""
function build_mc_vvvw_doe_fair_abs(pm::_PMD.AbstractExplicitNeutralIVRModel)
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
    objective_mc_fair_pg_abs(pm)
end

"""
	function solve_mc_vvvw_doe_equal(
		data::Union{Dict{String,<:Any},String},
		model_type::Type,
		solver;
		kwargs...
	)

Solve DOE quantification problem with volt-var/watt response and objective for mean absolute value differences 
"""
function solve_mc_vvvw_doe_equal(data::Union{Dict{String,<:Any},String},  solver; kwargs...)
    return _PMD.solve_mc_model(data, _PMD.IVRENPowerModel, solver, build_mc_vvvw_doe_equal; kwargs...)
end

"""
	function build_mc_vvvw_doe_equal(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for DOE quantification in current-voltage variable space with explicit neutrals and  with volt-var/watt response
    and objective for mean absolute value differences 
"""
function build_mc_vvvw_doe_equal(pm::_PMD.AbstractExplicitNeutralIVRModel)
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
    objective_mc_max_pg_equal(pm)
end


"""
	function solve_mc_vvvw_doe_log_fairness(
		data::Union{Dict{String,<:Any},String},
		model_type::Type,
		solver;
		kwargs...
	)

Solve DOE quantification problem with volt-var/watt response and objective for mean absolute value differences 
"""
function solve_mc_vvvw_doe_log_fairness(data::Union{Dict{String,<:Any},String},  solver; kwargs...)
    return _PMD.solve_mc_model(data, _PMD.IVRENPowerModel, solver, build_mc_vvvw_doe_log_fairness; kwargs...)
end

"""
	function build_mc_vvvw_doe_log_fairness(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for DOE quantification in current-voltage variable space with explicit neutrals and  with volt-var/watt response
    and objective for mean absolute value differences 
"""
function build_mc_vvvw_doe_log_fairness(pm::_PMD.AbstractExplicitNeutralIVRModel)
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
    objective_mc_max_log_fairness(pm)
end
