
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