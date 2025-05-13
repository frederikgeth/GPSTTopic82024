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
    # @show "vv  gen $id"
    if id == 1
        # do nothing for source bus
        @show  generator
    else
        if configuration==_PMD.WYE
            for (idx, p) in enumerate(phases)
                JuMP.@NLconstraint(pm.model, qg[idx] == vv_curve_pu(vpn[idx])*smax[idx])
                # @show "vv added for gen $id"
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
    phases = generator["connections"][1]

    if id == 1
        # do nothing for source bus
    else
        if configuration==_PMD.WYE
            for (idx, p) in enumerate(phases)
                JuMP.@NLconstraint(pm.model, pg[idx] <= vw_curve_pu(vpn[idx])*smax[idx])
            end
        else #Delta
            error("delta connections not supported")
        end
    end
end
