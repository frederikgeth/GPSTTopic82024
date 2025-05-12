using GPSTTopic82024
using PowerModelsDistribution
using Ipopt
using Plots
ipopt = Ipopt.Optimizer
## Main loop

file = "data/LV30_315bus/Master.dss"

# vscale = 1.07

for vscale in 0.98:0.01:1.07, loadscale in [1] #0.8:0.05:1.0
    # for vscale in 0.98:0.01:1.07, loadscale in [1] #0.8:0.05:1.0
    eng4w = parse_file(file, transformations=[transform_loops!,remove_all_bounds!])
    eng4w["settings"]["sbase_default"] = 1
    eng4w["voltage_source"]["source"]["rs"] *=0
    eng4w["voltage_source"]["source"]["xs"] *=0
    eng4w["voltage_source"]["source"]["vm"] *=vscale

    reduce_line_series!(eng4w)
    math4w = transform_data_model(eng4w, kron_reduce=false, phase_project=false)
    add_start_vrvi!(math4w)

    for (i,bus) in math4w["bus"]

        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
            bus["vmax"] = ones(4)*1.2
            bus["vmax"][end] = 0.05
            bus["vmin"] = zeros(4)
            # bus["grounded"] .=  0
        else
            @show bus
        end

    end

    for (g,gen) in math4w["gen"]
        gen["cost"] = 0.0
    end

    for (d,load) in math4w["load"]
        load["pd"] .*= loadscale
        load["qd"] .*= loadscale
    end

    function add_gens!(math4w)
        gen_counter = 2
        for (d, load) in math4w["load"]
            if mod(load["index"], 7) == 1
                # phases = 3
                phases = length(load["connections"])-1
                math4w["gen"]["$gen_counter"] = deepcopy(math4w["gen"]["1"])
                math4w["gen"]["$gen_counter"]["name"] = "$gen_counter"
                math4w["gen"]["$gen_counter"]["index"] = gen_counter
                math4w["gen"]["$gen_counter"]["cost"] = 1.0 #*math4w["gen"]["1"]["cost"]
                math4w["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
                math4w["gen"]["$gen_counter"]["pmax"] = 5.0*ones(phases)
                math4w["gen"]["$gen_counter"]["pmin"] = 0.0*ones(phases)
                math4w["gen"]["$gen_counter"]["qmax"] = 5.0*ones(phases)
                math4w["gen"]["$gen_counter"]["qmin"] = -5.0*ones(phases)
                math4w["gen"]["$gen_counter"]["connections"] = load["connections"]
                gen_counter = gen_counter + 1
            end
        end
    end
    add_gens!(math4w)


    res_comp = solve_mc_vvvw_doe_competitive(math4w, ipopt)
    @assert(res_comp["termination_status"]==LOCALLY_SOLVED || res_comp["termination_status"]==ALMOST_LOCALLY_SOLVED)
    res_comp_obj = round(res_comp["objective"], digits=2)
    pg_cost1 = [gen["pg_cost"] for (g,gen) in res_comp["solution"]["gen"] if g!="1"]

    res_ms = solve_mc_vvvw_doe_mse(math4w, ipopt)
    @assert(res_ms["termination_status"]==LOCALLY_SOLVED || res_ms["termination_status"]==ALMOST_LOCALLY_SOLVED)
    pg_cost2 = [gen["pg_cost"] for (g,gen) in res_ms["solution"]["gen"] if g!="1"]
    res_ms_obj = round(res_ms["objective"], digits=2)

    res_abs = solve_mc_vvvw_doe_abs(math4w, ipopt)
    @assert(res_abs["termination_status"]==LOCALLY_SOLVED || res_abs["termination_status"]==ALMOST_LOCALLY_SOLVED)
    pg_cost3 = [gen["pg_cost"] for (g,gen) in res_abs["solution"]["gen"] if g!="1"]
    res_abs_obj = round(res_abs["objective"], digits=2)

    res_eq_obj = 0
    res_eq = solve_mc_vvvw_doe_equal(math4w, ipopt)
    if res_eq["termination_status"]==LOCALLY_INFEASIBLE
        pg_cost4 = 0 .*pg_cost1
        res_eq_obj = 0
    else
        @assert(res_eq["termination_status"]==LOCALLY_SOLVED || res_eq["termination_status"]==ALMOST_LOCALLY_SOLVED)
        pg_cost4 = [gen["pg_cost"] for (g,gen) in res_eq["solution"]["gen"] if g!="1"]
        res_eq_obj = round(res_eq["objective"], digits=2)
    end

    bb = sortperm(pg_cost3)
    pg_cost1 = pg_cost1[bb]
    pg_cost2 = pg_cost2[bb]
    pg_cost3 = pg_cost3[bb]
    pg_cost4 = pg_cost4[bb]

    plot(pg_cost1, label="max competitive $res_comp_obj")
    plot!(pg_cost2, label="min deviation squared $res_ms_obj")
    plot!(pg_cost3, label="min absolute deviation $res_abs_obj ")
    plot!(pg_cost4, label="equal $res_eq_obj")
    xlabel!("PV system id (-)")

    ylabel!("Export DOE (kW)")
    ylims!(0,5.1)
    xlims!(0.5,length(pg_cost1)+0.5)
    title!("Export limits at reference voltage of $vscale pu with load at $loadscale")
    savefig("objective_comparison_vsource$vscale load$loadscale.pdf")

end
# v_mag = stack([hypot.(bus["vr"],bus["vi"]) for (b,bus) in res["solution"]["bus"]], dims=1)


# plot(v_mag, label=["a" "b" "c" "n"])
# plot!([0; length(res["solution"]["bus"])], [0.9; 0.9], label="vmin")
# plot!([0; length(res["solution"]["bus"])], [1.1; 1.1], label="vmax")
# ylabel!("V (pu)")
# xlabel!("bus id (-)")

# v_mag = [gen["vg_pn"] for (g,gen) in res["solution"]["gen"]]


# for (b,bus) in res["solution"]["bus"]
#     bus["vpn"] = vpn = abs.((bus["vr"][1:end-1] .+ im.*bus["vi"][1:end-1]) .- (bus["vr"][end] .+ im.*bus["vi"][end]))
#     bus["vm"] = vm = abs.(bus["vr"] .+ im.*bus["vi"])
# end


# for (g,gen) in math4w["gen"]
#     bus = gen["gen_bus"]
#     conn = gen["connections"][1:end-1]
#     @show conn
#     @show res["solution"]["bus"]["$bus"]["vpn"] #[conn]
#     @show res["solution"]["bus"]["$bus"]["vm"] #[conn]
#     @show res["solution"]["gen"]["$g"]["vg_pn"]
# end
