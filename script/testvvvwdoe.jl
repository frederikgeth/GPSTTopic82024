using GPSTTopic82024
using PowerModelsDistribution
using Ipopt
using Plots
using CSV
using DataFrames

# ipopt = Ipopt.Optimizer

ipopt = optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>1000, "tol" => 1e-6)

## Main loop
case = "LV9_258bus"
casepath = "data/$case"
file = "$casepath/Master.dss"
busdistancesdf = CSV.read("$casepath/busdistances.csv", DataFrame)
busdistances_dict = Dict(busdistancesdf.Bus[i] => busdistancesdf.busdistances[i] for i in 1:nrow(busdistancesdf))


vscale = 1.08
loadscale = 0.2

# for vscale in 1.06:0.01:1.10, loadscale in [0.2] #0.8:0.05:1.0
    # for vscale in 0.98:0.01:1.07, loadscale in [1] #0.8:0.05:1.0
    eng4w = parse_file(file, transformations=[transform_loops!,remove_all_bounds!])
    eng4w["settings"]["sbase_default"] = 1
    eng4w["voltage_source"]["source"]["rs"] *=0
    eng4w["voltage_source"]["source"]["xs"] *=0
    eng4w["voltage_source"]["source"]["vm"] *=vscale

    reduce_line_series!(eng4w)
    math4w = transform_data_model(eng4w, kron_reduce=false, phase_project=false)

    for (i,bus) in math4w["bus"]

        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            lb =0.9
            bus["vm_pair_lb"] = [(1, 4, lb);(2, 4, lb);(3, 4, lb)]
            ub = 1.1
            bus["vm_pair_ub"] = [(1, 4, ub);(2, 4, ub);(3, 4, ub)]
            bus["vmax"] = ones(4)*1.5
            # bus["vmax"][end] = 0.05
            bus["vmin"] = zeros(4)
            # bus["grounded"] .=  0
        else
            @show bus
        end

    end

    for (g,gen) in math4w["gen"]
        gen["cost"] = 0.00001
    end

    for (d,load) in math4w["load"]
        load["pd"] .*= loadscale
        load["qd"] .*= loadscale
    end

    function add_gens!(math4w)
        gen_counter = 2
        for (d, load) in math4w["load"]
            if mod(load["index"], 4) == 1
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
    add_start_vrvi!(math4w)


    res_comp = solve_mc_vvvw_doe_competitive(math4w, ipopt)
    @show res_comp["termination_status"]
    @assert(res_comp["termination_status"]==LOCALLY_SOLVED || res_comp["termination_status"]==ALMOST_LOCALLY_SOLVED)
    res_comp_obj = round(res_comp["objective"], digits=2)
    pg_cost1 = [gen["pg_cost"] for (g,gen) in res_comp["solution"]["gen"] if g!="1"]
    pg_ref_comp = res_comp["solution"]["gen"]["1"]["pg"]
    qg_ref_comp = res_comp["solution"]["gen"]["1"]["qg"]


    # res_ms = solve_mc_vvvw_doe_mse(math4w, ipopt)
    # @assert(res_ms["termination_status"]==LOCALLY_SOLVED || res_ms["termination_status"]==ALMOST_LOCALLY_SOLVED)
    # pg_cost2 = [gen["pg_cost"] for (g,gen) in res_ms["solution"]["gen"] if g!="1"]
    # res_ms_obj = round(res_ms["objective"], digits=2)
    # pg_ref_ms = res_ms["solution"]["gen"]["1"]["pg"]


    # res_abs = solve_mc_vvvw_doe_abs(math4w, ipopt)
    # @assert(res_abs["termination_status"]==LOCALLY_SOLVED || res_abs["termination_status"]==ALMOST_LOCALLY_SOLVED)
    # pg_cost3 = [gen["pg_cost"] for (g,gen) in res_abs["solution"]["gen"] if g!="1"]
    # res_abs_obj = round(res_abs["objective"], digits=2)
    # pg_ref_abs = res_abs["solution"]["gen"]["1"]["pg"]


    res_log = solve_mc_vvvw_doe_log_fairness(math4w, ipopt)
    @show res_log["termination_status"]
    @assert(res_log["termination_status"]==LOCALLY_SOLVED || res_log["termination_status"]==ALMOST_LOCALLY_SOLVED)
    pg_cost5 = [gen["pg_cost"] for (g,gen) in res_log["solution"]["gen"] if g!="1"]
    res_log_obj = round(res_log["objective"], digits=2)
    pg_ref_log = res_log["solution"]["gen"]["1"]["pg"]
    qg_ref_log = res_log["solution"]["gen"]["1"]["qg"]

    res_eq_obj = 0
    res_eq = solve_mc_vvvw_doe_equal(math4w, ipopt)
    if res_eq["termination_status"]==LOCALLY_INFEASIBLE
        pg_cost4 = 0 .*pg_cost1
        res_eq_obj = 0
        pg_ref_eq = NaN
        qg_ref_eq = NaN
    else
        @assert(res_eq["termination_status"]==LOCALLY_SOLVED || res_eq["termination_status"]==ALMOST_LOCALLY_SOLVED)
        pg_cost4 = [gen["pg_cost"] for (g,gen) in res_eq["solution"]["gen"] if g!="1"]
        res_eq_obj = round(res_eq["objective"], digits=2)
        pg_ref_eq = res_eq["solution"]["gen"]["1"]["pg"]
        qg_ref_eq = res_eq["solution"]["gen"]["1"]["qg"]

    end

    bb = sortperm(pg_cost5)
    pg_cost1a = pg_cost1[bb]
    # pg_cost2a = pg_cost2[bb]
    # pg_cost3a = pg_cost3[bb]
    pg_cost4a = pg_cost4[bb]
    pg_cost5a = pg_cost5[bb]

    plot(pg_cost1a, linestyle=:dash, label="max. competitive $res_comp_obj - net. P $(round(sum(pg_ref_comp),digits=1)), Q $(round(sum(qg_ref_comp),digits=1))")
    # plot!(pg_cost2a, linestyle=:dot,  label="min. deviation squared $res_ms_obj - net. cons. $(round(sum(pg_ref_ms),digits=2))")
    # plot!(pg_cost3a, linestyle=:dashdotdot, label="min. absolute deviation $res_abs_obj - net. cons. $(round(sum(pg_ref_abs),digits=2))")
    plot!(pg_cost4a, linstyle=:dashdot, label="equal $res_eq_obj  - net. P $(round(sum(pg_ref_eq),digits=1)), Q $(round(sum(qg_ref_eq),digits=1))")
    plot!(pg_cost5a, linstyle=:solid, label="log fairness  - net. P $(round(sum(pg_ref_log),digits=1)), Q $(round(sum(qg_ref_log),digits=1))", legend=:bottomright)

    xlabel!("PV system id (-)")

    ylabel!("Export DOE (kW)")
    ylims!(0,5.1)
    xlims!(0.5,length(pg_cost1)+0.5)
    title!("Voltage of $vscale pu, load at $loadscale")
    savefig("VVVW$vscale load$loadscale.pdf")

# end
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
