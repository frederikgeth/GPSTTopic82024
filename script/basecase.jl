using GPSTTopic82024
using PowerModelsDistribution
using Ipopt
using Plots
ipopt = Ipopt.Optimizer

## Main loop
case = "LV9_258bus"
casepath = "data/$case"
file = "$casepath/Master.dss"
busdistancesdf = CSV.read("$casepath/busdistances.csv", DataFrame)
busdistances_dict = Dict(busdistancesdf.Bus[i] => busdistancesdf.busdistances[i] for i in 1:nrow(busdistancesdf))


vscale = 1
loadscale = 0.
# for vscale in 0.98:0.01:1.07, loadscale in [1] #0.8:0.05:1.0
    # for vscale in 0.98:0.01:1.07, loadscale in [1] #0.8:0.05:1.0
    eng4w = parse_file(file, transformations=[transform_loops!,remove_all_bounds!])
    eng4w["conductor_ids"] = 1:4
    eng4w["settings"]["sbase_default"] = 1
    eng4w["voltage_source"]["source"]["rs"] =zeros(3,3)
    eng4w["voltage_source"]["source"]["xs"] =zeros(3,3)
    eng4w["voltage_source"]["source"]["vm"] *=vscale

    reduce_line_series!(eng4w)
    math4w = transform_data_model(eng4w, kron_reduce=false, phase_project=false)
    math4w["bus"]["190"]["grounded"] = Bool[0, 0, 0, 0]
    math4w["bus"]["190"]["vr_start"] = [1.0, -0.5, -0.5, 0.0]
    math4w["bus"]["190"]["vi_start"] = [0.0, -0.866025, 0.866025, 0.0]
    math4w["bus"]["190"]["va"] = [0.0, -2.0944, 2.0944, 0.0]
    math4w["bus"]["190"]["vm"] = [1.0, 1.0, 1.0, 0.0]
    math4w["bus"]["190"]["vmin"] = 0 .*[1.0, 1.0, 1.0, 1.0]
    math4w["bus"]["190"]["vmax"] = 2* [1.0, 1.0, 1.0, 1.0]
    math4w["bus"]["190"]["terminals"] = collect(1:4)
    
    
    add_start_vrvi!(math4w)

    for (i,bus) in math4w["bus"]
        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            # bus["vm_pair_lb"] = [(1, 4, 0.7);(2, 4, 0.7);(3, 4, 0.7)]
            # bus["vm_pair_ub"] = [(1, 4, 1.5);(2, 4, 1.5);(3, 4, 1.5)]
            bus["vmax"] = ones(4)*2
            # bus["vmax"][end] = 2
            # bus["vmin"] = zeros(4)
            # bus["grounded"] .=  0
        else
            @show bus
        end

    end

    for (g,gen) in math4w["gen"]
        math4w["gen"]["1"]["cost"] = [1, 1, 1]
        s = 100
        gen["pmin"] = -s*ones(3)
        gen["pmax"] = s*ones(3)
        gen["qmin"] = -s*ones(3)
        gen["qmax"] = s*ones(3)
        gen["connections"] = collect(1:4)
    end

    for (d,load) in math4w["load"]
        load["pd"] .*= loadscale
        load["qd"] .*= loadscale
    end

    # function add_gens!(math4w)
    #     gen_counter = 2
    #     for (d, load) in math4w["load"]
    #         if mod(load["index"], 200) == 1
    #             # phases = 3
    #             phases = length(load["connections"])-1
    #             math4w["gen"]["$gen_counter"] = deepcopy(math4w["gen"]["1"])
    #             math4w["gen"]["$gen_counter"]["name"] = "$gen_counter"
    #             math4w["gen"]["$gen_counter"]["index"] = gen_counter
    #             math4w["gen"]["$gen_counter"]["cost"] = 1.0 #*math4w["gen"]["1"]["cost"]
    #             math4w["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
    #             math4w["gen"]["$gen_counter"]["pmax"] = 5.0*ones(phases)
    #             math4w["gen"]["$gen_counter"]["pmin"] = 0.0*ones(phases)
    #             math4w["gen"]["$gen_counter"]["qmax"] = 5.0*ones(phases)
    #             math4w["gen"]["$gen_counter"]["qmin"] = -5.0*ones(phases)
    #             math4w["gen"]["$gen_counter"]["connections"] = load["connections"]
    #             gen_counter = gen_counter + 1
    #         end
    #     end
    # end
    # add_gens!(math4w)
    
    res_comp = solve_mc_vvvw_opf(math4w, ipopt)
    res = solve_mc_opf(math4w, IVRENPowerModel, ipopt)

    # res_comp = solve_mc_vvvw_opf(math4w, ipopt)
    # @assert(res_comp["termination_status"]==LOCALLY_SOLVED || res_comp["termination_status"]==ALMOST_LOCALLY_SOLVED)
    # res_comp_obj = round(res_comp["objective"], digits=2)
    # pg_cost1 = [gen["pg_cost"] for (g,gen) in res_comp["solution"]["gen"] if g!="1"]

res = res_comp
v_mag = stack([hypot.(bus["vr"][1:4],bus["vi"][1:4]) for (b,bus) in res["solution"]["bus"]], dims=1)
plot(v_mag, label=["a" "b" "c" "n"])
plot!([0; length(res["solution"]["bus"])], [0.9; 0.9], label="vmin")
plot!([0; length(res["solution"]["bus"])], [1.1; 1.1], label="vmax")
ylabel!("V (pu)")
xlabel!("bus id (-)")