using GPSTTopic82024
using PowerModelsDistribution
using Ipopt
ipopt = Ipopt.Optimizer

file = "data/LV30_315bus/Master.dss"

eng4w = parse_file(file, transformations=[transform_loops!,remove_all_bounds!])
eng4w["settings"]["sbase_default"] = 1
eng4w["voltage_source"]["source"]["rs"] *=0
eng4w["voltage_source"]["source"]["xs"] *=0

math4w = transform_data_model(eng4w, kron_reduce=false, phase_project=false)
add_start_vrvi!(math4w)

for (i,bus) in math4w["bus"]
    if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
        bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
        bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
        # bus["grounded"] .=  0
    end
end

for (g,gen) in math4w["gen"]
    gen["cost"] = 0.0
end

for (d,load) in math4w["load"]
    load["pd"] .*= 0.2
    load["qd"] .*= 0.2
end

function add_gens!(math4w)
    gen_counter = 2
    for (d, load) in math4w["load"]
        println("Load Index: $(load["index"]), Load Bus Internal: $(load["load_bus"])")
        if mod(load["index"], 4) == 1
            math4w["gen"]["$gen_counter"] = deepcopy(math4w["gen"]["1"])
            math4w["gen"]["$gen_counter"]["name"] = "$gen_counter"
            math4w["gen"]["$gen_counter"]["index"] = gen_counter
            math4w["gen"]["$gen_counter"]["cost"] = 1.0 #*math4w["gen"]["1"]["cost"]
            math4w["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
            math4w["gen"]["$gen_counter"]["pmax"] = 5*ones(3)
            math4w["gen"]["$gen_counter"]["pmin"] = 0.0*ones(3)
            math4w["gen"]["$gen_counter"]["connections"] = [1;2;3;4]
            gen_counter = gen_counter + 1
        end
    end
end
add_gens!(math4w)

math4w["gen"]["4"]["pmax"] = ones(3)
# pm4w = instantiate_mc_model(
#         math4w,
#         IVRENPowerModel,
#         GPSTTopic82024.build_mc_doe;
#         multinetwork=false,
#     )

res = GPSTTopic82024.solve_mc_doe_fair_pg_abs(math4w, ipopt)

pg_cost = [gen["pg_cost"] for (g,gen) in res["solution"]["gen"]]

v_mag = [hypot.(bus["vr"],bus["vi"]) for (b,bus) in res["solution"]["bus"]]

