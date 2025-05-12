using Distributions
using OpenDSSDirect
using PowerModelsDistribution
using Ipopt
using Plots
using GPSTTopic82024
# include("./../src/GPSTTopic82024.jl")
include("./../notebooks/PerPhasePlot.jl")
using GPSTTopic82024
file = "data/LV30_315bus/Master.dss"
# file = "C:\\Users\\Frederik Geth\\Documents\\GitHub\\GPSTTopic82024\\data\\LV30_315bus\\Master.dss"
simulatePvAgainstDoe(;file)
function simulatePvAgainstDoe(; file::String="data/LV30_315bus/Master.dss")
    
    pvBusses, passiveBusses = assignPvLoadBusses(file, 2)
    pvBusExportLims = runDoeSimulation(file, pvBusses)

    OpenDSSDirect.dss("""
	    clear
	    compile $(file)
	    closedi
	""")

    vLimsPu = [0.94, 1.1]
    voltageBase = 230.0

    genToBusNames = addPv(pvBusses, 17.0, 0.23, 17.0)
     
    samplePv(Uniform(0.5, 1.5))
    OpenDSSDirect.Solution.Solve()

    loadBusVoltages = collectAllLoadBusVMagAngle()
    checkVoltageLimits(loadBusVoltages, voltageBase, vLimsPu)

    lineCurrents = collectAllLineCurrentMagAngle()
    checkCurrentLimits(lineCurrents) 

    plot_3ph_load()

    assertDoeLims(genToBusNames, pvBusExportLims)

    OpenDSSDirect.Solution.Solve() 

    plot_3ph_load()
    loadBusVoltages = collectAllLoadBusVMagAngle()
    checkVoltageLimits(loadBusVoltages, voltageBase, vLimsPu)

    lineCurrents = collectAllLineCurrentMagAngle()
    checkCurrentLimits(lineCurrents) 
end

function runDoeSimulation(file::String, pvLoadBusses::Vector{String})

    eng4wModel = parse_file(file, transformations=[transform_loops!,remove_all_bounds!])
    eng4wModel["settings"]["sbase_default"] = 1
    math4wModel = transform_data_model(eng4wModel, kron_reduce=false, phase_project=false)
    add_start_vrvi!(math4wModel)

    for (_, bus) in math4wModel["bus"] # assign phase to neutral voltage bounds
        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
        end
    end

    for (_, gen) in math4wModel["gen"]
        gen["cost"] = 0.0
    end

    for (pvIndex, loadBus) in enumerate(pvLoadBusses)
        loadBusAndPhase = split(loadBus, '.')
        loadBusName = loadBusAndPhase[1]
        connections = parse.(Ref(Int64), loadBusAndPhase[2:end])
        threePhBinary = zeros(3)
        (conn -> if conn < 4 threePhBinary[conn] = 1 end).(connections) # set binary vector based on connections, omitting neutral: connections = [1,3,4] -> [1,0,1]

        interalPvBus = math4wModel["bus_lookup"][loadBusName]
        genOffs = pvIndex + 1

        math4wModel["gen"]["$genOffs"] = deepcopy(math4wModel["gen"]["1"])
        math4wModel["gen"]["$genOffs"]["name"] = "$genOffs"
        math4wModel["gen"]["$genOffs"]["index"] = genOffs
        math4wModel["gen"]["$genOffs"]["cost"] = 1.0 
        math4wModel["gen"]["$genOffs"]["gen_bus"] = interalPvBus
        math4wModel["gen"]["$genOffs"]["pmax"] = 5.0*ones(3) # TODO: tunable inputs
        math4wModel["gen"]["$genOffs"]["pmin"] = 0.0*ones(3)
        math4wModel["gen"]["$genOffs"]["connections"] = connections
    end
    
    math4wModel["gen"]["4"]["pmax"] = ones(3) # TODO: required for new testcase?
    math4wModel["gen"]["4"]["connections"] = [1;2;3;4]

    res = GPSTTopic82024.solve_mc_doe_fair_pg_abs(math4wModel, Ipopt.Optimizer)

    display(res["solution"]["gen"])
    dssLoadBusToLims = Dict()
    for (genName, genData) in res["solution"]["gen"]
        genNum = parse(Int64, genName)
        pvGen = genNum - 1
        if pvGen < 1
            continue
        end
        dssLoadBusName = pvLoadBusses[pvGen]
        dssLoadBusToLims[dssLoadBusName] = genData
    end
    
    return dssLoadBusToLims
    
    # pg_cost = [gen["pg_cost"] for (_, gen) in res["solution"]["gen"]]
    # v_mag = [hypot.(bus["vr"],bus["vi"]) for (_, bus) in res["solution"]["bus"]] TODO: useful voltage output format
end

function assignPvPluto(testNetwork, pvRate)
    pvBusses, passiveBusses = assignPvLoadBusses(testNetwork, pvRate)
	cd("../../")
	outputStr = "Assigned OpenDSS PV Busses: ($(length(pvBusses)))\n" 
	for i in 1:convert(Int, ceil(length(pvBusses)/5))
		startWindow = (i - 1) * 5 + 1
		endWindow = i * 5
		if endWindow > length(pvBusses)
		    outputStr *= "$(pvBusses[startWindow:end]) \n"
		else
			outputStr *= "$(pvBusses[startWindow:endWindow]) \n"
		end
	end
	outputStr *= "Assigned OpenDSS Passive Busses: ($(length(passiveBusses)))"
	for i in 1:convert(Int, ceil(length(passiveBusses)/5))
		startWindow = (i - 1) * 5 + 1
		endWindow = i * 5
		if endWindow > length(passiveBusses)
			outputStr *= "$(passiveBusses[startWindow:end])\n"
		else
			outputStr *= "$(passiveBusses[startWindow:endWindow])\n"
		end
	end
    return pvBusses, passiveBusses
end

function assignPvLoadBusses(file::String, loadPerPv::Int64)::Tuple{Vector{String}, Vector{String}}
    
    OpenDSSDirect.dss("""
	    clear
	    compile $(file)
	    closedi
	""")
    loadNames = OpenDSSDirect.Loads.AllNames()
    assignedPvLoads = loadNames[1:loadPerPv:end]
    
    pvLoadBusNames = Vector{String}(undef, length(assignedPvLoads))
    passiveLoadBusNames = Vector{String}(undef, length(loadNames) - length(assignedPvLoads))

    load = OpenDSSDirect.Loads.First()
    pvIndex = 1
    passiveIndex = 1
    while load > 0
        loadName = OpenDSSDirect.Loads.Name()
        loadBus = OpenDSSDirect.Properties.Value("bus1")
        if loadName in assignedPvLoads
            pvLoadBusNames[pvIndex] = loadBus
            pvIndex += 1
        else 
            passiveLoadBusNames[passiveIndex] = loadBus
            passiveIndex += 1
        end
        load = OpenDSSDirect.Loads.Next()
    end
    return pvLoadBusNames, passiveLoadBusNames
end

function addPv(pvRate::Int64)
    load = OpenDSSDirect.Loads.First()
    loadNum = 1
    pvSystems = ""
    connectedBusses = []
    while load > 0
        println("Load Num $loadNum")
        if load % pvRate == 0
            loadName = "$(OpenDSSDirect.Loads.Name())"
            pvName = "PV_$(loadName)"
            
            OpenDSSDirect.Circuit.SetActiveElement("Load.$loadName")
            connectedBus = OpenDSSDirect.Properties.Value("bus1")
            phases = OpenDSSDirect.Loads.Phases()
            kV = OpenDSSDirect.Loads.kV()
            kVA = OpenDSSDirect.Loads.kVABase()

            println("Assign PV for load num $loadNum, named $(OpenDSSDirect.Loads.Name())")
            
            pvSystems *= "New PVSystem.$(pvName) phases=$(phases) bus1=$(connectedBus) kV=$(kV) kVA=$(kVA) irrad=1 Pmpp=$(kVA*0.9) temperature=25\n"
            push!(connectedBusses, connectedBus) 
        end
        loadNum += 1
        load = OpenDSSDirect.Loads.Next()
    end
    OpenDSSDirect.dss(pvSystems)

    return connectedBusses
end

function addPv(loadBusNames::Vector{String}, kvaBase::Float64, kv::Float64, pmpp::Float64)::Dict{String, String}
    load = OpenDSSDirect.Loads.First()
    pvSystems = ""
    genToBusNames = Dict{String, String}()
    while load > 0
        loadBus = OpenDSSDirect.Properties.Value("bus1")
        if loadBus in loadBusNames
            loadName = OpenDSSDirect.Loads.Name()
            pvName = lowercase("PV_$(loadName)")
            
            phases = OpenDSSDirect.Loads.Phases()
            # kV = OpenDSSDirect.Loads.kV()
            # kVA = OpenDSSDirect.Loads.kVABase()

            pvSystems *= "New PVSystem.$(pvName) phases=$(phases) bus1=$(loadBus) kV=$(kv) kVA=$(kvaBase) irrad=1 Pmpp=$(pmpp) temperature=25\n"
            genToBusNames[pvName] = loadBus  
        end
        load = OpenDSSDirect.Loads.Next()
    end
    OpenDSSDirect.dss(pvSystems)

    return genToBusNames
end

function samplePv(irradianceRange::Uniform)
    pv = OpenDSSDirect.PVsystems.First()
    pvNum = 1
    while pv > 0
        irradianceSample = rand(irradianceRange)
        pvName = OpenDSSDirect.PVsystems.Name()
        OpenDSSDirect.PVsystems.Irradiance(irradianceSample)

        pvNum += 1
        pv = OpenDSSDirect.PVsystems.Next()
    end
end

function assertDoeLims(genToBusNames::Dict, busToDoe::Dict)

    pv = OpenDSSDirect.PVsystems.First()

    while pv > 0
        pvName = OpenDSSDirect.PVsystems.Name()
        if !haskey(genToBusNames, lowercase(pvName))
            pv = OpenDSSDirect.PVsystems.Next()
            continue
        end

        busName = genToBusNames[pvName]
        connectionLims = busToDoe[busName]["pg"]
        totalExportLim = sum(connectionLims) # this can only be done on the assumption the PV connections are the same across
        exportKw = OpenDSSDirect.PVsystems.kW()

        println("Export Comparison Bus $busName: $totalExportLim kW, Measured: $exportKw kW")

        if exportKw > totalExportLim
            pmpp = OpenDSSDirect.PVsystems.Pmpp()
            newIrrad = getIrradForDc(totalExportLim, pmpp)
            OpenDSSDirect.PVsystems.Irradiance(newIrrad)
            println("Export Correction Bus $busName, new irradiance assigned: $newIrrad")
        end
        pv = OpenDSSDirect.PVsystems.Next()
        
        # busConnIncN = parse.(Ref(Int64), split(busName, '.')[2:end])
        # busConn = filter(x -> (x != 4), busConnIncN)
    end
end

function getIrradForDc(desiredDc::Float64, pmpp::Float64; temp::Int64=25)
    if temp != 25
        println("In getting Irradiance for DC PV value getIrradForDc(), values for temp other than 25 not implemented")
    end # TODO implement for other temps
    return desiredDc/pmpp
end


function collectAllLoadBusVMagAngle()
    
    loadBusVoltages = Dict()
    load = OpenDSSDirect.Loads.First()

    loadBusNames = []
    loadNum = 1
    while load > 0
        loadName = OpenDSSDirect.Loads.Name()
        loadBus = OpenDSSDirect.Properties.Value("bus1")
        push!(loadBusNames, loadBus)

        loadNum += 1
        load = OpenDSSDirect.Loads.Next()
    end

    for loadBus in loadBusNames
        loadBusPhases = split(loadBus, ".")
        
        loadBusName = loadBusPhases[1]

        OpenDSSDirect.Circuit.SetActiveBus(convert(String, loadBusName))
        
        vmagAngle = OpenDSSDirect.Bus.VMagAngle()
        loadBusVoltages[loadBusName] = vmagAngle
    end
    return loadBusVoltages
end

function collectAllLineCurrentMagAngle()

    lineCurrents = Dict()
    line = OpenDSSDirect.Lines.First()
    while line > 0
        
        lineName = OpenDSSDirect.Lines.Name()
        currentMagAng = OpenDSSDirect.CktElement.CurrentsMagAng()
        lineCurrents[lineName] = currentMagAng

        line = OpenDSSDirect.Lines.Next()
    end
    return lineCurrents
end

function checkCurrentLimits(currentMagAngles::Dict)
    line = OpenDSSDirect.Lines.First()
    while line > 0
        
        lineName = OpenDSSDirect.Lines.Name()
        currentLim = OpenDSSDirect.Lines.NormAmps()
        currentMags = currentMagAngles[lineName][1:2:end]
        phaseBreaches = findall(currentMags .> currentLim)
        if length(phaseBreaches) > 0
            println("Line $lineName breached current limits. Phases $(phaseBreaches), currents $(currentMags)")
        end
        line = OpenDSSDirect.Lines.Next()
    end
end

function checkVoltageLimits(vmagAngles::Dict, ratedVoltage::Float64, lims::Vector{Float64}; puLims::Bool=true)
    lowerVLim = lims[1]
    upperVLim = lims[2]
    if puLims
        lowerVLim *= ratedVoltage
        upperVLim *= ratedVoltage
    end

    for (busname, vMagAngle) in vmagAngles
        voltageMags =  vMagAngle[1:2:(end-2)] # end - 2 as want to ignore neutral voltage for pu check

        indexes = findall((x -> x < lowerVLim || x > upperVLim).(voltageMags))
        if length(indexes) > 0
            println("Voltage Breach for Bus $busname, phases $(indexes), voltages: $(voltageMags).")
        end
    end
end

