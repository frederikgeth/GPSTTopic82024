using OpenDSSDirect
const ODSS = OpenDSSDirect
using DataFrames
using CSV
file = "Master.dss"
odssstring = """
clear
compile $file
solve
export voltages
"""

dss(odssstring)
voltages_file= "2bus_EXP_VOLTAGES.csv"

##
voltages_df = CSV.read(voltages_file, DataFrame)
voltages_df[!,:Bus] = lowercase.(voltages_df[!,:Bus])
voltages_df

vpn(pm,pa,nm,na) = hypot(pm*exp(im*pa/180*pi) - nm*exp(im*na/180*pi))

voltages_df[!,:Van] = vpn.(voltages_df[!,Symbol(" Magnitude1")], voltages_df[!,Symbol(" Angle1")], voltages_df[!,Symbol(" Magnitude4")], voltages_df[!,Symbol(" Angle4")])
voltages_df[!,:Vbn] = vpn.(voltages_df[!,Symbol(" Magnitude2")], voltages_df[!,Symbol(" Angle2")], voltages_df[!,Symbol(" Magnitude4")], voltages_df[!,Symbol(" Angle4")])
voltages_df[!,:Vcn] = vpn.(voltages_df[!,Symbol(" Magnitude3")], voltages_df[!,Symbol(" Angle3")], voltages_df[!,Symbol(" Magnitude4")], voltages_df[!,Symbol(" Angle4")])

##
using PowerModelsDistribution
using Ipopt
ipopt = Ipopt.Optimizer

vscale = 1
eng4w = parse_file(file, transformations=[transform_loops!,remove_all_bounds!])
# eng4w["conductor_ids"] = 1:4
eng4w["settings"]["sbase_default"] = 1
eng4w["voltage_source"]["source"]["rs"] *= 0  
eng4w["voltage_source"]["source"]["xs"] *= 0   
# eng4w["voltage_source"]["source"]["vm"] *=vscale

# reduce_line_series!(eng4w)
math4w = transform_data_model(eng4w, kron_reduce=false, phase_project=false)
# math4w["bus"]["190"]["grounded"] = Bool[0, 0, 0, 0]
# math4w["bus"]["190"]["vr_start"] = [1.0, -0.5, -0.5, 0.0]
# math4w["bus"]["190"]["vi_start"] = [0.0, -0.866025, 0.866025, 0.0]
# math4w["bus"]["190"]["va"] = [0.0, -2.0944, 2.0944, 0.0]
# math4w["bus"]["190"]["vm"] = [1.0, 1.0, 1.0, 0.0]
# math4w["bus"]["190"]["vmin"] = 0 .*[1.0, 1.0, 1.0, 1.0]
# math4w["bus"]["190"]["vmax"] = 2* [1.0, 1.0, 1.0, 1.0]
# math4w["bus"]["190"]["terminals"] = collect(1:4)
add_start_vrvi!(math4w)
res = solve_mc_opf(math4w, IVRENPowerModel, ipopt)
vr1 = res["solution"]["bus"]["1"]["vr"]
vi1 = res["solution"]["bus"]["1"]["vi"]
vr2 = res["solution"]["bus"]["2"]["vr"]
vi2 = res["solution"]["bus"]["2"]["vi"]
vr3 = res["solution"]["bus"]["3"]["vr"]
vi3 = res["solution"]["bus"]["3"]["vi"]

vm1 = hypot.(vr1, vi1).*400/sqrt(3)
vm2 = hypot.(vr2, vi2).*400/sqrt(3)
vm3 = hypot.(vr3, vi3).*400/sqrt(3)
va1 = angle.(vr1+im*vi1)*180/pi
va2 = angle.(vr2+im*vi2)*180/pi
va3 = angle.(vr3+im*vi3)*180/pi

# B2 voltages according to OpenDSS:
#  227.343, 231.458, 231.103, 3.55806, 
#      0.0,   119.9,  -120.0,    -0.8,  
@show transpose(vm2)
@show transpose(va2)