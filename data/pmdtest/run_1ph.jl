using OpenDSSDirect
const ODSS = OpenDSSDirect
using DataFrames
using CSV
file = "1ph.dss"
odssstring = """
clear
compile $file
solve
export voltages
"""

dss(odssstring)
voltages_file= "Test_EXP_VOLTAGES.csv"

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
eng4w["conductor_ids"] = 1:4
eng4w["settings"]["sbase_default"] = 1
eng4w["voltage_source"]["source"]["rs"] =zeros(3,3)
eng4w["voltage_source"]["source"]["xs"] =zeros(3,3)
eng4w["voltage_source"]["source"]["vm"] *=vscale

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

vm1 = hypot.(vr1, vi1).*400/sqrt(3)
vm2 = hypot.(vr2, vi2).*400/sqrt(3)
va1 = angle.(vr1+im*vi1)
va2 = angle.(vr2+im*vi2)

# B2 voltages according to OpenDSS without grounding on the load bus
# "va": {
#     "1": -0.10546374855838672
#     "2": -2.099194996742148,
#     "3": 2.0694633590902627,
#     "4": 0.4890829587275647,
#   },
#   "vm": {
#     "1": 193.70297223692225
#     "2": 236.99468209783933,
#     "3": 233.05331516232354,
#     "4": 43.40160528580862,


@show transpose(vm1)
@show transpose(va1)


vm1 = hypot.(vr1, vi1).*400/sqrt(3)
vm2 = hypot.(vr2, vi2).*400/sqrt(3)
va1 = angle.(vr1+im*vi1)*180/pi
va2 = angle.(vr2+im*vi2)*180/pi
# B2 voltages according to OpenDSS with 1 Ohm neutral grounding on the load bus
# "va": {
#     "1": -
#     "2": 
#     "3": 
#     "4": 
#   },
#   "vm": {
#     "1": 196.956
#     "2": 248.896,
#     "3": 222,
#     "4":  26.9617,

