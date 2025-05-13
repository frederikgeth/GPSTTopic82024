using OpenDSSDirect
const ODSS = OpenDSSDirect
using DataFrames
using CSV
using Plots
file = "Master.dss"
odssstring = """
clear
compile $file
solve
export voltages
"""

dss(odssstring)
voltages_file= "lv9_258bus_EXP_VOLTAGES.csv"

##
voltages_df = CSV.read(voltages_file, DataFrame)
voltages_df[!,:Bus] = lowercase.(voltages_df[!,:Bus])
voltages_df

vpn(pm,pa,nm,na) = hypot(pm*exp(im*pa/180*pi) - nm*exp(im*na/180*pi))

voltages_df[!,:Van] = vpn.(voltages_df[!,Symbol(" Magnitude1")], voltages_df[!,Symbol(" Angle1")], voltages_df[!,Symbol(" Magnitude4")], voltages_df[!,Symbol(" Angle4")])
voltages_df[!,:Vbn] = vpn.(voltages_df[!,Symbol(" Magnitude2")], voltages_df[!,Symbol(" Angle2")], voltages_df[!,Symbol(" Magnitude4")], voltages_df[!,Symbol(" Angle4")])
voltages_df[!,:Vcn] = vpn.(voltages_df[!,Symbol(" Magnitude3")], voltages_df[!,Symbol(" Angle3")], voltages_df[!,Symbol(" Magnitude4")], voltages_df[!,Symbol(" Angle4")])


##
busdistances = ODSS.Circuit.AllBusDistances()
busnames = ODSS.Circuit.AllBusNames()
busdistancesdf = DataFrame(Bus=busnames, busdistances=busdistances)
sort!(busdistancesdf, [:busdistances])
busdist_dict = Dict(busdistancesdf.Bus[i] => busdistancesdf.busdistances[i] for i in 1:nrow(busdistancesdf))
CSV.write("busdistances.csv", busdistancesdf)


##
voltages_df = innerjoin(voltages_df, busdistancesdf, on=:Bus)
sort!(voltages_df,[:busdistances])

##
scatter(voltages_df.busdistances, voltages_df[!,Symbol(" Magnitude1")], label="phase a")
scatter!(voltages_df.busdistances, voltages_df[!,Symbol(" Magnitude2")], label="phase b")
scatter!(voltages_df.busdistances, voltages_df[!,Symbol(" Magnitude3")], label="phase c")
xlabel!("Electrical distance from distribution tranformer (km)")
ylabel!("Voltage magnitude to ground (V)")
savefig("volt_distance_drop_ground.pdf")

##
scatter(voltages_df.busdistances, voltages_df.Van, label="phase a")
scatter!(voltages_df.busdistances, voltages_df.Vbn, label="phase b")
scatter!(voltages_df.busdistances, voltages_df.Vcn, label="phase c")
xlabel!("Electrical distance from distribution tranformer (km)")
ylabel!("Voltage magnitude to neutral (V)")
savefig("volt_distance_drop_neutral.pdf")

##
scatter(voltages_df.busdistances, voltages_df[!,Symbol(" Magnitude4")],label="neutral")
xlabel!("Electrical distance from distribution tranformer (km)")
ylabel!("Neutral voltage rise (V)")
savefig("volt_distance_drop_ground_neutral.pdf")


