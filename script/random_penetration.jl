using Random
using DataFrames
r = rand(Float64, 200)

r10 = r.<0.1
r20 = r.<0.2
r30 = r.<0.3
r40 = r.<0.4
r50 = r.<0.5
r60 = r.<0.6
r70 = r.<0.7
r80 = r.<0.8
r90 = r.<0.9

@show sum(r10)/2
@show sum(r20)/2
@show sum(r30)/2
@show sum(r40)/2
@show sum(r50)/2
@show sum(r60)/2
@show sum(r70)/2
@show sum(r80)/2
@show sum(r90)/2

df = DataFrame(r10=r10, r20=r20, r30=r30, r40=r40, r50=r50, r60=r60, r70=r70, r80=r80, r90=r90)
# CSV.write("penetration_samples.csv", df)

#
penetration_samples = CSV.read("data/penetration_samples.csv", DataFrame)
pen = "r".*string.(collect(10:10:90))
pen =  ["r10", "r20", "r30", "r40", "r50", "r60", "r70", "r80", "r90"]