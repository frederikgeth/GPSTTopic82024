using StatsFuns
using Plots
ϵ = 0.0001
r = -2:0.01:2.0
y(x,ϵ) = ϵ.*log.(1.0 .+ exp.(x./ϵ))
plot(r, y(r,ϵ))
plot!(r, ϵ.*log1pexp.(r./ϵ), linestyle=:dot)
#now you see that the function evalution is unreliable!


V_vw = [195; 253; 260; 276]./230
P_vw = [100; 100; 20; 20]./100

V_vv = [195; 207; 220; 240; 258; 276]./230
Q_vv = [44; 44;   0;   0;  -60;  -60]./100

vw = plot(V_vw, P_vw)
title!("Volt-Watt characteristic")
xlabel!("Voltage (V pu)")
ylabel!("Active power (W pu)")

vv = plot(V_vv, Q_vv)
title!("Volt-var characteristic")
xlabel!("Voltage (V pu)")
ylabel!("Reactive power (var pu)")


function rectifier(x,y,a;type="smooth", ϵ=10^-4)
    if type=="nonsmooth"
        return f = vals -> a*max(0, vals - x) + y
    elseif type=="smooth"
        return f = vals -> a*ϵ*log1pexp((vals-x)/ϵ) + y
    end
end

ϵ=1
type = "smooth"
vwplot = plot(legend=:bottomleft,ylims=(0,100))
for ϵ in [3,2,1,0.5,0.01]
    relu = rectifier(253,100,-80/7;type=type, ϵ=ϵ)
    relu2 =rectifier(260,0,+80/7;type=type, ϵ=ϵ)
    vw_curve(x) = relu(x) + relu2(x)
    rr =185:0.1:276
    plot!(rr,vw_curve.(rr), label="epsilon=$ϵ")
end
title!("Volt-Watt characteristic")
xlabel!("Voltage (V)")
ylabel!("Active power (W % pu)")
xticks!([207, 230, 253, 260])
yticks!([100,20,0])
savefig("voltwatt.pdf")


vwplotpu = plot(legend=:bottomleft,ylims=(0,1))
relupu = rectifier(253/230,1,-80/7*(230/100);type=type, ϵ=ϵ)
relu2pu =rectifier(260/230,0,+80/7*(230/100);type=type, ϵ=ϵ)
vw_curvepu(x) = relupu(x) + relu2pu(x)
rr =(185:0.1:276)./230
plot!(rr,vw_curvepu.(rr), label="epsilon=$ϵ")
title!("Volt-Watt characteristic")
xlabel!("Voltage (V pu)")
ylabel!("Active power (W pu)")
xticks!([0.9, 1.0, 1.1, 1.13])
yticks!([1,0.2,0])
savefig("voltwattpu.pdf")

##
vvplot = plot(legend=:bottomleft)
for ϵ in [3,2,1,0.5,0.01]
    r1 = rectifier(207,44,-44/13;type=type, ϵ=ϵ)
    r2 = rectifier(220,0,+44/13;type=type, ϵ=ϵ)
    r3 = rectifier(240,0,-60/18;type=type, ϵ=ϵ)
    r4 = rectifier(258,0,+60/18;type=type, ϵ=ϵ)
    vv_curve(x) = r1(x) + r2(x) + r3(x) + r4(x)
    rr =185:0.1:276
    plot!(rr,vv_curve.(rr), label="epsilon=$ϵ")
end
title!("Volt-var characteristic")
xlabel!("Voltage (V)")
ylabel!("Reactive power (var % pu)")
xticks!([207,220,230, 240, 258])
yticks!([44,0,-60])
savefig("voltvar.pdf")


plot(vvplot, vwplot,   layout=(2,1), size=(600,800))
savefig("voltvarwatt.pdf")
##
ϵ=0.001
vvplotpu = plot(legend=:bottomleft)
r1pu = rectifier(207/230,44/100,-44/13*(230/100);type=type, ϵ=ϵ)
r2pu = rectifier(220/230,0,+44/13*(230/100);type=type, ϵ=ϵ)
r3pu = rectifier(240/230,0,-60/18*(230/100);type=type, ϵ=ϵ)
r4pu = rectifier(258/230,0,+60/18*(230/100);type=type, ϵ=ϵ)
vv_curve_pu(x) = r1pu(x) + r2pu(x) + r3pu(x) + r4pu(x)
rr =(185:0.1:280)./230
plot!(rr,vv_curve_pu.(rr), label="epsilon=$ϵ")
title!("Volt-var characteristic")
xlabel!("Voltage (V pu)")
ylabel!("Reactive power (var pu)")
xticks!([0.9,0.96,1.0, 1.04, 1.12])
yticks!([44,0,-60]./100)
savefig("voltvarpu.pdf")

##
using JuMP
using Ipopt

m = JuMP.Model(Ipopt.Optimizer)
register(m, :vv_curve_pu, 1, vv_curve_pu; autodiff = true)
@variable(m, 0.9 <=x<=1.04)
@variable(m, y11)
@NLconstraint(m, vv_curve_pu(x-3) == y11)
@objective(m, Min, y11)
optimize!(m)
@show value.(x)
@show value.(y11)
termination_status(m)
scatter!([value.(x)],[value.(y11)],markersize=5 )
##
m2 = JuMP.Model(Ipopt.Optimizer)
@variable(m2, x2>=xmin)
@variable(m2, y22>=0)
@NLconstraint(m2, ϵ*log(1+exp(x2/ϵ)) == y22)
@objective(m2, Min, y22)
optimize!(m2)
@show value.(x2)
@show value.(y22)
termination_status(m2)

