function rectifier(x,y,a;type="smooth", ϵ=0.0001)
    if type=="nonsmooth"
        return f = vals -> a*max(0, vals - x) + y
    elseif type=="smooth"
        #we use a numerically safe version from StatsFuns.jl, this avoid unnecessary over/underflow in the evaluation of the exp/log
        return f = vals -> a*ϵ*StatsFuns.log1pexp((vals-x)/ϵ) + y
    end
end

```
Inputs in PU voltage
outputs in pu power relative to apparent power rating
```
function voltvar_handle(;ϵ=0.0001, type="smooth")
    # V_vv = [195; 207; 220; 240; 258; 276]./230
    # Q_vv = [44; 44;   0;   0;  -60;  -60]./100

    r1pu = rectifier(207/230,44/100,-44/13*(230/100);type=type, ϵ=ϵ)
    r2pu = rectifier(220/230,0,+44/13*(230/100);type=type, ϵ=ϵ)
    r3pu = rectifier(240/230,0,-60/18*(230/100);type=type, ϵ=ϵ)
    r4pu = rectifier(258/230,0,+60/18*(230/100);type=type, ϵ=ϵ)
    vv_curve_pu(x) = r1pu(x) + r2pu(x) + r3pu(x) + r4pu(x)
end

```
Inputs in PU voltage
outputs in pu power relative to apparent power rating
```
function voltwatt_handle(;ϵ=0.0001, type="smooth")
    # V_vw = [195; 253; 260; 276]./230
    # P_vw = [100; 100; 20; 20]./100

    relupu = rectifier(253/230,1,-80/7*(230/100);type=type, ϵ=ϵ)
    relu2pu =rectifier(260/230,0,+80/7*(230/100);type=type, ϵ=ϵ)
    vw_curvepu(x) = relupu(x) + relu2pu(x)
end