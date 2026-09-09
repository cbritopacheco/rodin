import numpy as np
# (t, |qIn|) straight from the run log, cycle 2
log = """0.781 0.00233662
0.800 0.00219964
0.821 0.00206554
0.841 0.00195194
0.861 0.00184995
0.881 0.00175785
0.901 0.00167438
0.921 0.00162667
0.941 0.00162908
0.961 0.00167091
0.981 0.00173929
1.001 0.00181993
1.021 0.00189828
1.041 0.00196079
1.061 0.00199612
1.081 0.00199595
1.101 0.00195544
1.121 0.00187352
1.141 0.00178060
1.161 0.00169713
1.181 0.00162175
1.201 0.00158288
1.221 0.00169830
1.241 0.00194746
1.261 0.00227955
1.281 0.00263529
1.301 0.00295538
1.321 0.00318900
1.341 0.00329931
1.361 0.00326622
1.381 0.00308694
1.401 0.00286853
1.421 0.00268027
1.441 0.00251558
1.461 0.00241603
1.481 0.00244301
1.501 0.00254599
1.521 0.00265979
1.541 0.00272175
1.561 0.00268627
1.581 0.00254126"""
d = np.array([[float(x) for x in l.split()] for l in log.strip().split("\n")])
t, q = d[:,0], d[:,1]
phi = t - 0.8

def load(f):
    return np.array([[float(x) for x in l.split()]
                     for l in open(f) if l.strip() and not l.startswith('#')])
import sys
base = sys.argv[1]
pv = load(base+'/presion_PV_SR.dat'); mv = load(base+'/presion_MV_SR.dat')
tw, dpw = pv[:,0], pv[:,1]-mv[:,1]
dp = np.interp(phi, tw, dpw)

# fit  L dq/dt + R q = dp(t)
dqdt = np.gradient(q, phi)
A = np.column_stack([dqdt, q])
(L, R), *_ = np.linalg.lstsq(A, dp, rcond=None)
pred = L*dqdt + R*q
ss = 1 - np.sum((dp-pred)**2)/np.sum((dp-dp.mean())**2)
print("fit of  L dq/dt + R q = dp(t)  to the logged flow")
print("   L = %.0f Pa s^2/m^2      R = %.0f Pa s/m^2      R^2 = %.3f"%(L, R, ss))
print("   inertial time constant  tau = L/R = %.2f s   vs cycle T = 0.80 s   -> %.1f cycles to coast down"%(L/R, (L/R)/0.8))
print()
closed = np.interp(phi, tw, (dpw<=1e-9).astype(float)) > 0.5
# per-sample width, so the masked integral does not bridge the gaps between
# the closed windows
w = np.gradient(phi)
vol_closed = np.sum(q[closed]*w[closed])
vol_tot = np.sum(q*w)
print("flow during the CLOSED windows")
print("   |q| range while closed : %.5f - %.5f m^2/s   (never reaches zero)"%(q[closed].min(), q[closed].max()))
print("   |q| range while open   : %.5f - %.5f m^2/s"%(q[~closed].min(), q[~closed].max()))
print("   ratio of means         : %.2f  (closed / open)"%(q[closed].mean()/q[~closed].mean()))
print("   valve shut %.0f%% of the sampled cycle (32%% exactly, from the .dat),\n   and %.0f%% of the transported volume still crosses the MV while it is shut"
      %(100*closed.mean(), 100*vol_closed/vol_tot))
print()
i_dp, i_q = np.argmax(dp), np.argmax(q)
print("phase lag: dp peaks at phi = %.3f s, q peaks at phi = %.3f s  -> lag %.3f s = %.0f deg"
      %(phi[i_dp], phi[i_q], phi[i_q]-phi[i_dp], 360*(phi[i_q]-phi[i_dp])/0.8))
