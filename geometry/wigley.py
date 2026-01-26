import rhinoscriptsyntax as rs
# -----------
# PARAMETERS
# -----------
L = 3.0
B = 0.4
T = 0.3
# Sampling
nx   = 16
nz   = 12
npts = 41
#
rs.EnableRedraw(False)
#
x_list = [ (-L/2.0) + L * i/(nx-1.0) for i in range(nx) ]
z_list = [ (-T) + T * i/(nz-1.0) for i in range(nz) ]
# ---------
# Stations
# ---------
station_curves = []
for xi in x_list:
    pts = []
    for k in range(npts):
        zk = (-T) + T * k/(npts-1.0)
        # Wigley half-breadth (one side)
        yk = (B/2.0) * (1.0 - (2.0*xi/L)**2) * (1.0 - (zk/T)**2)
        if yk < 0: yk = 0.0
        pts.append((xi, yk, zk))
    crv = rs.AddInterpCurve(pts, degree=3)
    if crv: station_curves.append(crv)
# -----------
# Waterlines
# -----------
waterline_curves = []
for zi in z_list:
    pts = []
    for k in range(npts):
        xk = (-L/2.0) + L * k/(npts-1.0)
        yk = (B/2.0) * (1.0 - (2.0*xk/L)**2) * (1.0 - (zi/T)**2)
        if yk < 0: yk = 0.0
        pts.append((xk, yk, zi))
    crv = rs.AddInterpCurve(pts, degree=3)
    if crv: waterline_curves.append(crv)
# ----------
# Half Hull
# ----------
all_crvs = station_curves + waterline_curves
side_srf = rs.AddNetworkSrf(all_crvs)
#
for cid in all_crvs:
    rs.DeleteObject(cid)
#
if not side_srf:
    rs.EnableRedraw(True)
    print("Network surface creation failed. Try increasing nx/nz/npts.")
    raise SystemExit
#--------
# Mirror
#--------
p0 = (0,0,0)
p1 = (1,0,0)
p2 = (0,0,1)
#
mir_srf = rs.MirrorObject(side_srf, p0, p1, copy=True)
#
rs.EnableRedraw(True)
