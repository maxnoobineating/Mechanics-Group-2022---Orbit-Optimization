from orb_control import*

# Data definition:
# Optimization Arguments:
# each of:
# mu    0 : Gravitaional constant
# aE    1 : Earth's semi-major axis
# eE    2 : Earth's eccentricity
# aM    3 : Mars' ''
# eM    4 : ''
# periM 5 : Mars' perigee measured from absolute anamoly (Earth)
# Isp   6 : Specific impulse
# fM    7 : maximum Specific thrust (relative to initial mass)
# DV    8 : deltaV budget
#
# args = (mu, aE, eE, aM, eM, periM, Isp, fM, DV)

# Orbit element:
# One of:
# Oi = (ai, ei, f_1i, f_2i), 0 <= i <= n-1
# OE = (aE, eE, None, f2E)  // f2E = ths[0], variable
# OM = (aM, eM, f1M, None)  // f1M = sum(ths), variable
#
# Oi = len(4) tuple

# Orbit elements list:
# list of Oi
#
# os = (n+1) list of (len(4) tuple) // including OE, OM

# Optimization variables:
# vars =
# list of:
# r, 0 <= i < n
# ci, n <= i < 2n-1
# th, 2n, 2n+1
# where vars[2n] = Earth departure anomaly, vars[2n+1] = Mars arrival anomaly
#
# vars = 2*n+1 array

############
# not for efficiency
# Plot
# xy, options -> plot
def makeplot(xy
            , color='k'
            , line=''
            , linewidth=.5
            , marker=''
            , markersize=.1
            , label=None):
    return {'xy': xy
        , 'color': color
        , 'line': line
        , 'linewidth': linewidth
        , 'marker': marker
        , 'markersize': markersize
        , 'label': label}

# plot it
# plot = {'xy': xy, 'color': 'b', 'line': None, 'marker': 'o'}
def xyplot(*plots):
    lgs = []
    lbs = []
    for p in plots:
        lg = plt.plot(p['xy'][0]
                    , p['xy'][1]
                    , linestyle = p['line']
                    , marker = p['marker']
                    , color = p['color']
                    , markersize=p['markersize']
                    , linewidth=p['linewidth'])
        if p['label']!=None:
            lgs.append(lg)
            lbs.append(p['label'])
    hds = tuple((lg[0] for lg in lgs))
    tuplbs = tuple((lb for lb in lbs))
    plt.legend(handles=hds
            , labels=tuplbs)

# plot datas and their decorators
# decs = list of functions of decorators containing plt.plot()
#       , will pass in axs data to execute.
# dec = xys, axs, ex(tra) -> None (plot!)
# xys = list of datas
def plotall(plots=[], decs=[], **exs):
    fig, axs = plt.subplots(1, 1)
    xyplot(*plots)
    for dec in decs:
        dec(plots, axs, exs)
    axs.axis('equal')
    fig.tight_layout()
    plt.title(exs.get('title', ''))
    plt.show()

#
# Plot Decorators
#
factor = [0.8]
red_black5 = ['#a91e25', '#7c1822', '#510a12', '#420a0d', '#0e0e10']
orange_brown5 = ['#c75d05', '#ff760c', '#ff9549', '#f1cb8c', '#915634']
light_darkBrown5 = ['#b58181', '#a36767', '#814d4d', '#683d3d', '#562b2b']
light_darkBlue5 = ['#5766bd', '#4b54a0', '#3b468a', '#374081', '#2e3979']
light_darkCyan5 = ['#4ef1ef', '#2acaea', '#34bdc6', '#02a9b9', '#0095a4']
light_darkGreen = ['#93a47d', '#93b47d', '#93c47d', '#93d47d', '#93e47d']

def sumf(lst, func=lambda x:x):
    sum = 0
    for ele in lst:
        sum += func(ele)
    return sum

def axis_scale(axs):
    return max(sumf(axs.get_xlim(), abs), sumf(axs.get_ylim(), abs))

# Axis
def paxis(xys, axs):
    xlim = list(axs.get_xlim())
    ylim = list(axs.get_ylim())

    xlim[0] = min(-1, xlim[0])
    xlim[1] = max(1, xlim[1])
    ylim[0] = min(-1, ylim[0])
    ylim[1] = max(1, ylim[1])

    x0 = mat([xlim[0], 0])*factor[0]
    x1 = mat([(xlim[1]-xlim[0]), 0])*factor[0]
    y0 = mat([0, ylim[0]])*factor[0]
    y1 = mat([0, (ylim[1]-ylim[0])])*factor[0]

    plt.quiver(*[x0, y0]
            , x1, y1
            , color='k'
            , width=0.003
            , angles='xy'
            , scale_units='xy'
            , scale=1)

# vector dot to dot:
def pvec_dtd(dot1, dot2):
    dot1 = mat(dot1)
    dot2 = mat(dot2)
    plt.quiver(*dot1
            , *(dot2-dot1)
            , color='k'
            , width=0.003
            , angles='xy'
            , scale_units='xy'
            , scale=1)

# draw velocity vector
def pthstvec(plots, axs, exs):
    vars = exs['vars']
    args = exs['args']
    xys = to_xys(vars)
    vels = velvecs(vars, args)
    maxvel = max(map(lambda v: norm(v), vels))
    maxaxis = axis_scale(axs)
    print(vels)
    print(maxaxis)
    print(maxvel)
    plt.quiver(*xys
            , *(-vels.T)
            , color='k'
            , width=0.001
            , angles='xy'
            , scale_units='xy'
            , scale=maxvel/maxaxis*20)

def dot(arr1, arr2):
    arr1 = mat(arr1)
    arr2 = mat(arr2)
    assert arr1.ndim == 1 and arr2.ndim == 1, "fun.dot: input not vectors"
    assert arr1.size == arr2.size, "fun.dot: vec dimension mismatch"
    return sum(arr1*arr2)

def norm(arr):
    arr = mat(arr)
    assert arr.ndim == 1, "fun.norm: input not vector"
    return sqrt(dot(arr, arr))

cross = lambda a, b: a[0]*b[1]-a[1]*b[0]

# counter-clockwise
def rotate(v, th):
    R = np.vstack((
        [cos(th), -sin(th)],
        [sin(th), cos(th)]))
    return sum((R*mat(v)).T)

def phaseAngle(arr):
    arr = mat(arr)
    angle = acos(arr[0]/norm(arr))
    return angle if arr[1]>=0 else -angle

def polarToCart(arr_pol):
    arr_pol = mat(arr_pol)
    arr_cart = []
    for pol in arr_pol.T:
        arr_cart.append([pol[0]*cos(pol[1]), pol[0]*sin(pol[1])])
    return (mat(arr_cart)).T

# generate set of ellipse point in xy
# elps (section) -> xys
def ellipse_2d(elps, fidel=128):
    peri = elps['peri']
    e = elps['e']
    assert e >= 0 and e < 1, "fun.ellipse_2d: invalid eccentricity"
    f1 = elps['fs'][0]
    f2 = elps['fs'][1]
    r_ths = []
    dth = (f2-f1)/fidel
    for n in range(fidel+1):
        r_ths.append((ellipseEq(fth=n*dth+f1, elps=elps), n*dth+f1+peri))
    return polarToCart(mat(r_ths).T)

# fth, elps -> r
# fth:individual anamoly
def ellipseEq(fth, elps):
    e = elps['e']
    assert e >= 0 and e < 1, "fun.ellipse_2d: invalid eccentricity"
    a = elps['a']
    return a*(1-e**2)/(1+e*cos(fth))

# vars to r-theta
def to_rths(vars):
    n = (len(vars)-1)//2
    rs = vars[:n]   # n
    thE = vars[-2]
    thM = vars[-1]  # n
    dth = (thM-thE)/(n-1)
    return np.vstack((rs, [thE+dth*i for i in range(n)]))

def to_xys(vars):
    return polarToCart(to_rths(vars))

# Orbit arc
# vars -> list of {'a': ai, 'e':ei, 'peri': perigeeAngle, 'fs': (f1i, f2i)}
def orbArcs(vars, args):
    n = (len(vars)-1)//2
    os = to_EosM(vars, args, n)[1:n] # n-1
    # os = list of (ai, ei, f1i, f2i)
    thE = vars[-2]
    thM = vars[-1]
    dth = (thM-thE)/(n-1)
    arcs = []
    for i, o, c in zip(range(n-1), os, vars[n:2*n-1]):
        absperi = thE+dth*i-o[2]
        arcs.append({'a': o[0]
                    ,'e': o[1]
                    , 'peri': absperi
                    , 'fs': (o[2], o[3])})
#         print("c: {0:.3f}\
# , a: {1:.3f}\
# , e: {2:.3f}\
# , peri: {3:.3f}\
# , fs: ({4:.3f}, {5:.3f})\
# , th: {6:.3f}".format(c
#             , o[0]
#             , o[1]
#             , absperi
#             , o[2]
#             , o[3]
#             , vars[-2]+dth*i))
    return arcs


# Orbit arc to points
# list of {'a': ai, 'e': ei, 'peri': perigeeAngle, 'fs': (f1i, f2i)} -> xy
def arcsampling(arcs):
    xy = np.empty(shape=(0, 2))
    for arc in arcs:
        xy = np.append(xy, ellipse_2d(arc, fidel=128).T, axis=0)
    return xy.T

# Smooth Orbit sampling wrapper
# vars -> xy
def smoothxy(vars, args):
    return arcsampling(orbArcs(vars, args))


# interpolating an orbit arc (mainly for ci) #
# concatenate two consecutive vars section
# vars1 of n, vars2 of m -> vars of n+m-1
def con_vars(vars1, vars2):
    v1 = list(vars1)
    v2 = list(vars2)
    n1 = (len(v1)-1)//2
    n2 = (len(v2)-1)//2
    assert v1[n1-1] == v2[0], "fun.varscon: erroneous vars sections concatenation"
    return mat(v1[:n1]+v2[1:n2]+v1[n1:-2]+v2[n2:-2]+[None, None], np.float64)

# interp helper (sectioned vars: ri, rj, dth, ci), n -> vars of n
# !!! duped from EosM
def arc_interp(ri, rj, ci, dth, nfld=1):
    if nfld == 0:
        return mat([ri, rj, ci, None, None], np.float64)
    ki = sqrt(rj**2+ri**2-2*rj*ri*cos(dth))
    xi = (ri**2-rj**2)/2/ki
    yi = rj*ri*sin(dth)/ki
    ui = (rj-ri)/2*sqrt(1+2*ci**2/rj/ri/(1-cos(dth)))
    ai = (sqrt(ci**2+(ui+ki/2)**2)+sqrt(ci**2+(ui-ki/2)**2)+ri+rj)/4
    ei = sqrt((xi-ui)**2+(yi-ci)**2)/2/ai

    F1 = mat([xi, yi])
    F2 = mat([ui, ci])
    P1 = mat([-ki/2, 0])
    P2 = mat([ki/2, 0])
    r1i = P1 - F1
    r1j = P2 - F1

    F1F2 = F1 - F2
    f1i = acos(dot(F1F2, r1i)/norm(F1F2)/norm(r1i))
    if cross(F1F2, r1i)<0:
        f1i = 2*pi - f1i
    f2i = f1i + dth # f2i can exceed 2pi
    r_itp = ellipseEq(f1i+dth/2, {'a': ai, 'e': ei})
    ritp = rotate(r1i, dth/2)/norm(r1i)*r_itp

    F2P1 = F2 - P1
    F2P2 = F2 - P2
    pc1 = ritp - r1i
    pc2 = r1j - ritp

    ci1 = cross(pc1, F2P1)/norm(pc1)
    ci2 = cross(pc2, F2P2)/norm(pc2)

    vars1 = arc_interp(ri, r_itp, ci1, dth/2, nfld-1)
    vars2 = arc_interp(r_itp, rj, ci2, dth/2, nfld-1)
    return con_vars(vars1, vars2)

# vars -> list of (ri, rj, ci, dth)
def decomp_vars(vars):
    n = (len(vars)-1)//2
    thE = vars[-2]
    thM = vars[-1]
    dth = (thM-thE)/(n-1)
    rs = iter(vars[:n])
    cs = iter(vars[n:2*n-1])
    ri = next(rs)
    var_lst = []
    for i in range(n-1):
        rj = next(rs)
        ci = next(cs)
        var_lst.append((ri, rj, ci, dth))
        ri = rj
    return var_lst

# vars, nfld -> vars of n*(2**nfld)
def orbit_interp(vars, nfld):
    var_lst = decomp_vars(vars)
    vs = iter(var_lst)
    vars_itped = arc_interp(*next(vs), nfld=nfld)
    for var in vs:
        vars_itped = con_vars(vars_itped, arc_interp(*var, nfld=nfld))
    vars_itped[-2] = vars[-2]
    vars_itped[-1] = vars[-1]
    return vars_itped

# vars -> list of vec(x, y)
def velvecs(vars, args):
    fid = (len(vars)-1)//2
    os = to_EosM(vars, args, fid)
    urs = mat(list(map(lambda v: v/norm(v), to_xys(vars).T)))
    unrs = mat(list(map(lambda v: rotate(v, pi/2), urs)))
    vecs = []
    mu = args[0]
    for i in range(fid):
        ai = os[i][0]
        ei = os[i][1]
        f2i = os[i][3]
        hi = sqrt(mu*ai*(1-ei**2))
        aj = os[i+1][0]
        ej = os[i+1][1]
        f1j = os[i+1][2]
        hj = sqrt(mu*aj*(1-ej**2))

        vec = mu*((ej*sin(f1j)/hj-ei*sin(f2i)/hi)*urs[i]\
                +((1+ej*cos(f1j))/hj-(1+ei*cos(f2i))/hi)*unrs[i])
        vecs.append(vec)
    return mat(vecs)
###################

#---------------------------------------#
# !!! Resource Heavy, optimize these !!!

# from Optimization Variable to Orbit Elements:
@njit(cache=True)
def to_EosM(vars, args, fid):
    rs = vars[0:fid] # n
    cs = vars[fid:2*fid-1] # n-1
    thE = vars[-2]
    thM = vars[-1]
    dth = (thM-thE)/(fid-1)
    os = [(args[1], args[2], 0.0, thE)]
    osappend = os.append
    ri = rs[0]

    for i in range(fid-1):
        rj = rs[i+1]
        ci = cs[i]
        ki = sqrt(rj**2+ri**2-2*rj*ri*cos(dth))
        # =>
        xi = (ri**2-rj**2)/2/ki
        yi = rj*ri*sin(dth)/ki
        ui = (rj-ri)/2*sqrt(1+2*ci**2/rj/ri/(1-cos(dth)))
        # =>
        ai = (sqrt(ci**2+(ui+ki/2)**2)+sqrt(ci**2+(ui-ki/2)**2)+ri+rj)/4
        ei = sqrt((xi-ui)**2+(yi-ci)**2)/2/ai

        r_x = -xi-ki/2
        r_y = -yi
        axis_x = xi-ui
        axis_y = yi-ci
        f1i = acos((axis_x*r_x+axis_y*r_y)\
            /sqrt(axis_x**2+axis_y**2)\
            /sqrt(r_x**2+r_y**2))
        if axis_x*r_y-axis_y*r_x < 0:
            f1i = 2*pi - f1i
        f2i = f1i + dth # f2i can exceed 2pi
        osappend((ai, ei, f1i, f2i))
        ri = rj

    osappend((args[3], args[4], thM-args[5], 0.0))
    # = (aM, eM, f1M, None)
    return os


# os, args -> float
@njit(parallel=True, cache=True)
def sumdV(os, args, fid):
    mu = args[0]
    ai = os[0][0]
    ei = os[0][1]
    f2i = os[0][3]
    sumdv = 0.0
    ios = iter(os)
    next(ios) # skip first
    for O in ios:
        aj = O[0]
        ej = O[1]
        f1j = O[2]
        li = sqrt(ai*(1-ei**2))
        lj = sqrt(aj*(1-ej**2))
        sumdv += sqrt(mu)*sqrt((ej*sin(f1j)/lj-ei*sin(f2i)/li)**2 \
                                + ((1+ej*cos(f1j))/lj-(1+ei*cos(f2i))/li)**2)
        f2i = O[3]
        ai = aj
        ei = ej
    return sumdv

# os, args -> float list
@njit(cache=True)
def dVs(os, args, fid):
    mu = args[0]
    ai = os[0][0]
    ei = os[0][1]
    f2i = os[0][3]
    dvs = []
    dvsappend = dvs.append
    ios = iter(os)
    next(ios) # skip first
    for O in ios:
        aj = O[0]
        ej = O[1]
        f1j = O[2]
        li = sqrt(ai*(1-ei**2))
        lj = sqrt(aj*(1-ej**2))
        dvsappend(sqrt(mu)*sqrt((ej*sin(f1j)/lj-ei*sin(f2i)/li)**2 \
                                + ((1+ej*cos(f1j))/lj-(1+ei*cos(f2i))/li)**2))
        f2i = O[3]
        ai = aj
        ei = ej
    return dvs

# DeltaT

# os, args -> float
@njit(parallel=True, cache=True)
def sumdT(os, args, fid):
    mu = args[0]
    sumdt = 0.0
    for O in os[1:fid]:
        ai = O[0]
        ei = O[1]
        f1i = O[2]
        f2i = O[3]
        E1i = 2*pi*(f1i//(2*pi))+acos((ei+cos(f1i))/(1+ei*cos(f1i)))\
            if f1i%(2*pi)<=pi\
            else 2*pi*((f1i//(2*pi))+1)-acos((ei+cos(f1i))/(1+ei*cos(f1i)))
        E2i = 2*pi*(f2i//(2*pi))+acos((ei+cos(f2i))/(1+ei*cos(f2i)))\
            if f2i%(2*pi)<=pi\
            else 2*pi*((f2i//(2*pi))+1)-acos((ei+cos(f2i))/(1+ei*cos(f2i)))
        sumdt += abs(sqrt(ai**3/mu)*(E2i-ei*sin(E2i))-sqrt(ai**3/mu)*(E1i-ei*sin(E1i)))
    return sumdt

# os, args -> float list
@njit(cache=True)
def dTs(os, args, fid):

    mu = args[0]
    dts = []
    dtsappend = dts.append
    for O in os[1:fid]:
        ai = O[0]
        ei = O[1]
        f1i = O[2]
        f2i = O[3]
        E1i = 2*pi*(f1i//(2*pi))+acos((ei+cos(f1i))/(1+ei*cos(f1i)))\
            if f1i%(2*pi)<=pi\
            else 2*pi*((f1i//(2*pi))+1)-acos((ei+cos(f1i))/(1+ei*cos(f1i)))
        E2i = 2*pi*(f2i//(2*pi))+acos((ei+cos(f2i))/(1+ei*cos(f2i)))\
            if f2i%(2*pi)<=pi\
            else 2*pi*((f2i//(2*pi))+1)-acos((ei+cos(f2i))/(1+ei*cos(f2i)))
        dtsappend(abs(sqrt(ai**3/mu)*(E2i-ei*sin(E2i))-sqrt(ai**3/mu)*(E1i-ei*sin(E1i))))
    return dts

#----------------------------#

# Minimum Energy (sum deltaV)
# @njit(parallel=True)
TX = lambda x: "0.0" if x==0\
                    else "{0:.3f}e{1}".format(x/10**(log(x, 10)//1), int(log(x, 10)//1))

# print("r: {0}, {1} |".format(TX(vars[0]), TX(vars[1])), end='')
# print("c: {0} |".format(("-" if vars[2]<0 else "")+TX(abs(vars[2]))), end='')
# print("thE: {0:.4f}, thM: {1:.4f}\r".format(vars[-2], vars[-1]), end='')
def objA(vars, args, fid):
    return sumdV(to_EosM(vars, args, fid), args, fid)

# (Option B) Minimum Time (sum dT)
@njit(parallel=True, cache=True)
def objB(vars, args, fid):
    # print("# Obj B")
    return sumdT(to_EosM(vars, args, fid), args, fid)

# R = njit(lambda th, e, a: a*(1-e**2)/(1+e*cos(th)))

# OE boundary constraint
@njit(parallel=True, cache=True)
def consOE(vars, args, fid):
    # print("# consOEM")
    # thM = sum(islice(vars, fid, 2*fid)) - args[5] # phase_sum-periM
    # rM = vars[fid-1]
    # aM = args[3]
    # eM = args[4]
    return (args[3]*(1-args[4]**2)/(1+args[4]*cos(vars[-1]-args[5])) - vars[fid-1])/vars[fid-1]

# OM boundary constraint
@njit(parallel=True, cache=True)
def consOM(vars, args, fid):
    # print("# consOEM")
    # rE = vars[0]
    # thE = vars[fid]
    # aE = args[1]
    # eE = args[2]
    return (args[1]*(1-args[2]**2)/(1+args[2]*cos(vars[-2])) - vars[0])/vars[0]


# !!! This constraint is not strict enough!!
# coasting phase will allows some unusually high deltaV node!
# Thrust limit constraint
# DM = njit(lambda dv, Isp: exp(-dv/Isp))
# vars, args -> float list (of n)
@njit(parallel=True, cache=True)
def consThrslim(vars, args, fid):
    # print("# ConsThrust")
    os = to_EosM(vars, args, fid)
    dvs = iter(dVs(os, args, fid))
    dts = dTs(os, args, fid)
    Isp = args[6]
    dv = next(dvs)
    multi = exp(-dv/Isp)
    fs = [Ratio*dv*multi/dts[0]] # exp(-dv/Isp)/dt
    fsappend = fs.append

    for dt in dts:
        dv = next(dvs)
        multi *= exp(-dv/Isp)
        fsappend(Ratio*dv*multi/dt)
    return fs


# (Option B) DeltaV budget constraint
@njit(parallel=True, cache=True)
def consdeltaV(vars, args, fid):
    # print("# ConsDeltaV")
    return (args[8]-sumdV(to_EosM(vars, args, fid), args, fid))*Ratio

# # Callback ( xk[, opt_obj] -> bool, opt_obj for 'trust_constr')
# def callbackme(xk):
#     print("callback!")
#     global count
#     count+=1
#     if count%20 == 0:
#         for x in xk:
#             print(x)
#     return false


#
# elpsM = Earth_Mars.elpsM
# elpsE = Earth_Mars.elpsE
# fid = Earth_Mars.fid
# verb = Earth_Mars.verb
# maxiter = Earth_Mars.maxiter
# xtol = Earth_Mars.xtol
# gtol = Earth_Mars.gtol
# periE = Earth_Mars.periE
# periM = Earth_Mars.periM
#
#
# th0E = 0
# th0M = pi
#
# r0E = ellipseEq(elps=elpsE, fth=th0E-periE)
# r0M = ellipseEq(elps=elpsM, fth=th0M-periM)
#
# rs0 = np.linspace(r0E, r0M, fid)
# ths0 = np.linspace((pi+periM)/(fid-1), (pi+periM)/(fid-1), fid)
# ths0 = [th0E, th0M]
# cs0 = np.linspace(0.0, 0.0, fid-1)
# vars0 = np.concatenate((rs0, cs0, ths0), axis=0)
#
# args = EM['VASIMR'].args
# pvars0 = makeplot(smoothxy(vars0, args))
# plots = [pvars0
#     , makeplot(ellipse_2d(fidel=128, elps=Earth_Mars.elpsE)\
#                 , line=':'\
#                 , label=None\
#                 , linewidth=.2)\
#     , makeplot(ellipse_2d(fidel=128, elps=Earth_Mars.elpsM)\
#                 , line=':'\
#                 , label=None\
#                 , linewidth=.2)]
# decs = [pthstvec]
#
#
# plotall(plots, decs, vars=vars0, args=EM['VASIMR'].args)
