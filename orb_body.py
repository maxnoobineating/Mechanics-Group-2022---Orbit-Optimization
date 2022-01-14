from orb_functions import *
# parameter unload...
# dunno why directly insert as argument will not immediately evaluate

elpsM = EM['VSMIR'].elpsM
elpsE = EM['VSMIR'].elpsE
args = EM['VSMIR'].args
fid = EM['VSMIR'].fid
verb = EM['VSMIR'].verb
maxiter = EM['VSMIR'].maxiter
xtol = EM['VSMIR'].xtol
gtol = EM['VSMIR'].gtol
fM = EM['VSMIR'].fM
ubs = EM['VSMIR'].ubs
lbs = EM['VSMIR'].lbs
periE = EM['VSMIR'].periE
periM = EM['VSMIR'].periM

#-------------------------------------------#
# "Initial" initial guess
th0E = 0
th0M = pi

r0E = ellipseEq(elps=elpsE, fth=th0E-periE)
r0M = ellipseEq(elps=elpsM, fth=th0M-periM)

rs0 = np.linspace(r0E, r0M, fid)
# ths0 = np.linspace((pi+periM)/(fid-1), (pi+periM)/(fid-1), fid)
ths0 = [th0E, th0M]
cs0 = np.linspace(0.0, 0.0, fid-1)
vars0 = np.concatenate((rs0, cs0, ths0), axis=0)
vars0 = orbit_interp(np.array([147.0997134 , 216.89938903
                                , 0.        , 0.
                                , 3.14159265], np.float64), EM['VSMIR'].nfld)
os0 = to_EosM(vars0, args, fid)
dvs0 = dVs(os0, args, fid)
dts0 = dTs(os0, args, fid)

#
# xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
#     , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
#     , makeplot(smoothxy(vars0), color='b'))

# xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
#     , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
#     , makeplot(smoothxy(vars0, args), color='b'))
#-------------------------------------------#



# function repackaging
objectiveA = \
    lambda vars, *args: objA(vars, args, fid=fid)

objectiveB = \
    lambda vars, *args: objB(vars, args, fid=fid)

consOMbnd_eq = \
    lambda vars, *args: consOM(vars, args, fid=fid)

consOEbnd_eq = \
    lambda vars, *args: consOE(vars, args, fid=fid)

consdeltaV_OptB_ineq = \
    lambda vars, *args: consdeltaV(vars, args, fid=fid)

# constraint packaging

con1 = {'type': 'eq', 'fun': consOEbnd_eq, 'args': args}
con2 = {'type': 'eq', 'fun': consOMbnd_eq, 'args': args}
conB = {'type': 'ineq', 'fun': consdeltaV_OptB_ineq, 'args': args}

bnds = opt.Bounds(lb=lbs, ub=ubs, keep_feasible=False)

nlc_lb = np.linspace(0, 0, fid)
nlc_ub = np.linspace(Ratio*fM, Ratio*fM, fid)

nlc = opt.NonlinearConstraint(\
    njit(lambda vars: consThrslim(vars, args, fid))\
    , nlc_lb\
    , nlc_ub)

# Dictionary constraint 和 NonlinearConstraint Object 可以混用!
# Option A
cons1 = [con1, con2, nlc]

with cProfile.Profile() as pr1:
    sol1 = minimize(objectiveA
        , vars0.copy()
        , method='trust-constr'
        , bounds=bnds
        , constraints=cons1
        , args=args
        , options={'xtol': xtol
                , 'maxiter': maxiter
                , 'verbose': verb})

vars1 = sol1.x.copy()
rths1 = rths(vars1)
xys1 = xys(vars1)
smxys1 = smoothxy(vars1, args)
os1 = to_EosM(vars1, args, fid=fid)
dvs1 = dVs(os1, args, fid=fid)
dts1 = dTs(os1, args, fid=fid)

xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
    , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
    , makeplot(smxys1, color='b'))

xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
    , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
    , makeplot(xys1, color='b'))

stats = pstats.Stats(pr1)
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats(.01)



#---- Basinhopping ----#
cons1 = [con1, con2, nlc]

min_keywords = {"method": 'trust-constr'
    , "bounds": bnds
    , "constraints": cons1
    , "args": args
    , "options": {'xtol' : xtol
                , 'maxiter': maxiter
                , 'verbose': verb}}
with cProfile.Profile() as pr1b:
    sol1b = basinhopping(objectiveA
                        , vars0.copy()
                        , T = 50
                        , niter = 20
                        , minimizer_kwargs=min_keywords
                        , niter_success=5)
# niter_success
# Stop the run if the global minimum candidate remains the same
# for this number of iterations.

# T: temperature
# The “temperature” parameter for the accept or reject criterion.
# Higher “temperatures” mean that larger jumps in function value
# will be accepted. For best results T should be comparable to
# the separation (in function value) between local minima.

vars1b = sol1b.x.copy()
rths1b = rths(vars1b)
xys1b = xys(vars1b)
smxys1b = smoothxy(vars1b, args)
os1b = to_EosM(vars1b, args, fid=fid)
dvs1b = dVs(os1b, args, fid=fid)
dts1b = dTs(os1b, args, fid=fid)

xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
    , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
    , makeplot(smxys1b, color='b'))

xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
    , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
    , makeplot(xys1b, color='b'))

stats = pstats.Stats(pr1b)
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats(.01)

#--------------------------------#


# Option B
cons2 = [con1, con2, conB, nlc]

with cProfile.Profile() as pr2:
    sol2 = minimize(objectiveB
        , sol1.x.copy()
        , method='trust-constr'
        , bounds=bnds
        , constraints=cons2
        , args=args
        , options={'xtol' : xtol
                , 'maxiter': maxiter
                , 'verbose': verb})

vars2 = sol2.x.copy()
rths2 = rths(vars2)
xys2 = xys(vars2)
smxys2 = smoothxy(vars2, args)
os2 = to_EosM(vars2, args, fid=fid)
dvs2 = dVs(os2, args, fid=fid)
dts2 = dTs(os2, args, fid=fid)

xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
    , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
    , makeplot(smxys2, color='b'))

xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
    , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
    , makeplot(xys2, color='r'))

stats = pstats.Stats(pr2)
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats(.01)

#---- Basinhopping ----#
min_keywords = {"method": 'trust-constr'
    , "bounds": bnds
    , "constraints": cons2
    , "args": args
    , "options": {'xtol' : xtol
                , 'maxiter': maxiter
                , 'verbose': verb}}
with cProfile.Profile() as pr2b:
    sol2b = basinhopping(objectiveB
                        , vars0.copy()
                        , T = 50
                        , niter = 20
                        , minimizer_kwargs=min_keywords
                        , niter_success=5)
# niter_success
# Stop the run if the global minimum candidate remains the same
# for this number of iterations.

# T: temperature
# The “temperature” parameter for the accept or reject criterion.
# Higher “temperatures” mean that larger jumps in function value
# will be accepted. For best results T should be comparable to
# the separation (in function value) between local minima.

vars2b = sol2b.x.copy()
rths2b = rths(vars2b)
xys2b = xys(vars2b)
smxys2b = smoothxy(vars2b, args)
os2b = to_EosM(vars2b, args, fid=fid)
dvs2b = dVs(os2b, args, fid=fid)
dts2b = dTs(os2b, args, fid=fid)

xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
    , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
    , makeplot(smxys2b, color='b'))

xyplot(makeplot(ellipse_2d(fidel=128, elps=elpsE), line='--')
    , makeplot(ellipse_2d(fidel=128, elps=elpsM), line='--')
    , makeplot(xys2b, color='b'))

stats = pstats.Stats(pr2b)
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats(.01)
