from orb_functions import *
# parameter unload...
# dunno why directly insert as argument will not immediately evaluate

def main(str):
    vars = {}
    pr = {}

    elpsM = EM[str].elpsM
    elpsE = EM[str].elpsE
    args = EM[str].args
    fid = EM[str].fid
    verb = EM[str].verb
    maxiter = EM[str].maxiter
    xtol = EM[str].xtol
    gtol = EM[str].gtol
    fM = EM[str].fM
    ubs = EM[str].ubs
    lbs = EM[str].lbs
    periE = EM[str].periE
    periM = EM[str].periM

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
    # vars['init'] = np.concatenate((rs0, cs0, ths0), axis=0)
    # Hohmann transfer inital guess
    vars['init'] = orbit_interp(np.array([147.0997134 , 216.89938903
                                        , 0.        , 0.
                                        , 3.14159265], np.float64)
                                , EM[str].nfld)

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
    #---- Basinhopping ----#
    # Option A
    consA = [con1, con2, nlc]

    min_keywords = {"method": 'trust-constr'
        , "bounds": bnds
        , "constraints": consA
        , "args": args
        , "options": {'xtol' : xtol
                    , 'maxiter': maxiter
                    , 'verbose': verb}}

    with cProfile.Profile() as pr['A']:
        solA = basinhopping(objectiveA
                            , vars['init'].copy()
                            , T = 50
                            , niter = 20
                            , minimizer_kwargs=min_keywords
                            , niter_success=5)
    vars['A'] = solA.x.copy()


    #--------------------------------#
    #---- Basinhopping ----#
    # Option B
    consB = [con1, con2, conB, nlc]

    min_keywords = {"method": 'trust-constr'
        , "bounds": bnds
        , "constraints": consB
        , "args": args
        , "options": {'xtol' : xtol
                    , 'maxiter': maxiter
                    , 'verbose': verb}}

    with cProfile.Profile() as pr['B']:
        solB = basinhopping(objectiveB
                            , vars['A'].copy()
                            , T = 50
                            , niter = 20
                            , minimizer_kwargs=min_keywords
                            , niter_success=5)
    vars['B'] = solB.x.copy()

    return (vars, pr)

vars = {}
pr = {}

vars['RD-0410'], pr['RD-0410'] = main('RD-0410')
vars['VASIMR'], pr['VASIMR'] = main('VASIMR')

print("VASIMR Option A Statistics:")
stats = pstats.Stats(pr['VASIMR']['A'])
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats(.01)
print("VASIMR Option B Statistics:")
stats = pstats.Stats(pr['VASIMR']['B'])
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats(.01)
print("RD-0410 Option A Statistics:")
stats = pstats.Stats(pr['RD-0410']['A'])
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats(.01)
print("RD-0410 Option B Statistics:")
stats = pstats.Stats(pr['RD-0410']['B'])
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats(.01)

osRA = to_EosM(vars['RD-0410']['A'], EM['RD-0410'].args, fid=Earth_Mars.fid)
osRB = to_EosM(vars['RD-0410']['B'], EM['RD-0410'].args, fid=Earth_Mars.fid)
osVA = to_EosM(vars['VASIMR']['A'], EM['VASIMR'].args, fid=Earth_Mars.fid)
osVB = to_EosM(vars['VASIMR']['B'], EM['VASIMR'].args, fid=Earth_Mars.fid)

dvsRA = dVs(osRA, EM['RD-0410'].args, fid=Earth_Mars.fid)
dvsRB = dVs(osRB, EM['RD-0410'].args, fid=Earth_Mars.fid)
dvsVA = dVs(osVA, EM['VASIMR'].args, fid=Earth_Mars.fid)
dvsVB = dVs(osVB, EM['VASIMR'].args, fid=Earth_Mars.fid)

dtsRA = dTs(osRA, EM['RD-0410'].args, fid=Earth_Mars.fid)
dtsRB = dTs(osRB, EM['RD-0410'].args, fid=Earth_Mars.fid)
dtsVA = dTs(osVA, EM['VASIMR'].args, fid=Earth_Mars.fid)
dtsVB = dTs(osVB, EM['VASIMR'].args, fid=Earth_Mars.fid)

payloadRatioRA = exp(-sum(dvsRA)/EM['RD-0410'].Isp)
payloadRatioVA = exp(-sum(dvsVA)/EM['VASIMR'].Isp)

xyplot(makeplot(ellipse_2d(fidel=128, elps=Earth_Mars.elpsE)\
            , line=':'\
            , label=None\
            , linewidth=.2)\
    , makeplot(ellipse_2d(fidel=128, elps=Earth_Mars.elpsM)\
            , line=':'\
            , label=None\
            , linewidth=.2)\
    , makeplot(smoothxy(vars['RD-0410']['A'], EM['RD-0410'].args)\
            , color='y'\
            , label='RD-0410 A'\
            , line=':'
            , linewidth=1.0)\
    , makeplot(smoothxy(vars['RD-0410']['B'], EM['RD-0410'].args)\
            , color='r'\
            , label='RD-0410 B'\
            , line=':'
            , linewidth=1.0)\
    , makeplot(smoothxy(vars['VASIMR']['A'], EM['VASIMR'].args)\
            , color='b'\
            , label='VASIMR A'\
            , line=':'
            , linewidth=1.0)\
    , makeplot(smoothxy(vars['VASIMR']['B'], EM['VASIMR'].args)\
            , color='m'\
            , label='VASIMR B'\
            , line=':'
            , linewidth=1.0))
