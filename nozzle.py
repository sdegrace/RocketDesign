import math
from utils import *


class Nozzle:

    def __init__(self, runiv, mweight, altin, lconv, athroat, aconv, arat, azrat, ptin, pconv, psin, ttin, tref, tconv, gamopt, gam0):

        self.gam0 = gam0
        self.gamopt = gamopt
        self.tconv = tconv
        self.tref = tref
        self.ttin = ttin
        self.psin = psin
        self.pconv = pconv
        self.ptin = ptin
        self.azrat = azrat
        self.aconv = aconv
        self.arat = arat
        self.lconv = lconv
        self.athroat = athroat
        self.altin = altin
        self.mweight = mweight
        self.runiv = runiv

    def compute(self):
        g0 = 32.2
        rgas = self.runiv * g0 / self.mweight
        alt = self.altin / self.lconv
        ath = self.athroat / self.aconv
        rthrt = math.sqrt(ath / 3.1415926)
        aexit = self.arat * ath
        rexit = math.sqrt(aexit / 3.1415926)
        azero = self.azrat * ath
        rzero = math.sqrt(azero / 3.1415926)
        pt = self.ptin / self.pconv

        pamb = self.psin / self.pconv
        tt = (self.ttin + self.tref) / self.tconv
        if self.gamopt == 1:
            gamtfac = calc_gamma(tt, self.gamopt)
            gamma = self.gam0 * gamtfac / 1.4

        gm1 = gamma - 1.0
        fac1 = gm1 / gamma

        counter = 0
        machth = 1.0  # assume flow is choked
        aircor = weightflow_per_area_given_mach(1.0, rgas, gamma)
        mexit = get_mach(2, (aircor / self.arat), rgas, gamma)
        psup = get_pressure_ratio_isentropic(mexit, gamma) * pt
        mflow = aircor * ath * (pt / 14.7) / math.sqrt(tt / 518.)

        if pamb <= psup:
            mode = 0
            trat = get_temp_ratio_isentropic(mexit, gamma)
            uex = mexit * math.sqrt(gamma * rgas * trat * tt)
            pexit = psup

        # over expanded nozzle
        if pamb > psup:
            # find exit pressure at which normal shock leaves the nozzle
            mnsup = mexit
            psub = psup * normal_shock_static_pressure_ratio(mnsup, gamma)

            # slightly overexpanded .. no shock in nozzle
            if pamb <= psub:
                mode = 1
                pexit = psup
                trat = get_temp_ratio_isentropic(mexit, gamma)
                uex = mexit * math.sqrt(gamma * rgas * trat * tt)

            # highly overexpanded .. normal shock in nozzle
            if pamb > psub:
                mode = 2
                pexit = pamb
                anso = aexit
                mnsup = mexit
                mnsub = normal_shock_mach_after(mnsup, gamma)
                ptrat = normal_shock_total_pressure_ratio(mnsup, gamma)
                pso = get_pressure_ratio_isentropic(mnsub, gamma) * ptrat * pt
                ansn = anso - 1.
                while (abs(pexit - pso) > .001) and (counter < 20):
                    counter += 1
                    mnsup = get_mach(2, (aircor * ath / ansn), rgas, gamma)
                    mnsub = normal_shock_mach_after(mnsup, gamma)
                    ptrat = normal_shock_total_pressure_ratio(mnsup, gamma)
                    mexit = get_mach(0, (aircor / self.arat / ptrat), rgas, gamma)
                    psn = get_pressure_ratio_isentropic(mexit, gamma) * ptrat * pt
                    deriv = (psn - pso) / (ansn - anso)
                    pso = psn
                    anso = ansn
                    ansn = anso + (pexit - pso) / deriv

                ans = anso
                rns = math.sqrt(ans / 3.1415926)
                trat = get_temp_ratio_isentropic(mexit, gamma)
                uex = mexit * math.sqrt(gamma * rgas * trat * tt)
        if pamb > .0001:
            npr = pt / pamb
        else:
            npr = 1000.

        fgros = mflow * uex / g0 + (pexit - pamb) * aexit

        return fgros


# lunits = 0
# aconv = 1.  #/ *area    sq    inches * /
# lconv = 1.  #/ *length    feet * /
# fconv = 1.0  #/ *pounds * /
# pconv = 1.0  #/ *lb / sq in * /
# tref = 459.7  #/ *zero    rankine * /
# tconv = 1.0  #/ *degrees    F * /
# mconv1 = 1.0  #/ *airflow    rate    lbs / sec * /
#
# athmx = 300.
# athmn = .1
# aexmx = 100.
# aexmn = 1.
# azmx = 10.
# azmn = 1.1
# ptmx = 3000.
# ptmn = 1.
# pemx = 15.
# pemn = 0.0
# ttmx = 6500.
# ttmn = 500.
# altmn = 0.0
# altmx = 100000.
# altin = 0.0
# alt = 0.0
# flomode = 0
#
# mwtab = mweight = 16.0
# rgas = 1716. #/ *air - ft2 / sec2 R * /
# molopt = 1
# fuelold = fuelopt = 4
# oxopt = 1
# ttab = tcomb = ttin = 5870.
# gamopt = 1
# gam0 = 1.32
# gamma = 1.22
# tcomopt = 1
# ofopt = 0
# ofrat = 8.0
#
# athroat = 150.0
# arat = 50.0
# azero = 600.0
# azrat = 4.0
# ptin = 2000.
# psin = 14.7
#
# mode = 0
# ans = 0.0
# fact = .2
# xt = 40
# yt = 0
# sldloc = 40
# lngth = .5
# lngmx = (athroat / 144.) * 20.0
# lngmn = (athroat / 144.) * .01
#
# runiv = 1545.
#
# print(Nozzle().compute(runiv, mweight, altin, lconv, athroat, aconv, arat, azrat, ptin, pconv, psin, ttin, tref, tconv, gamopt,
#                 gam0))