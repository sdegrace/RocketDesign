import math
from enum import Enum
from utils import *

math.pi


# ShockMode = Enum("ShockMode", "under_expanded over_expanded shock_in_nozzle")
class ShockMode(Enum):
    under_expanded = 0
    slight_overexpanded = 1
    shock_in_nozzle = 2
    # shock_in_nozzle = 3


class Nozzle:

    def __init__(self, molecular_weight_exhaust, area_throat, area_ratio_exit, area_ratio_plenum,
                 total_pressure_chamber, temperature_total, gamopt,
                 specific_heat_ratio_std):

        self.thrust_gross = None
        self.specific_heat_ratio_std = specific_heat_ratio_std
        self.gamopt = gamopt
        self.temperature_total = temperature_total
        self.total_pressure_chamber = total_pressure_chamber
        self.area_ratio_plenum = area_ratio_plenum
        self.area_ratio_exit = area_ratio_exit
        self.area_throat = area_throat
        self.molecular_weight_exhaust = molecular_weight_exhaust

        self.rgas = R_u * g_0 / self.molecular_weight_exhaust
        self.radius_throat = math.sqrt(self.area_throat / math.pi)
        self.area_exit = self.area_ratio_exit * self.area_throat
        self.radius_exit = math.sqrt(self.area_exit / math.pi)
        self.area_plenum = self.area_ratio_plenum * self.area_throat
        self.radius_plenum = math.sqrt(self.area_plenum / math.pi)
        # self.pt = self.total_pressure_chamber / kP_to_psi

        # self.tt = (self.temperature_total + temperature_reference) / kelvin_per_rankine

    def compute(self, pressure_free_stream):

        # pamb = pressure_free_stream / kP_to_psi

        if self.gamopt == 1:
            gamtfac = calc_specific_heat_ratio(self.temperature_total, self.gamopt)
            specific_heat_ratio = self.specific_heat_ratio_std * gamtfac / 1.4

        gm1 = specific_heat_ratio - 1.0
        fac1 = gm1 / specific_heat_ratio

        counter = 0
        machth = 1.0  # assume flow is choked
        aircor = weightflow_per_area_given_mach(1.0, self.rgas, specific_heat_ratio)
        print(aircor)
        mach_exit = get_mach(2, (aircor / self.area_ratio_exit), self.rgas, specific_heat_ratio)
        psup = get_pressure_ratio_isentropic(mach_exit, specific_heat_ratio) * self.total_pressure_chamber
        mflow = aircor * self.area_throat * (self.total_pressure_chamber / p_0) / math.sqrt(
            self.temperature_total / T_0)

        if pressure_free_stream <= psup:
            mode = ShockMode.under_expanded
            temp_ratio = get_temp_ratio_isentropic(mach_exit, specific_heat_ratio)
            uex = mach_exit * math.sqrt(specific_heat_ratio * self.rgas * temp_ratio * self.temperature_total)
            pressure_exit = psup

        # over expanded nozzle
        elif pressure_free_stream > psup:
            # find exit pressure at which normal shock leaves the nozzle
            mach_number_supersonic = mach_exit
            psub = psup * normal_shock_static_pressure_ratio(mach_number_supersonic, specific_heat_ratio)

            # slightly overexpanded - no shock in nozzle
            if pressure_free_stream <= psub:
                mode = ShockMode.slight_overexpanded
                pressure_exit = psup
                temp_ratio = get_temp_ratio_isentropic(mach_exit, specific_heat_ratio)
                uex = mach_exit * math.sqrt(specific_heat_ratio * self.rgas * temp_ratio * self.temperature_total)

            # highly overexpanded - normal shock in nozzle
            if pressure_free_stream > psub:
                mode = ShockMode.shock_in_nozzle
                pressure_exit = pressure_free_stream
                anso = self.area_exit
                mach_number_supersonic = mach_exit
                mach_number_subsonic = normal_shock_mach_after(mach_number_supersonic, specific_heat_ratio)
                total_pressure_ratio = normal_shock_total_pressure_ratio(mach_number_supersonic, specific_heat_ratio)
                pso = get_pressure_ratio_isentropic(mach_number_subsonic, specific_heat_ratio) * total_pressure_ratio * self.total_pressure_chamber
                ansn = anso - 1.
                while (abs(pressure_exit - pso) > .001) and (counter < 20):
                    counter += 1
                    mach_number_supersonic = get_mach(2, (aircor * self.area_throat / ansn), self.rgas, specific_heat_ratio)
                    mach_number_subsonic = normal_shock_mach_after(mach_number_supersonic, specific_heat_ratio)
                    total_pressure_ratio = normal_shock_total_pressure_ratio(mach_number_supersonic, specific_heat_ratio)
                    mach_exit = get_mach(0, (aircor / self.area_ratio_exit / total_pressure_ratio), self.rgas, specific_heat_ratio)
                    psn = get_pressure_ratio_isentropic(mach_exit,
                                                        specific_heat_ratio) * total_pressure_ratio * self.total_pressure_chamber
                    deriv = (psn - pso) / (ansn - anso)
                    pso = psn
                    anso = ansn
                    ansn = anso + (pressure_exit - pso) / deriv

                ans = anso
                rns = math.sqrt(ans / math.pi)
                temp_ratio = get_temp_ratio_isentropic(mach_exit, specific_heat_ratio)
                uex = mach_exit * math.sqrt(specific_heat_ratio * self.rgas * temp_ratio * self.temperature_total)
        elif pressure_free_stream > .0001:
            self.nozzle_pressure_ratio = self.total_pressure_chamber / pressure_free_stream
        else:
            self.nozzle_pressure_ratio = 1000.

        self.thrust_gross = mflow * uex / g_0 + (pressure_exit - pressure_free_stream) * self.area_exit

        return self.thrust_gross


lunits = 0
sq_cm_to_sq_in = 1.  # / *area    sq    inches * /
# feet_to_meters = 1.  # / *length    feet * /
fconv = 1.0  # / *pounds * /
kP_to_psi = 1.0  # / *lb / sq in * /
# temperature_reference = 459.7  # / *zero    rankine * /
kelvin_per_rankine = 1.0  # / *degrees    F * /
mconv1 = 1.0  # / *airflow    rate    lbs / sec * /

athmx = 300.
athmn = .1
aexmx = 100.
aexmn = 1.
azmx = 10.
azmn = 1.1
ptmx = 3000.
ptmn = 1.
pemx = 15.
pemn = 0.0
ttmx = 6500.
ttmn = 500.
alt_min = 0.0
altmx = 100000.
altitude = 0.0
alt = 0.0
flomode = 0

mwtab = molecular_weight_exhaust = 16.0
# rgas = 1716.  # / *air - ft2 / sec2 R * /
molopt = 1
fuelold = fuelopt = 4
oxopt = 1
ttab = tcomb = temperature_total = 5870. + T_0
gamopt = 1
specific_heat_ratio_std = 1.32
specific_heat_ratio = 1.22
tcomopt = 1
ofopt = 0
ofrat = 8.0

area_throat = 150.0
area_ratio_exit = 50.0
azero = 600.0
area_ratio_plenum = 4.0
total_pressure_chamber = 2000.
pressure_free_stream = 14.7

mode = 0
ans = 0.0
fact = .2
xt = 40
yt = 0
sldloc = 40
lngth = .5
lngmx = (area_throat / 144.) * 20.0
lngmn = (area_throat / 144.) * .01

aconv = 6.4516;      #               /* sq cm */
lconv = .3048;       #              /* meters */
fconv = 4.448 ;      #              /* newtons */
pconv = 6.891 ;      #        /* kilo-pascals */
tref = 273.1 ;       #        /* zero kelvin */
tconv = 0.555555 ;   #         /* degrees C */
mconv1 = .4536 ;     #          /* kg/sec */

# athroat = aths * aconv ;
# athmx = athms * aconv ;
# ptin = pts * pconv ;
# ttin = tts * tconv ;
# ttmx = ttmxs * tconv ;
# tcomb = ttcos * tconv ;
# ttab = ttabs * tconv ;
# psin = pss * pconv ;
# ptmx = ptmxs * pconv ;
# altin = altns * lconv ;
# altmx = altmxs * lconv ;
# athmn = .01 * aconv ;
# ptmn = 1. * pconv ;
# pemx = 15. * pconv ;
# pemn = 0.0 * pconv ;
# ttmn = 500. * tconv ;

molecular_weight_exhaust =  16
area_throat =               967.74/10000
area_ratio_exit =           49.906
area_ratio_plenum =         3.9925
total_pressure_chamber =    5690.07*1000#13782.0*1000
temperature_total =         3261.108 + T_0
# gamopt =                    1
specific_heat_ratio_std =   1.32
pressure_free_stream =      101.297*1000

print('F', Nozzle(molecular_weight_exhaust, area_throat, area_ratio_exit, area_ratio_plenum, total_pressure_chamber,
             temperature_total, gamopt, specific_heat_ratio_std).compute(pressure_free_stream))

