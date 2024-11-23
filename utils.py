import math

def get_temp_ratio_isentropic(machin, gam):
    # /* Utility to get the isentropic temperature ratio given the mach number*/       double number,gm1,mach1s;
    mach1s = machin * machin
    gm1 = gam - 1.0
    number = 1.0 / (1.0 + .5 * gm1 * mach1s)

    return number


def normal_shock_total_pressure_ratio(machin, gam):
    # // NACA 1135 - normal shock relation pt ratio - eq 97

    msq = machin * machin
    gm1 = gam - 1.0
    gp1 = gam + 1.0

    fac2 = (2.0 * gam * msq - gm1) / gp1
    fac1 = (gp1 * msq) / (gm1 * msq + 2.0)
    number = (math.pow(fac1, (gam / gm1))) * (math.pow((1.0 / fac2), (1.0 / gm1)))

    return number


def normal_shock_mach_after(machin, gam):
    # // NACA 1135 - normal shock relation  mach - eq 96... but more complicated?

    msq = machin * machin
    gm1 = gam - 1.0
    gp1 = gam + 1.0

    fac2 = (2.0 * gam * msq - gm1) / gp1
    fac1 = (gp1 * msq) / (gm1 * msq + 2.0)
    number = math.sqrt(msq / (fac2 * fac1))

    return number


def weightflow_per_area_given_mach(mach, gascon, gam):
    # /* Utility to get the corrected weightflow per area given the Mach number */

    fac2 = (gam + 1.0) / (2.0 * (gam - 1.0))
    fac1 = math.pow((1.0 + .5 * (gam - 1.0) * mach * mach), fac2)
    number = 20.78 * math.sqrt(gam / gascon) * mach / fac1

    return number


def calc_gamma(temp, opt):
    # Utility to get gamma as a function of temp
    a = -7.6942651e-13
    b = 1.3764661e-08
    c = -7.8185709e-05
    d = 1.436914
    if opt == 0:
        number = 1.4

    else:
        number = a * temp * temp * temp + b * temp * temp + c * temp + d

    return number

def get_pressure_ratio_isentropic(machin, gam):
    # /* Utility to get the isentropic pressure ratio given the mach number */
    mach1s = machin * machin
    gm1 = gam - 1.0
    fac1 = 1.0 + .5 * gm1 * mach1s
    number = math.pow(1.0 / fac1, gam / gm1)

    return number

def normal_shock_static_pressure_ratio(machin, gam):
    # // NACA 1135 - normal shock relation ps ratio - eq. 93

    msq = machin * machin
    gm1 = gam - 1.0
    gp1 = gam + 1.0

    fac2 = (2.0 * gam * msq - gm1) / gp1
    number = fac2

    return number

def get_mach(sub, corair, gascon, gam):
    # /* Utility to get the Mach number given the corrected airflow per area */
    # /* iterate for mach number */

    a = (gam - 1) / 2.0
    b = -(gam + 1.0) / (2.0 * (gam - 1.0))
    k = 20.78 * math.sqrt(gam / gascon)

    chokair = weightflow_per_area_given_mach(1.0, gascon, gam)
    if corair > chokair:
        number = 1.0
        return number

    else:
        if sub == 1:
            macho = 1.0  # /* sonic */
        else:
            if sub == 2:
                macho = 1.703  # /* supersonic */
            else:
                macho = .25  # /* subsonic */
            airo = weightflow_per_area_given_mach(macho, gascon, gam)  # /* initial guess */
            iter = 1
            machn = macho - .2
            while abs(corair - airo) > .0001 and iter < 50:
                airn = weightflow_per_area_given_mach(machn, gascon, gam)
                deriv = (airn - airo) / (machn - macho)
                airo = airn
                macho = machn
                machn = macho + (corair - airo) / deriv
                iter += 1
        number = macho

    return number