


"""
p_t = Total Pressure
T_t = Total Temperature
P_0 = Free Stream Pressure
gamma = Specific Heat Ratio
R = Gas Constant
A = Area
 - A* = Combustion Chamber throat area
 - A_e = Nozzle exit area
r = Density
V = Velocity
M = Mach
a = speed of sound

Mass flow Rate = m_dot = r * V * A

m_dot = Mass Flow Rate = (A* * p_t)/sqrt(T_t)*sqrt(gamma/R)*((gamma+1)/2)**-((gamma+1)/(2*(gamma-1))
Exit Mach = M_e = A_e/A* = ((gamma + 1)/2)**-((gamma + 1)/(2*(gamma-1))*(1+((gamma-1)/2)*M_e**2)**((gamma+1)/(2*(gamma-1)))/M_e
"""

from math import sqrt, tan, atan

g0 = 9.81


# gamma = 1
# M_e = 1
#
# ((gamma + 1)/2)**-((gamma + 1)/(2*(gamma-1)))*(1+((gamma-1)/2)*M_e**2)**((gamma+1)/(2*(gamma-1)))/M_e

# class EquationSet:
#
#     def __init__(self, equations):
#         self.equations = equations
#
#     def solve(self):


class Isentropic:

    @staticmethod
    def speed_of_sound(specific_heat_ratio, specific_gas_constant, temperature):
        gam, R, T = specific_heat_ratio, specific_gas_constant, temperature
        return sqrt(gam * R * T)

    @staticmethod
    def velcoity_from_mach(speed_of_sound, mach_number):
        a, M = speed_of_sound, mach_number
        return a * M

    @staticmethod
    def mach_number(velocity, speed_of_sound):
        return velocity / speed_of_sound

    @staticmethod
    def speed_of_sound_given_pressure(specific_heat_ratio, pressure, density):
        gam, p, rho = specific_heat_ratio, pressure, density
        return sqrt(gam * p / rho)

    @staticmethod
    def dynamic_pressure_given_density(density, velocity):
        rho, v = density, velocity
        return (1 / 2) * rho * v ** 2

    @staticmethod
    def dynamic_pressure_given_pressure(specific_heat_ratio, pressure, mach):
        gam, p, M = specific_heat_ratio, pressure, mach
        return (gam / 2) * p * M ** 2

    @staticmethod
    def total_pressure_ratio(specific_heat_ratio, mach):
        gam, M = specific_heat_ratio, mach
        return (1 + ((gam - 1) / 2) * M ** 2) ** (-gam / (gam - 1))

    @staticmethod
    def total_temperature_ratio(specific_heat_ratio, mach):
        gam, M = specific_heat_ratio, mach
        return 1 / (1 + ((gam - 1) / 2) * M ** 2)

    @staticmethod
    def total_density_ratio(specific_heat_ratio, mach):
        gam, M = specific_heat_ratio, mach
        return (1 + ((gam - 1) / 2) * M ** 2) ** (-1 / (gam - 1))

    @staticmethod
    def area_ratio(specific_heat_ratio, mach):
        gam, M = specific_heat_ratio, mach
        exp = -(gam + 1) / (2 * (gam - 1))
        p1 = ((gam + 1) / 2) ** exp
        p2 = ((1 + ((gam - 1) / 2) * M ** 2) ** exp) / M
        return p1 * p2

    # @staticmethod
    # def prandtl_meyer(specific_heat_ratio, mach):


################
# Fundamentals #
################

def h(internal_energy, pressure, volume):
    u, p, v = internal_energy, pressure, volume
    return u + p * v


def enthalpy(enthalpy_of_formation_std, sensible_enthalpy_ref, sensible_enthalpy_std):
    h_bar_fo, h_bar, h_bar_o = enthalpy_of_formation_std, sensible_enthalpy_ref, sensible_enthalpy_std
    return h_bar_fo + (h_bar - h_bar_o)


def specific_heat_ratio(specific_heat_constant_pressre, specific_heat_constant_volume):
    c_p, c_v = specific_heat_constant_pressre, specific_heat_constant_volume
    return c_p / c_v


def gas_constant(specific_heat_constant_pressre, specific_heat_constant_volume):
    c_p, c_v = specific_heat_constant_pressre, specific_heat_constant_volume
    return c_p - c_v


# TODO: Eq 5

def equation_of_state(pressure=None, volume=None, moles=None, specific_gas_constant=None, temperature=None):
    p, v, n, R, T = pressure, volume, moles, specific_gas_constant, temperature
    degrees_of_freedom = sum(parameter is None for parameter in [p, n, R, T])
    assert degrees_of_freedom <= 1, "Cannot solve for multiple parameters"
    if degrees_of_freedom == 0:
        return p * v == n * R * T
    elif p is None:
        return (n * R * T) / v
    elif v is None:
        return (n * R * T) / p
    elif n is None:
        return (p * v) / (R * T)
    elif R is None:
        return (p * v) / (n * T)
    elif T is None:
        return (p * v) / (n * R)
    else:
        raise Exception('wut')


def massflow_basic(area, velocity, specific_volume):
    A, v, V = area, velocity, specific_volume
    return A * v / V


##########
# THRUST #
##########

def thrust(massflow, velocity_exit, pressure_exit, pressure_ambient, area_exit):
    m_dot, v_2, p_2, p_3, A_2 = massflow, velocity_exit, pressure_exit, pressure_ambient, area_exit
    return m_dot * v_2 + (p_2 - p_3) * A_2


def v(k, R, T_0, p_2, p_0):
    p1 = k - 1
    return sqrt(2 * k / p1 * R * T_0 * (1 - (p_2 / p_0) ** (p1 / k)))


def c_star(p_1, A_t, m_dot):
    return p_1 * A_t / m_dot


########################
# ISENTROPIC RELATIONS #
########################

def temp_ratio_pressure(p_x, p_y, k):
    return (p_x / p_y) ** ((k - 1) / k)


def temp_ratio_volume(V_y, V_x, k):
    return (V_y / V_x) ** (k - 1)


##################
# MACH RELATIONS #
##################

def a(k, R, T):
    return sqrt(k * R * T)


def M(v, a):
    return v / a


def area_ratio(M_x, M_y, k):
    kplus = k + 1
    kminus = k - 1
    p1 = (kminus / 2)

    return (M_x / M_y) * sqrt(((1 + p1 * M_y ** 2) / (1 + p1 * M_x ** 2)) ** (kplus / kminus))


#########################
# STAGNATION CONDITIONS #
#########################

def T0(T, k, M):
    return T * (1 + (k - 1) / 2 * M ** 2)


def p0(p, k, M):
    kminus = k - 1
    return p * (1 + kminus / 2 * M ** 2) ** (k / kminus)


def rho0(rho, k, M):
    kminus = k - 1
    return rho * (1 + kminus / 2 * M ** 2) ** (1 / kminus)


################################
# CRITICAL / THROAT CONDITIONS #
################################

def Tt(T1, k):
    return 2 * T1 / (k + 1)


def pt(p1, k):
    return p1 * (2 / (k + 1)) ** (k / (k - 1))


def Vt(V1, k):
    return V1 * ((k + 1) / 2) ** (1 / (k - 1))


def m_dot_crit(A_t, p1, k, R, T1):
    kplus = k + 1
    num = k * sqrt((2 / kplus) ** (kplus / (k - 1)))
    den = sqrt(k * R * T1)
    return A_t * p1 * (num / den)


#################
# NOZZLE RATIOS #
#################

def velocity_ratio(k, p_x, p1):
    kminus = k - 1
    return sqrt((k + 1) / kminus * (1 - (p_x / p1) ** (kminus / k)))


def area_ratio(k, p_x, p1):
    kplus = k + 1
    kminus = k - 1
    prat = p_x / p1
    return (kplus / 2) ** (1 / kminus) * prat ** (1 / k) * sqrt(kplus / kminus * (1 - prat ** (kminus / k)))


######################
# CHAMBER PROPERTIES #
######################

def L_star(V_c, A_t):
    return V_c / A_t


def t_s(V_c, m_dot, V_1):
    return V_c / (m_dot * V_1)


'''
The conservation of mass (continuity) tells us that the mass flow rate mdot through a tube is a constant and equal to the product of the density r, velocity V, and flow area A:

Eq #1:

mdot = r * V * A

Considering the mass flow rate equation, it appears that for a given area and a fixed density, we could increase the
mass flow rate indefinitely by simply increasing the velocity. In real fluids, however, the density does not remain 
fixed as the velocity increases because of compressibility effects. We have to account for the change in density to 
determine the mass flow rate at higher velocities. If we start with the mass flow rate equation given above and use the 
isentropic flow relations and the equation of state, we can derive a compressible form of the mass flow rate equation.
'''


def thrust(free_stream_massflow, exit_massflow, exit_velocity, free_stream_velocity, exit_area, exit_pressure,
           free_stream_pressure):
    return ((exit_massflow * exit_velocity) -  # massflow of exit stream
            (free_stream_massflow * free_stream_velocity) +  # massflow of free stream
            (
                    exit_pressure - free_stream_pressure) * exit_area)  # "Suction" force across the nozzle exit due to pressure differential


'''
F=m_dot*V_e+A_e(p_e−p_0)
F = (m dot * V)e – (m dot * V)0 + (pe – p0) * Ae

force is change in momentum over change in time. Engine provides thrust by accelerating a flow stream
flow stream has mass per time, so the mass flow rate already contains the time dependence (mass/time). We can express 
the change in momentum across the propulsion device as the change in the mass flow rate times the velocity.

 If mass doesn't change force
is entire dependent on velocity. So force is proportional to the change in stream velocity.
Rockets have initial stream velocity of zero, so are 100% dependent on exit velocity.

fgros = mflow * uex / g0 + (pexit - pamb) * aexit ;
'''


def thrust_rocket(massflow, exit_velocity, exit_area, exit_pressure, free_stream_pressure):
    return thrust(0, massflow, exit_velocity,
                  0, exit_area, exit_pressure,
                  free_stream_pressure)  # massflow * exit_velocity + exit_area * (exit_pressure - free_stream_pressure)


def mass_flow_rate(mach_number, specific_heat_ratio, specific_gas_constant, temperature, total_pressure, area):
    M, gam, R, T, p_t, A = mach_number, specific_heat_ratio, specific_gas_constant, temperature, total_pressure, area

    a = sqrt(gam * R * T)  # speed of sound
    V = M * a
    M = V / a  # Mach number
    # T/Tt = (1 + .5 * (gam -1) * M^2) ^-1
    # p = p_t * (T / T_t) ** (gam / (gam - 1))
    ratio_T2T_t = 1 / (1 + .5 * (gam - 1) * M ** 2)
    p = p_t * ratio_T2T_t ** (gam / (gam - 1))
    r = p / (R * T)  # Equation of State solved for density

    mdot = r * V * A
    return mdot


# def equivalent_velocity(free_stream_massflow, exit_massflow, exit_velocity, free_stream_velocity, exit_area,
#                         exit_pressure, free_stream_pressure):
#     return thrust(free_stream_massflow, exit_massflow, exit_velocity, free_stream_velocity, exit_area, exit_pressure,
#                   free_stream_pressure
#                   )/(exit_massflow*g0)

def isp(massflow, exit_velocity, exit_area, exit_pressure, free_stream_pressure):
    return thrust_rocket(massflow, exit_velocity, exit_area, exit_pressure, free_stream_pressure) / (massflow * g0)


def equation_of_state(pressure=None, volume=None, moles=None, specific_gas_constant=None, temperature=None):
    p, v, n, R, T = pressure, volume, moles, specific_gas_constant, temperature
    degrees_of_freedom = sum(parameter is None for parameter in [p, n, R, T])
    assert degrees_of_freedom <= 1, "Cannot solve for multiple parameters"
    if degrees_of_freedom == 0:
        return p * v == n * R * T
    elif p is None:
        return (n * R * T) / v
    elif v is None:
        return (n * R * T) / p
    elif n is None:
        return (p * v) / (R * T)
    elif R is None:
        return (p * v) / (n * T)
    elif T is None:
        return (p * v) / (n * R)
    else:
        raise Exception('wut')


# print(mass_flow_rate(2, 1.4, 1716, 100, 2000, 200))
#
# M = 2
# gam = 1.4
# mgam = gam - 1
# pgam = gam + 1
# R = 1716
# T = 100
# pt = 2116
# Tt = 518
# A = 1.001
# val1 = A * (pt / (sqrt(Tt)))
# val2 = sqrt(gam/R)
# val3 = 1 + (0.5 * mgam * M**2)
# val4 = pgam / (2*mgam)
# mfr = val1 * val2 * M * (val3**(-1*val4))
# print(
# (A * pt/sqrt(Tt)) * sqrt(gam/R) * M * (1 + .5 * (gam-1) * M**2 )**-((gam + 1)/(gam - 1)/2)
# )
# print(val1, val2, val3, val4, mfr)

# print(thrust_rocket(643, 3743, 48387, 19.86, 101.297))
massflow, exit_velocity, exit_area, exit_pressure, free_stream_pressure = 643, 3743, 48387 / 10000, 19.86 * 1000, 101.297 * 1000
print((massflow * exit_velocity) + (exit_area * (exit_pressure - free_stream_pressure)))
print(thrust_rocket(massflow, exit_velocity, exit_area, exit_pressure, free_stream_pressure))
print(thrust_rocket(massflow, exit_velocity, exit_area, exit_pressure, free_stream_pressure) / (massflow * g0))
T = 3534.15
M = 0.0288
gam = 1.4

A0 = sqrt(gam * 8.314 * T / M)
Umax = A0 * (2 / (gam - 1)) ** 0.5
print(A0, Umax, Umax / g0)
