import math
from functools import reduce


# from enum import Enum


class Fluid:

    def __init__(self, molar_mass: float, cp: float, cv: float):
        self.molar_mass = molar_mass  # g/mol
        self.cp = cp
        self.cv = cv
        self.gamma = cp / cv
        self.R_sp = cp - cv


class FluidMix(Fluid):
    def __init__(self, constituents: list[Fluid], molar_fractions: list[float]):
        self.constituents = constituents
        self.molar_fractions = molar_fractions
        self.molar_mass = sum(f * c.molar_mass for f, c in zip(molar_fractions, constituents))
        self.cp = sum(f * c.cp for f, c in zip(molar_fractions, constituents))
        self.cv = sum(f * c.cv for f, c in zip(molar_fractions, constituents))
        super().__init__(self.molar_mass, self.cp, self.cv)
        # self.gamma = sum(f * c.gamma for f, c in zip(molar_fractions, constituents))
        # self.R_sp = sum(f * c.R_sp for f, c in zip(molar_fractions, constituents))


class Source:

    def __init__(self, temperature: float, pressure: float, fluid_mix: FluidMix):
        self.temperature = temperature
        self.pressure = pressure
        self.fluid = fluid_mix
        self.density = pressure * fluid_mix.molar_mass / (fluid_mix.R_sp * temperature)


class Flow:

    def __init__(self, fluid: Fluid, area_in: float = None, massflow: float = None,
                 temperature: float = None, pressure: float = None, area_exit: float = None):
        self.massflow = massflow
        self.fluid = fluid
        self.temperature = temperature
        self.pressure = pressure
        self.density = pressure * fluid.molar_mass / (fluid.R_sp * temperature)
        self.area_in = area_in
        if area_exit is not None:
            self.area_exit = area_exit
        else:
            self.area_exit = area_in

    @property
    def volume_flow(self):
        return self.massflow * self.fluid.gamma

    @property
    def speed_of_sound(self):
        return math.sqrt(self.fluid.gamma * self.fluid.R_sp * self.temperature)

    @property
    def mach_in(self):
        return self.flow_velocity_in / self.speed_of_sound

    @property
    def mach_out(self):
        return self.flow_velocity_out / self.speed_of_sound

    @property
    def temperature_total_in(self):
        return self.temperature * (1 + self.mach_in ** 2 * (self.fluid.gamma - 1) / 2)

    @property
    def temperature_total_out(self):
        return self.temperature * (1 + self.mach_out ** 2 * (self.fluid.gamma - 1) / 2)

    @property
    def flow_velocity_in(self):
        r = (self.pressure / (self.fluid.R_sp * self.temperature))
        return self.massflow / (self.area_in * r)

    @property
    def flow_velocity_out(self):
        r = (self.pressure / (self.fluid.R_sp * self.temperature))
        return self.massflow / (self.area_exit * r)

    def combine_flows(self, other, mixing_area):

        f1 = self.fluid
        f2 = other.fluid

        total_flow = self.massflow + other.massflow
        mixture_temp = (self.temperature * self.massflow * self.fluid.cv + other.temperature * other.massflow * other.fluid.cv) / (
                    total_flow * (self.fluid.cv + other.fluid.cv))
        mixture_pres = (self.pressure * self.massflow + other.pressure * other.massflow) / total_flow
        if f1 == f2:
            return Flow(f1, mixing_area, massflow=total_flow, temperature=mixture_temp, pressure=mixture_pres)

        constituents1 = [f1] if type(f1) == Fluid else f1.constituents
        constituents2 = [f1] if type(f1) == Fluid else f1.constituents

        moles1 = [self.massflow / f1.molar_mass] \
            if type(f1) == Fluid else \
            [self.massflow / c.molar_mass for c in f1.constituents]
        moles2 = [other.massflow / f2.molar_mass] \
            if type(f2) == Fluid else \
            [other.massflow / c.molar_mass for c in f2.constituents]

        molar_sum_1 = sum(moles1)
        molar_sum_2 = sum(moles2)

        moles = moles1 + moles2
        molar_sum = sum(moles)

        molar_fractions = [m/molar_sum for m in moles]

        fluid = FluidMix(constituents1 + constituents2, molar_fractions)

        return Flow(fluid, mixing_area, massflow=total_flow, temperature=mixture_temp, pressure=mixture_pres)


class Component:
    """
    =================================================================================
    Generalized Energy Equation:
    (p1/gamma) + z1 + (v1**2/(2g)) + h_A - h_R -h_L = (p2/gamma) + z_2 + (v2**2/(2g))
    =================================================================================
    p = pressure
    z = elevation
    v = velocity
    h_A = energy added
    h_R = energy removed
    h_L = energy lost

    ====================
    Continuity Equation:
    Q = vA
    ====================

    Q = Volume Flow Rate
    A = cross sectional area

    ==================================
    Friction Losses - Darcy's Equation
    h_L = f * (L/D) * (v**2/(2g))
    ==================================

    f = friction factor
    L = pipe length
    D = pipe diameter

    =====================
    Minor Losses
    h_L = K * (v**2/(2g))
    p_L -  K * (rho * (v**2))/2
    =====================

    K = Resistive Coefficient
    """

    def __init__(self, resistive_coefficient: float, flow_in: Flow, area_in: float, area_exit: float = None):
        self.resistive_coefficient = resistive_coefficient
        self.area = area_in
        if area_exit is None:
            self.area_exit = area_in  # TODO: This currently does nothing
        self.flow_in = flow_in
        self.dp = resistive_coefficient * flow_in.density * flow_in.flow_velocity_out ** 2 / 2
        self.flow_out = Flow(self.flow_in.fluid, self.area_exit, massflow=flow_in.massflow, temperature=flow_in.temperature,
                             pressure=flow_in.pressure - self.dp)


class Valve(Component):

    def __init__(self, setting: float, max_setting: float, resistive_coefficient: float, flow_in: Flow, area_in: float,
                 area_exit: float = None):
        super().__init__(resistive_coefficient, flow_in, area_in, area_exit)
        self.setting = setting
        self.max_setting = max_setting

    @property
    def flow(self):
        return self.setting * self.max_setting


class Turbine:
    # Shaft Power = Mass Flow Rate x Specific Heat x Temperature Drop
    # Shaft Power = eff * c_p * Tt_in * (1 - TPR ** ((gamma - 1)/gamma) )

    def __init__(self, efficiency: float, flow_in: Flow, pressure_ratio: float, inlet_area: float, outlet_area: float):
        self.efficiency = efficiency
        self.flow_in = flow_in
        self.pressure_ratio = pressure_ratio
        self.power = efficiency * flow_in.fluid.cp * flow_in.temperature_total_out * (
                1 - pressure_ratio ** ((flow_in.fluid.gamma - 1) / flow_in.fluid.gamma))
        self.temperature_ratio = pressure_ratio ** (1 - 1 / flow_in.fluid.gamma)
        self.flow_out = Flow(flow_in.fluid, massflow=flow_in.massflow,
                             temperature=flow_in.temperature * self.temperature_ratio,
                             pressure=flow_in.pressure * pressure_ratio)


class Compressor:

    def __init__(self, efficiency: float, flow_in: Flow, pressure_ratio: float, inlet_area: float):
        self.efficiency = efficiency
        self.flow_in = flow_in
        self.pressure_ratio = pressure_ratio
        self.work = (flow_in.fluid.cp * flow_in.temperature_total_out * pressure_ratio ** (
                (flow_in.fluid.gamma - 1) / flow_in.fluid.gamma) - 1) / efficiency
        self.power = self.work * flow_in.massflow
        self.temperature_ratio = pressure_ratio ** (1 - 1 / flow_in.fluid.gamma)
        self.flow_out = Flow(flow_in.fluid, massflow=flow_in.massflow,
                             temperature=flow_in.temperature * self.temperature_ratio,
                             pressure=flow_in.pressure * pressure_ratio)


class Plenum:

    def __init__(self, in_flows: list[Flow], area_exit=None):
        self.in_flows = in_flows
        self.area_exit = area_exit
        in_flow_area = sum(f.area_exit for f in self.in_flows)
        self.flow_out = reduce(lambda x, y: x.combine_flows(y, in_flow_area), in_flows)

    @property
    def area_throat(self):
        return self.area_exit


class Expander:

    def __init__(self, flow_in: Flow, area_in: float, area_exit: float = None):
        self.flow_in = flow_in
        self.flow_out = flow_out
        self.area_in = area_in
        self.area_exit = area_exit

    @property
    def area_throat(self):
        return self.area_in


class Nozzle:

    def __init__(self, in_flows: list[Flow], area_throat: float, area_exit: float):
        self.in_flows = in_flows
        # self.flow_out = flow_out
        self.area_throat = area_throat
        self.area_exit = area_exit

        self.plenum = Plenum(in_flows, area_throat)
        self.bell = Expander(self.plenum.flow_out, )


if __name__ == '__main__':
    hydrogen = Fluid(molar_mass=2.02, cp=14.29, cv=10.16)
    oxygen = Fluid(molar_mass=32.0, cp=29.38, cv=21.0)
    mix = FluidMix([hydrogen, oxygen], [2 / 3, 1 / 3])
    source = Source(500, 600, mix)
    flow = Flow(mix, massflow=80, temperature=500, pressure=600)
    valve = Valve(0.01, 1.0, 0.001, flow, area=.1)
    compressor = Compressor(.9, flow, 5, .1)
    print()
