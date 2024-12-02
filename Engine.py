import math


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

    def __init__(self, massflow: float, fluid: Fluid, temperature: float, pressure: float, area_in: float,
                 area_out: float):
        self.massflow = massflow
        self.fluid = fluid
        self.temperature = temperature
        self.pressure = pressure
        self.density = pressure * fluid.molar_mass / (fluid.R_sp * temperature)
        self.area_in = area_in
        self.area_out = area_out

    @property
    def volume_flow(self):
        return self.massflow * self.fluid.gamma

    @property
    def speed_of_sound(self):
        return math.sqrt(self.fluid.gamma * self.fluid.R_sp * self.temperature)

    @property
    def mach_in(self):
        return self.flow_velocity(self.area_in) / self.speed_of_sound

    @property
    def mach_out(self):
        return self.flow_velocity(self.area_out) / self.speed_of_sound

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
        return self.massflow / (self.area_out * r)


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

    def __init__(self, resistive_coefficient: float, flow_in: Flow, area_in: float, area_out: float = None):
        self.resistive_coefficient = resistive_coefficient
        self.area = area_in
        if area_out is None:
            self.area_out = area_in  # TODO: This currently does nothing
        self.flow_in = flow_in
        self.dp = resistive_coefficient * flow_in.density * flow_in.flow_velocity(area_in) ** 2 / 2
        self.flow_out = Flow(flow_in.massflow, flow_in.fluid, flow_in.temperature, flow_in.pressure - self.dp)


class Valve(Component):

    def __init__(self, setting: float, max_setting: float, resistive_coefficient: float, flow_in: Flow, area_in: float,
                 area_out: float = None):
        super().__init__(resistive_coefficient, flow_in, area_in, area_out)
        self.setting = setting
        self.max_setting = max_setting

    @property
    def flow(self):
        return self.flow_fraction * self.max_setting


class Turbine:
    # Shaft Power = Mass Flow Rate x Specific Heat x Temperature Drop
    # Shaft Pawer = eff * c_p * Tt_in * (1 - TPR ** ((gamma - 1)/gamma) )

    def __init__(self, efficiency: float, flow_in: Flow, pressure_ratio: float, inlet_area: float, outlet_area: float):
        self.efficiency = efficiency
        self.flow_in = flow_in
        self.pressure_ratio = pressure_ratio
        self.power = efficiency * flow_in.fluid.cp * flow_in.temperature_total(area=inlet_area) * (
                1 - pressure_ratio ** ((flow_in.fluid.gamma - 1) / flow_in.fluid.gamma))
        self.temperature_ratio = pressure_ratio ** (1 - 1 / flow_in.fluid.gamma)
        self.flow_out = Flow(flow_in.massflow, flow_in.fluid,
                             flow_in.temperature * self.temperature_ratio,
                             flow_in.pressure * pressure_ratio)


class Compressor:

    def __init__(self, efficiency: float, flow_in: Flow, pressure_ratio: float, inlet_area: float):
        self.efficiency = efficiency
        self.flow_in = flow_in
        self.pressure_ratio = pressure_ratio
        self.work = (flow_in.fluid.cp * flow_in.temperature_total(area=inlet_area) * pressure_ratio ** (
                (flow_in.fluid.gamma - 1) / flow_in.fluid.gamma) - 1) / efficiency
        self.power = self.work * flow_in.massflow
        self.temperature_ratio = pressure_ratio ** (1 - 1 / flow_in.fluid.gamma)
        self.flow_out = Flow(flow_in.massflow, flow_in.fluid,
                             flow_in.temperature * self.temperature_ratio,
                             flow_in.pressure * pressure_ratio)


class Plenum:


if __name__ == '__main__':
    hydrogen = Fluid(molar_mass=2.02, cp=14.29, cv=10.16)
    oxygen = Fluid(molar_mass=32.0, cp=29.38, cv=21.0)
    mix = FluidMix([hydrogen, oxygen], [2 / 3, 1 / 3])
    source = Source(500, 600, mix)
    flow = Flow(80, mix, 500, 600)
    valve = Valve(0.01, 1.0, 0.001, flow, area=.1)
    compressor = Compressor(.9, flow, 5, .1)
    print()
