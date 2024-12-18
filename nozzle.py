from enum import Enum
from utils import *


class ShockMode(Enum):
    under_expanded = 0
    slight_overexpanded = 1
    shock_in_nozzle = 2


class Nozzle:

    def __init__(self, molecular_weight_exhaust, area_throat, area_ratio_exit, area_ratio_plenum,
                 total_pressure_chamber, temperature_total,
                 specific_heat_ratio_std):

        self.mach_exit = None
        self.pressure_exit = None
        self.velocity_exhaust = None
        self.nozzle_pressure_ratio = None
        self.specific_heat_ratio = None
        self.massflow = None
        self.thrust_gross = None
        self.mode = None
        self.specific_heat_ratio_std = specific_heat_ratio_std
        self.temperature_total = temperature_total
        self.total_pressure_chamber = total_pressure_chamber
        # self.area_ratio_plenum = area_ratio_plenum
        self.area_ratio_exit = area_ratio_exit
        self.area_throat = area_throat
        self.molecular_weight_exhaust = molecular_weight_exhaust

        self.R_gas = R_u * g_0 / self.molecular_weight_exhaust
        self.radius_throat = math.sqrt(self.area_throat / math.pi)
        self.area_exit = self.area_ratio_exit * self.area_throat
        self.radius_exit = math.sqrt(self.area_exit / math.pi)
        # self.area_plenum = self.area_ratio_plenum * self.area_throat
        # self.radius_plenum = math.sqrt(self.area_plenum / math.pi)
        self.specific_heat_ratio = self.specific_heat_ratio_std * calc_specific_heat_ratio(self.temperature_total) / 1.4
        # self.pt = self.total_pressure_chamber / kP_to_psi

        # self.tt = (self.temperature_total + temperature_reference) / kelvin_per_rankine

    def compute(self, pressure_free_stream):
        counter = 0
        corrected_flow_per_area = weightflow_per_area_given_mach(1.0, self.R_gas, self.specific_heat_ratio)
        self.mach_exit = get_mach(2, (corrected_flow_per_area / self.area_ratio_exit), self.R_gas, self.specific_heat_ratio)
        psup = get_pressure_ratio_isentropic(self.mach_exit, self.specific_heat_ratio) * self.total_pressure_chamber
        self.massflow = corrected_flow_per_area * self.area_throat * (self.total_pressure_chamber / p_0) / math.sqrt(
            self.temperature_total / T_0)

        if pressure_free_stream <= psup:
            self.mode = ShockMode.under_expanded
            temp_ratio = get_temp_ratio_isentropic(self.mach_exit, self.specific_heat_ratio)
            self.velocity_exhaust = self.mach_exit * math.sqrt(
                self.specific_heat_ratio * self.R_gas * temp_ratio * self.temperature_total)
            self.pressure_exit = psup

        # over expanded nozzle
        elif pressure_free_stream > psup:
            # find exit pressure at which normal shock leaves the nozzle
            mach_number_supersonic = self.mach_exit
            psub = psup * normal_shock_static_pressure_ratio(mach_number_supersonic, self.specific_heat_ratio)

            # slightly overexpanded - no shock in nozzle
            if pressure_free_stream <= psub:
                self.mode = ShockMode.slight_overexpanded
                self.pressure_exit = psup
                temp_ratio = get_temp_ratio_isentropic(self.mach_exit, self.specific_heat_ratio)
                self.velocity_exhaust = self.mach_exit * math.sqrt(
                    self.specific_heat_ratio * self.R_gas * temp_ratio * self.temperature_total)

            # highly overexpanded - normal shock in nozzle
            if pressure_free_stream > psub:
                self.mode = ShockMode.shock_in_nozzle
                self.pressure_exit = pressure_free_stream
                anso = self.area_exit
                mach_number_supersonic = self.mach_exit
                mach_number_subsonic = normal_shock_mach_after(mach_number_supersonic, self.specific_heat_ratio)
                total_pressure_ratio = normal_shock_total_pressure_ratio(mach_number_supersonic,
                                                                         self.specific_heat_ratio)
                pso = get_pressure_ratio_isentropic(mach_number_subsonic,
                                                    self.specific_heat_ratio) * total_pressure_ratio * self.total_pressure_chamber
                ansn = anso - 1.
                while (abs(self.pressure_exit - pso) > .001) and (counter < 20):
                    counter += 1
                    mach_number_supersonic = get_mach(2, (corrected_flow_per_area * self.area_throat / ansn), self.R_gas,
                                                      self.specific_heat_ratio)
                    mach_number_subsonic = normal_shock_mach_after(mach_number_supersonic, self.specific_heat_ratio)
                    total_pressure_ratio = normal_shock_total_pressure_ratio(mach_number_supersonic,
                                                                             self.specific_heat_ratio)
                    self.mach_exit = get_mach(0, (corrected_flow_per_area / self.area_ratio_exit / total_pressure_ratio), self.R_gas,
                                              self.specific_heat_ratio)
                    psn = get_pressure_ratio_isentropic(self.mach_exit,
                                                        self.specific_heat_ratio) * total_pressure_ratio * self.total_pressure_chamber
                    deriv = (psn - pso) / (ansn - anso)
                    pso = psn
                    anso = ansn
                    ansn = anso + (self.pressure_exit - pso) / deriv

                # ans = anso
                # rns = math.sqrt(ans / math.pi)
                temp_ratio = get_temp_ratio_isentropic(self.mach_exit, self.specific_heat_ratio)
                self.velocity_exhaust = self.mach_exit * math.sqrt(
                    self.specific_heat_ratio * self.R_gas * temp_ratio * self.temperature_total)
        if pressure_free_stream > .0001:
            self.nozzle_pressure_ratio = self.total_pressure_chamber / pressure_free_stream
        else:
            self.nozzle_pressure_ratio = 1000.

        self.thrust_gross = self.massflow * self.velocity_exhaust + (
                self.pressure_exit - pressure_free_stream) * self.area_exit

        return self.thrust_gross


if __name__ == '__main__':
    mwe = 16
    A_th = 967.74 / 10000
    Arat_e = 49.906
    Arat_pl = 3.9925
    P_t_chamber = 13782.0 * 1000
    T_tot = 3261.108 + T_0
    gamma_std = 1.32
    P_fs = 101.297 * 1000

    throats = [100, 1000, 5000]
    exhaust_ratios = [10, 50, 200]
    print('      \t\t\t 10\t\t\t\t\t  50 \t\t\t\t 200')
    for t in throats:
        print(f'throat {t}:   ', end='')
        for er in exhaust_ratios:
            print(Nozzle(mwe, t, er, Arat_pl, P_t_chamber, T_tot, gamma_std).compute(P_fs)/1000000, end=' ')
        print()

    # print('F', Nozzle(mwe, A_th, Arat_e, Arat_pl, P_t_chamber, T_tot, gamma_std).compute(P_fs))
