"""Modèles thermodynamiques : Raoult's Law (avec Antoine) et Peng-Robinson.

Unités internes (SI) : T en K, P en Pa, enthalpies massiques en kJ/kg.
Les constantes du composé sont lues dans composes.json :
  M (g/mol), rho_liq (kg/m3), antoine {A,B,C} (P en mmHg, T en °C),
  cp_liq (kJ/kg/K), dH_vap (kJ/kg), cp_vap (kJ/kg/K),
  Tc (K), Pc (bar), omega, cp_ig = [a,b,c,d] (J/mol/K, cp = a + bT + cT² + dT³).
"""
import json
import math
import os

import numpy as np
from scipy.optimize import brentq

R = 8.314462618      # J/mol/K
MMHG = 133.322368    # Pa
T_REF = 273.15       # K, température de référence des enthalpies


def charger_composes(chemin=None):
    chemin = chemin or os.path.join(os.path.dirname(os.path.abspath(__file__)), "composes.json")
    with open(chemin, encoding="utf-8") as f:
        return json.load(f)


class Raoult:
    """Liquide idéal + gaz parfait : P = P_sat(T), P_sat par Antoine."""
    nom = "Raoult's Law"

    def __init__(self, c):
        self.c = c
        self.A, self.B, self.C = (c["antoine"][k] for k in ("A", "B", "C"))

    def Psat(self, T):
        return 10 ** (self.A - self.B / (self.C + T - 273.15)) * MMHG

    def Tsat(self, P):
        denom = self.A - math.log10(P / MMHG)
        if denom <= 0:
            raise ValueError("Pression hors du domaine de l'équation d'Antoine.")
        return self.B / denom - self.C + 273.15

    def h(self, T, P, phase):
        """Enthalpie massique (kJ/kg), référence : liquide à 0 °C."""
        c = self.c
        if phase == "L":
            return c["cp_liq"] * (T - T_REF)
        Tb = self.Tsat(P)
        return c["cp_liq"] * (Tb - T_REF) + c["dH_vap"] + c["cp_vap"] * (T - Tb)

    def Z_vap(self, T, P):
        return 1.0


class PengRobinson:
    """Équation d'état cubique de Peng-Robinson (corps pur)."""
    nom = "Peng-Robinson"

    def __init__(self, c):
        self.c = c
        self.M = c["M"] / 1000.0
        self.Tc = c["Tc"]
        self.Pc = c["Pc"] * 1e5
        w = c["omega"]
        self.ac = 0.45724 * R**2 * self.Tc**2 / self.Pc
        self.b = 0.07780 * R * self.Tc / self.Pc
        self.kappa = 0.37464 + 1.54226 * w - 0.26992 * w**2

    def _a(self, T):
        alpha = (1 + self.kappa * (1 - math.sqrt(T / self.Tc))) ** 2
        a = self.ac * alpha
        dadT = -self.ac * self.kappa * math.sqrt(alpha) / math.sqrt(T * self.Tc)
        return a, dadT

    def _racines(self, T, P):
        a, _ = self._a(T)
        A = a * P / (R * T) ** 2
        B = self.b * P / (R * T)
        coeffs = [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]
        reels = sorted(r.real for r in np.roots(coeffs) if abs(r.imag) < 1e-9 and r.real > B)
        return reels, A, B

    @staticmethod
    def _lnphi(Z, A, B):
        s2 = math.sqrt(2)
        return (Z - 1 - math.log(Z - B)
                - A / (2 * s2 * B) * math.log((Z + (1 + s2) * B) / (Z + (1 - s2) * B)))

    def _Z(self, T, P, phase):
        reels, A, B = self._racines(T, P)
        if not reels:
            raise ValueError("Aucune racine de Z trouvée.")
        return (reels[0] if phase == "L" else reels[-1]), A, B

    def Tsat(self, P):
        """T telle que ln(phi_L) = ln(phi_V) (égalité des fugacités)."""
        def f(T):
            reels, A, B = self._racines(T, P)
            if len(reels) < 2:
                return None
            return self._lnphi(reels[0], A, B) - self._lnphi(reels[-1], A, B)

        T_prec = f_prec = None
        for T in np.linspace(0.40 * self.Tc, 0.999 * self.Tc, 400):
            fT = f(T)
            if fT is None:
                f_prec = None
                continue
            if f_prec is not None and f_prec * fT < 0:
                return brentq(f, T_prec, T, xtol=1e-9)
            T_prec, f_prec = T, fT
        raise ValueError("T_sat introuvable (P hors du domaine liquide-vapeur).")

    def Z_vap(self, T, P):
        return self._Z(T, P, "V")[0]

    def h(self, T, P, phase):
        """h = h_gaz parfait + h_résiduelle (kJ/kg), référence : gaz parfait à 0 °C."""
        Z, A, B = self._Z(T, P, phase)
        a, dadT = self._a(T)
        s2 = math.sqrt(2)
        h_res = (R * T * (Z - 1)
                 + (T * dadT - a) / (2 * s2 * self.b)
                 * math.log((Z + (1 + s2) * B) / (Z + (1 - s2) * B)))
        ca, cb, cc, cd = self.c["cp_ig"]
        h_ig = (ca * (T - T_REF) + cb / 2 * (T**2 - T_REF**2)
                + cc / 3 * (T**3 - T_REF**3) + cd / 4 * (T**4 - T_REF**4))
        return (h_ig + h_res) / self.M / 1000.0


MODELES = {Raoult.nom: Raoult, PengRobinson.nom: PengRobinson}
