"""Flux, séparateur, heater et courbe d'ébullition."""
from scipy.optimize import brentq

from thermo import R

# Conversions d'unités vers le SI
UNITES_T = {"K": (lambda v: v, lambda k: k),
            "°C": (lambda v: v + 273.15, lambda k: k - 273.15)}
UNITES_P = {"Pa": 1.0, "kPa": 1e3, "bar": 1e5, "atm": 101325.0}
UNITES_DEBIT = {
    "massique": {"kg/h": 1 / 3600, "kg/s": 1.0},
    "molaire": {"mol/s": 1.0, "kmol/h": 1000 / 3600},
    "volumique": {"m³/h": 1 / 3600, "m³/s": 1.0, "L/h": 1e-3 / 3600},
}


def debit_massique_SI(c, type_debit, valeur, unite):
    """Convertit le débit saisi (un seul champ) en kg/s."""
    f = UNITES_DEBIT[type_debit][unite]
    if type_debit == "massique":
        return valeur * f
    if type_debit == "molaire":
        return valeur * f * c["M"] / 1000
    return valeur * f * c["rho_liq"]          # volumique (liquide)


def debits_phase(modele, c, m, phase, T, P):
    """Les trois débits d'une phase : kg/h, mol/s, m³/h."""
    n = m / (c["M"] / 1000)
    if phase == "L":
        V = m / c["rho_liq"]
    else:
        V = n * modele.Z_vap(T, P) * R * T / P
    return {"kg/h": m * 3600, "mol/s": n, "m³/h": V * 3600}


def phase_flux(modele, T, P):
    return "L" if T < modele.Tsat(P) else "V"


def separer(m, noms, pourcentages):
    if len(set(noms)) != len(noms) or any(not n for n in noms):
        raise ValueError("Les noms des flux doivent être non vides et différents.")
    if any(p < 0 for p in pourcentages):
        raise ValueError("Un pourcentage est négatif.")
    s = sum(pourcentages)
    if abs(s - 100) > 1e-6:
        raise ValueError(f"La somme des pourcentages vaut {s:.4g} % (il faut 100 %).")
    return {n: m * p / 100 for n, p in zip(noms, pourcentages)}


def chauffer(modele, c, T_in, P_in, m, dP, spec, valeur):
    """Heater. spec : 'x' (fraction vapeur), 'T' (T_out en K) ou 'Q' (kW).
    m en kg/s. Retourne un dictionnaire de résultats."""
    P_out = P_in - dP
    if P_out <= 0:
        raise ValueError("P_out doit être positive (ΔP trop grand).")
    Tb = modele.Tsat(P_out)
    h_in = modele.h(T_in, P_in, phase_flux(modele, T_in, P_in))
    hL = modele.h(Tb, P_out, "L")
    hV = modele.h(Tb, P_out, "V")

    if spec == "x":
        x = valeur
        if not 0 <= x <= 1:
            raise ValueError("La fraction vapeur doit être comprise entre 0 et 1.")
        T_out, h_out = Tb, (1 - x) * hL + x * hV
    elif spec == "T":
        T_out = valeur
        if abs(T_out - Tb) < 1e-6:
            raise ValueError("T_out = T_eb : la fraction vapeur est indéterminée, "
                             "utilisez la spécification « Fraction vapeur ».")
        x = 0.0 if T_out < Tb else 1.0
        h_out = modele.h(T_out, P_out, "L" if x == 0 else "V")
    elif spec == "Q":
        h_out = h_in + valeur / m
        try:
            if h_out <= hL:
                x, T_out = 0.0, brentq(lambda T: modele.h(T, P_out, "L") - h_out, 100.0, Tb)
            elif h_out <= hV:
                x, T_out = (h_out - hL) / (hV - hL), Tb
            else:
                x, T_out = 1.0, brentq(lambda T: modele.h(T, P_out, "V") - h_out, Tb, Tb + 800)
        except ValueError:
            raise ValueError("Q hors du domaine de calcul.")
    else:
        raise ValueError("Spécification inconnue.")

    Q = m * (h_out - h_in)                       # kW
    res = {"P_out": P_out, "T_eb": Tb, "T_out": T_out, "x": x, "Q_kW": Q}
    res["vapeur"] = debits_phase(modele, c, x * m, "V", T_out, P_out) if x > 0 else None
    res["liquide"] = debits_phase(modele, c, (1 - x) * m, "L", T_out, P_out) if x < 1 else None
    return res


def courbe_ebullition(modele, P_min, P_max, n):
    """Retourne (liste de P en Pa, liste de T_eb en K)."""
    Ps = [P_min + (P_max - P_min) * i / (n - 1) for i in range(n)]
    return Ps, [modele.Tsat(P) for P in Ps]
