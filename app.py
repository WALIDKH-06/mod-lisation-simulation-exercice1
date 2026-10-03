"""Interface : flux d'entrée -> séparateur -> heater (style DWSIM).
Lancer :  python app.py   (nécessite numpy, scipy, matplotlib)
"""
import tkinter as tk
from tkinter import ttk, messagebox

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import thermo
import operations as op

SPECS = {"Fraction vapeur x (0 à 1)": "x",
         "Température de sortie": "T",
         "Chaleur Q (kW)": "Q"}


def nombre(texte):
    return float(texte.strip().replace(",", "."))


class Application(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Simulation : flux → séparateur → heater")
        self.composes = thermo.charger_composes()
        self.lignes = []                      # lignes (nom, %) du séparateur
        self._construire()
        self._creer_lignes()
        self._maj_unites_debit()
        self._maj_spec()

    # ------------------------------------------------------------ interface
    def _construire(self):
        gauche = ttk.Frame(self, padding=8)
        gauche.grid(row=0, column=0, sticky="n")
        droite = ttk.Frame(self, padding=8)
        droite.grid(row=0, column=1, sticky="nsew")
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        # 0. Modèle et composé
        f0 = ttk.LabelFrame(gauche, text="0. Modèle thermodynamique et composé", padding=6)
        f0.grid(sticky="ew", pady=3)
        self.modele = ttk.Combobox(f0, values=list(thermo.MODELES), state="readonly", width=18)
        self.modele.current(0)
        self.compose = ttk.Combobox(f0, values=list(self.composes), state="readonly", width=18)
        self.compose.current(0)
        ttk.Label(f0, text="Modèle").grid(row=0, column=0, sticky="w")
        self.modele.grid(row=0, column=1, padx=4)
        ttk.Label(f0, text="Composé").grid(row=1, column=0, sticky="w")
        self.compose.grid(row=1, column=1, padx=4)
        self.compose.bind("<<ComboboxSelected>>", lambda e: self._maj_flux())

        # 1. Flux d'entrée
        f1 = ttk.LabelFrame(gauche, text="1. Flux d'entrée F", padding=6)
        f1.grid(sticky="ew", pady=3)
        self.T = tk.StringVar(value="10")
        self.uT = ttk.Combobox(f1, values=list(op.UNITES_T), state="readonly", width=6)
        self.uT.current(1)
        self.P = tk.StringVar(value="1")
        self.uP = ttk.Combobox(f1, values=list(op.UNITES_P), state="readonly", width=6)
        self.uP.current(3)
        ttk.Label(f1, text="Température").grid(row=0, column=0, sticky="w")
        ttk.Entry(f1, textvariable=self.T, width=10).grid(row=0, column=1)
        self.uT.grid(row=0, column=2)
        ttk.Label(f1, text="Pression").grid(row=1, column=0, sticky="w")
        ttk.Entry(f1, textvariable=self.P, width=10).grid(row=1, column=1)
        self.uP.grid(row=1, column=2)

        self.type_debit = ttk.Combobox(f1, values=list(op.UNITES_DEBIT), state="readonly", width=10)
        self.type_debit.current(0)
        self.type_debit.bind("<<ComboboxSelected>>", lambda e: self._maj_unites_debit())
        self.D = tk.StringVar(value="30")
        self.uD = ttk.Combobox(f1, state="readonly", width=6)
        self.uD.bind("<<ComboboxSelected>>", lambda e: self._maj_flux())
        ttk.Label(f1, text="Débit (un seul champ)").grid(row=2, column=0, columnspan=3, sticky="w", pady=(6, 0))
        self.type_debit.grid(row=3, column=0, sticky="w")
        ttk.Entry(f1, textvariable=self.D, width=10).grid(row=3, column=1)
        self.uD.grid(row=3, column=2)
        self.affichage_flux = ttk.Label(f1, text="", foreground="#555")
        self.affichage_flux.grid(row=4, column=0, columnspan=3, sticky="w", pady=(4, 0))
        for var in (self.T, self.P, self.D):
            var.trace_add("write", lambda *a: self._maj_flux())
        for cb in (self.uT, self.uP):
            cb.bind("<<ComboboxSelected>>", lambda e: self._maj_flux())

        # 2. Séparateur
        f2 = ttk.LabelFrame(gauche, text="2. Séparateur", padding=6)
        f2.grid(sticky="ew", pady=3)
        ttk.Label(f2, text="Nombre de flux").grid(row=0, column=0, sticky="w")
        self.N = tk.IntVar(value=2)
        ttk.Spinbox(f2, from_=1, to=6, textvariable=self.N, width=4,
                    command=self._creer_lignes).grid(row=0, column=1, sticky="w")
        ttk.Button(f2, text="Créer", command=self._creer_lignes).grid(row=0, column=2, padx=4)
        self.zone_lignes = ttk.Frame(f2)
        self.zone_lignes.grid(row=1, column=0, columnspan=3, pady=4)
        ttk.Label(f2, text="Flux envoyé au heater").grid(row=2, column=0, columnspan=2, sticky="w")
        self.vers_heater = ttk.Combobox(f2, state="readonly", width=8)
        self.vers_heater.grid(row=2, column=2)

        # 3. Heater
        f3 = ttk.LabelFrame(gauche, text="3. Heater", padding=6)
        f3.grid(sticky="ew", pady=3)
        ttk.Label(f3, text="ΔP (même unité que P)").grid(row=0, column=0, sticky="w")
        self.dP = tk.StringVar(value="0")
        ttk.Entry(f3, textvariable=self.dP, width=10).grid(row=0, column=1)
        ttk.Label(f3, text="Spécification").grid(row=1, column=0, sticky="w")
        self.spec = ttk.Combobox(f3, values=list(SPECS), state="readonly", width=24)
        self.spec.current(0)
        self.spec.bind("<<ComboboxSelected>>", lambda e: self._maj_spec())
        self.spec.grid(row=1, column=1, columnspan=2, padx=4)
        self.valeur = tk.StringVar(value="1")
        ttk.Label(f3, text="Valeur").grid(row=2, column=0, sticky="w")
        ttk.Entry(f3, textvariable=self.valeur, width=10).grid(row=2, column=1)
        self.unite_spec = ttk.Label(f3, text="")
        self.unite_spec.grid(row=2, column=2, sticky="w")

        # Courbe
        f4 = ttk.LabelFrame(gauche, text="Courbe T_eb = f(P)", padding=6)
        f4.grid(sticky="ew", pady=3)
        self.pmin, self.pmax, self.npts = (tk.StringVar(value=v) for v in ("0.1", "3", "25"))
        for i, (txt, var) in enumerate((("P min (atm)", self.pmin), ("P max (atm)", self.pmax),
                                        ("Nb de points", self.npts))):
            ttk.Label(f4, text=txt).grid(row=i, column=0, sticky="w")
            ttk.Entry(f4, textvariable=var, width=8).grid(row=i, column=1)

        boutons = ttk.Frame(gauche)
        boutons.grid(pady=6)
        ttk.Button(boutons, text="Calculer", command=self.calculer).pack(side="left", padx=4)
        ttk.Button(boutons, text="Tracer la courbe", command=self.tracer).pack(side="left", padx=4)

        # Résultats
        ttk.Label(droite, text="Résultats").pack(anchor="w")
        self.sortie = tk.Text(droite, width=62, height=36, font=("Consolas", 10), state="disabled")
        self.sortie.pack(fill="both", expand=True)

    def _creer_lignes(self):
        for w in self.zone_lignes.winfo_children():
            w.destroy()
        n = max(1, min(6, self.N.get()))
        defauts = [("f1", "40"), ("f2", "60")] if n == 2 else \
                  [(f"f{i+1}", f"{100 / n:g}") for i in range(n)]
        self.lignes = []
        ttk.Label(self.zone_lignes, text="Nom").grid(row=0, column=0)
        ttk.Label(self.zone_lignes, text="%").grid(row=0, column=1)
        for i, (nom, pct) in enumerate(defauts):
            vn, vp = tk.StringVar(value=nom), tk.StringVar(value=pct)
            ttk.Entry(self.zone_lignes, textvariable=vn, width=10).grid(row=i + 1, column=0, padx=2)
            ttk.Entry(self.zone_lignes, textvariable=vp, width=8).grid(row=i + 1, column=1, padx=2)
            vn.trace_add("write", lambda *a: self._maj_noms())
            self.lignes.append((vn, vp))
        self._maj_noms(choisir_dernier=True)

    def _maj_noms(self, choisir_dernier=False):
        noms = [vn.get() for vn, _ in self.lignes]
        courant = self.vers_heater.get()
        self.vers_heater["values"] = noms
        if choisir_dernier or courant not in noms:
            self.vers_heater.set(noms[-1] if noms else "")

    def _maj_unites_debit(self):
        unites = list(op.UNITES_DEBIT[self.type_debit.get()])
        self.uD["values"] = unites
        self.uD.current(0)
        self._maj_flux()

    def _maj_spec(self):
        s = SPECS[self.spec.get()]
        self.unite_spec.config(text={"x": "(0 à 1)", "T": self.uT.get(), "Q": "kW"}[s])
        self.valeur.set({"x": "1", "T": "", "Q": ""}[s])

    # ------------------------------------------------------------ lecture
    def _c(self):
        return self.composes[self.compose.get()]

    def _modele(self):
        return thermo.MODELES[self.modele.get()](self._c())

    def _flux_entree(self):
        T = op.UNITES_T[self.uT.get()][0](nombre(self.T.get()))
        P = nombre(self.P.get()) * op.UNITES_P[self.uP.get()]
        m = op.debit_massique_SI(self._c(), self.type_debit.get(), nombre(self.D.get()), self.uD.get())
        return T, P, m

    def _maj_flux(self):
        try:
            _, _, m = self._flux_entree()
            d = op.debits_phase(None, self._c(), m, "L", 0, 1)
            self.affichage_flux.config(
                text=f"= {d['kg/h']:.4g} kg/h | {d['mol/s']:.4g} mol/s | {d['m³/h']:.4g} m³/h")
        except Exception:
            self.affichage_flux.config(text="")

    # ------------------------------------------------------------ calcul
    def calculer(self):
        try:
            c, modele = self._c(), self._modele()
            T, P, m = self._flux_entree()
            noms = [vn.get().strip() for vn, _ in self.lignes]
            pcts = [nombre(vp.get()) for _, vp in self.lignes]
            flux = op.separer(m, noms, pcts)
            cible = self.vers_heater.get()
            if cible not in flux:
                raise ValueError("Choisissez le flux envoyé au heater.")
            dP = nombre(self.dP.get()) * op.UNITES_P[self.uP.get()]
            s = SPECS[self.spec.get()]
            v = nombre(self.valeur.get())
            if s == "T":
                v = op.UNITES_T[self.uT.get()][0](v)
            res = op.chauffer(modele, c, T, P, flux[cible], dP, s, v)
            phase = op.phase_flux(modele, T, P)
        except Exception as e:
            messagebox.showerror("Erreur", str(e))
            return

        L = [f"Modèle : {modele.nom} | Composé : {self.compose.get()}", "",
             "=== Flux F ===",
             f"T = {T - 273.15:.2f} °C ({T:.2f} K) | P = {P:.0f} Pa | "
             f"phase : {'liquide' if phase == 'L' else 'vapeur'}"]
        d = op.debits_phase(modele, c, m, phase, T, P)
        L.append(f"Débit : {d['kg/h']:.4g} kg/h | {d['mol/s']:.4g} mol/s | {d['m³/h']:.4g} m³/h")
        L += ["", "=== Séparateur (T, P, composition inchangées) ==="]
        for n, mi in flux.items():
            di = op.debits_phase(modele, c, mi, phase, T, P)
            L.append(f"{n}: {di['kg/h']:.4g} kg/h | {di['mol/s']:.4g} mol/s | {di['m³/h']:.4g} m³/h")
        L += ["", f"=== Heater (alimenté par {cible}) ===",
              f"P_out = {res['P_out']:.0f} Pa",
              f"Température d'ébullition T_eb = {res['T_eb'] - 273.15:.2f} °C ({res['T_eb']:.2f} K)",
              f"T_out = {res['T_out'] - 273.15:.2f} °C | fraction vapeur x = {res['x']:.4g}",
              f"Q = {res['Q_kW']:.4g} kW"]
        for cle, titre in (("vapeur", "Vapeur"), ("liquide", "Liquide")):
            dd = res[cle]
            if dd:
                L.append(f"{titre} : {dd['kg/h']:.4g} kg/h | {dd['mol/s']:.4g} mol/s | {dd['m³/h']:.4g} m³/h")
        self.sortie.config(state="normal")
        self.sortie.delete("1.0", "end")
        self.sortie.insert("1.0", "\n".join(L))
        self.sortie.config(state="disabled")

    def tracer(self):
        try:
            modele = self._modele()
            n = int(nombre(self.npts.get()))
            if n < 2:
                raise ValueError("Il faut au moins 2 points.")
            pmin, pmax = nombre(self.pmin.get()), nombre(self.pmax.get())
            if not 0 < pmin < pmax:
                raise ValueError("Il faut 0 < P min < P max.")
            Ps, Ts = op.courbe_ebullition(modele, pmin * 101325, pmax * 101325, n)
        except Exception as e:
            messagebox.showerror("Erreur", str(e))
            return
        fen = tk.Toplevel(self)
        fen.title("Courbe d'ébullition")
        fig = Figure(figsize=(6, 4.2))
        ax = fig.add_subplot(111)
        ax.plot([p / 101325 for p in Ps], [t - 273.15 for t in Ts], "o-", ms=3)
        ax.set_xlabel("Pression (atm)")
        ax.set_ylabel("T d'ébullition (°C)")
        ax.set_title(f"{self.compose.get()} : T_eb = f(P), {modele.nom}")
        ax.grid(True)
        FigureCanvasTkAgg(fig, master=fen).get_tk_widget().pack(fill="both", expand=True)


if __name__ == "__main__":
    Application().mainloop()
