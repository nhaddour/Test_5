# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 07:51:32 2023

@author: naouf
"""

import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

t_exp = None
Dg_exp = None

def fenton_reaction(t, y, I, V, F, *ki):
    Fe2p, H2O2, Fe3p, OHm, Hp, OHpt, HO2pt, O2ptm = y
    V_m3 = V / 1000

    rates = [
        ki[0] * Fe2p * H2O2,
        ki[1] * Fe3p * H2O2,
        ki[2] * H2O2 * OHpt,
        ki[3] * HO2pt,
        ki[4] * O2ptm * H2O2,
        ki[5] * OHpt * Fe2p,
        ki[6] * HO2pt * Fe2p,
        ki[7] * HO2pt * Fe3p,
        ki[8] * O2ptm * Fe2p,
        ki[9] * O2ptm * Fe3p,
        ki[10] * OHpt**2,
        ki[11] * HO2pt**2,
        ki[12] * O2ptm * Hp,
        ki[13] * OHpt * HO2pt,
        ki[14] * OHpt * O2ptm,
        ki[15] * HO2pt * O2ptm,
        ki[16] * HO2pt * H2O2,
        ki[17] * O2ptm * H2O2,
    ]

    dFe2p_dt = (
        I / (2 * F * V_m3)
        - rates[0]
        + rates[1]
        - rates[5]
        - rates[6]
        + rates[7]
        - rates[8]
        + rates[9]
    )
    dFe3p_dt = (
        rates[0]
        - rates[1]
        + rates[5]
        + rates[6]
        - rates[7]
        + rates[8]
        - rates[9]
    )
    dH2O2_dt = (
        -rates[0]
        - rates[1]
        - rates[2]
        + rates[6]
        + rates[8]
        + rates[10]
        + rates[11]
        + rates[15]
        - rates[16]
        - rates[17]
    )
    dOHm_dt = (
        rates[0]
        + rates[5]
        + rates[6]
        + 2 * rates[8]
        + rates[14]
        + rates[15]
        + rates[17]
    )
    dHp_dt = rates[1] + rates[3] - rates[4] - rates[7] - rates[12]
    dOHpt_dt = (
        rates[0]
        - rates[2]
        - rates[10]
        - rates[13]
        - rates[14]
        + rates[16]
        + rates[17]
    )
    dHO2pt_dt = (
        rates[1]
        + rates[2]
        - rates[3]
        + rates[4]
        - rates[6]
        - rates[7]
        - rates[11]
        + rates[12]
        - rates[13]
        - rates[15]
        - rates[16]
    )
    dO2ptm_dt = (
        rates[3]
        - rates[4]
        - rates[8]
        - rates[9]
        - rates[12]
        - rates[14]
        - rates[15]
        - rates[17]
    )

    dydt = [
        dFe2p_dt,
        dH2O2_dt,
        dFe3p_dt,
        dOHm_dt,
        dHp_dt,
        dOHpt_dt,
        dHO2pt_dt,
        dO2ptm_dt,
    ]
    return dydt


def degradation_model(t, ka, kb, kc, I, V, H2O2_0, ki, pH):
    
    # Constants
    F = 9.64853e4 # Faraday Constant in C.mol-1
    
    # Concentrations initiales
    Fe2p_0 = 0.0
    Fe3p_0 = 0.0
    OHm_0 = 10**(-pH)
    Hp_0 = 10**(-14) / OHm_0
    OHpt_0 = 0.0
    HO2pt_0 = 0.0
    O2ptm_0 = 0.0
    
    # Les valeurs des constantes cinétiques ki
    ki = np.array([
        6.3e-2, # 1
        2.0e-6, # 2
        3.3e4,  # 3
        1.58e5, # 4
        1.0e7,  # 5
        3.2e5,  # 6
        1.2e3,  # 7
        3.6e2,  # 8
        1e4,    # 9
        5e4,    # 10
        5.2e6,  # 11
        8.3e2,  # 12
        1e7,    # 13
        7.1e6,  # 14
        1.01e7, # 15
        9.7e4,  # 16
        5e-4,   # 17
        1.3e-4  # 18
    ])  # in mol.m-3.s-1

    # Conditions initiales
    y0 = [Fe2p_0, H2O2_0, Fe3p_0, pH, OHpt_0, HO2pt_0, O2ptm_0]

    # Résoudre le système d'équations différentielles
    sol = solve_ivp(fenton_reaction, (t[0], t[-1]), y0, args=(I, V, ki), t_eval=t, method='Radau', rtol=1e-6, max_step=1e-2)

    # Extraire les concentrations des radicaux à partir de la solution
    OHpt = sol.y[5]
    HO2pt = sol.y[6]
    O2ptm = sol.y[7]

    # Calculer la dégradation du VM
    Dg = ((ka * OHpt) + (kb * HO2pt) + (kc * O2ptm))

    return Dg


def import_excel():
    global t_exp, Dg_exp
    file_path = filedialog.askopenfilename(filetypes=[("Fichiers Excel", "*.xlsx;*.xls")])
    if file_path:
        try:
            data = pd.read_excel(file_path)
            t_exp = data['t'].to_numpy()
            Dg_exp = data['Dg'].to_numpy()
        except KeyError as e:
            messagebox.showerror("Erreur", f"Colonne manquante dans le fichier Excel : {e}")

def analyze():
    global t_exp, Dg_exp
    if t_exp is None or Dg_exp is None:
        messagebox.showerror("Erreur", "Veuillez importer un fichier Excel avant de lancer l'analyse.")
        return

    initial_guess = (
        0.1, 0.1, 0.1,  # ka, kb, kc
        float(entry_I.get()),  # I
        float(entry_V.get()),  # V
        float(entry_H2O2_0.get()),  # H2O2_0
        float(entry_pH.get()),  # pH
        9.64853e4, #F
    )
    params, pcov = curve_fit(degradation_model, t_exp, Dg_exp, p0=initial_guess)
    print("Nombre de paramètres retournés par curve_fit:", len(params))
    ka_fit = params[0]
    kb_fit = params[1]
    kc_fit = params[2]
    VM_0_fit = params[3]
    I_fit = params[4]
    V_fit = params[5]
    H2O2_0_fit = params[6]
    pH_fit = params[7]
    

    # Calcul des prédictions du modèle
    t_model = np.linspace(t_exp.min(), t_exp.max(), 100)
    Dg_model = degradation_model(t_model, ka_fit, kb_fit, kc_fit, VM_0_fit, I_fit, V_fit, H2O2_0_fit, pH_fit)

    # Tracer les résultats
    ax.clear()
    ax.plot(t_exp, Dg_exp, 'o', label='Données expérimentales')
    ax.plot(t_model, Dg_model, '-', label='Modèle ajusté')
    # Ajoutez ici le code pour calculer et afficher les incertitudes sur le graphe
    ax.set_xlabel('Temps')
    ax.set_ylabel('Dégradation du VM')
    ax.set_title('Dégradation du vert de malachite en fonction de la concentration des radicaux')
    ax.legend()
    canvas.draw()

    label_ka_value.config(text=f"{ka_fit:.4f}")
    label_kb_value.config(text=f"{kb_fit:.4f}")
    label_kc_value.config(text=f"{kc_fit:.4f}")

root = tk.Tk()
root.title("Analyse de dégradation du colorant vert de malachite")

# Importer le fichier Excel
button_import = tk.Button(root, text="Importer le fichier Excel", command=import_excel)
button_import.pack(padx=10, pady=10)

# Concentration initiale de VM
label_VM_0 = tk.Label(root, text="Concentration initiale de VM (VM_0):")
label_VM_0.pack()
entry_VM_0 = tk.Entry(root)
entry_VM_0.pack(padx=10, pady=10)

# Valeur de I
label_I = tk.Label(root, text="Intensité du courant (I):")
label_I.pack()
entry_I = tk.Entry(root)
entry_I.pack(padx=10, pady=10)

# Valeur de V
label_V = tk.Label(root, text="Volume de la solution (V):")
label_V.pack()
entry_V = tk.Entry(root)
entry_V.pack(padx=10, pady=10)

# Concentration initiale de H2O2
label_H2O2_0 = tk.Label(root, text="Concentration initiale de H2O2:")
label_H2O2_0.pack()
entry_H2O2_0 = tk.Entry(root)
entry_H2O2_0.pack(padx=10, pady=10)

# Valeur du pH
label_pH = tk.Label(root, text="Valeur du pH:")
label_pH.pack()
entry_pH = tk.Entry(root)
entry_pH.pack(padx=10, pady=10)

# Lancer l'analyse
button_analyze = tk.Button(root, text="Lancer l'analyse", command=analyze)
button_analyze.pack(padx=10, pady=10)

# Créer le graphe
fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Afficher les constantes cinétiques
label_ka = tk.Label(root, text="Constante cinétique (ka):")
label_ka.pack()
label_ka_value = tk.Label(root, text="")
label_ka_value.pack()

label_kb = tk.Label(root, text="Constante cinétique (kb):")
label_kb.pack()
label_kb_value = tk.Label(root, text="")
label_kb_value.pack()

label_kc = tk.Label(root, text="Constante cinétique (kc):")
label_kc.pack()
label_kc_value = tk.Label(root, text="")
label_kc_value.pack()

root.mainloop()


