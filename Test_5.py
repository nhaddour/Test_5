# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 07:51:32 2023

@author: naouf
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import tkinter as tk
from tkinter import *
from tkinter import filedialog, messagebox
import matplotlib
matplotlib.use('tkagg')
import pandas as pd



def fenton_reaction(t, y, I, V, F, ki):
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


def run_simulation():
    I = float(entry_I.get())
    V = float(entry_V.get())
    H2O2 = float(entry_H2O2.get())
    pH = float(entry_pH.get())

    # Constants
    F = 9.64853e4 # Faraday Constant in C.mol-1

    # System parameters
    pH = pH
    I = I
    V = V

    # Kinetic constants
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

    # Initial concentrations
    y0 = np.array([
        0,                  # Fe2p_0
        2e-3,               # H2O2_0
        0,                  # Fe3p_0
        10**(pH-14),        # OHm_0
        10**(-pH),          # Hp_0
        0,                  # OHpt_0
        0,                  # HO2pt_0
        0                   # O2ptm_0
    ])


    # Calculer les concentrations initiales et autres paramètres
    # ...

    t_start = 0
    t_end = 1800
    num_points = 10000

    # Résoudre le système d'équations différentielles
    sol = solve_ivp(
        lambda t, y: fenton_reaction(t, y, I, V, F, ki),
        (t_start, t_end),
        y0,
        method="Radau",
        t_eval=np.linspace(t_start, t_end, num_points),
    )

 
class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.master.title("Importer fichier Excel et Fitting")

        self.lbl_a = tk.Label(self.master, text="a:")
        self.lbl_a.grid(row=0, column=0)
        self.entry_a = tk.Entry(self.master, state="readonly")
        self.entry_a.grid(row=0, column=1)

        self.lbl_b = tk.Label(self.master, text="b:")
        self.lbl_b.grid(row=1, column=0)
        self.entry_b = tk.Entry(self.master, state="readonly")
        self.entry_b.grid(row=1, column=1)

        self.lbl_c = tk.Label(self.master, text="c:")
        self.lbl_c.grid(row=2, column=0)
        self.entry_c = tk.Entry(self.master, state="readonly")
        self.entry_c.grid(row=2, column=1)

        self.canvas = tk.Canvas(self.master, width=1000, height=1000)
        self.canvas.grid(row=3, columnspan=2)

        self.btn_import = tk.Button(self.master, text="Importer fichier Excel", command=self.import_excel)
        self.btn_import.grid(row=4, columnspan=2)

    def import_excel(self):
        file_path = filedialog.askopenfilename()
        df = pd.read_excel(file_path)

        t = df['temps']
        Dg = df['degradation']

        try:
            popt, pcov = curve_fit(self.model, t, Dg)
        except Exception as e:
            messagebox.showerror("Erreur de Fitting", str(e))
            return

        a, b, c = popt
        a_err, b_err, c_err = np.sqrt(np.diag(pcov))

        self.entry_a.config(state="normal")
        self.entry_a.delete(0, tk.END)
        self.entry_a.insert(0, f"{a:.3f} +/- {a_err:.3f}")
        self.entry_a.config(state="readonly")

        self.entry_b.config(state="normal")
        self.entry_b.delete(0, tk.END)
        self.entry_b.insert(0, f"{b:.3f} +/- {b_err:.3f}")
        self.entry_b.config(state="readonly")

        self.entry_c.config(state="normal")
        self.entry_c.delete(0, tk.END)
        self.entry_c.insert(0, f"{c:.3f} +/- {c_err:.3f}")
        self.entry_c.config(state="readonly")
        
        self.plot_graph(t, Dg, popt, pcov)

    def model(self, y, a, b, c):
        return (a * y[5]) + (b * y[6])+ (c * y[7])

    def plot_graph(self, t, Dg, popt, pcov):
       self.canvas.delete("all")
    
       plt.figure(figsize=(10,10))
       plt.plot(t, Dg, 'bo', label="Données")
    
       t_fit = np.linspace(t.min(), t.max(), 100)
       try:
        Dg_fit = self.model(t_fit, *popt)
        Dg_upper = self.model(t_fit, *(popt + np.sqrt(np.diag(pcov))))
        Dg_lower = self.model(t_fit, *(popt - np.sqrt(np.diag(pcov))))
       except Exception as e:
            messagebox.showerror("Erreur de génération de courbe", str(e))
            return

       plt.plot(t_fit, Dg_fit, 'r-', label="Modèle")
    
       if isinstance(Dg_lower, np.ndarray):
            t_stack = np.hstack((t_fit, t_fit[::-1]))
            Dg_stack = np.hstack((Dg_upper, Dg_lower[::-1]))
            plt.fill(t_stack, Dg_stack, facecolor='gray', alpha=0.2, label="Incertitude")
    
       plt.xlabel('temps')
       plt.ylabel('degradation')
       plt.title('Graphe Dg=f(t)')
       plt.legend()
    
       plt.savefig('graph.png')

       self.photo = tk.PhotoImage(file='graph.png')
       self.canvas.create_image(0, 0, anchor='nw', image=self.photo)
       
root = tk.Tk()
root.title("Importer fichier Excel et Fitting")
app = Application(master=root)

tk.Label(root, text="I (mA)").grid(row=0, column=0)
entry_I = Entry(root)
entry_I.grid(row=0, column=1)

tk.Label(root, text="V (mL)").grid(row=0, column=2)
entry_V = Entry(root)
entry_V.grid(row=0, column=3)

tk.Label(root, text="H2O2 (M)").grid(row=0, column=4)
entry_H2O2 = Entry(root)
entry_H2O2.grid(row=0, column=5)

tk.Label(root, text="pH").grid(row=0, column=6)
entry_pH = Entry(root)
entry_pH.grid(row=0, column=7)

tk.Label(root, text="a:").grid(row=1, column=0)
entry_a = tk.Entry(root, state="readonly")
entry_a.grid(row=1, column=1)

tk.Label(root, text="b:").grid(row=1, column=2)
entry_b = tk.Entry(root, state="readonly")
entry_b.grid(row=1, column=3)

tk.Label(root, text="c:").grid(row=1, column=4)
entry_c = tk.Entry(root, state="readonly")
entry_c.grid(row=1, column=5)

tk.Button(root, text="Importer fichier Excel", command=app.import_excel).grid(row=2, columnspan=2)
tk.Button(root, text="Run Simulation", command=run_simulation).grid(row=3, columnspan=2)

canvas = tk.Canvas(root, width=1000, height=1000)
canvas.grid(row=10, columnspan=2)

root.mainloop()


