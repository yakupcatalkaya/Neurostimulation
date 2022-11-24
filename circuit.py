# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 18:59:16 2022

@author: yakupcatalkaya
"""

from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np


class Item():
    """
    
    All electrical properties of any material that can be seen in the model
    such as potential(v), current(i), conductance(g), resistance(r), and its 
    name to be able to identify later.
    
    """
    def __init__(self, resistance=None, name=None):
        global freq, time
        
        self.resistance = resistance        
        self.name = name
        self.input_v = 0
        self.potential = 1
        self.period = 1 / freq
        self.time = time
        


class Resistor(Item):
    """
    
    g (conductance = constant/R) of a voltage-gated ion channels, axial resistance,
    and membrane resistance
    
    """
    def __init__(self, name, memb_pot=-65, current=0, leak=0,
                 start_point=None, conductance=None, time=None):
        Item.__init__(self)
        self.name = name.lower().strip()
        self.memb_pot = memb_pot
        self.current = current
        self.conductance = conductance
        self.const_cond = conductance
        self.leak_potential = leak
        self.time = time
        self.m = Vector(self.name + "_m")
        self.h = Vector(self.name + "_h")
        self.n = Vector(self.name + "_n")
        self.cond_record = Vector(self.name + "_conductance")
        self.i_record = Vector(self.name + "_current")


    def update_variables(self):
        for index, item in enumerate(self.time):
            if self.name == "g_na":
                self.cond_record.record(self.const_cond *self.h.val[index] * self.m.val[index] ** 3)
            elif self.name == "g_k":
                self.cond_record.record(self.const_cond * self.n.val[index] ** 4)
            elif self.name == "g_cl":
                self.cond_record.record(self.const_cond)
            else:
                raise SystemError("Look at here.")
            # self.i_record.val = (self.cond_record.val * self.memb_pot)


    def inform(self):
        print("\n", self.name, "\n", "-" * 20)
        print(" h:", self.h, "\n", "m:", self.m, "\n", "n:", self.n)
        print(" g:", self.conductance, "\n", "i:", self.current, "\n")
        print(" R:", self.resistance)
        print(" Leak:", self.leak_potential, "\n", "Output:", self.output_v, "\n")
            
            

class Capacitor(Item):
    """
    
    membrance capacitance that is caused by lipid layer that is simly an leaky 
    insulator which is typically 1uF/cm2
    
    """
    def __init__(self, init, name):
        self.value = init
        self.name = name.lower().strip()
        
        
    def update(self, new_value):
        self.value = new_value



class VoltageSource(Item):
    """
    
    Ek, Ena ,Ecl values are resting membrane potentials for Potassium (K+1), 
    Sodium(Na+1), Chlor (Cl-1) ions respectively which can be found by 
    Ex = [outward movement of X] - [inward movement of X] 
    Ek = -75 mV       Ena = +55 mV       Ecl = -69 mV        
    
    """
    def __init__(self, name, init_voltage=0):
        Item.__init__(self)
        self.name = name.lower().strip()
        self.potential = init_voltage
        self.val = init_voltage



class CurrentSource(Item):
    """
    
    Current source will be inserted with glass micropippet so that ions can 
    flow from extracellular fluid into axon's cytoplasmic fluid, especially Na+
    ion will pass so that membrane poteintal can be increased to reach a point
    that will start action potential.
    
    """
    def __init__(self, initial, time, name, value):
        self.time = time
        self.initial = initial
        self.value = value
        self.name = name.lower().strip()
        self.values = np.zeros((len(self.time),))
        for idx, item in enumerate(self.values): 
            self.values[idx] = self.get_value(self.time[idx])

            
    def get_value(self, t):
        if 10.0 < t < 10.5:
            return self.value
        return self.initial


    def update(self, at_t=0, new=None):
        if new==None:
            self.value = self.values[at_t]
        else:
            self.value = new
        return self.value



class VoltagePoint(Item):
    """
    
    Voltage point object to be able change voltage value. 
    
    """
    def __init__(self, value, name):
        self.value = value
        self.name = name.lower().strip()
    
    
    def update_v(self, new_value):
        self.value = new_value
    
    
    def inform(self):
        print(self.name + " value is " + str(self.value))



class Vector():
    """
    
    Vector for recording values.
    
    """    
    def __init__(self, name):
        global freq, start_time, stop_time, time
        self.name = name.lower().strip()
        self.update_index = 0
        self.time = time
        self.val = np.zeros((len(self.time),), dtype=object)

    def record(self, new_value):
        self.val[self.update_index] = new_value
        self.update_index += 1
        
        
    def get_values(self):
        return self.val
    
    
    def plotter(self, labelx="", labely="", title="", label="", limx=None):
        plt.figure()
        plt.plot(self.time[200:], self.val[200:], color='blue', label=label)
        plt.ylabel(labely)
        plt.xlabel(labelx)
        plt.title(title)
        plt.xlim(limx)
        plt.savefig(title)
        plt.show()
    

def initial(u=0.0):
    m_alpha = (0.1*(u-25))/(1-np.exp(-(u-25)/10))
    m_beta = 4*np.exp(-u/18)

    h_alpha = 0.07*np.exp(-u/20)
    h_beta = 1/(1+np.exp(-(u-30)/10))
    
    n_alpha = (0.01*(u-10))/(1-np.exp(-(u-10)/10))
    n_beta = 0.125*np.exp(-u/80)
    
    n_inf = n_alpha / (n_alpha + n_beta)
    h_inf = h_alpha / (h_alpha + h_beta)
    m_inf = m_alpha / (m_alpha + m_beta)
    
    return np.array([u, n_inf, h_inf, m_inf])


def ode_solver(y, t):
    u, n, m, h = y[0], y[1], y[2], y[3]
    dy = np.zeros((4,))
    
    gk = (g_k.const_cond / memb_capacitance.value) * np.power(n, 4.0)
    gna = (g_na.const_cond / memb_capacitance.value) * np.power(m, 3.0) * h
    gcl = (g_cl.const_cond / memb_capacitance.value)
    
    m_alpha = (0.1 * (25.0 - u)) / (np.exp(2.5 - (0.1 * u)) - 1.0)
    m_beta = 4.0 * np.exp(-u / 18.0)
    h_alpha = 0.07 * np.exp(-u / 20.0)
    h_beta = 1.0 / (np.exp(3.0 - (0.1 * u)) + 1.0)
    n_alpha = (0.01 * (10.0 - u)) / (np.exp(1.0 - (0.1 * u)) - 1.0)
    n_beta = 0.125 * np.exp(-u / 80.0)
    
    dv = stimuli.get_value(t) / memb_capacitance.value
    dv -= gk * (u-e_k.val) + gna * (u-e_na.val) + gcl * (u-e_cl.val)
    dm = m_alpha * (1 - m) - m_beta * m
    dn = n_alpha * (1 - n) - n_beta * n
    dh = h_alpha * (1 - h) - h_beta * h

    dy[0], dy[1], dy[2], dy[3]= dv, dn, dm, dh

    return dy


freq = int(10e3)
start_time = 0.0
v_intracellular_rest = -65
stop_time = 20.0
np.random.seed(1000)
time = np.linspace(start_time, stop_time, freq)

v_extracellular = VoltagePoint(0, "v_extracellular")
v_intracellular = VoltagePoint(v_intracellular_rest, "v_intracellular")
v_memb_potential = VoltagePoint(v_intracellular_rest,"v_memb")

e_k = VoltageSource("e_k", +12)
e_na = VoltageSource("e_na", -115)
e_cl = VoltageSource("e_cl", -10.613)

memb_capacitance = Capacitor(1.0, "Membrance_Capacitance")
stimuli = CurrentSource(0, time, "Stimuli_current", 200)

g_k = Resistor(name="g_k", start_point=v_extracellular, 
                leak=e_k.potential, conductance=36.0, time=time)
g_na = Resistor(name="g_na", start_point=v_extracellular, 
                leak=e_na.potential, conductance=120.0, time=time)
g_cl = Resistor(name="g_cl", start_point=v_extracellular, 
                leak=e_cl.potential, conductance=0.3, time=time)
   
v_memb_record = Vector("v_membrane")
            
solution = odeint(ode_solver, initial(), time).T
new_v_memb, m, n, h = solution[0], solution[1], solution[2], solution[3]
v_memb_record.val = new_v_memb + v_intracellular_rest

v_memb_record.plotter(labelx="time(ms)", labely="Membrance potential (mV)", 
                       title="Action Potential Graph")

m_record = Vector("m")
h_record = Vector("h")
n_record = Vector("n")
m_record.val, h_record.val, n_record.val = m, h, n

for item in [g_k, g_na, g_cl]:
    item.m.val, item.n.val, item.h.val= m, n, h

stimuli_record = Vector("stimuli")
stimuli_record.val = stimuli.values
stimuli_record.plotter(labelx="time(ms)", labely="Current (uA/cm2)", title="Stimuli graph")

g_k.update_variables()
g_na.update_variables()
g_cl.update_variables()

g_k.cond_record.plotter(labelx="time(ms)", labely="Conductance(uS)", title="g_k graph")
g_na.cond_record.plotter(labelx="time(ms)", labely="Conductance(uS)", title="g_na graph")
g_cl.cond_record.plotter(labelx="time(ms)", labely="Conductance(uS)", title="g_cl graph")

m_record.plotter(labelx="time(ms)", labely="m", title="m values")
h_record.plotter(labelx="time(ms)", labely="h", title="h values")
n_record.plotter(labelx="time(ms)", labely="n", title="n values")
