# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 18:59:16 2022

@author: yakupcatalkaya
"""

# import math
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
        global freq
        
        self.resistance = resistance        
        self.name = name
        self.input_v = 0
        self.potential = 1
        self.period = 1 / freq

    
    def update_values(self, value_names, new_values):
        for value_name, new_value in value_names, new_values:
            setattr(self, value_name, new_value)
        


class Resistor(Item):
    """
    
    g (conductance = constant/R) of a voltage-gated ion channels, axial resistance,
    and membrane resistance
    
    """
    def __init__(self, name, m=0, n=0, h=1, memb_pot=-65, current=0, leak=0,
                 start_point=None, conductance=None):
        Item.__init__(self)
        self.name = name.lower().strip()
        self.m_init, self.n_init, self.h_init = m, n, h
        self.m, self.n, self.h = m, n, h
        self.memb_pot = memb_pot
        self.u = memb_pot
        self.current = current
        self.conductance = conductance
        self.constant_conductance = conductance
        self.leak_potential = leak
        self.input_v = start_point.value
        self.output_v = self.input_v + self.potential + self.leak_potential
        
        
    def ode_for_mhn(self, y0, t, mhn):
        global v_intracellular_rest
        
        self.u = self.memb_pot - v_intracellular_rest
        
        if mhn == "m":
            self.m_alpha = (0.1*(self.u-25))/(1-np.exp(-(self.u-25)/10))
            self.m_beta = 4*np.exp(-self.u/18)
            self.dm = self.m_alpha*(1-self.m)-self.m_beta*self.m
            self.m_tau = 1 / (self.m_alpha + self.m_beta)
            self.m_zero = y0
            return self.dm

        elif mhn == "n":
            self.n_alpha = (0.01*(self.u-10))/(1-np.exp(-(self.u-10)/10))
            self.n_beta = 0.125*np.exp(-self.u/80)
            self.dn = self.n_alpha*(1-self.n)-self.n_beta*self.n
            self.n_tau = 1 / (self.n_alpha + self.n_beta)
            self.n_zero = y0
            return self.dn
            
        elif mhn == "h":
            self.h_alpha = 0.07*np.exp(-self.u/20)
            self.h_beta = 1/(1+np.exp(-(self.u-30)/10))
            self.dh = self.h_alpha*(1-self.h)-self.h_beta*self.h
            self.h_tau = 1 / (self.h_alpha + self.h_beta)
            self.h_zero = y0
            return self.dh
        

    def update_variables(self, v_intracellular, v_extracellular, at_t,
                         n_zero=0.01, m_zero=0.001, h_zero=0.9):
        global time
        
        self.memb_pot = (v_intracellular.value - v_extracellular.value)
        self.m_ode = odeint(self.ode_for_mhn, y0=m_zero, t=time, args=("m",))
        self.h_ode = odeint(self.ode_for_mhn, y0=h_zero, t=time, args=("h",))
        self.n_ode = odeint(self.ode_for_mhn, y0=n_zero, t=time, args=("n",))
        
        self.m_inf = self.m_ode[-1]
        self.h_inf = self.h_ode[-1]
        self.n_inf = self.n_ode[-1]
        
        self.m = self.m_zero-((self.m_zero-self.m_inf)*(1-np.exp(-at_t/self.m_tau)))
        self.h = self.h_zero-((self.h_zero-self.h_inf)*(1-np.exp(-at_t/self.h_tau)))
        self.n = self.n_zero-((self.n_zero-self.n_inf)*(1-np.exp(-at_t/self.n_tau)))
        
        if self.name == "g_na":
            self.conductance = self.constant_conductance *self.h * self.m ** 3
            self.current = self.conductance * (self.u - self.leak_potential)
            
        elif self.name == "g_k":
            self.conductance = self.constant_conductance * self.n ** 4
            self.current = self.conductance * (self.u - self.leak_potential)
            
        elif self.name == "g_cl":
            self.conductance = self.constant_conductance
            self.current = self.conductance * (self.u - self.leak_potential)
            
        else:
            raise SystemError("Look at here.")
        
    
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



class CurrentSource(Item):
    """
    
    Current source will be inserted with glass micropippet so that ions can 
    flow from extracellular fluid into axon's cytoplasmic fluid, especially Na+
    ion will pass so that membrane poteintal can be increased to reach a point
    that will start action potential.
    
    """
    def __init__(self, init, name):
        self.value = init
        self.name = name.lower().strip()
    
    
    def update(self, new_value):
        self.value = new_value



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
        self.val = np.zeros(int((stop_time - start_time) * freq))
        self.update_index = 0
        self.time = time


    def record(self, new_value):
        self.val[self.update_index] = new_value
        self.update_index += 1
        
        
    def get_values(self):
        return self.val
    
    
    def plotter(self, labelx="", labely="", title="", label="", limx=None):
        plt.figure()
        plt.plot(self.time, self.val, color='blue', label=label)
        plt.ylabel(labely)
        plt.xlabel(labelx)
        plt.title(title)
        plt.xlim(limx)
        plt.show()
    


def ode_for_vm(y0, t, v_in_nominator, v_in_denominator, stimuli_i, memb_capacitance):
    dvdt = ((v_in_nominator / v_in_denominator) + stimuli_i) / memb_capacitance
    print(dvdt)
    return dvdt


freq = 40e3
start_time = 0
period = 1 / freq
repeat = 1500
v_intracellular_rest = -65
stop_time = repeat * period
time = np.arange(start_time, stop_time, period, dtype=np.float32)

v_extracellular = VoltagePoint(0, "v_extracellular")
v_intracellular = VoltagePoint(v_intracellular_rest, "v_intracellular")

e_k = VoltageSource("e_k", -75)
e_na = VoltageSource("e_na", +55)
e_cl = VoltageSource("e_cl", -69)

memb_capacitance = Capacitor(1, "Membrance_Capacitance")

stimuli = CurrentSource(0, "Stimuli_current")

g_k = Resistor(name="g_k", start_point=v_extracellular, 
                leak=e_k.potential, conductance=35)
g_na = Resistor(name="g_na", start_point=v_extracellular, 
                leak=e_na.potential, conductance=40)
g_cl = Resistor(name="g_cl", start_point=v_extracellular, 
                leak=e_cl.potential, conductance=0.3)
   
v_intra_record = Vector("v_intra")
gk_record = Vector("gk_conductance") 
gna_record = Vector("gna_conductance")
gcl_record = Vector("gcl_conductance")
gk_m_record = Vector("gk_m")
gk_n_record = Vector("gk_n")
gk_h_record = Vector("gk_h")

for rep in range(repeat):
    v_in_nominator = 0
    v_in_denominator = 0
    if rep == 1000: stimuli.update(1)
    if rep == 1200: stimuli.update(0)
    
    for item in [g_k, g_na, g_cl]: 
        item.update_variables(v_intracellular, v_extracellular, at_t=rep, 
                              n_zero=item.n, m_zero=item.m, h_zero=item.h)
        v_in_nominator += item.conductance * item.leak_potential
        v_in_denominator += item.conductance
        
    # ode_for_vm(v_in_nominator, v_in_denominator, stimuli.value, memb_capacitance.value)
    
    new_v_intracellulars = odeint(ode_for_vm, v_intracellular.value, time, 
                                 args=(v_in_nominator, v_in_denominator, 
                                       stimuli.value, memb_capacitance.value))
    new_v_intracellular = new_v_intracellulars[0]
    v_intracellular.update_v(new_v_intracellular)
    v_intra_record.record(v_intracellular.value)
    gk_record.record(g_k.conductance)
    gna_record.record(g_na.conductance)
    gcl_record.record(g_cl.conductance)
    gk_m_record.record(g_k.m)
    gk_n_record.record(g_k.n)
    gk_h_record.record(g_k.h)


v_intra_record.plotter(labelx="time(s)", labely="Membrance potentiak(mV)", 
                       title="Action Potential Graph", limx=(990 * period, 1100 * period))

gk_record.plotter(labelx="time(s)", labely="Conductance(uS)", title="g_k graph")
gna_record.plotter(labelx="time(s)", labely="Conductance(uS)", title="g_na graph")
gcl_record.plotter(labelx="time(s)", labely="Conductance(uS)", title="g_cl graph")

gk_m_record.plotter(labelx="time(s)", labely="m", title="m values")
gk_n_record.plotter(labelx="time(s)", labely="n", title="n values")
gk_h_record.plotter(labelx="time(s)", labely="h", title="h values")

# plt.plot(v_intra_record.val[990:1010])
