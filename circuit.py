# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 18:59:16 2022

@author: yakupcatalkaya
"""

import math



class Item():
    """
    All electrical properties of any material that can be seen in the model
    such as potential(v), current(i), conductance(g), resistance(r), and its 
    name to be able to identify later.
    
    """
    def __init__(self, resistance=None, name=None):
        
        self.resistance = resistance        
        self.name = name
        self.input_v = 0
        self.potential = 1

    
    def update_values(self, value_names, new_values):
        for value_name, new_value in enumerate(value_names, new_values):
            setattr(self, value_name, new_value)
        


class Resistor(Item):
    """
    g (conductance = constant/R) of a voltage-gated ion channels, axial resistance,
    and membrane resistance
    
    """
    def __init__(self, name, m=0, n=0, h=0, memb_pot=-65, current=1e-12, leak=0,
                 start_point=None):
        Item.__init__(self)
        self.name = name.lower().strip()
        self.m = m 
        self.n = n
        self.h = h 
        self.memb_pot = memb_pot
        self.current = current
        self.leak_potential = leak
        self.input_v = start_point

    def update_m_n_h_g(self, v_intracellular, v_extracellular):
        self.memb_pot = (v_intracellular.value - v_extracellular.value)
        
        self.n_alpha = (0.02*(self.memb_pot-25))/(1-math.exp(-(self.memb_pot-25)/9))
        self.m_alpha = (0.182*(self.memb_pot+35))/(1-math.exp(-(self.memb_pot+35)/9))
        self.h_alpha = 1/(1+math.exp(-(self.memb_pot-62)/6))
        
        self.n_beta = (-0.002*(self.memb_pot-25))/(1-math.exp((self.memb_pot-25)/9))
        self.m_beta = (-0.124*(self.memb_pot+35))/(1-math.exp((self.memb_pot+35)/9))
        self.h_beta = (4*math.exp((self.memb_pot+90)/12))/(1+math.exp(-(self.memb_pot-62)/6))
        
        self.dm = self.m_alpha*(1-self.m)-self.m_beta*self.m
        self.dn = self.n_alpha*(1-self.n)-self.n_beta*self.n
        self.dh = self.h_alpha*(1-self.h)-self.h_beta*self.h
    
        self.m += self.dm
        self.n += self.dn
        self.h += self.dh
        
        if self.name == "g_na":
            self.conductance = self.current / (self.potential * self.h * (self.m ** 3))
            self.resistance = 1 / (self.conductance * self.h * (self.m ** 3))
            self.ina = self.potential / self.resistance
            
        elif self.name == "g_k":
            self.conductance = self.current / (self.potential * (self.n ** 4))
            self.resistance = 1 / (self.conductance * (self.n ** 4))
            self.ik = self.potential / self.resistance
            
        elif self.name == "g_cl":
            self.conductance = self.current / self.potential
            self.resistance = 1 / self.conductance
            self.icl = self.potential / self.resistance
        else:
            raise SystemError("Look at here.")
            self.conductance = 1 / self.resistance
            self.i = self.potential / self.resistance
            
            

class Capacitor(Item):
    """
    
    membrance capacitance that is caused by lipid layer that is simly an leaky 
    insulator which is typically 1uF/cm2
    
    """
    def __init__(self, init, dc):
        pass



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
    def __init__(self, init, di):
        pass
    


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

    
v_extracellular = VoltagePoint(0, "v_extracellular")
v_intracellular = VoltagePoint(-65, "v_intracellular")
potential_k = VoltagePoint(-75, "potential_k")
potential_na = VoltagePoint(+55, "potential_na")
potential_cl = VoltagePoint(-69, "potential_cl")

g_k = Resistor(name="g_k", start_point=v_extracellular, leak = potential_k.value)
g_na = Resistor(name="g_na", start_point=v_extracellular, leak = potential_na.value)
g_cl = Resistor(name="g_cl", start_point=v_extracellular, leak = potential_cl.value)
   
e_k = VoltageSource(name="e_k", init_voltage=-75)
e_na = VoltageSource(name="e_na", init_voltage=+55)
e_cl = VoltageSource(name="e_cl", init_voltage=-69)

for repeat in range(50):
    for item in [g_k, g_na, g_cl]: item.update_m_n_h_g(v_intracellular, v_extracellular)
    v_intracellular.update_v((e_k.output_v + e_na.output_v + e_cl.output_v) / 3)
    v_intracellular.inform()
    v_extracellular.inform()




# def main():
#     return True

# if __name__ == "__main__":
#     main()
