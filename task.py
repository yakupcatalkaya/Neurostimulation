# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 09:37:13 2022

@author: yakupcatalkaya
"""

import model as m

my_neurons = {}
membrane = m.CreateNeuron()
my_neurons["membrane"] = membrane

stimuli = m.CreateStimuli(stimuli_type="IClamp", stim_delay=1, duration=10.0, amplitude=2.0,
                          position=my_neurons[list(my_neurons.keys())[0]].axon(0.0))

m.h.topology()
frequency = 40  # kHz
m.h.dt = 1 / frequency
t = m.h.Vector().record(m.h._ref_t) 


m.start_simulator(f_init=-65, duration=25)
my_neurons[list(my_neurons.keys())[0]].inform()

membrane.plot_recording_potential(t=t, label="Potential for " + list(my_neurons.keys())[0],
                                  axon_pos=0.0)
membrane.plot_recording_potential(t=t, label="Potential for " + list(my_neurons.keys())[0],
                                  axon_pos=0.5)
membrane.plot_recording_potential(t=t, label="Potential for " + list(my_neurons.keys())[0],
                                  axon_pos=1.0)
m.plot_spikes(my_neurons)
m.plotter(t.as_numpy(), stimuli._stimuli_current.as_numpy(), labelx="time (ms)", 
          labely="Current (nA)", title="IClamp Current Plot", label="IClamp")

