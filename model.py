# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 11:09:02 2022

@author: yakupcatalkaya
"""

from neuron import h, gui
from neuron.units import ms, mV
import pickle
import matplotlib.pyplot as plt
# %matplotlib qt


class CreateStimuli:
    def __init__(self, stimuli_type="IClamp", neuron=None, stim_num=1, stim_start=9,
                 stim_delay=1, duration=10.0, amplitude=2.0, stimuli_w=0.04,
                 asyn_e=0, asyn_i=0, tau=0, onset=0, gmax=0, position=None):
        self._position = position
        if stimuli_type.lower()=="NetStim".lower():
            self._stimuli = h.NetStim(self._position)             
            self._stimuli.number = stim_num
            self._stimuli.start = stim_start
            self._nc = h.NetCon(self._stimuli, neuron)
            self._nc.delay = stim_delay
            self._nc.weight[0] = stimuli_w
            
        elif stimuli_type.lower()=="IClamp".lower():
            self._stimuli = h.IClamp(self._position)
            self._stimuli.delay = stim_delay
            self._stimuli.dur = duration
            self._stimuli.amp = amplitude
            
        elif stimuli_type.lower()=="AlphaSynapse".lower():
            self._stimuli = h.AlphaSynapse(self._position)
            # self._stimuli.alpha() = alpha
            self._stimuli.e = asyn_e
            self._stimuli.tau = tau
            self._stimuli.onset = onset
            self._stimuli.gmax = gmax
            self._stimuli.i = asyn_i

            

class ConnectSynapses:
    def __init__(self, neuron_one, neuron_two, neuron_one_synapse_pos=1.0,
                 neuron_two_synapse_pos=0.5, connection_delay=5.0, 
                 connection_weight=0.01, connection_threshold=10):

        self._ncstim = h.NetCon(neuron_one.axon(neuron_one_synapse_pos)._ref_v,
                                neuron_two.syn, sec=neuron_one.axon)            
        self._ncstim.delay = connection_delay
        self._ncstim.weight[0] = connection_weight
        neuron_one.ncl.append(self._ncstim)
        # self._ncstim.threshold = connection_threshold
        
    

class CreateNeuron:
    def __init__(self,neuron_type="N", Ra=100, cm=1, soma_L=12.6, soma_diam=12.6,
                 dend_L=200, dend_diam=1, nseg=11, axon_L=400, synapse_position=0.5, 
                 axon_diam=1, soma_pos=0.5, dend_pos=0.5, axon_pos=0.5):
        
        self._neuron_type = neuron_type
        self._Ra = Ra
        self._cm = cm
        self._soma_L = soma_L 
        self._soma_diam = soma_diam
        self._nseg = nseg
        self.ncl = []
        self.soma, self.dend, self.axon = None, None, None
        self._create_soma()
        self._create_dendrit(dend_L, dend_diam)
        self._create_axon(axon_L, axon_diam)
        self.all = self.soma.wholetree()
        self._create_synapse(section=self.dend(synapse_position))
        self._setup_biophysics()
        self._soma_v = h.Vector().record(self.soma(soma_pos)._ref_v)
        self._dend_v = h.Vector().record(self.dend(dend_pos)._ref_v)
        self._axon_v = h.Vector().record(self.axon(axon_pos)._ref_v)
        self._syn_i = h.Vector().record(self.syn._ref_i)
        self._spike_detector = h.NetCon(self.soma(soma_pos)._ref_v, None, sec=self.soma)
        self._spike_times = h.Vector()
        self._spike_detector.record(self._spike_times)
        
        
    def _create_soma(self):
        self.soma = h.Section(name="soma", cell=self)
        self.soma.L = self._soma_L
        self.soma.diam = self._soma_diam
        self.soma.nseg = self._nseg
        self.soma.cm = self._cm
        self.insert_hh(section=self.soma, gnabar=0.12, gkbar=0.036,
                       gl=0.0003, el=-54.3)


    def _create_dendrit(self, dend_L, dend_diam, name="dend", position=1.0):
        self.dend = h.Section(name=name, cell=self)
        self.dend.L = dend_L
        self.dend.nseg = self._nseg
        self.dend.diam = dend_diam
        self.dend.connect(self.soma(position))
        self.insert_pas(section=self.dend)
        self.insert_hh(section=self.dend, gnabar=0.012, gkbar=0.036,
                        gl=0.0003, el= -64)
        

    def _create_axon(self, axon_L, axon_diam, name="axon", position=0.0):
        self.axon = h.Section(name=name, cell=self)
        self.axon.L = axon_L
        self.axon.nseg = self._nseg
        self.axon.diam = axon_diam
        self.axon.connect(self.soma(position))
        self.insert_pas(section=self.axon)
        self.insert_hh(section=self.axon, gnabar=0.12, gkbar=0.036,
                        gl=0.0003, el=-54.3)
        
        
    def _create_synapse(self, section, synapse_position=1.0, tau=2, e=0):
        self.syn = h.ExpSyn(section)
        self.syn.tau = tau * ms
        self.syn.e = e
        
        
    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = self._Ra  
            sec.cm = self._cm  
            
        
    def insert_pas(self, section, g=0.001, e=-65):
        section.insert("pas") 
        for seg in section:  
            seg.pas.g = g         
            seg.pas.e = e  
    
    
    def insert_hh(self, section, gnabar=0.12, gkbar=0.036, gl=0.0003, el=-54.3):
        section.insert("hh")
        for seg in section:
            seg.hh.gnabar = gnabar  
            seg.hh.gkbar = gkbar  
            seg.hh.gl = gl 
            seg.hh.el = el
    
    
    def insert_extracellular(self, section, xraxial=0, xc=0, xg=0, e=0):
        section.insert("extracellular")
        for seg in section:
            seg.extracellular.xraxial = xraxial
            seg.extracellular.xc = xc
            seg.extracellular.xg = xg
            seg.extracellular.e = e
            
    
    def plot_recording_potential(self, t, soma_pos=0.5, dend_pos=0.5,
                                 axon_pos=0.5, label=""):
        soma_v, dend_v = self._soma_v.as_numpy(), self._dend_v.as_numpy() 
        axon_v, t = self._axon_v.as_numpy(), t.as_numpy()
        plt.figure()
        plt.plot(t, soma_v, color='blue', label='soma(' + str(soma_pos) + ')')
        plt.plot(t, dend_v, color='red', label='dend(' + str(dend_pos) + ')')
        plt.plot(t, axon_v, color='green', label='axon(' + str(axon_pos) + ')')
        plt.plot([t[0], t[-1]], [self.syn.e, self.syn.e], 
                          label='syn reversal',color='black', linestyle=':')
        plt.legend()
        plt.ylabel('potential (mV)')
        plt.xlabel('time (ms)')
        plt.title(label)
        plt.show()
        return soma_v
    
        
    def plot_recording_current(self, t, synapse_name=""):
        t, syn_i = t.as_numpy(), self._syn_i.as_numpy()
        label='Synaptic Current of ' + synapse_name
        plt.figure()
        plt.plot(t, syn_i, color='blue', label="Synapse current")
        plt.legend()
        plt.ylabel(h.units('ExpSyn.i'))
        plt.xlabel('time (ms)')
        plt.title(label)
        plt.show()        
        
        
    def inform(self):
        for sec in self.all:
            infos = sec.psection()
            for info in infos:
                print(info, ":==:", infos[info])
            print("\n")
    
    
    def __repr__(self):
        return "Neuron {} : ".format(self._neuron_type)


def pickle_it(files=[], names=None, instruction="wb"):
    if instruction == "wb":
        for idx, file in enumerate(files):
            pickle_out = open(names[idx] + ".pickle", instruction)
            pickle.dump(file,pickle_out)
            pickle_out.close()
    else:
        for name in names:
            pickle_out = open(name + ".pickle", instruction)
            file = pickle.load(pickle_out)
            files.append(file)
            pickle_out.close()
        return files    


def start_simulator(f_init=-65, duration=50):
    h.finitialize(f_init * mV)
    h.continuerun(duration)
    

def plot_spikes(neurons):
    plt.figure()
    colors = list({u'c':(0.0, 0.75, 0.75), u'b':(0.0, 0.0, 1.0),
              u'g':(0.0, 0.5, 0.0), u'y':(0.75, 0.75, 0), u'k':(0.0, 0.0, 0.0),
              u'r':(1.0, 0.0, 0.0), u'm':(0.75, 0, 0.75)}.values())
    for neuron_index, neuron_key in enumerate(neurons):
        neuron = neurons[neuron_key]
        plt.vlines(neuron._spike_times.as_numpy(), neuron_index + 0.5, 
                   neuron_index + 1.5, colors = colors[neuron_index], 
                   label=neuron_key)
    plt.xlabel("time (ms)")
    plt.ylabel("Neurons")
    plt.title("Spike Times of Neurons")
    plt.legend()
    plt.show()    
    
    
def main():
    # h.load_file("stdrun.hoc")
    my_neurons = {}
    
    
    for times in range(5):
        neuron = CreateNeuron(neuron_type = str(times))
        neuron._setup_biophysics()
        my_neurons["neuron" + str(times)] = neuron
    
    stimuli = CreateStimuli(stimuli_type="IClamp", 
                            position=my_neurons[list(my_neurons.keys())[0]].dend(0.5))
    
    for key_index, neuron_one_key in enumerate(list(my_neurons.keys())):
        if key_index < len(list(my_neurons.keys())) - 1:
            neuron_two_key = list(my_neurons.keys())[key_index + 1]
            neuron_two = my_neurons[neuron_two_key]
            neuron_one = my_neurons[neuron_one_key]
            ConnectSynapses(neuron_one = neuron_one, neuron_two = neuron_two)
            neuron_one.inform()
    
    h.topology()
    frequency = 40  # kHz
    h.dt = 1 / frequency
    t = h.Vector().record(h._ref_t)   
        
    start_simulator()
    
    for neuron_key in my_neurons:
        neuron = my_neurons[neuron_key]
        neuron.plot_recording_potential(t=t, label="Potential for " + neuron_key)
        neuron.plot_recording_current(t=t, synapse_name=neuron_key)
        
    plot_spikes(my_neurons)


if __name__ == "__main__":
    main()
