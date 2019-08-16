# NC-Network-of-LIF-Neurons-Coupled-by-Chemical-Synapses
This program generates a network of LIF neurons which are connected through chemical synapses. It solves LIF equation as follows:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\tau\dot{v}=-(v-v_{rest})+I_{ext}-I_{syn}" title="LIF Formula" />

where &tau; is the time constant of neuron, v is the membrane potential of the postsynaptic neuron, v<sub>rest</sub> is the
resting state membrane potential, I_ext is the external input, which is equivalent to a sensory stimulus or any
input sourced from the activity of population of neurons nearby, and I<sub>syn</sub> is accumulative synaptic inputs 
arrived from all the pre-synaptic neurons to the given post-synaptic one, which is calculated through

<img src="https://latex.codecogs.com/svg.latex?\Large&space;I_{syn,j}(t)=\sum_{i}g_{ij}*S_{syn,ij}(t)(v(t)-E_{syn,ij})}." />

here, i is associated index of pre synaptic neurons, v is the postsynaptic potential, and E<sub>syn</sub> is the synaptic 
reversal potential. Whether a synapse is 
excitatory or inhibitory is determined by E<sub>syn</sub>. For inhibitory synapses, E<sub>syn</sub> equals to 
E_syn_inh (here it is -80), and for excitatory 
ones the value is set to zero.
g is maximal conductance of the synapse, and S represents the fraction of bound receptors. Its kinetics are described by 
the following equation,

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\dot{S}(t)=\alpha&&N(1-S(t))&&-\beta&&S(t)&&,"/>

 where N is equivalent to the concentration of transmitters released into the synaptic cleft following the arrival of an action potential
 at the pre-synaptic terminal. Since the concentration of neurotransmitters in the cleft rises and falls very rapidly,
it is assumed that N occurs as a rectangular pulse. &alpha; and &beta; are the
rate constants for transmitter binding to and unbinding from post-synaptic receptor, respectively. &alpha; and &beta; and duration of N 
characterize the dynamics of post-synaptic potential (EPSP & IPSP). They directly specify the shape of S. 
[(Destexhe, 1993)](https://www.researchgate.net/profile/Terrence_Sejnowski/publication/220499797_An_Efficient_Method_for_Computing_Synaptic_Conductances_Based_on_a_Kinetic_Model_of_Receptor_Binding/links/54a4afee0cf267bdb90679ca/An-Efficient-Method-for-Computing-Synaptic-Conductances-Based-on-a-Kinetic-Model-of-Receptor-Binding.pdf)


This program returns 4 graphs. Figure(1) depicts the sub-threshold
dynamics of neurons' membrane potential.  Figure(2) represents the dynamics of N. It also shows the dynamics of 
S<sub>12</sub> and S<sub>21</sub>, which is meant to demonstrate the typical profile of S signal,
which in turn determine the input I<sub>syn</sub>. You can produce EPSP and IPSP
corresponding to NMDA, GABA<sub>A</sub>, GABA<sub>B</sub>, etc. by setting relevant values for 
&alpha; and &beta; and therefore for the shape the S profile.
Figure(3) demonstrates the accumulative synaptic inputs from pre-synaptic 
neurons arrived to the post-synaptic ones. Note that negative I<sub>syn</sub> corresponds 
to EPSP and positive I<sub>syn</sub> corresponds to IPSP (Consider the negative 
sign before I<sub>syn</sub> in the LIF equation). Figure(4) is raster plot of the neuronal activity in the network.
