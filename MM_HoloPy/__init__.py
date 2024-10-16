""" Multi-Map Holographic Analyais in Python (MM_HoloPy))

***Introduction:
MM_HoloPy is a Python package that is designed to study the developed Multi-Map Holographic
Metrology for the Fred Young Sub-mm Telescope which is a coma-corrected Crossed-Dragone type
antenna and consists of two large reflectors. This data analysis method is also available for 
holography analysis of other Dual-reflector antennas, such as Gregory type antenna. 
Conventional one-beam holographic experiment also can use this software for reflector surface analysis. 

The package offers the way to simulate a practical holographic measurement system, such as
random noise of the cross-correlating receiver, defocus errors in the optics, atmospheric fluct-
ations, telescope scanning errors, and inaccurancy of reference receiver positions. 


The package consists of two main parts:
1. EM simulations for FYST, Two-step Kirchoff-Fresnel intergration, Two-step Physical Optics (or A-PO).
   a. An efficient and accurate physical optics analysis method, 'Two-step' Kirchhoff-Fresnel
   intergration method, is implemented in this package to predict the diffraction beam of the FYST
   telescope in far and near field. The method is based on scalar diffraction theory, therefore, the
   cx-polarization beam pattern of the telescope cannot be computed.

   b. The package offers the way of using full physical optics anlaysis. And the same idea of breaking PO 
   integration into two steps are also used, if the optical component has two focal plane or two beamwaist
   plane (Quasi-Optics)

   C. Hope GO, GTD and PTD analysis methods will be included in future.

   For holographic analysis, only the scalar diffraction techniques are used to produce the observed beams.

2. Inference data anlaysis software for Multi-map holography
   a. package for different random noise and systemtic noise simulations.

"""

import sys

__minimum_python_version__ = "3.7"
