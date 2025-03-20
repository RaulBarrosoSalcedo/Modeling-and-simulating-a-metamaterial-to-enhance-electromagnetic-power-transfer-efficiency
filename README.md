## Modeling-and-simulating-a-metamaterial-to-enhance-electromagnetic-power-transfer-efficiency
This is the repository of the codes related to the IEEE Latam Journal Publication entitled " Modeling and simulating a metamaterial to enhance electromagnetic power transfer efficiency"

# Analitycally calculated metamaterial parameters and format of the paper Figures.

The file "Codes of Figures" is the Latex code that allowed to authors to construct the annalitical calculation of the metamaterial parameters which were presented in the paper in the Figs. 7, 8, 9, 11 and 13. This code performed in the "tikz" latex environment includes each one of the utilized and programmed annalytical formulation that allows the calculus of each parameter, being this step properly commented in the code.

# Electric permitivitty estimated by means of a numerical Finite Elements estimation

The file named "Epsilonrcalculation.aedt" in HFSS format can be utilized to estimate the value of epsilon in presence of lossless materials in the squeme of the Fig. 6.a. The existent capacitance Cm between the PEC surfaces must be calculated in presence of the conductor cube of edge "b". Then, the existent capacitance C0 between the PEC surfaces must be calculated in absence of the conductor cube. Finally, the relative epsilon er is estimated as the ratio er = Cm/C0. Hence, Fig. 7 can be constructed variating the ratio b/a for each calculation.

# Magnetic permeabilitty estimated by menas of a numerical Finite Elements estimation

The file named "Murcalculation.aedt" in HFSS format can be utilized to estimate the value of mu in presence of lossless materials in a squeme that is equivalent of the Fig. 6.b. For the simulations, such squeme was performed with an electric field instead a magentic field, but with a perfect magnetic cube (PMC) conductor in the cell, side PMC surfaces and upper and lower PEC surfaces that yields border conditions that restrict the E behaviour that are equievalent to a electric perfect conductor cube produces for H in Fig.6b. Then, in the simualtions, the existent capacitance Cm´ between the PEC surfaces must be calculated in presence of the PMC cube of edge "b". Then, the existent capacitance C0´ between the PEC surfaces must be calculated in absence of the PMC conductor cube. Finally, the relative mu given by mur is estimated as the ratio mur = Cm´/C0´. Hence, Fig. 8 can be constructed variating the ratio b/a for each calculation.

# Attenuation constant estimated by means of a numerical FDTD estimation

The file "alfa.m" in Matlab format can be utilized to estimate the attenueation constant of the Fig. 10 scheme, using a typical FDTD transient excitation. The user must be able to extract that parameter from the efield-standing-wave-ratio pattern which is calculated with this code. Follow the instructions included in the comments of the code. The final results must be the ones presented in Fig. 11.

# Power efficiency of the proposed system estimated by means of a numerical FDTD estimation

The file "MatchingSystem.m" in Matlab format can be utilized to estimate the fields of the Fig. 12 scheme, using a typical FDTD transient excitation. After that, it is necessary to run the Efficiency_vs_f.m file too, in the same folder if it is desired to reproduce the power efficiency results presented in Fig. 13.



