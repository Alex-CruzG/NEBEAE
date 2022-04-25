# NEBEAE
 Matlab implementation of Nonlinear Extended Blind End-member and Abundance Extraction (NEBEAE) algorithm. <br><br>
NEBEAE implements a blind hyperspectral unmixing based on the multilinear mixing model [(MMM)](10.1109/TGRS.2015.2453915).<br>

In the problem formulation, we include a normalization step in the hyperspectral measurements for the end-members and abundances to improve robustness. 
The blind unmixing process can be separated into three estimation subproblems for each component in the model, which are solved by a cyclic coordinate descent algorithm and quadratic constrained optimizations. Each problem is mathematically formulated and derived to construct the overall nonlinear iterative unmixing technique. <br>

We evaluated our proposal with synthetic and experimental datasets from the remote sensing literature (Cuprite and Urban datasets). <br>

The performance is compared with two state-of-the-art unmixing methods based on MMM: (i) Multilinear Mixing Model for Nonlinear Spectral Unmixing [(MMMNSU)](https://doi.org/10.1109/TGRS.2015.2453915) (end-members initialized by [VCA](https://doi.org/10.1109/TGRS.2005.844293)), and (ii) Unsupervised Nonlinear Spectral Unmixing Based on MMM [(UNSUBMMM)](https://doi.org/10.1109/TGRS.2017.2693366).
 
 
 
 
 | [Paper]()  <br>
 Nonlinear Extended Blind End-member and Abundance Extraction for Hyperspectral Images <br>
 [Daniel U. Campos-Delgado]()<sup>1,2</sup>, 
 [Inés A. Cruz-Guerrero]()<sup>2</sup> 
 [Juan N. Mendoza-Chavarría]()<sup>2</sup>, 
 [Aldo R. Mejía-Rodríguez]()<sup>2</sup>, 
 [Samuel Ortega]()<sup>3,4</sup>,
 [Himar Fabelo]()<sup>4</sup>,
 [Gustavo M. Callico]()<sup>4</sup> <br>
 <sup>1</sup>Optical Communication Research Institute (IICO), Autonomous University of San Luis Potosí, Av. Karakorum 1470, 78210, S.L.P., México<br>
<sup>2</sup>Faculty of Science, Autonomous University of San Luis Potosí, Av. Parque Chapultepec 1570, 78290, S.L.P., Mexico<br>
<sup>3</sup>Norwegian Institute of Food Fisheries and Aquaculture Research (NOFIMA), 9019 Tromsø, Norway<br>
<sup>4</sup>nstitute for Applied Microelectronics (IUMA), University of Las Palmas de Gran Canaria, E35017 Las Palmas de Gran Canaria, Spain<br>
In Signal Processing 2022.


 
