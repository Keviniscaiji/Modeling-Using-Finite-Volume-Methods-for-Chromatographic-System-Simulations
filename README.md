# Modeling-Using-Finite-Volume-Methods-for-Chromatographic-System-Simulations

**Group 8**  
_Members:_ Kaiwen Gong, Sze Yan How, Weilai Li, Mengjiao Wu, Yuang Zhang


## Abstract

In this project, we modeled the chromatography process using a numerical scheme developed by Medi and Amanullah. Our work extended previous studies by applying this numerical scheme to a wider range of parameters and incorporating a more complex Freundlich isotherm model, which better predicts gas adsorption on heterogeneous surfaces than the simpler Henry's law isotherm.

We employed the finite volume method with Van Leer flux limiters to simulate the chromatographic adsorption process. This choice was made for its precision in accurately simulating sharp concentration fronts without oscillations. Initially, we established a basic upwind finite volume scheme with Henry's law isotherm and subsequently extended it to include the Freundlich model to investigate the impact of various parameters, such as bed length and isotherm model parameters.

Our study also involved validating our results against existing literature, demonstrating the effectiveness of the Van Leer scheme in capturing sharp adsorption fronts. Key findings of our project include the impact of bed length and the parameter 'v' in the Freundlich model on the saturation time and adsorption capacity. Overall, our project provides valuable insights into the behavior of chromatographic adsorption processes and enhances the understanding of adsorption dynamics in heterogeneous systems.

## Project Description

This project aims to model the chromatography process using the finite volume scheme developed by Medi and Amanullah. Unlike previous studies that solely focused on numerical scheme accuracy and computational cost, this work applies the numerical scheme to a broader parameter spectrum and incorporates a more complex Freundlich isotherm model, which is widely adopted for predicting gas adsorption processes on heterogeneous surfaces. The primary goal is to gain insights into the behavior of realistic adsorbents using the Freundlich isotherm model and to assess the impact of various parameters on the adsorption process.

## References

1. **Nimibofa Ayawei, Augustus Newton Ebelegi, and Wankasi Bonbebe.** "Modelling and Interpretation of Adsorption Isotherms." _Journal of Chemistry_, 2017: 1-11. [DOI: 10.1155/2017/3039817](https://doi.org/10.1155/2017/3039817).

2. **R. Haghpanah, A. Majumder, R. Nilam, A. Rajendran, S. Farooq, I. A. Karimi, and M. Amanullah.** "Multiobjective Optimization of a Four-Step Adsorption Process for Postcombustion Co2 Capture via Finite Volume Simulation." _Industrial & Engineering Chemistry Research_, 52, no. 11 (2013): 4249–4265. [DOI: 10.1021/ie302658y](https://doi.org/10.1021/ie302658y).

3. **S. Javeed, S. Qamar, A. Seidel-Morgenstern, and G. Warnecke.** "Efficient and Accurate Numerical Simulation of Nonlinear Chromatographic Processes." _Computers & Chemical Engineering_, 35, no. 11 (2011): 2294–2305. [DOI: 10.1016/j.compchemeng.2010.10.002](https://doi.org/10.1016/j.compchemeng.2010.10.002).

4. **B. Medi and M. Amanullah.** "Application of a Finite-Volume Method in the Simulation of Chromatographic Systems: Effects of Flux Limiters." _Industrial & Engineering Chemistry Research_, 50, no. 3 (2010): 1739–1748. [DOI: 10.1021/ie100617c](https://doi.org/10.1021/ie100617c).

5. **N. S. Wilkins, A. Rajendran, and S. Farooq.** "Dynamic Column Breakthrough Experiments for Measurement of Adsorption Equilibrium and Kinetics." _Adsorption_, 27, no. 3 (2020): 397–422. [DOI: 10.1007/s10450-020-00269-6](https://doi.org/10.1007/s10450-020-00269-6).

6. **Hardiljeet K Boparai, Meera Joseph, and Denis M O’Carroll.** "Kinetics and Thermodynamics of Cadmium Ion Removal by Adsorption onto Nano Zerovalent Iron Particles." _Journal of Hazardous Materials_, 186, no. 1 (2011): 458–465.

7. **Nick D Hutson and Ralph T Yang.** "Theoretical Basis for the Dubinin-Radushkevitch (DR) Adsorption Isotherm Equation." _Adsorption_, 3 (1997): 189–195.

