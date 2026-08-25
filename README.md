# W Boson Mass Sensitivity from Muon Transverse Momentum

Monte Carlo simulation study of the statistical sensitivity to the W-boson mass using the transverse-momentum distribution of muons produced in

\[
W \rightarrow \mu\nu.
\]

The project generates W-boson decays, constructs the muon \(p_T\) spectrum for different W-boson masses, and estimates the corresponding statistical sensitivity using the Fisher information.

---

## Overview

The goal of the project is to study how much information on the W-boson mass \(M_W\) is contained in the muon transverse-momentum distribution.

The simulation follows the chain

\[
p_T^W
\longrightarrow
W
\longrightarrow
W\rightarrow\mu\nu
\longrightarrow
p_T^\mu
\longrightarrow
f(p_T^\mu|M_W)
\longrightarrow
I(M_W)
\longrightarrow
\sigma_{M_W}.
\]

The W-boson transverse momentum distribution can either be sampled from an input ROOT histogram or fixed to a given value for dedicated kinematic studies.

The Fisher information is evaluated numerically from two simulated spectra corresponding to slightly different W-boson masses.

---

## Physics

For a given W-boson mass \(M_W\), the decay

\[
W\rightarrow\mu\nu
\]

is simulated in the W rest frame and subsequently transformed into the laboratory frame.

The main observable is the transverse momentum of the muon,

\[
p_T^\mu.
\]

A probability distribution

\[
f(p_T^\mu;M_W)
\]

is obtained from the simulated events.

The sensitivity to the W-boson mass is quantified through the Fisher information

\[
I(M_W)
=
\int
\frac{1}{f(p_T^\mu;M_W)}
\left(
\frac{\partial f(p_T^\mu;M_W)}
{\partial M_W}
\right)^2
dp_T^\mu.
\]

The derivative with respect to the mass is approximated using a finite difference,

\[
\frac{\partial f}{\partial M_W}
\simeq
\frac{
f(p_T^\mu;M_W+\Delta M)
-
f(p_T^\mu;M_W)
}
{\Delta M}.
\]

For \(N\) independent events, the Cramér–Rao bound gives

\[
\sigma_{M_W}
\geq
\frac{1}{\sqrt{N I(M_W)}}.
\]

The current implementation focuses on the statistical and kinematic aspects of the problem. Detector effects, backgrounds, higher-order corrections and the full experimental likelihood are not included.

---

## Project structure

```text
bosoneW/
├── include/
│   ├── FisherInformation.hpp
│   ├── Generator.hpp
│   ├── HistUtils.hpp
│   ├── SpectrumBuilder.hpp
│   └── WDecaySampler.hpp
│
├── src/
│   ├── FisherInformation.cpp
│   ├── Generator.cpp
│   ├── HistUtils.cpp
│   ├── SpectrumBuilder.cpp
│   └── WDecaySampler.cpp
│
├── input/
│   └── parameters.txt
│
├── results/
│   └── ...
│
├── distribution_pTW.root
├── main.cpp
├── CMakeLists.txt
└── README.md