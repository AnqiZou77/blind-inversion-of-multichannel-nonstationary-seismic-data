# Description

Implementation of  (Powered by [Matlab](www.matlab.com)).

The algorithm are designed to solve blind Inversion of Multichannel Nonstationary Seismic Data for Acoustic Impedance(AI) and wavelet (Powered by [Matlab](www.matlab.com)). Its key features are:

- the algorithem is for blind inversion which inveted the AI and wavelet simultaneously  ;
- for AI inversion, a new kind of weighted isotropic total variation (WITV) regularization scheme is proposed to preserve the texture of the data and suppress the noise ;
- the alternating-direction method of multipliers (ADMM) iterative scheme is used for the optimization method which helps efficiently refine the inverse procedure.

Developer: Zou Anqi (zouanqi@mail.iggcas.ac.cn)

## Reference

Zou, A.Q, Wang, Y.F, Wang, D.H, 2022, Blind Inversion of Multichannel Nonstationary Seismic Data for Acoustic Impedance and Wavelet: Pure and Applied Geophysics, 179: 2147-2166.

## Usage

Using matlab IDE, run the default `training.py` or `testing.py`.
Run Demo4nonstatAIinversion.m for AI inversion when wavelet is known.
Run Demo4nonstatWaveEst.m for wavelet estimation when AI is known.
Run Demo4NSBD.m for blind inversion for AI and wavelet.

the parameter


### dataset

The test dataset is not included in this repo. User can put into their on sythetic or real seismic data. If you need the test data, please contact with me zouanqi@mail.iggcas.ac.cn.

