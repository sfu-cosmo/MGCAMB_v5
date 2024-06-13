
# [ReACT v2: ACT-io et ReACT-io](reactions/data/doggie.jpeg)

## This repository contains version 2 of the [ReACT](https://github.com/nebblu/ReACT) code. 

## Table of contents: 

1. [Introduction](https://github.com/nebblu/ACTio-ReACTio#1-introduction)
2. [Requirements](https://github.com/nebblu/ACTio-ReACTio#2-requirements)
3. [Installation](https://github.com/nebblu/ACTio-ReACTio#3-installation)
4. [Running ReACT](https://github.com/nebblu/ACTio-ReACTio#4-running-react)
5. [Models of gravity and dark energy](https://github.com/nebblu/ACTio-ReACTio#5-available-models)
6. [Adding in new models](https://github.com/nebblu/ACTio-ReACTio#6-adding-in-new-models)
7. [Parameter description](https://github.com/nebblu/ACTio-ReACTio#7-parameters)
8. [Citation](https://github.com/nebblu/ACTio-ReACTio#8-citation)
9. [Notes](https://github.com/nebblu/ACTio-ReACTio#9-notes)
10. [What's new?](https://github.com/nebblu/ACTio-ReACTio#10-what-is-new)


## 1. Introduction

ReACT is a halo model and standard perturbation theory code based on the software packages [Copter](http://mwhite.berkeley.edu/Copter/) (0905.0479) and MG-Copter (1606.02520). It allows the efficient computation of many large scale structure observables for a wide class of gravity and dark energy models. 

The package is made up of 5 folders: 

* `reactions` : contains all the c++ source code for the calculations. Specifically, these can be found in `reactions/src`. 
* `pyreact` : contains the wrapper functions calling the c++ code from python. 
* `cosmosis` : contains the cosmosis module for ReACT. 
* `notebooks` : example python notebooks for ReACT computations plus two Mathematica notebooks performing translations between covariant theories and modifications to the perturbative Poisson equation  and comparisons of different parametrisations of the nonlinear Poisson equation.


ReACT can perform the follwing calculations for general theories beyond LCDM

* Spherical collapse  (1812.05594): `reactions/src/SCOL.cpp`.

* Halo model power spectrum (1812.05594):  `reactions/src/HALO.cpp`.

* Real and redshift space LSS 2 point statistics (1606.02520): `reactions/src/SPT.cpp`.

* Numerical perturbation theory kernels up to 3rd order (1606.02168): `reactions/src/SpericalFunctions.cpp`.

* Real space bispectrum at 1-loop level** (1808.01120): `reactions/src/BSPT.cpp`.

* Exact perturbation theory kernels up to 4th order** (1808.01120): `reactions/src/BSPTN.cpp`.

* Beyond LCDM modifications to background, linear, 2nd, 3rd order Poisson equation and fully nonlinear Poisson equation: `reactions/src/BeyondLCDM.cpp`.

** for LCDM, nDGP and f(R) gravity only. 


## 2. Requirements

### C++ Compiler and automake
ReACT is written in C++, so you'll need a modern C++ compiler.

### GSL
The ReACT extension makes use of a number of [gsl](http://www.gnu.org/software/gsl/) packages such as odeiv2 package which is part of the gsl library.  You will need to have a version of gsl installed (gsl 2. or later). 

### SUNDIALS

For Spherical Collapse module (`reactions/src/SCOL.cpp`) you will also need the [SUNDIALS](https://computing.llnl.gov/projects/sundials) package version 4.1.0 (also tested to work with version 5.0) **Note** that later versions of sundials may produce errors as they update their functions so please work with version 4.1.0. 

Get these packages using your package manager of choice (e.g., `homebrew` on mac OS).

One should also make sure that sundials is on your library path, for example 

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/sundials/install_dir/lib64:${LD_LIBRARY_PATH} 
```

### Python
A recent version of python should be installed. I have worked with Python 3.6.8 without issue. 

## 3. Installation 

ReACT can be installed by following the proceeding instructions. If you experience any issues, please have a look at the [Issues tab of version 1](https://github.com/nebblu/ReACT/issues?q=is%3Aissue+is%3Aclosed). Leave an issue in this repository if you experience any problems and we'll get back to you as soon as possible, or directly contact ben.bose@ed.ac.uk. 

We have also included a Docker file to run ReACT within a container. This is explained below. 

### General installation
The Python interface can then be installed with
```
$ python setup.py install
```
or 
```
$ pip install .
```
If you want to work on the library, install it in "developer" mode:
```
$ python setup.py develop
```
or
```
$ pip install -e .
```

### Linux specific installation 

**In all the following change the paths accordingly**

1. Make sure the following are installed: python, sundials, g++, gsl. The versions I've tested with are:

* Python 3.6.8
* g++ (GCC) 8.5.0
* sundials-4.1.0

2. Clone ACTio-ReACTio:

```
$ git clone https://github.com/nebblu/ACTio-ReACTio.git
```

3. Add in sundials directory to pyreact/Makefile :

**Example:** 

```
LDFLAGS += -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -L/home/bose/sundials/instdir/lib64
```
```
CPPFLAGS += -I/home/bose/sundials/instdir/include
```

4. Export sundials library to LD_LIBRARY_PATH:

```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/sundials/instdir/lib64:${LD_LIBRARY_PATH}
```

5. Install react

```
$ python3 setup.py develop --user
```

6. Check that everything is working. Go to the reactions/examples/ directory and try to run one of the example files. For example:

```
$ g++ -I/home/bose/react_tutorial/ACTio-ReACTio/reactions/include -L/home/bose/react_tutorial/ACTio-ReACTio/reactions/lib -lcopter -lgsl -lstdc++ halo_ps.cpp -o test
```

```
$ ./test
```

Example output:

```
$ 0 1.000000e-03 4.039658e+02 4.213977e+02 9.959108e-01 ...
```

**Note** you may also need to export the copter library path to run the example files: 

```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/react_tutorial/ACTio-ReACTio/reactions/lib:${LD_LIBRARY_PATH}
```


### Docker

We have also included a Docker file, `Dockerfile'. Docker allows you to run the code within what is called a container which is a dedicated environment built according to specifications given in the Docker file. To run ReACT in a Docker container you may follow these steps:

1. Install [docker](https://www.docker.com/) and open the app to get the docker daemon running
2. Create an empty folder on you laptop/machine
3. Copy the Dockerfile from this github-repo
4. Go to the folder with the Dockerfile you've created and run

```
docker build -t mybuild . 
```

5. Once this builds, you can jump into the container (which has ReACT installed) using the following command

```
docker run -v path-to-the-folder-with-the-Dockerfile:/home -i -t mybuild
```

This will also automatically take all the local files in `path-to-the-folder-with-the-Dockerfile` (e.g. `/Users/bbose/Desktop/myfolder` on your machine) to the `/home` folder inside the container. Anything placed in the home folder from within the container will then automatically show up locally. This lets you transfer ReACT output produced within the container to the local system. 

Alternatively, if you would like to use jupyter notebooks, then you should assign a port to your container in order to allow the connection from a host browser
```
docker run -v path-to-the-folder-with-the-Dockerfile:/home -it -p 8888:8888 mybuild
```

I've tested this works on the tests and example files within the `reactions/examples` and `reactions/tests` directories within the container. You can use the following command to run them

```
g++ -I/ACTio-ReACTio/reactions/include -L/ACTio-ReACTio/reactions/lib  spt.cpp -lgsl -lcopter
```
```
./a.out
```

To run jupyter notebooks you can jump to the 'ACTio-ReACTio/notebooks'-folder and run the following command:
```
jupyter notebook --ip 0.0.0.0 --allow-root
```
You can access the notebooks through your desktops browser on http://localhost:8888 , i.e. copy-paste the second automatically generated hyperlink.  

To exit the container simply use the exit command 

```
exit 
```

You can later come back to your image in the individual container by opening the Docker app -> going to Containers, where you will see a list of all your containers -> choose and run a container with the 'mybuild'-image you used before (e.g. charming_faraday: mybuild RUNNING PORT:8888) -> open the command line (CLI-option in the Docker app when you hover over a container). Alternatively, there is a Docker extension in VS-code where you can (re-)start your container and see the files there all at ones, which makes the navigation easier.

## 4. Running ReACT

### Python
The pyreact wrapper allows the C++ code (the native language of ReACT) to be called by Python code. Example jupyter notebooks that demonstrate the usage of ReACT can be found in the `notebooks` folder. Specifically we have included : 

* pyreact_demo.ipynb : Demonstrates the basic halo model reaction [2005.12184](https://arxiv.org/abs/2005.12184).
* pyreact_demo-eft.ipynb : Demonstrates the implementation of the EFTofDE with a screened and unscreened case [2210.01094](https://arxiv.org/abs/2210.01094). 
* pyreact_demo-fofr.ipynb : Demonstrates the implementation of various parametrisations for beyond LCDM physics against the exact solution in Hu-Sawicki f(R) [2210.01094](https://arxiv.org/abs/2210.01094).
* pyreact_demo-neutrinos.ipynb : Demonstrates the implementation of the massive neutrinos halo model reaction [2105.12114](https://arxiv.org/abs/2105.12114).
* pyreact_demo-rsd.ipynb : Demonstrates the implementation of the redshift space power spectrum multipoles [1606.02520](https://arxiv.org/abs/1606.02520).
* pyreact_demo-bspt.ipynb : Demonstrates the implementation of the real space bispectrum [1808.01120](https://arxiv.org/abs/1808.01120).


###  C++
One can also run ReACT and MG-Copter in their native C++. Again, a number of example output C++ scripts have been included in `reactions/examples` as well as a number of cosmologies in `reactions/examples/transfers`. Specifically we have included : 

* actio_reactio_phen.cpp : Example code to output the reaction and halo spectra for EFTofDE linear & phenomenological nonlinear modifications (see [this section](https://github.com/nebblu/ACTio-ReACTio#5-available-models)).
* actio_reactio_ppf.cpp : Example code to output the reaction and halo spectra for EFTofDE linear & nPPF parametrised nonlinear modifications (see [this section](https://github.com/nebblu/ACTio-ReACTio#5-available-models)).
* reaction_mnu.cpp : Example code to output the reaction and halo spectra for beyond LCDM physics & massive neutrinos.
* halo_ps.cpp : Example code to output the halo model powerspectrum.
* spt_rsd.cpp : Example code to output the 1-loop powerspectrum in real and redshift space.
* spt.cpp : Example code to output the 1-loop powerspectrum in real space.
* bs.cpp : Example code to output the real space bispectrum.


We can compile these examples with a command similar to : 

```
$ gcc -I/Users/bbose/Desktop/ACTio-ReACTio/reactions/include -L/Users/bbose/Desktop/ACTio-ReACTio/reactions/lib -lcopter -lgsl -lstdc++ bs.cpp -o test
```

Then just run 

```
$ ./test 
```

All example files require the specification of a transfer function or linear power spectrum. The given examples take transfer function inputs. You can specify the transfer function in two possible ways: 

1. Directly as produced using a Boltzmann solver (e.g. [MGCAMB](https://github.com/sfu-cosmo/MGCAMB)) as in reaction_mnu.cpp. **Note** that this file should be headerless. 
2. Otherwise, you must create your own cosmology file as in `reactions/examples/transfers' within which you should specify the associated z=0 LCDM transfer function file with the following two column format: {k[h/Mpc], T(k)} with the transfer function normalised to 1 at small k. This is then assumed by the code to rescale the power spectrum to the target, beyond LCDM cosmology  using internally computed LCDM and non-standard growth factors.  

The internal flag **modcamb** tells ReACT whether or not to treat the input transfer function/power spectrum as specifying a z=0 LCDM cosmology (false) or specifying a transfer function/ power spectrum at the target z and in the target beyond LCDM theory (true) (as produced by EFTCAMB, MGCAMB etc). We must set modcamb = false for option 2. Option 1 can have modcamb assume both false or true values depending on whether or not the specified function is LCDM z=0 or not. 

## 5. Available Models

In Pyreact we currently have the following models and model parameters 

### Specific theories:

1. gr : general relativity. No extra parameters.  
2. f(r) : [Hu-Sawicki f(R)](https://arxiv.org/abs/0705.1158). **extpars[0]** = $f_{R0}$. 
3. dgp : normal branch of [DGP gravity](https://arxiv.org/abs/hep-th/0005016). **extpars[0]** = $\Omega_{rc}$. 
4. quintessence : Quintessence. **extpars[0]** = $w_0$.
5. cpl : [CPL evolving dark energy](https://arxiv.org/abs/gr-qc/0009008) $w = w_0 + (1-a)w_a$ . **extpars[0-1]**=  $w_0,w_a$.
6. ds : [Dark Scattering with CPL background](https://arxiv.org/abs/1605.05623). **extpars[0-2]** = $w_0,w_a,\xi*h$. Note the factor of h included in $\xi$!

### Parametrised theory (all assume no 2nd and 3rd order modifications to the Poisson equation )
These models assume the [effective field theory of dark energy](https://arxiv.org/abs/1907.03150v2) (EFTofDE) at the background and linear level. We have

* **extpars[0-4]** = $\alpha_{K0} , \alpha_{B0} , \alpha_{M0} , \alpha_{T0} , M_0$ , where the `0` indicates the value today and $a M_0=\log(M^2/m_P^2)$ where $m_P^2$ is the Planck mass. 

Models assuming the small scale (large k) limit in linear modification (see for example Eq.26 of [this reference](https://arxiv.org/abs/1606.05339))

7. eftppf :  EFTofDE with a [post parametrised friedmannian (nPPF)](https://arxiv.org/abs/1608.00522) $G_{eff,non-linear}$ in spherical collapse equations. **extpars[5-12]** = $p_1,...,p_8$. 
8. eftus :  EFTofDE without screening, i.e. $G_{eff, non-linear}$ = $G_{eff,linear}$.  
9. eftss : EFTofDE with superscreening , i.e. $G_{eff, non-linear}$ = $G_{N}$.  

Full k-dependence in linear modification: 

10. fulleftppf : same as case 7. 
11. fulleftus : same as case 8
12. fullefterf : EFTofDE with a phenomenological $G_{eff,non-linear}$ in spherical collapse equations. **extpars[5-8]** = $p_1,...,p_4$. See [this paper](https://arxiv.org/abs/2210.01094) for details.

**Note** for models 10-12 we use an additional parameter extpars[19] which multiplies the EFTofDE background function c(a). We use this to set $c(a) = 0$ manually if needed. See [this paper](https://arxiv.org/abs/2210.01094) for more details.

**Note** For EFTofDE models (7-12), the scale factor dependence, 1st and 2nd scale factor derivatives of the alpha functions must be specified in `reactions/src/BeyondLCDM.cpp` - see `alphai_eft`, `dalphai_eft` and `ddalphai_eft` functions respectively. The default time dependencies are

* $\alpha_i(a)=a\alpha_{i,0}$
* $M^2/m^2_P=e^{M_0*a}$ 

13. minimal : CPL for background (**extpars[0-1] = {w0,wa}**), growth index for linear perturbations (**extpars[2] = $\gamma$**), phenomenological $G_{eff,non-linear}$ in spherical collapse (**extpars[3-6]** = $p_1,...,p_4$). 

In the C++ code the model is specified as the integer value of each model in the last argument of the functions , e.g. for the 1-loop real space power spectrum in f(R) gravity we would specify the following functional call 

```
PLOOPn2(1, k, pars, extpars, err, 2)
```

where the arguments are: 

* 1 specifies a call to the 1-loop matter-matter calculation. 
*  k is the wave mode. 
*  pars hold base cosmological parameters. 
*  extpars[0] = $|f_{R0}|$. 
*  err is the absolute error on the 1-loop integrals.
*  2 specifies Hu-Sawicki f(R) gravity 

See `reactions/src/examples` for more parameter and function details.  

**Note** that all cases except 4,5,6 assume a LCDM background expansion. 


## 6. Adding in new models

One can add in new models by going to the source code, `reactions/src/BeyondLCDM.cpp` and adding in a new `case n:` in the relevant functions. The array extpars stores the theory parameters. The default size of this array is 20 and can be increased by changing the maxpars parameter in `reactions/src/BeyondLCDM.h`. 

The functions one should consider when adding in a new models are as follows: 

### A) Functions to edit for new models 

* `HAg` : Normalised Hubble expansion $\frac{H}{H_0}$
* `HA1g` : Normalised scale factor derivative of Hubble  $aH \frac{dH}{da} \frac{1}{H_0^2}$
* `myfricF` : Euler equation friction term (only present in Dark Scattering)
* `mu` : linear modification to the Poisson equation
* `gamma2` : 2nd order modification to the Poisson equation
* `gamma3` : 3rd order modification to the Poisson equation
* `mymgF` : nonlinear modification to the Poisson equation [Needed only for spherical collapse SCOL.cpp, not for 1-loop calculations]
* `WEFF` : Effective dark energy fluid contribution to virial theorem [Needed only for spherical collapse SCOL.cpp, not for 1-loop calculations]


### B) Functions to edit for EFTofDE parametrisations (cases 7-12):

* `riccibackgroundp` : scale factor derivative of background Ricci scalar in FRLW
* `alphai_eft` : scale factor dependence of parameter : $\alpha_i(a)$
* `dalphai_eft` : scale factor derivative :  $\frac{d \alpha_i(a)}{da}$
* `ddalphai_eft` : 2nd scale factor derivative : $\frac{d^2 \alpha_i(a)}{da^2} $
* `HA3g` : normalised 2nd scale factor derivative of Hubble : $\frac{d^2 H(a)}{da^2} \frac{1}{H_0}$


### C) Functions to edit for custom background expansion: 
**Note** to use these functions in your new model, you should initialise a spline with H(a) - run hubble_init : see SpecialFunctions.cpp

* `bespokehub` : Normalised Hubble expansion: $\frac{H}{H_0}$ (this can involve a solution to some ODE)
* `bespokehubd` : $aH \frac{dH}{da} \frac{1}{H_0^2}$
* `bespokehubdd` : $\frac{d^2 H(a)}{da^2} \frac{1}{H_0}$

Once specified, you can use the following splines (initialised with hubble_init in SpecialFunctions class) in HAg and HA1g (and HA3g if EFTofDE is being considered)
* `myhubble` : $\frac{H}{H_0}$
* `myhubbled` : $aH \frac{dH}{da} \frac{1}{H_0^2}$
* `myhubbledd` : $\frac{d^2 H(a)}{da^2}\frac{1}{H_0}$


Once you have added in the required functions you need to re-compile the source code by going to the `reactions` directory and running 

```
make && make install
```

or for a full re-installation go to the top directory and run

```
python setup.py clean && python setup.py develop
```

If you want to add in this model to the Pyreact wrapper, you will also need to add in the nth model in `pyreact/react.py` to the appropriate function you want to call in python. 

## 7. Parameters

### Python parameters 

**General parameters** 

* **model** Selects which theoretical model is to be applied (see section: [Models of gravity and dark energy](https://github.com/nebblu/ACTio-ReACTio#5-available-models)). 

* **mass_loop** is the number of mass bins to be sampled over the range  5<Log10[M]<20, M in solar masses. Default is 30. 

* **is_transfer** Is the input Pk a transfer function or power spectrum? Default is linear power spectrum (false value). 

* **h,n_s,omega_m,omega_b** are the base LCDM cosmological parameters assumed in producing the linear power spectrum or transfer function input.


**Specific to the compute_reaction function : basic halo model reaction function for specific theories** 

* **sigma_8** is the input power spectrum sigma8 value at z=0.

* **z** array of output redshifts in ascending order (note that z<2.5 to ensure code stability).

* **k** array of wave mode outputs in ascending order. 

* **Pk** Linear z=0 LCDM power spectrum array. 

* **fR0,Omega_rc,w,wa** values of f(R), DGP, Quintessence or CPL model parameters. 


**Specific to the compute_reaction_ext function - halo model reaction function for more general theories** 

* **extpars** array of theory parameters (see section: [Models of gravity and dark energy](https://github.com/nebblu/ACTio-ReACTio#5-available-models)). 

**Specific to the compute_reaction_nu_ext function - for massive neutrino cosmologies** 

* **omega_nu** massive neutrino density fraction at z=0. 

* **As** is the primordial spectrum amplitude. Needed to rescale the various input transfer functions

* **pscale** is the pivot scalar.

* **Tm,Tcb** are transfer function arrays of z x k for the total matter and total CDM+baryons in the target cosmology. 

* **Tcblcdm** is the LCDM CDM+baryon transfer function array of z x k.


### C++ parameters 

**modcamb** tells ReACT whether or not to treat the input transfer function or linear power spectrum file as specified at the target theory and redshift  - true value - or as specified for LCDM z=0 which the code will then rescale using internally computed growth - false value. Default is false as in the original version of the code. 

**modg**: Tells ReACT to manually set $k_\star$ and $\mathcal{E}$ to LCDM values (1e-6 and 1. resp.). This is needed because of sensitivity of $\mathcal{E}$ to the ratio of 1-halo terms which may not be exactly equal at large scales for different cosmologies even when modified gravity isn't present. 

**Note** the cosmoSIS module has not yet been extended to include massive neutrinos. 

## 8. Citation

When using ReACT in a publication, please acknowledge the code by citing the relevant papers from the following:

[arXiv:1812.05594](https://arxiv.org/abs/1812.05594) : "On the road to percent accuracy I: non-linear reaction of the matter power spectrum to dark energy and modified gravity"

[arXiv:1909.02561](https://arxiv.org/abs/1909.02561) : "On the road to percent accuracy III: non-linear reaction of the matter power spectrum to massive neutrinos"

[arXiv:2005.12184](https://arxiv.org/abs/2005.12184) : "On the road to per-cent accuracy IV:  ReACT -- computing the non-linear power spectrum beyond LCDM"

[arXiv:2105.12114](https://arxiv.org/abs/2105.12114) : "On the road to percent accuracy V: the non-linear power spectrum beyond LCDM with massive neutrinos and baryonic feedback"

[arXiv:2111.13598](https://arxiv.org/abs/2111.13598)  : "On the road to per cent accuracy VI: the non-linear power spectrum for interacting dark energy with baryonic feedback and massive neutrinos
"

[arXiv:2210.01094](https://arxiv.org/abs/2210.01094) : "Fast and accurate predictions of the nonlinear matter power spectrum for general models of Dark Energy and Modified Gravity"

[arXiv:1606.02520](https://arxiv.org/abs/1606.02520) : "A Perturbative Approach to the Redshift Space Power Spectrum: Beyond the Standard Model"

[arXiv:1808.01120](https://arxiv.org/abs/1808.01120) : "The one-loop matter bispectrum as a probe of gravity and dark energy"



Respective bibtex entries:

```
 @article{Cataneo:2018cic,
   author = "Cataneo, Matteo and Lombriser, Lucas and Heymans, Catherine and Mead, Alexander and Barreira, Alexandre and Bose, Sownak and Li, Baojiu",
    title = "{On the road to percent accuracy: non-linear reaction of the matter power spectrum to dark energy and modified gravity}",
     eprint = "1812.05594",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1093/mnras/stz1836",
    journal = "Mon. Not. Roy. Astron. Soc.",
    volume = "488",
    number = "2",
    pages = "2121--2142",
    year = "2019"
}
```

```
@article{Cataneo:2019fjp,
    author = "Cataneo, Matteo and Emberson, J. D. and Inman, Derek and Harnois-Deraps, Joachim and Heymans, Catherine",
    title = "{On the road to per cent accuracy \textendash{} III. Non-linear reaction of the matter power spectrum to massive neutrinos}",
    eprint = "1909.02561",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1093/mnras/stz3189",
    journal = "Mon. Not. Roy. Astron. Soc.",
    volume = "491",
    number = "3",
    pages = "3101--3107",
    year = "2020"
}
```



```
@article{Bose:2020wch,
    author = "Bose, Benjamin and Cataneo, Matteo and Tröster, Tilman and Xia, Qianli and Heymans, Catherine and Lombriser, Lucas",
    title = "{On the road to per-cent accuracy IV: ReACT -- computing the non-linear power spectrum beyond $\Lambda$CDM}",
    eprint = "2005.12184",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    month = "5",
    year = "2020"
}
```

```
@article{Bose:2021mkz,
    author = "Bose, Benjamin and Wright, Bill S. and Cataneo, Matteo and Pourtsidou, Alkistis and Giocoli, Carlo and Lombriser, Lucas and McCarthy, Ian G. and Baldi, Marco and Pfeifer, Simon and Xia, Qianli",
    title = "{On the road to percent accuracy V: the non-linear power spectrum beyond $\Lambda$CDM with massive neutrinos and baryonic feedback}",
    eprint = "2105.12114",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    month = "5",
    year = "2021"
}
```

```
@article{Carrilho:2021rqo,
    author = "Carrilho, Pedro and Carrion, Karim and Bose, Benjamin and Pourtsidou, Alkistis and Hidalgo, Juan Carlos and Lombriser, Lucas and Baldi, Marco",
    title = "{On the road to~per\,cent accuracy VI: the non-linear power spectrum for interacting dark energy with baryonic feedback and massive neutrinos}",
    eprint = "2111.13598",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1093/mnras/stac641",
    journal = "Mon. Not. Roy. Astron. Soc.",
    volume = "512",
    number = "3",
    pages = "3691--3702",
    year = "2022"
}
```

```
@article{Bose:2022vwi,
    author = "Bose, B. and Tsedrik, M. and Kennedy, J. and Lombriser, L. and Pourtsidou, A. and Taylor, A.",
    title = "{Fast and accurate predictions of the nonlinear matter power spectrum for general models of Dark Energy and Modified Gravity}",
    eprint = "2210.01094",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1093/mnras/stac3783",
    month = "10",
    year = "2022"
}
```

```
@article{Bose:2016qun,
    author = "Bose, Benjamin and Koyama, Kazuya",
    title = "{A Perturbative Approach to the Redshift Space Power Spectrum: Beyond the Standard Model}",
    eprint = "1606.02520",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1088/1475-7516/2016/08/032",
    journal = "JCAP",
    volume = "08",
    pages = "032",
    year = "2016"
}
```

```
@article{Bose:2018zpk,
    author = "Bose, Benjamin and Taruya, Atsushi",
    title = "{The one-loop matter bispectrum as a probe of gravity and dark energy}",
    eprint = "1808.01120",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    reportNumber = "YITP-18-88",
    doi = "10.1088/1475-7516/2018/10/019",
    journal = "JCAP",
    volume = "10",
    pages = "019",
    year = "2018"
}
```

## 9. Notes 

### Relation between $\xi$ in CLASS and ReACT 
The parameter $\xi_{class}$ in the [Dark Scattering version of CLASS](https://github.com/PedroCarrilho/class_public/tree/IDE_DS) and the value of $\xi_{react}$ in ReACT can be related by the following formula: 

$$  \xi_{react} = h   \xi_{class}  \frac{\Omega_{cdm}}{\Omega_m} / (1. + 0.194407  h  (1. - \Omega_m) * (1. + w_0) \xi_{class} (1.-\frac{\Omega_{cdm}}{\Omega_m})) $$

The relation accounts for the fact that the CLASS code only couples Dark Energy to CDM while ReACT couples to CDM + baryons. Note we have also included the factor of h needed to relate ReACT's $\xi$ to the one used in [the original paper](https://arxiv.org/abs/1412.1080).

### Parameter ranges (Updated: 22/03/21)
* To optimise root finding within the spherical collapse portion of the code, the maximum redshift that one can solve the reaction for currently is z=2.5. 
* There are some current issues in the wCDM part of the code. Namely for very particular values of w0 and wa in the CPL evolving dark energy case, the spherical collapse library cannot solve the virial theorem. We advise sticking to the ranges 
-1.3<w0<-0.7 and -1.5<wa<0.6 to avoid these issues. 
* Spherical collapse will not solve for values of sigma_8(z=0)< 0.55 or 1.4<sigma_8(z=0).  
* We have tested the massive neutrino code for m_nu=0.48eV (Omega_nu = 0.0105) without issue in the absence of modified gravity or evolving dark energy. 
* Spherical collapse may encounter issues for high fr0 values in combination with high Omega_nu. We have tested the code for Log[fr0]=-4 with m_nu=0.3eV with success but do not recommend higher values than these in combination.  

### Miscellaneous (Updated: 28/05/20)
* If errors in spherical collapse are experienced for non-f(R) theories, try setting yenvf=0 in the scol_init function in reactions/src/HALO.cpp.
* Note if using the stand-alone version of ReACT, the reaction may have deviations away from unity of the order of ~0.1-0.3% for k<1e-3. Pyreact automatically sets it to unity at such large scales. 


## 10. What is new ?

We have implemented the following to v.2:

* Added in EFTofDE linear and background parametrisations for model independent predictions. 
* Implemented two parametrisations of the nonlinear Poisson equation modification : a Parametrised Post Friedmannian approach and a phenomenological approach based on the error function. See reactions/examples .
* Cleaned up directories. 
* Added heavy commenting in source code and notebooks and re ordered some source code for clarity. 
* Merged [Dark Scattering model](https://github.com/PedroCarrilho/ReACT/tree/react_with_interact_baryons). 
* Split off beyond LCDM functions to BeyondLCDM.cpp in src directory for easily adding in new models. 
* Upgraded all beyond LCDM functions to use an n-dimensional array for parameters allowing for arbitrary number of theory parameters to be used. 
* Upgraded redshift space functions to use an n-dimensional array for rsd and bias parameters.
* Optimised RSD multipole computation. 
* Created python wrapper for RSD multipoles. 
* Created new example python notebooks for rsd multipoles and parametrised-beyond LCDM calculations. 
* New option in compute_react_ext and compute_react_nu_ext (compute_pseudo) that computes the non-linear pseudo spectrum using Takahashi et al Halofit's prediction internally.
* Added in a new function in SpecialFunctions.cpp (hubble_init) to initialise a spline of a Hubble function H(a) calculated using the bespokehub function in BeyondLCDM.cpp. This is useful if there is not an analytic form for the Hubble function.  
* Added in Hu-Sawicki forms for $\alpha$ EFTofDE basis in BeyondLCDM.cpp. 
