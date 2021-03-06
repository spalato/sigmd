# sigmd.py - Multidimensional spectroscopy signals analysis

This module defines tools helpful to the study and isolation of 2D signals.
The complete expressions for the detected non-linear signals can easily contain
over 100 terms. This makes bookkeeping unwieldy, despite the expressions being
relatively simple sums and products.

This module uses symbolic mathematics (`sympy`) to keep track of all the
terms in the detected signals. This allows the study of experimental signal 
retrieval for multiple pulse experiments up to arbitrary order.

The module is meant as a toolkit for scripted work (ie: in a jupyter notebook).
Representative notebooks are supplied as examples of signal isolation using
phase matching, phase cycling and chopping. These notebooks aim to reproduce
results from the literature and techniques common to the multidimensional
spectroscopy community.
 
This toolkit is not comprehensive: it aims at helping with the study of phase
cycling schemes. This of course not the only way to isolate signals: rotating
frames, frequency filtering and other physical arguments are not covered
(although some cases shouldn't be too hard.). 
If you want to expand it or integrate this into a larger
project, let us know!

There may be some difficulty rendering equations in github markdown, sorry for any eyesores!
Equations rendered using: Math to image VSCode extension. 
In case of malfunction, direct your criticism at: https://github.com/contact
 
# Installing
The module requires python and sympy. To view and edit the notebooks, you will
need Jupyter.
The quickest way to install all these things is by installing the Anaconda
distribution (https://www.continuum.io/downloads). 

If you prefer a manual installation, have a look at:

- Python (though, if you need this, we strongly recommend Anaconda.): https://www.python.org/downloads/
- Sympy: http://docs.sympy.org/latest/index.html
- Jupyter: http://jupyter.readthedocs.io/en/latest/install.html

Once these things are setup, you can get the project by git clone.

# How it works
This module defines a `Signal` object, which describes a given non-linear response
field, and some helper functions that automate their generation (ex: all 3rd
order signals). 

For a given detection scheme, some signals need to be neglected. For example,
the pump-probe geometry will only detect signals emitted in the <!-- $k_3$ --> <img src="https://render.githubusercontent.com/render/math?math=k_3"> direction.
This is implemented by `filters`, which can be supplied during the generation of
the fields or applied afterwards.

Once this is done, the detected fields need to be added (ie: <!-- $E_3$ --> <img src="https://render.githubusercontent.com/render/math?math=E_3"> in pump-probe
geometry) and the whole expression abs-squared and further manipulated (neglect
signal-signal interference). The complete mess is then phase cycled (using
substitution), to yield the detected signals. This part can be done 
entirely using sympy functions.

The notebooks show this complete process in action; the pump-probe notebook is
a good place to start.

## Signals
Signals are defined in `signals.py`. A `Signal` object is a sympy symbol
representing a non-linear response. The prefered way to create them is
by using the helper functions: `signal(k)` for a single signal and 
`signal_for_order(n, m)`, for multiple signals.

```python
>>> signal([-1, 2, 3])
```
returns: <!-- $\chi_{-k_1+k_2+k_3} A_{2} A_{3} \overline{A_{1}}$ --> <img src="https://render.githubusercontent.com/render/math?math=%5Cchi_%7B-k_1%2Bk_2%2Bk_3%7D%20A_%7B2%7D%20A_%7B3%7D%20%5Coverline%7BA_%7B1%7D%7D">


The signal <!-- $\chi$ --> <img src="https://render.githubusercontent.com/render/math?math=%5Cchi"> is indexed by the light-matter interactions
that gave rise to it: `signal([-1, 2, 3])` generates the signal resulting from the
interactions with <!-- $-k_1+k_2+k_3$ --> <img src="https://render.githubusercontent.com/render/math?math=-k_1%2Bk_2%2Bk_3">. It is thus a sum of Feynmann
pathways. The signals are multiplied by complex prefactors <!-- $A_i$ --> <img src="https://render.githubusercontent.com/render/math?math=A_i"> used for phase cycling and chopping. 

Not that this is **not** the coherence transfer pathways <!-- $\vec{\alpha}$ --> <img src="https://render.githubusercontent.com/render/math?math=%5Cvec%7B%5Calpha%7D"> used in [our paper][Seiler JCP 2017].  This notation allows the separation of signals with
identical phases, such as the linear absorption <!-- $k_3$ --> <img src="https://render.githubusercontent.com/render/math?math=k_3">, the pump-probe signal <!-- $k_1-k_1+k_3$ --> <img src="https://render.githubusercontent.com/render/math?math=k_1-k_1%2Bk_3"> and the transient grating <!-- $k_3-k_3+k_3$ --> <img src="https://render.githubusercontent.com/render/math?math=k_3-k_3%2Bk_3">, which all have <!-- $\alpha=(0,0,1)$ --> <img src="https://render.githubusercontent.com/render/math?math=%5Calpha%3D(0%2C0%2C1)">.

To generate multiple signals at once, use the `signals_for_order(n, m)`
function. This function will return the sum of all signals of order `n` from `m`
pulses. This function takes two optional arguments:

1. `strict` for strict ordering (ex: pulse 1 never follows pulse 2) (default: True)
2. `filters` to filter out some pulse orders. They are detailed in the next section.

```python
>>> signals_for_order(3,3, filters=pump_probe)
```
<!-- $$
\begin{aligned}
&\chi_{+k_1-k_1+k_3} A_{1} A_{3} \overline{A_{1}} + \chi_{+k_1-k_2+k_3} A_{1} A_{3} \overline{A_{2}} + \chi_{+k_2-k_2+k_3} A_{2} A_{3} \overline{A_{2}} \\ &+ \chi_{+k_3+k_3-k_3} A_{3}^{2} \overline{A_{3}} + \chi_{+k_3-k_3+k_3} A_{3}^{2} \overline{A_{3}} + \chi_{-k_1+k_1+k_3} A_{1} A_{3} \overline{A_{1}} \\ &+ \chi_{-k_1+k_2+k_3} A_{2} A_{3} \overline{A_{1}} + \chi_{-k_2+k_2+k_3} A_{2} A_{3} \overline{A_{2}} + \chi_{-k_3+k_3+k_3} A_{3}^{2} \overline{A_{3}}
\end{aligned}
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0D%0A%26%5Cchi_%7B%2Bk_1-k_1%2Bk_3%7D%20A_%7B1%7D%20A_%7B3%7D%20%5Coverline%7BA_%7B1%7D%7D%20%2B%20%5Cchi_%7B%2Bk_1-k_2%2Bk_3%7D%20A_%7B1%7D%20A_%7B3%7D%20%5Coverline%7BA_%7B2%7D%7D%20%2B%20%5Cchi_%7B%2Bk_2-k_2%2Bk_3%7D%20A_%7B2%7D%20A_%7B3%7D%20%5Coverline%7BA_%7B2%7D%7D%20%5C%5C%20%26%2B%20%5Cchi_%7B%2Bk_3%2Bk_3-k_3%7D%20A_%7B3%7D%5E%7B2%7D%20%5Coverline%7BA_%7B3%7D%7D%20%2B%20%5Cchi_%7B%2Bk_3-k_3%2Bk_3%7D%20A_%7B3%7D%5E%7B2%7D%20%5Coverline%7BA_%7B3%7D%7D%20%2B%20%5Cchi_%7B-k_1%2Bk_1%2Bk_3%7D%20A_%7B1%7D%20A_%7B3%7D%20%5Coverline%7BA_%7B1%7D%7D%20%5C%5C%20%26%2B%20%5Cchi_%7B-k_1%2Bk_2%2Bk_3%7D%20A_%7B2%7D%20A_%7B3%7D%20%5Coverline%7BA_%7B1%7D%7D%20%2B%20%5Cchi_%7B-k_2%2Bk_2%2Bk_3%7D%20A_%7B2%7D%20A_%7B3%7D%20%5Coverline%7BA_%7B2%7D%7D%20%2B%20%5Cchi_%7B-k_3%2Bk_3%2Bk_3%7D%20A_%7B3%7D%5E%7B2%7D%20%5Coverline%7BA_%7B3%7D%7D%0D%0A%5Cend%7Baligned%7D%0D"></div>


Notebooks provide complete examples.

## Filters

Filters help in eliminating contributions to the detected signal. They can be
used to enforce experimental phase matching constraints or to neglect weak
contributions. Useful filters are collected in the `filters` module. We 
distinguish two kinds of filters: index filters and
term filters.

Index filters are applied during signal generation by `signals_for_order`.
The filters are applied to the indices before the signal itself is created.
Filtering out elements earlier in the process makes everything simpler. An
example of this filter is `emit_along`, which lets you select signals which
phase match a given wavevector (ex: <!-- $k_3$ --> <img src="https://render.githubusercontent.com/render/math?math=k_3">). Theses filters take a list of
indices (ex: `[-1, 2, 3]`) and return true if they should be kept.

Term filters are used to eliminate terms in the sum of contributions. They
operate on the symbolic maths. This requires working with the internal sympy
structures. As such, the definitions of these filters tends to be more
complicated. There are two strategies to achieve the desired effect.

In the first (simpler) strategy, you convert the symbol to it's string
representation 
(using `repr`) then apply string matching logic. For example, you can neglect
signal-signal interference terms by keeping only the terms whose string
representation contains `'E_'`, indicating an exciting field.

The second technique revolves around the underlying sympy expression tree. This
requires understanding of the sympy internals, but allows more powerful
specifications, such as filtering by coherence transfer. For more details, see:
http://docs.sympy.org/dev/tutorial/manipulation.html .


# Notebooks
The supplied notebooks act simultaneously as examples, benchmarks and reference.
They are located in `notebooks/`.

The most important notebooks are:

- `signals.ipynb` Demonstrates the basic properties of the signals.
- `pump_probe.ipynb` Examples of signal isolation using phase matching, phase
cycling and chopping. This notebook is hopefully well documented and serves as
an interactive tutorial of sorts.
- `collinear_3p.ipynb` Collinear detection for 3 pulses. Demonstrates the phase
 cycling from [Seiler JCP 2017].
 
Other notebooks: None yet!

Whishlist:

- Tan 3 step PP (or merge in PP)?
- boxcars
- action based.

## Extensions
The following notebook extensions will make for a more pleasurable experience
 (see https://github.com/ipython-contrib/jupyter_contrib_nbextensions):

- Collapsible Headings
- Table of Contents.


# Future

I have since moved on to other problems, but feel free to run with it.


On the horizon was:

Action detection (maybe already works), Feymann/Liouville Pathways expansion,
your contribution.

# requirements

check if python3 works. (I think it does)

# Links

[Seiler JCP 2017]: http://aip.scitation.org/doi/10.1063/1.4990500
[Rendering of equations]: https://github.com/TeamMeow/vscode-math-to-image

Seiler JCP 2017: http://aip.scitation.org/doi/10.1063/1.4990500

Rendering of equations: https://github.com/TeamMeow/vscode-math-to-image
