# NeverMindTheMolluscs

This is an attempt to pythonize the BASIC programs included with
[The Algorithmic Beauty of Sea Shells](https://doi.org/10.1007/978-3-662-03617-4) by
[Hans Meinhardt](https://www.eb.tuebingen.mpg.de/emeriti/hans-meinhardt).

Although an e-book is still [available from Springer](https://www.springer.com/gb/book/9783540921417), the simulations consisting of over 20K lines of BASIC code are not included with it. Second-hand copies including the original CD or even diskette can still be found.

Here, the activator-inhibitor diffusion simulations are implemented in vectorized Numpy.

Each example in the book has a reference number. To display a simulation in a Jupyter notebook:

```python
from IPython import display
from mollusc.examples import MolluscExample

lioconcha_h = MolluscExample(88)
# The red, greem anf blue channels correspond to substances 0, 1 and 5:
display.Image(data=lioconcha_h.render(substances=(0,1,5)))
```

### Simulation 88: [Lioconcha Hieroglyphica](https://en.wikipedia.org/wiki/Lioconcha_hieroglyphica)

![Lioconcha Hieroglyphica](/img/lioconcha_hieroglyphica.png)
*Lioconcha Hieroglyphica*

![Lioconcha Hieroglyphica in FreeBASIC](/img/sp_88.png)
*Original FreeBASIC*

![Lioconcha Hieroglyphica in Python](/img/nmtm_88.png)
*Python*


