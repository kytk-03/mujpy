Syntax examples (to be removed)
+++++++++++++++++++++++++++++++
useful examples of reStructuredText directives

Local fields in MuSR
--------------------
Inline: in :py:mod:`mujpy`, only a scalar coupling between :math:`\mathbf{B_{\mathrm{cont}}}` and 
:math:`\mathbf{m}_e` is allowed, proportional to :math:`\vert \psi_s (\mathbf{r}_\mu) \vert ^2`.

It is **very badly** approximated by considering

.. math::

   \mathbf{B_{\mathrm{dip}}^\prime} = \frac{\mu_0}{4 \pi} \sum _{i=1} ^N \left( -\frac{\mathbf{m}}{r^3 _i} + \frac{3 (\mathbf{m}_i \cdot \mathbf{r}_i)\mathbf{r}_i }{r^5 _i} \right)

As we said :mod:`mujpy` can only handle 1-k magnetic structures.

Example of reference
--------------------
Find below an example notebook that activates the Gui interface. An alternative way to familiarise with mujpy is to run directly the **here comes the sun** :ref:`intro` first and then come back here for a more systematic introduction. (fake bit from muesr for rst syntax) 

.. _intro: ./Introduction.html

.. code-block:: python

    >>> from mujpy.mujpy import mugui
    

(other fake bit for syntax)
A quick way to do that is using the helper function :py:func:`~muesr.utilities.ms.mago_add` from
:py:mod:`~muesr.utilities.ms`. 

link to www page `Python <http://www.python.org/>`_

example table 

     Idx | Sel | Desc. 
      0  |     | No title
      1  |  *  | FM m//c

by 
Jun Sugiyama `et al.` [Sugiyama2011]_.


.. [Sugiyama2011] Jun Sugiyama `et al.`, Phys. Rev. B 84, 054430 (2011)
