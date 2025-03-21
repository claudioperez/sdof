
``sdof``
========

Lightning-fast integration of single degree-of-freedom systems.

.. container::

   |Latest PyPI version|

--------------

This package solves scalar differential equations of the form

.. math::


   m \ddot{u} + c \dot{u} + k u = f(t)

Integration is carried out using a Generalized - :math:`\alpha`
integrator that is implemented under the hood in highly optimized
multi-threaded C code.

Generalized - :math:`\alpha` is an implicit method that allows for high
frequency energy dissipation and second order accuracy. With the right
selection of parameters, the method can be specialized to the
Hibert-Hughes-Taylor (HHT), or Newmark families of integration schemes.

.. raw:: html

   <hr />


   
**Version**: |version|

**Download documentation**:
`Historical versions of documentation <https://numpy.org/doc/>`_
   
**Useful links**:
`Installation <https://numpy.org/install/>`_ |
`Source Repository <https://github.com/numpy/numpy>`_ |
`Issue Tracker <https://github.com/numpy/numpy/issues>`_ |
`Q&A Support <https://numpy.org/gethelp/>`_ |
`Mailing List <https://mail.python.org/mailman/listinfo/numpy-discussion>`_


.. grid:: 2

    .. grid-item-card::
        :img-top: _static/index-images/getting_started.svg

        Getting Started
        ^^^^^^^^^^^^^^^

        New to NumPy? Check out the Absolute Beginner's Guide. It contains an
        introduction to NumPy's main concepts and links to additional tutorials.

        +++

        .. button-ref:: user/absolute_beginners
            :expand:
            :color: secondary
            :click-parent:

            To the absolute beginner's guide

    .. grid-item-card::
        :img-top: _static/index-images/user_guide.svg

        User Guide
        ^^^^^^^^^^

        The user guide provides in-depth information on the
        key concepts of NumPy with useful background information and explanation.

        +++

        .. button-ref:: user
            :expand:
            :color: secondary
            :click-parent:

            To the user guide

    .. grid-item-card::
        :img-top: _static/index-images/api.svg

        Library Reference
        ^^^^^^^^^^^^^^^^^

        The reference guide contains a detailed description of the functions,
        modules, and objects included in NumPy. The reference describes how the
        methods work and which parameters can be used. It assumes that you have an
        understanding of the key concepts.

        +++

        .. button-ref:: library
            :expand:
            :color: secondary
            :click-parent:

            To the reference guide

    .. grid-item-card::
        :img-top: _static/index-images/contributor.svg

        Examples
        ^^^^^^^^

        Want to add to the codebase? Can help add translation or a flowchart to the
        documentation? The contributing guidelines will guide you through the
        process of improving NumPy.

        +++

        .. button-ref:: examples
            :expand:
            :color: secondary
            :click-parent:

            To the contributor's guide



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   library/index.rst
   examples/index.rst

..
   theory/index


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`


.. |Latest PyPI version| image:: https://img.shields.io/pypi/v/sdof?logo=pypi&style=for-the-badge
   :target: https://pypi.python.org/pypi/sdof

