Bessel Functions
================

These are Bessel and related functions.

Bessel Functions
----------------

These are the solutions to the Bessel differential equation,

.. math::

    y'' z^2 + y' z + \left(z^2 - \nu^2\right)y = 0


Complex Only
^^^^^^^^^^^^

These functions have real counterparts defined in the C standard, but their complex versions were left out.

.. c:function:: double _Complex cj0(double _Complex x)
.. c:function:: float _Complex cj0f(float _Complex x)
.. c:function:: long double _Complex cj0l(long double _Complex x)

    :param x: Argument of the Bessel function.
    :return: The value of the Bessel function of the first kind.

    Finds the value of :math:`J_0(x)`. Calls :code:`cjn(0, x)`.

.. c:function:: double _Complex cj1(double _Complex x)
.. c:function:: float _Complex cj1f(float _Complex x)
.. c:function:: long double _Complex cj1l(long double _Complex x)

    :param x: Argument of the Bessel function.
    :return: The value of the Bessel function of the first kind.

    Finds the value of :math:`J_1(x)`. Calls :code:`cjn(1, x)`.

.. c:function:: double _Complex cjn(int n, double _Complex x)
.. c:function:: float _Complex cjnf(int n, float _Complex x)
.. c:function:: long double _Complex cjnl(int n, long double _Complex x)

    :param n: The parameter of the Bessel function.
    :param x: The argument of the Bessel function.
    :return: The value of the Bessel function of the first kind.

    Finds the value of :math:`J_n(x)`. If :code:`x` is zero, then the value returned is zero. Otherwise, the following formula is used.

    .. math::

        \frac{x^n}{2^n n!} \sum_{k = 0}^\infty \frac{(-1)^k z^{2k}}{4^k (n + 1)_k k!}

    This was chosen over the standard definition due to the :math:`4^k` in the denominator, which should improve convergence.


.. c:function:: double _Complex cy0(double _Complex x)
.. c:function:: float _Complex cy0f(float _Complex x)
.. c:function:: long double _Complex cy0l(long double _Complex x)

    :param x: Argument of the Bessel function.
    :return: The value of the Bessel function of the second kind.

    Finds the value of :math:`Y_0(x)`. Calls :code:`cyn(0, x)`.

.. c:function:: double _Complex cy1(double _Complex x)
.. c:function:: float _Complex cy1f(float _Complex x)
.. c:function:: long double _Complex cy1l(long double _Complex x)

    :param x: Argument of the Bessel function.
    :return: The value of the Bessel function of the second kind.

    Finds the value of :math:`Y_1(x)`. Calls :code:`cyn(1, x)`.

.. c:function:: double _Complex cyn(int n, double _Complex x)
.. c:function:: float _Complex cynf(int n, float _Complex x)
.. c:function:: long double _Complex cynl(int n, long double _Complex x)

    :param n: The parameter of the Bessel function.
    :param x: The argument of the Bessel function.
    :return: The value of the Bessel function of the second kind.

    Finds the value of :math:`Y_n(x)`. If :code:`x` is zero, then the value returned is :code:`NAN` to represent complex infinity. This is computed using a limit definition, since there does not seem to be a very good way to represent this function as a sum.

Real and Complex
^^^^^^^^^^^^^^^^

These are definitions of functions for real and complex arguments.

.. c:function:: double jnu(double nu, double x)
.. c:function:: float jnuf(float nu, float x)
.. c:function:: long double jnul(long double nu, long double x)
.. c:function:: double _Complex cjnu(double _Complex nu, double _Complex x)
.. c:function:: float _Complex cjnuf(float _Complex nu, float _Complex x)
.. c:function:: long double _Complex cjnul(long double _Complex nu, long double _Complex x)

    :param nu: The parameter of the Bessel function.
    :param x: The argument of the Bessel function.
    :return: The value of the Bessel function of the first kind.

    Finds the value of :math:`J_\nu(x)`. This is computed using the following sum.

    .. math::

        \frac{x^\nu}{2^\nu \Gamma(\nu + 1)} \sum_{k = 0}^\infty \frac{(-1)^k z^{2k}}{4^k (\nu + 1)_k k!}

    This was chosen over the standard definition due to the :math:`4^k` in the denominator, which should improve convergence. In some cases, this sum is skipped. These are

    .. math::

        J_\nu(0) = 0; \Re(\nu) > 0 \vee \nu \in \mathbb{Z}

    .. math::

        J_\nu(0) = \infty; \Re(\nu) < 0 \wedge \nu \not\in \mathbb{Z}

    .. math::

        J_\nu(0) = \mathrm{NaN}; \Re(\nu) = 0 \wedge \nu \ne 0
        

.. c:function:: double ynu(double nu, double x)
.. c:function:: float ynuf(float nu, float x)
.. c:function:: long double ynul(long double nu, long double x)
.. c:function:: double _Complex cynu(double _Complex nu, double _Complex x)
.. c:function:: float _Complex cynuf(float _Complex nu, float _Complex x)
.. c:function:: long double _Complex cynul(long double _Complex nu, long double _Complex x)

    :param nu: The parameter of the Bessel function.
    :param x: The argument of the Bessel function.
    :return: The value of the Bessel function of the second kind.

    Finds the value of :math:`Y_\nu(x)`. If :math:`\nu` is an integer, then this uses the definition of :code:`yn` for the appropriate type. Otherwise, it uses the following formula.

    .. math::

        Y_\nu(z) = \csc(\pi\nu) \left(\cos(\pi\nu) J_\nu(z) - J_{-\nu}(z)\right)

Modified Bessel Functions
-------------------------

These are the solutions to the modified Bessel differential equation,

.. math::

    y'' z^2 + y' z - \left(z^2 + \nu^2\right)y = 0

.. c:function:: double i0(double x)
.. c:function:: float i0f(float x)
.. c:function:: long double i0l(long double x)
.. c:function:: double _Complex ci0(double _Complex x)
.. c:function:: float _Complex ci0f(float _Complex x)
.. c:function:: long double _Complex ci0l(long double _Complex x)

    :param x: The argument to the Bessel function.
    :return: The value of the modified Bessel function of the first kind.

    Calculates :math:`I_0(x)`. Uses the following sum.

    .. math::

        I_0(x) = \sum_{k = 0}^\infty \frac{z^{2k}}{4^k (k!)^2}

    This skips the sum when :math:`x = 0`, since the value is simply 0.

.. c:function:: double i1(double x)
.. c:function:: float i1f(float x)
.. c:function:: long double i1l(long double x)
.. c:function:: double _Complex ci1(double _Complex x)
.. c:function:: float _Complex ci1f(float _Complex x)
.. c:function:: long double _Complex ci1l(long double _Complex x)

    :param x: The argument to the Bessel function.
    :return: The value of the modified Bessel function of the first kind.

    Calculates :math:`I_1(x)`. Uses the following sum.

    .. math::

        I_n(x) = \frac{z}{2} \sum_{k = 0}^\infty \frac{z^{2k}}{4^k (k + 1)! k!}

    This skips the sum when :math:`x = 0`, since the value is simply 0.

.. c:function:: double in(int n, double x)
.. c:function:: float inf(int n, float x)
.. c:function:: long double inl(int n, long double x)
.. c:function:: double _Complex cin(int n, double _Complex x)
.. c:function:: float _Complex cinf(int n, float _Complex x)
.. c:function:: long double _Complex cinl(int n, long double _Complex x)

    :param n: The parameter of the Bessel function.
    :param x: The argument of the Bessel function.
    :return: The value of the modified Bessel function of the first kind.

    Calculates :math:`I_n(x)`. Uses the following sum.

    .. math::

        I_n(x) = \frac{z^n}{2^{n} n!} \sum_{k = 0}^\infty \frac{z^{2k}}{4^k (n + 1)_k k!}

.. c:function:: double inu(double nu, double x)
.. c:function:: float inuf(float nu, float x)
.. c:function:: long double inul(long double nu, long double x)
.. c:function:: double _Complex cinu(double _Complex nu, double _Complex x)
.. c:function:: float _Complex cinuf(float _Complex nu, float _Complex x)
.. c:function:: long double _Complex cinul(long double _Complex nu, long double _Complex x)

    :param nu: The parameter of the Bessel function.
    :param x: The argument of the Bessel function.
    :return: The value of the modified Bessel function of the first kind.

    Calculates :math:`I_\nu(x)`. Uses the following sum.

    .. math::

        I_\nu(x) = \frac{z^\nu}{2^{\nu} \Gamma(\nu + 1)} \sum_{k = 0}^\infty \frac{z^{2k}}{4^k (\nu + 1)_k k!}

    The function returns some special values without computing the sum. These are

    .. math::

        I_\nu(0) = 0; \Re(\nu) > 0 \vee \nu \in \mathbb{Z}

    .. math::

        I_\nu(0) = \infty; \Re(\nu) < 0 \wedge \nu \not\in \mathbb{Z}

    .. math::

        I_\nu(0) = \mathrm{NaN}; \Re(\nu) = 0 \wedge \nu \ne 0


.. c:function:: double k0(double x)
.. c:function:: float k0f(float x)
.. c:function:: long double k0l(long double x)
.. c:function:: double _Complex ck0(double _Complex x)
.. c:function:: float _Complex ck0f(float _Complex x)
.. c:function:: long double _Complex ck0l(long double _Complex x)

    :param x: The argument to the Bessel function.
    :return: The value of the modified Bessel function of the second kind.

    Calculates :math:`K_0(x)`. Calls :code:`kn` for the appropriate arguments.

.. c:function:: double k1(double x)
.. c:function:: float k1f(float x)
.. c:function:: long double k1l(long double x)
.. c:function:: double _Complex ck1(double _Complex x)
.. c:function:: float _Complex ck1f(float _Complex x)
.. c:function:: long double _Complex ck1l(long double _Complex x)

    :param x: The argument to the Bessel function.
    :return: The value of the modified Bessel function of the second kind.

    Calculates :math:`K_1(x)`. Calls :code:`kn` for the appropriate arguments.

.. c:function:: double kn(int n, double x)
.. c:function:: float knf(int n, float x)
.. c:function:: long double knl(int n, long double x)
.. c:function:: double _Complex ckn(int n, double _Complex x)
.. c:function:: float _Complex cknf(int n, float _Complex x)
.. c:function:: long double _Complex cknl(int n, long double _Complex x)

    :param n: The parameter of the Bessel function.
    :param x: The argument of the Bessel function.
    :return: The value of the modified Bessel function of the second kind.

    Calculates :math:`K_n(x)`. Uses the following definition.

    .. math::

        K_n(x) = -\frac{\pi}{2} i^n Y_n(ix) + (-1)^n \left(i\ln(ix) - \ln(x)\right) I_n(x)

   For the real functions, if this ends up being an imaginary number, :code:`NAN` is returned instead.

.. c:function:: double knu(double nu, double x)
.. c:function:: float knuf(float nu, float x)
.. c:function:: long double knul(long double nu, long double x)
.. c:function:: double _Complex cknu(double _Complex nu, double _Complex x)
.. c:function:: float _Complex cknuf(float _Complex nu, float _Complex x)
.. c:function:: long double _Complex cknul(long double _Complex nu, long double _Complex x)

    :param nu: The parameter of the Bessel function.
    :param x: The argument of the Bessel function.
    :return: The value of the modified Bessel function of the second kind.

    Calculates :math:`K_\nu(x)`. Uses the following definition for non-integer :math:`\nu`. For integer :math:`\nu`, the appropriate :c:func:`kn` is called.

    .. math::

        K_\nu(x) = \frac{\pi}{2}\csc(\pi\nu)\left(I_{-\nu}(x) - I_{\nu}(x)\right)

Airy Functions
--------------

These are the solutions to the Airy differential equation, and their derivatives.

.. math::

    y'' - zy = 0

.. c:function:: double ai(double x)
.. c:function:: float aif(float x)
.. c:function:: long double ail(long double x)
.. c:function:: double _Complex cai(double _Complex x)
.. c:function:: float _Complex caif(float _Complex x)
.. c:function:: long double _Complex cail(long double _Complex x)

    :param x: The argument to the Airy function.
    :return: The value of the Airy function of the first kind.

    Returns the value of :math:`\mathrm{Ai}(x)`.

.. c:function:: double bi(double x)
.. c:function:: float bif(float x)
.. c:function:: long double bil(long double x)
.. c:function:: double _Complex cbi(double _Complex x)
.. c:function:: float _Complex cbif(float _Complex x)
.. c:function:: long double _Complex cbil(long double _Complex x)

    :param x: The argument to the Airy function.
    :return: The value of the Airy function of the second kind.

    Returns the value of :math:`\mathrm{Bi}(x)`.

.. c:function:: double aip(double x)
.. c:function:: float aipf(float x)
.. c:function:: long double aipl(long double x)
.. c:function:: double _Complex caip(double _Complex x)
.. c:function:: float _Complex caipf(float _Complex x)
.. c:function:: long double _Complex caipl(long double _Complex x)

    :param x: The argument to the Airy function.
    :return: The value of the first derivative of the Airy function of the first kind.

    Returns the value of :math:`\mathrm{Ai}'(x)`.

.. c:function:: double bip(double x)
.. c:function:: float bipf(float x)
.. c:function:: long double bipl(long double x)
.. c:function:: double _Complex cbip(double _Complex x)
.. c:function:: float _Complex cbipf(float _Complex x)
.. c:function:: long double _Complex cbipl(long double _Complex x)

    :param x: The argument to the Airy function.
    :return: The value of the first derivative of the Airy function of the second kind.

    Returns the value of :math:`\mathrm{Bi}'(x)`.

Struve Functions
----------------

These are the solutions to the Struve differential equation.

.. math::

    y'' z^2 + y' z + \left(z^2 - \nu^2\right)y = \frac{4}{\sqrt{\pi} \Gamma\left(\nu + \frac{1}{2}\right)} \left(\frac{z}{2}\right)^{\nu + 1}

.. c:function:: double struveh0(double x)
.. c:function:: float struveh0f(float x)
.. c:function:: long double struveh0l(long double x)
.. c:function:: double _Complex cstruveh0(double _Complex x)
.. c:function:: float _Complex cstruveh0f(float _Complex x)
.. c:function:: long double _Complex cstruveh0l(long double _Complex x)

    :param x: The argument to the Struve H function.
    :return: The value of the Struve H function.

    Calculates :math:`H_0(x)`.

.. c:function:: double struveh1(double x)
.. c:function:: float struveh1f(float x)
.. c:function:: long double struveh1l(long double x)
.. c:function:: double _Complex cstruveh1(double _Complex x)
.. c:function:: float _Complex cstruveh1f(float _Complex x)
.. c:function:: long double _Complex cstruveh1l(long double _Complex x)

    :param x: The argument to the Struve H function.
    :return: The value of the Struve H function.

    Calculates :math:`H_1(x)`.

.. c:function:: double struvehn(int n, double x)
.. c:function:: float struvehnf(int n, float x)
.. c:function:: long double struvehnl(int n, long double x)
.. c:function:: double _Complex cstruvehn(int n, double _Complex x)
.. c:function:: float _Complex cstruvehnf(int n, float _Complex x)
.. c:function:: long double _Complex cstruvehnl(int n, long double _Complex x)

    :param n: The parameter of the Struve H function.
    :param x: The argument to the Struve H function.
    :return: The value of the Struve H function.

    Calculates :math:`H_n(x)`.

.. c:function:: double struvehnu(double nu, double x)
.. c:function:: float struvehnuf(float nu, float x)
.. c:function:: long double struvehnul(long double nu, long double x)
.. c:function:: double _Complex cstruvehnu(double _Complex nu, double _Complex x)
.. c:function:: float _Complex cstruvehnuf(float _Complex nu, float _Complex x)
.. c:function:: long double _Complex cstruvehnul(long double _Complex nu, long double _Complex x)

    :param nu: The parameter of the Struve H function.
    :param x: The argument to the Struve H function.
    :return: The value of the Struve H function.

    Calculates :math:`H_\nu(x)`.

.. c:function:: double struvel0(double x)
.. c:function:: float struvel0f(float x)
.. c:function:: long double struvel0l(long double x)
.. c:function:: double _Complex cstruvel0(double _Complex x)
.. c:function:: float _Complex cstruvel0f(float _Complex x)
.. c:function:: long double _Complex cstruvel0l(long double _Complex x)

    :param x: The argument to the Struve L function.
    :return: The value of the Struve L function.

    Calculates :math:`L_0(x)`.

.. c:function:: double struvel1(double x)
.. c:function:: float struvel1f(float x)
.. c:function:: long double struvel1l(long double x)
.. c:function:: double _Complex cstruvel1(double _Complex x)
.. c:function:: float _Complex cstruvel1f(float _Complex x)
.. c:function:: long double _Complex cstruvel1l(long double _Complex x)

    :param x: The argument to the Struve L function.
    :return: The value of the Struve L function.

    Calculates :math:`L_1(x)`.

.. c:function:: double struveln(int n, double x)
.. c:function:: float struvelnf(int n, float x)
.. c:function:: long double struvelnl(int n, long double x)
.. c:function:: double _Complex cstruveln(int n, double _Complex x)
.. c:function:: float _Complex cstruvelnf(int n, float _Complex x)
.. c:function:: long double _Complex cstruvelnl(int n, long double _Complex x)

    :param n: The parameter of the Struve L function.
    :param x: The argument to the Struve L function.
    :return: The value of the Struve L function.

    Calculates :math:`L_n(x)`.

.. c:function:: double struvelnu(double nu, double x)
.. c:function:: float struvelnuf(float nu, float x)
.. c:function:: long double struvelnul(long double nu, long double x)
.. c:function:: double _Complex cstruvelnu(double _Complex nu, double _Complex x)
.. c:function:: float _Complex cstruvelnuf(float _Complex nu, float _Complex x)
.. c:function:: long double _Complex cstruvelnul(long double _Complex nu, long double _Complex x)

    :param nu: The parameter of the Struve L function.
    :param x: The argument to the Struve L function.
    :return: The value of the Struve L function.

    Calculates :math:`L_\nu(x)`.

Spherical Bessel Functions
--------------------------

These are the solutions to the spherical Bessel differential equation.

.. math::

    y'' z^2 + 2y' z + \left(z^2 - \nu(\nu + 1)\right)y = 0

.. c:function:: double spj0(double x)
.. c:function:: float spj0f(float x)
.. c:function:: long double spj0l(long double x)
.. c:function:: double _Complex cspj0(double _Complex x)
.. c:function:: float _Complex cspj0f(float _Complex x)
.. c:function:: long double _Complex cspj0l(long double _Complex x)

    :param x: The argument to the Bessel function.
    :return: The value of the spherical Bessel function of the first kind.

    Calculates :math:`j_0(x)`.

.. c:function:: double spj1(double x)
.. c:function:: float spj1f(float x)
.. c:function:: long double spj1l(long double x)
.. c:function:: double _Complex cspj1(double _Complex x)
.. c:function:: float _Complex cspj1f(float _Complex x)
.. c:function:: long double _Complex cspj1l(long double _Complex x)

    :param x: The argument to the Bessel function.
    :return: The value of the spherical Bessel function of the first kind.

    Calculates :math:`j_1(x)`.

.. c:function:: double spjn(int n, double x)
.. c:function:: float spjnf(int n, float x)
.. c:function:: long double spjnl(int n, long double x)
.. c:function:: double _Complex cspjn(int n, double _Complex x)
.. c:function:: float _Complex cspjnf(int n, float _Complex x)
.. c:function:: long double _Complex cspjnl(int n, long double _Complex x)

    :param n: The parameter to the Bessel function.
    :param x: The argument to the Bessel function.
    :return: The value of the spherical Bessel function of the first kind.

    Calculates :math:`j_n(x)`.

.. c:function:: double spjnu(double nu, double x)
.. c:function:: float spjnuf(float nu, float x)
.. c:function:: long double spjnul(long double nu, long double x)
.. c:function:: double _Complex cspjnu(double _Complex nu, double _Complex x)
.. c:function:: float _Complex cspjnuf(float _Complex nu, float _Complex x)
.. c:function:: long double _Complex cspjnul(long double _Complex nu, long double _Complex x)

    :param nu: The parameter to the Bessel function.
    :param x: The argument to the Bessel function.
    :return: The value of the spherical Bessel function of the first kind.

    Calculates :math:`j_\nu(x)`.

.. c:function:: double spy0(double x)
.. c:function:: float spy0f(float x)
.. c:function:: long double spy0l(long double x)
.. c:function:: double _Complex cspy0(double _Complex x)
.. c:function:: float _Complex cspy0f(float _Complex x)
.. c:function:: long double _Complex cspy0l(long double _Complex x)

    :param x: The argument to the Bessel function.
    :return: The value of the spherical Bessel function of the second kind.

    Calculates :math:`y_0(x)`.

.. c:function:: double spy1(double x)
.. c:function:: float spy1f(float x)
.. c:function:: long double spy1l(long double x)
.. c:function:: double _Complex cspy1(double _Complex x)
.. c:function:: float _Complex cspy1f(float _Complex x)
.. c:function:: long double _Complex cspy1l(long double _Complex x)

    :param x: The argument to the Bessel function.
    :return: The value of the spherical Bessel function of the second kind.

    Calculates :math:`y_1(x)`.

.. c:function:: double spyn(int n, double x)
.. c:function:: float spynf(int n, float x)
.. c:function:: long double spynl(int n, long double x)
.. c:function:: double _Complex cspyn(int n, double _Complex x)
.. c:function:: float _Complex cspynf(int n, float _Complex x)
.. c:function:: long double _Complex cspynl(int n, long double _Complex x)

    :param n: The parameter to the Bessel function.
    :param x: The argument to the Bessel function.
    :return: The value of the spherical Bessel function of the second kind.

    Calculates :math:`y_n(x)`.

.. c:function:: double spynu(double nu, double x)
.. c:function:: float spynuf(float nu, float x)
.. c:function:: long double spynul(long double nu, long double x)
.. c:function:: double _Complex cspynu(double _Complex nu, double _Complex x)
.. c:function:: float _Complex cspynuf(float _Complex nu, float _Complex x)
.. c:function:: long double _Complex cspynul(long double _Complex nu, long double _Complex x)

    :param nu: The parameter to the Bessel function.
    :param x: The argument to the Bessel function.
    :return: The value of the spherical Bessel function of the second kind.

    Calculates :math:`y_\nu(x)`.
