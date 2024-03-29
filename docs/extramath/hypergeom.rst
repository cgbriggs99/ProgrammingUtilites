Hypergeometric Functions
========================

These are hypergeometric and related functions.

Hypergeometric Functions
------------------------

These may be considered analytic continuations of the sum

.. math::

    _pF_q(a_1, a_2, \ldots, a_p; b_1, b_2, \ldots , b_q; z) = \sum_{k = 0}^\infty \frac{\prod_{i = 1}^p (a_i)_k}{\prod_{i = 1}^q (b_i)_q} z^k

In the case of :math:`U(a, b, z)`, it can be considered one of the solutions to the following differential equation, where the other is :math:`_1F_1(a; b; z)`.

.. math::

    zy'' + (b - z)y' - ay = 0


.. c:function:: double hyper1f1(double a, double b, double z)
.. c:function:: float hyper1f1f(float a, float b, float z)
.. c:function:: long double hyper1f1l(long double a, long double b, long double z)
.. c:function:: double _Complex chyper1f1(double _Complex a, double _Complex b, double _Complex z)
.. c:function:: float _Complex chyper1f1f(float _Complex a, float _Complex b, float _Complex z)
.. c:function:: long double _Complex chyper1f1l(long double _Complex a, long double _Complex b, long double _Complex z)

    :param a: The first parameter to the hypergeometric series.
    :param b: The second parameter to the hypergeometric series.
    :param z: The value to compute at.
    :return: The value of the given hypergeometric series.

    Computes :math:`_1F_1(a; b; z)`.

.. c:function:: double hyperu(double a, double b, double z)
.. c:function:: float hyperuf(float a, float b, float z)
.. c:function:: long double hyperul(long double a, long double b, long double z)
.. c:function:: double _Complex chyperu(double _Complex a, double _Complex b, double _Complex z)
.. c:function:: float _Complex chyperuf(float _Complex a, float _Complex b, float _Complex z)
.. c:function:: long double _Complex chyperul(long double _Complex a, long double _Complex b, long double _Complex z)

    :param a: The fisrt parameter to the hypergeometric series.
    :param b: The second parameter to the hypergeometric series.
    :param z: The value to compute at.
    :return: The value of the given hypergeometric series.

    Computes :math:`U(a, b, z)`.

.. c:function:: double hyper2f1(double a, double b, double c, double z)
.. c:function:: float hyper2f1f(float a, float b, float c, float z)
.. c:function:: long double hyper2f1l(long double a, long double b, long double c, long double z)
.. c:function:: double _Complex chyper2f1(double _Complex a, double _Complex b, double _Complex c, double _Complex z)
.. c:function:: float _Complex chyper2f1f(float _Complex a, float _Complex b, float _Complex c, float _Complex z)
.. c:function:: long double _Complex chyper2f1l(long double _Complex a, long double _Complex b, long double _Complex c, long double _Complex z)

    :param a: The first numerator argument.
    :param b: The second numerator argument.
    :param c: The denominator argument.
    :param z: The position to evaluate the function at.
    :return: The value of the hypergeometric function.

    Calculates :math:`_2F_1(a, b; c; z)`.

.. c:function:: double hyperpfq(unsigned int p, unsigned int q, const double *a, const double *b, double z)
.. c:function:: float hyperpfqf(unsigned int p, unsigned int q, const float *a, const float *b, float z)
.. c:function:: long double hyperpfql(unsigned int p, unsigned int q, const long double *a, const long double *b, long double z)
.. c:function:: double _Complex chyperpfq(unsigned int p, unsigned int q, const double _Complex *a, const double _Complex *b, double _Complex z)
.. c:function:: float _Complex chyperpfqf(unsigned int p, unsigned int q, const float _Complex *a, const float _Complex *b, float _Complex z)
.. c:function:: long double _Complex chyperpfql(unsigned int p, unsigned int q, const long double _Complex *a, const long double _Complex *b, long double _Complex z)

    :param p: The number of numerator arguments.
    :param q: The number of denominator arguments.
    :param a: The numerator arguments.
    :param b: The denominator arguments.
    :param z: The value to compute at.
    :return: The value of the general hypergeometric function.

    Computes :math:`_pF_q\left(a_1, a_2, \ldots , a_p; b_1, b_2, \ldots b_q; z\right)`.

Whittaker Functions
-------------------

These are the solutions to the differential equation,

.. math::

    y'' - \frac{\left(z^2 - 4az + 4b^2 - 1\right)}{4z^2} y = 0

.. c:function:: double whittakerm(double a, double b, double z)
.. c:function:: float whittakermf(float a, float b, float z)
.. c:function:: long double whittakerml(long double a, long double b, long double z)
.. c:function:: double _Complex cwhittakerm(double _Complex a, double _Complex b, double _Complex z)
.. c:function:: float _Complex cwhittakermf(float _Complex a, float _Complex b, float _Complex z)
.. c:function:: long double _Complex cwhittakerml(long double _Complex a, long double _Complex b, long double _Complex z)

    :param a: The first parameter to the Whittaker function.
    :param b: The second parameter to the Whittaker function.
    :param z: The position to compute at.
    :return: The value of the Whittaker M function.

    Computes :math:`M_{a,b}(z)`.

.. c:function:: double whittakerw(double a, double b, double z)
.. c:function:: float whittakerwf(float a, float b, float z)
.. c:function:: long double whittakerwl(long double a, long double b, long double z)
.. c:function:: double _Complex cwhittakerw(double _Complex a, double _Complex b, double _Complex z)
.. c:function:: float _Complex cwhittakerwf(float _Complex a, float _Complex b, float _Complex z)
.. c:function:: long double _Complex cwhittakerwl(long double _Complex a, long double _Complex b, long double _Complex z)

    :param a: The first parameter to the Whittaker function.
    :param b: The second parameter to the Whittaker function.
    :param z: The position to compute at.
    :return: The value of the Whittaker W function.

    Computes :math:`W_{a,b}(z)`.

Meijer Function
---------------

A rather complicated function with several parameters.

.. c:function:: double meijerg(unsigned int p, unsigned int q, unsigned int m, unsigned int n, const double *a, const double *b, double z)
.. c:function:: float meijergf(unsigned int p, unsigned int q, unsigned int m, unsigned int n, const float *a, const float *b, float z)
.. c:function:: long double meijergl(unsigned int p, unsigned int q, unsigned int m, unsigned int n, const long double *a, const long double *b, long double z)
.. c:function:: double _Complex cmeijerg(unsigned int p, unsigned int q, unsigned int m, unsigned int n, const double _Complex *a, const double _Complex *b, double _Complex z)
.. c:function:: float _Complex cmeijergf(unsigned int p, unsigned int q, unsigned int m, unsigned int n, const float _Complex *a, const float _Complex *b, float _Complex z)
.. c:function:: long double _Complex cmeijergl(unsigned int p, unsigned int q, unsigned int m, unsigned int n, const long double _Complex *a, const long double _Complex *b, long double _Complex z)

    :param p: Length of the a's.
    :param q: Length of the b's.
    :param m: Index for splitting the b's.
    :param n: Index for splitting the a's.
    :param a: First set of values.
    :param b: Second set of values.
    :param z: The value to compute the function at.
    :return: The Meijer G function.

    Computes the following:

    .. math::

        G_{p, q}^{m, n}\left(z \left|\begin{matrix}
                  a_1, a_2, \ldots , a_p\\
                  b_1, b_2, \ldots , b_q
                  \end{matrix}\right.\right)
            
Appel Function
--------------

This is a generalization of the hypergeometric function, defined as the following.

.. math::

    F_1\left(a; b_1, b_2; c; z_1, z_2\right) = \sum_{k = 0}^\infty \sum_{l = 0}^\infty \frac{(a)_{k + l} (b_1)_k (b_2)_l z_1^k z_2^l}{(c)_{k + l} k! l!}

.. c:function:: double apellf1(double a, double b1, double b2, double c, double z1, double z2)
.. c:function:: float apellf1f(float a, float b1, float b2, float c, float z1, float z2)
.. c:function:: long double apellf1l(long double a, long double b1, long double b2, long double c, long double z1, long double z2)
.. c:function:: double _Complex capellf1(double _Complex a, double _Complex b1, double _Complex b2, double _Complex c, double _Complex z1, double _Complex z2)
.. c:function:: float _Complex capellf1f(float _Complex a, float _Complex b1, float _Complex b2, float _Complex c, float _Complex z1, float _Complex z2)
.. c:function:: long double _Complex capellf1l(long double _Complex a, long double _Complex b1, long double _Complex b2, long double _Complex c, long double _Complex z1, long double _Complex z2)

    :params a, b1, b2, c: See the equation above.
    :params z1, z2: The positions to evaluate this function at.
    :return: The Appel F1 function.

    Computes :math:`F_1(a; b_1, b_2; c; z_1, z_2)`.
