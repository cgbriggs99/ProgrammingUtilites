Integrated Functions
====================

These are functions defined by integrals.

Gamma Functions
---------------

These tend to be integrals of the form

.. math::

    \int_0^\infty t^{z-1} e^{-t} dt

Complex Only
^^^^^^^^^^^^

.. c:function:: double _Complex clgamma(double _Complex z)
.. c:function:: float _Complex clgammaf(float _Complex z)
.. c:function:: long double _Complex clgammal(long double _Complex z)

    :param z: The argument to the log gamma function.
    :return: The log-gamma function.

    Computes the logarithmic gamma function for complex arguments. It is implemented as :code:`clog(ctgamma(z))`, though mathematically this is wrong. A better implementation may be added in the future.

.. c:function:: double _Complex ctgamma(double _Complex z)
.. c:function:: float _Complex ctgammaf(float _Complex z)
.. c:function:: long double _Complex ctgammal(long double _Complex z)

    :param z: The argument to the gamma function.
    :return: The gamma function.

    Computes the gamma function for complex arguments. Uses the Lanczos approximation with :math:`g = 8` and :math:`n = 12`. Values taken from Wikipedia. May have values computed at compile time for different parameters. This function is defined by

    .. math::

        \Gamma(z) = \int_0^\infty t^{z - 1} e^{-t} dz

    for real :math:`z`. 

Real and Complex
^^^^^^^^^^^^^^^^

.. c:function:: double inGamma(double a, double z)
.. c:function:: float inGammaf(float a, float z)
.. c:function:: long double inGammal(long double a, long double z)
.. c:function:: double _Complex cinGamma(double _Complex a, double _Complex z)
.. c:function:: float _Complex cinGammaf(float _Complex a, float _Complex z)
.. c:function:: long double _Complex cinGammal(long double _Complex a, long double _Complex z)

    :param a: The exponent for the incomplete gamma function.
    :param z: The lower bound for the incomplete gamma function.

    Computes the upper incomplete gamma function. This is defined by the following integral.

    .. math::

        \Gamma(a, z) = \int_z^\infty t^{a - 1} e^{-t} dt

    for real :math:`a, z`. This is implemented as the following series.

    .. math::

        \Gamma(a, z) = \Gamma(a) - z^a \sum_{k = 0}^\infty \frac{(-z)^k}{(a + k) k!}

.. c:function:: double ingamma(double a, double z)
.. c:function:: float ingammaf(float a, float z)
.. c:function:: long double ingammal(long double a, long double z)
.. c:function:: double _Complex cingamma(double _Complex a, double _Complex z)
.. c:function:: float _Complex cingammaf(float _Complex a, float _Complex z)
.. c:function:: long double _Complex cingammal(long double _Complex a, long double _Complex z)

    :param a: The exponent for the incomplete gamma function.
    :param z: The lower bound for the incomplete gamma function.

    Computes the lower incomplete gamma function. This is defined by the following integral.

    .. math::

        \gamma(a, z) = \int_0^z t^{a - 1} e^{-t} dt

    for real :math:`a, z`. This is implemented as the following series.

    .. math::

        \gamma(a, z) = z^a \sum_{k = 0}^\infty \frac{(-z)^k}{(a + k) k!}

.. c:function:: double ingamma2(double a, double z1, double z2)
.. c:function:: float ingamma2f(float a, float z1, float z2)
.. c:function:: long double ingamma2l(long double a, long double z1, long double z2)
.. c:function:: double _Complex cingamma2(double _Complex a, double _Complex z1, double _Complex z2)
.. c:function:: float _Complex cingamma2f(float _Complex a, float _Complex z1, float _Complex z2)
.. c:function:: long double _Complex cingamma2l(long double _Complex a, long double _Complex z1, long double _Complex z2)

    :param a: The exponent for the incomplete gamma function.
    :param z1: The lower bound for the incomplete gamma function.
    :param z2: The upper bound for the incomplete gamma function.
    :return: The incomplete gamma function between two bounds.

    Computes the generalized incomplete gamma function. This is defined by the following integral.

    .. math::

        \gamma(a, z) = \int_{z_1}^{z_2} t^{a - 1} e^{-t} dt

    for real :math:`a, z_1, z_2`. This is implemented by subtracting lower incomplete gamma functions.

.. c:function:: double subfact(double z)
.. c:function:: float subfactf(float z)
.. c:function:: long double subfactl(long double z)
.. c:function:: double _Complex csubfact(double _Complex z)
.. c:function:: float _Complex csubfactf(float _Complex z)
.. c:function:: long double _Complex csubfactl(long double _Complex z)

    :param z: The number for the subfactorial.
    :return: The subfactorial.

    Computes :math:`!z`, the subfactorial. This is defined and computed as

    .. math::

        !z = \frac{\Gamma(z + 1, -1)}{e}

.. c:function:: double inGammareg(double a, double z)
.. c:function:: float inGammaregf(float a, float z)
.. c:function:: long double inGammaregl(long double a, long double z)
.. c:function:: double _Complex cinGammareg(double _Complex a, double _Complex z)
.. c:function:: float _Complex cinGammaregf(float _Complex a, float _Complex z)
.. c:function:: long double _Complex cinGammaregl(long double _Complex a, long double _Complex z)

    :param a: The exponent to the gamma function.
    :param z: The lower bound for the gamma function.
    :return: The regularized incomplete gamma function.

    Computes the regularized upper incomplete gamma function. This is defined and implemented as

    .. math::

        Q(a, z) = \frac{\Gamma(a, z)}{\Gamma(a)}

.. c:function:: double ingammareg(double a, double z)
.. c:function:: float ingammaregf(float a, float z)
.. c:function:: long double ingammaregl(long double a, long double z)
.. c:function:: double _Complex cingammareg(double _Complex a, double _Complex z)
.. c:function:: float _Complex cingammaregf(float _Complex a, float _Complex z)
.. c:function:: long double _Complex cingammaregl(long double _Complex a, long double _Complex z)

    :param a: The exponent to the gamma function.
    :param z: The lower bound for the gamma function.
    :return: The regularized incomplete gamma function.

    Computes the regularized lower incomplete gamma function. This is defined and implemented as

    .. math::

        q(a, z) = \frac{\gamma(a, z)}{\Gamma(a)}

.. c:function:: double ingamma2reg(double a, double z1, double z2)
.. c:function:: float ingamma2regf(float a, float z1, float z2)
.. c:function:: long double ingamma2regl(long double a, long double z1, long double z2)
.. c:function:: double _Complex cingamma2reg(double _Complex a, double _Complex z1, double _Complex z2)
.. c:function:: float _Complex cingamma2regf(float _Complex a, float _Complex z1, float _Complex z2)
.. c:function:: long double _Complex cingamma2regl(long double _Complex a, long double _Complex z1, long double _Complex z2)

    :param a: The exponent for the incomplete gamma function.
    :param z1: The lower bound for the incomplete gamma function.
    :param z2: The upper bound for the incomplete gamma function.
    :return: The regularized incomplete gamma function between two bounds.

    Computes the regularized incomplete gamma function. This is defined and implemented as

    .. math::

        Q(a, z_1, z_2) = \frac{\Gamma(a, z_1, z_2)}{\Gamma(a)}

.. c:function:: double polygamma(double z)
.. c:function:: float polygammaf(float z)
.. c:function:: long double polygammal(long double z)
.. c:function:: double _Complex cpolygamma(double _Complex z)
.. c:function:: float _Complex cpolygammaf(float _Complex z)
.. c:function:: long double _Complex cpolygammal(long double _Complex z)

    :param z: The argument to the polygamma function.
    :return: The value of the polygamma function.

    Polygamma functions are defined as derivatives of the logarithmic gamma function. This function computes :math:`\psi^(0)(z) = \frac{d}{dz} \ln\Gamma(z)`.

.. c:function:: double polygamman(int n, double z)
.. c:function:: float polygammanf(int n, float z)
.. c:function:: long double polygammanl(int n, long double z)
.. c:function:: double _Complex cpolygamman(int n, double _Complex z)
.. c:function:: float _Complex cpolygammanf(int n, float _Complex z)
.. c:function:: long double _Complex cpolygammanl(int n, long double _Complex z)

    :param n: The derivative number.
    :param z: The argument to the polygamma function.
    :return: The value of the polygamma function.

    Polygamma functions are defined as derivatives of the logarithmic gamma function. This function computes :math:`\psi^(n)(z) = \frac{d^{n + 1}}{dz^{n + 1}} \ln\Gamma(z)`.

.. c:function:: double polygammanu(double nu, double z)
.. c:function:: float polygammanuf(float nu, float z)
.. c:function:: long double polygammanul(long double nu, long double z)
.. c:function:: double _Complex cpolygammanu(double _Complex nu, double _Complex z)
.. c:function:: float _Complex cpolygammanuf(float _Complex nu, float _Complex z)
.. c:function:: long double _Complex cpolygammanul(long double _Complex nu, long double _Complex z)

    :param nu: The generalized derivative number.
    :param z: The argument to the polygamma function.
    :return: The value of the polygamma function.

    Polygamma functions are defined as derivatives of the logarithmic gamma function. This function computes :math:`\psi^(\nu)(z)`, which is a continuation from the natural numbers to the complex numbers.

Beta Functions
--------------

.. c:function:: double beta(double a, double b)
.. c:function:: float betaf(float a, float b)
.. c:function:: long double betal(long double a, long double b)
.. c:function:: double _Complex cbeta(double _Complex a, double _Complex b)
.. c:function:: float _Complex cbetaf(float _Complex a, float _Complex b)
.. c:function:: long double _Complex cbetal(long double _Complex a, long double _Complex b)

    :param a,b: The exponents for the beta function.
    :return: The beta function on the parameters.

    Computes the beta function. This is defined as

    .. math::

        B(a, b) = \int_0^1 t^{a - 1} (t - 1)^{b - 1} dt = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a + b)}

    It is computed this way, but using log-gammas to reduce the chances of infinities.

.. c:function:: double betainc(double z, double a, double b)
.. c:function:: float betaincf(float z, float a, float b)
.. c:function:: long double betaincl(long double z, long double a, long double b)
.. c:function:: double _Complex cbetainc(double _Complex z, double _Complex a, double _Complex b)
.. c:function:: float _Complex cbetaincf(float _Complex z, float _Complex a, float _Complex b)
.. c:function:: long double _Complex cbetaincl(long double _Complex z, long double _Complex a, long double _Complex b)

    :param z: The upper bound to the beta integral.
    :param a,b: The exponents for the beta function.
    :return: The incomplete beta function on the parameters.

    Computes the incomplete beta function. This is defined as the following integral.

    .. math::

        B_z(a, b) = \int_0^z t^{a - 1} (t - 1)^{b - 1} dt


.. c:function:: double betareg(double z, double a, double b)
.. c:function:: float betaregf(float z, float a, float b)
.. c:function:: long double betaregl(long double z, long double a, long double b)
.. c:function:: double _Complex cbetareg(double _Complex z, double _Complex a, double _Complex b)
.. c:function:: float _Complex cbetaregf(float _Complex z, float _Complex a, float _Complex b)
.. c:function:: long double _Complex cbetaregl(long double _Complex z, long double _Complex a, long double _Complex b)

    :param z: The upper bound to the beta integral.
    :param a,b: The exponents for the beta function.
    :return: The incomplete beta function on the parameters.

    Computes the regularized incomplete beta function. This is defined as the following

    .. math::

        I_z(a, b) = \frac{B_z(a, b)}{B(a, b)}

Error Functions
---------------

Important in statistics due to their relation to normal distributions. The complex form is generally usefull outside of statistics due to Gaussian distributions being important basis functions for machine learning and quantum chemistry.

.. c:function:: double _Complex cerf(double _Complex x)
.. c:function:: float _Complex cerff(float _Complex x)
.. c:function:: long double _Complex cerfl(long double _Complex x)

    :param x: The distance from the average.
    :return: The error function.

    Computes the error function for complex arguments.

.. c:function:: double ierf(double x)
.. c:function:: float ierff(float x)
.. c:function:: long double ierfl(long double x)
.. c:function:: double _Complex cierf(double _Complex x)
.. c:function:: float _Complex cierff(float _Complex x)
.. c:function:: long double _Complex cierfl(long double _Complex x)

    :param x: The proportion of the population.
    :return: The inverse error function of the proportion.

    Computes the inverse of the error function.

.. c:function:: double ierfc(double x)
.. c:function:: float ierfcf(float x)
.. c:function:: long double ierfcl(long double x)
.. c:function:: double _Complex cierfc(double _Complex x)
.. c:function:: float _Complex cierfcf(float _Complex x)
.. c:function:: long double _Complex cierfcl(long double _Complex x)

    :param x: The proportion of the population.
    :return: The inverse error function complement of the proportion.

    Computes the inverse of the error function complement.

Other Integrals
---------------

These are various integrated functions.

Fresnel Integrals
^^^^^^^^^^^^^^^^^

.. c:function:: double fresnels(double x)
.. c:function:: float fresnelsf(float x)
.. c:function:: long double fresnelsl(long double x)
.. c:function:: double _Complex cfresnels(double _Complex x)
.. c:function:: float _Complex cfresnelsf(float _Complex x)
.. c:function:: long double _Complex cfresnelsl(long double _Complex x)

    :param x: The parameter for the integral.
    :return: The Fresnel sine integral.

    Computes the Fresnel sine integral :math:`S(x)`. This is defined as

    .. math::

        S(x) = \int_0^x \sin\left(\frac{\pi t^2}{2}\right) dt

    It is implemented as the sum

    .. math::

        S(x) = x^3 \sum_{k = 0}^\infty \frac{2^{-2k - 1}\pi^{2k + 1} \left(-z^4\right)^k}{(4 k + 3) (2k + 1)!}

.. c:function:: double fresnelc(double x)
.. c:function:: float fresnelcf(float x)
.. c:function:: long double fresnelcl(long double x)
.. c:function:: double _Complex cfresnelc(double _Complex x)
.. c:function:: float _Complex cfresnelcf(float _Complex x)
.. c:function:: long double _Complex cfresnelcl(long double _Complex x)

    :param x: The parameter for the integral.
    :return: The Fresnel cosine integral.

    Computes the Fresnel cosine integral :math:`C(x)`. This is defined as

    .. math::

        C(x) = \int_0^x \cos\left(\frac{\pi t^2}{2}\right) dt

    It is implemented as the sum

    .. math::

        C(x) = x^3 \sum_{k = 0}^\infty \frac{2^{-2k}\pi^{2k} \left(-z^4\right)^k}{(4 k + 1) (2k)!}

Exponential-type Integrals
^^^^^^^^^^^^^^^^^^^^^^^^^^

Not all of these are purely exponential, and some are trigonometric.

.. c:function:: double integrale(double nu, double x)
.. c:function:: float integralef(float nu, float x)
.. c:function:: long double integralel(long double nu, long double x)
.. c:function:: double _Complex cintegrale(double _Complex nu, double _Complex x)
.. c:function:: float _Complex cintegralef(float _Complex nu, float _Complex x)
.. c:function:: long double _Complex cintegralel(long double _Complex nu, long double _Complex x)

    :param nu: The exponent for the integral.
    :param x: The multiplier for the exponent.
    :return: The E integral.

    Computes :math:`E_\nu(x)`. This is the integral

    .. math::

        E_\nu(x) = \int_1^\infty \frac{e^{-zt}}{t^\nu} dt

.. c:function:: double integralei(double x)
.. c:function:: float integraleif(float x)
.. c:function:: long double integraleil(long double x)
.. c:function:: double _Complex cintegralei(double _Complex x)
.. c:function:: float _Complex cintegraleif(float _Complex x)
.. c:function:: long double _Complex cintegraleil(long double _Complex x)

    :param x: The multiplier for the exponent.
    :return: The Ei integral.

    Computes :math:`Ei(x)` which is a rather complicated integral.

.. c:function:: double logint(double x)
.. c:function:: float logintf(float x)
.. c:function:: long double logintl(long double x)
.. c:function:: double _Complex clogint(double _Complex x)
.. c:function:: float _Complex clogintf(float _Complex x)
.. c:function:: long double _Complex clogintl(long double _Complex x)

    :param x: The argument to the function.
    :return: The logarithmic integral.

    Computes the logarithmic integral :math:`li(x) = \int_0^x \frac{1}{\ln(t)} dt`.

.. c:function:: double sinint(double x)
.. c:function:: float sinintf(float x)
.. c:function:: long double sinintl(long double x)
.. c:function:: double _Complex csinint(double _Complex x)
.. c:function:: float _Complex csinintf(float _Complex x)
.. c:function:: long double _Complex csinintl(long double _Complex x)

    :param x: The argument to the function.
    :return: The sine integral.

    Computes the sine integral :math:`Si(x) = \int_0^x \frac{\sin(t)}{t} dt`.

.. c:function:: double cosint(double x)
.. c:function:: float cosintf(float x)
.. c:function:: long double cosintl(long double x)
.. c:function:: double _Complex ccosint(double _Complex x)
.. c:function:: float _Complex ccosintf(float _Complex x)
.. c:function:: long double _Complex ccosintl(long double _Complex x)

    :param x: The argument to the function.
    :return: The cosine integral.

    Computes the cosine integral :math:`Ci(x) = \int_0^x \frac{\cos(t)}{t} dt`.

.. c:function:: double sinhint(double x)
.. c:function:: float sinhintf(float x)
.. c:function:: long double sinhintl(long double x)
.. c:function:: double _Complex csinhint(double _Complex x)
.. c:function:: float _Complex csinhintf(float _Complex x)
.. c:function:: long double _Complex csinhintl(long double _Complex x)

    :param x: The argument to the function.
    :return: The hyperbolic sine integral.

    Computes the hyberbolic sine integral :math:`Shi(x) = \int_0^x \frac{\sinh(t)}{t} dt`.

.. c:function:: double coshint(double x)
.. c:function:: float coshintf(float x)
.. c:function:: long double coshintl(long double x)
.. c:function:: double _Complex ccoshint(double _Complex x)
.. c:function:: float _Complex ccoshintf(float _Complex x)
.. c:function:: long double _Complex ccoshintl(long double _Complex x)

    :param x: The argument to the function.
    :return: The hyperbolic cosine integral.

    Computes the hyperbolic cosine integral :math:`Chi(x) = \int_0^x \frac{\cosh(t)}{t} dt`.
