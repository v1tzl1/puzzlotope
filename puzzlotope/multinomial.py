"""

This file is taken from the sympy.ntheory.multinomial module from
http://docs.sympy.org/latest/_modules/sympy/ntheory/multinomial.html

The following license applies according to
http://docs.sympy.org/latest/aboutus.html#license


Copyright (c) 2006-2017 SymPy Development Team

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of SymPy nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

--------------------------------------------------------------------------------

Patches taken from the Diofant project (https://github.com/diofant/diofant)
are licensed as:

Copyright (c) 2006-2017 SymPy Development Team,
              2013-2017 Sergey B Kirpichev

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of Diofant, or SymPy nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

"""

def binomial_coefficients(n):
    """Return a dictionary containing pairs :math:`{(k1,k2) : C_kn}` where
    :math:`C_kn` are binomial coefficients and :math:`n=k1+k2`.
    Examples
    ========

    >>> from sympy.ntheory import binomial_coefficients
    >>> binomial_coefficients(9)
    {(0, 9): 1, (1, 8): 9, (2, 7): 36, (3, 6): 84,
     (4, 5): 126, (5, 4): 126, (6, 3): 84, (7, 2): 36, (8, 1): 9, (9, 0): 1}

    See Also
    ========

    binomial_coefficients_list, multinomial_coefficients
    """
    d = {(0, n): 1, (n, 0): 1}
    a = 1
    for k in range(1, n//2 + 1):
        a = (a * (n - k + 1))//k
        d[k, n - k] = d[n - k, k] = a
    return d



def binomial_coefficients_list(n):
    """ Return a list of binomial coefficients as rows of the Pascal's
    triangle.

    Examples
    ========

    >>> from sympy.ntheory import binomial_coefficients_list
    >>> binomial_coefficients_list(9)
    [1, 9, 36, 84, 126, 126, 84, 36, 9, 1]

    See Also
    ========

    binomial_coefficients, multinomial_coefficients
    """
    d = [1] * (n + 1)
    a = 1
    for k in range(1, n//2 + 1):
        a = (a * (n - k + 1))//k
        d[k] = d[n - k] = a
    return d



def multinomial_coefficients0(m, n, _tuple=tuple, _zip=zip):
    """Return a dictionary containing pairs ``{(k1,k2,..,km) : C_kn}``
    where ``C_kn`` are multinomial coefficients such that
    ``n=k1+k2+..+km``.

    For example:

    >>> from sympy import multinomial_coefficients
    >>> multinomial_coefficients(2, 5) # indirect doctest
    {(0, 5): 1, (1, 4): 5, (2, 3): 10, (3, 2): 10, (4, 1): 5, (5, 0): 1}

    The algorithm is based on the following result:

       Consider a polynomial and its ``n``-th exponent::

         P(x) = sum_{i=0}^m p_i x^i
         P(x)^n = sum_{k=0}^{m n} a(n,k) x^k

       The coefficients ``a(n,k)`` can be computed using the
       J.C.P. Miller Pure Recurrence [see D.E.Knuth, Seminumerical
       Algorithms, The art of Computer Programming v.2, Addison
       Wesley, Reading, 1981;]::

         a(n,k) = 1/(k p_0) sum_{i=1}^m p_i ((n+1)i-k) a(n,k-i),

       where ``a(n,0) = p_0^n``.
    """

    if not m:
        if n:
            return {}
        return {(): 1}
    if m == 2:
        return binomial_coefficients(n)
    symbols = [(0,)*i + (1,) + (0,)*(m - i - 1) for i in range(m)]
    s0 = symbols[0]
    p0 = [_tuple(aa - bb for aa, bb in _zip(s, s0)) for s in symbols]
    r = {_tuple(aa*n for aa in s0): 1}
    l = [0] * (n*(m - 1) + 1)
    l[0] = r.items()
    for k in range(1, n*(m - 1) + 1):
        d = defaultdict(int)
        for i in range(1, min(m, k + 1)):
            nn = (n + 1)*i - k
            if not nn:
                continue
            t = p0[i]
            for t2, c2 in l[k - i]:
                tt = _tuple([aa + bb for aa, bb in _zip(t2, t)])
                d[tt] += nn*c2
                if not d[tt]:
                    del d[tt]
        r1 = [(t, c//k) for (t, c) in d.items()]
        l[k] = r1
        r.update(r1)
    return r


def multinomial_coefficients(m, n):
    r"""Return a dictionary containing pairs ``{(k1,k2,..,km) : C_kn}``
    where ``C_kn`` are multinomial coefficients such that
    ``n=k1+k2+..+km``.

    For example:

    >>> from sympy.ntheory import multinomial_coefficients
    >>> multinomial_coefficients(2, 5) # indirect doctest
    {(0, 5): 1, (1, 4): 5, (2, 3): 10, (3, 2): 10, (4, 1): 5, (5, 0): 1}

    The algorithm is based on the following result:

    .. math::
        \binom{n}{k_1, \ldots, k_m} =
        \frac{k_1 + 1}{n - k_1} \sum_{i=2}^m \binom{n}{k_1 + 1, \ldots, k_i - 1, \ldots}

    Code contributed to Sage by Yann Laigle-Chapuy, copied with permission
    of the author.

    See Also
    ========

    binomial_coefficients_list, binomial_coefficients
    """
    if not m:
        if n:
            return {}
        return {(): 1}
    if m == 2:
        return binomial_coefficients(n)
    if m >= 2*n and n > 1:
        return dict(multinomial_coefficients_iterator(m, n))
    t = [n] + [0] * (m - 1)
    r = {tuple(t): 1}
    if n:
        j = 0  # j will be the leftmost nonzero position
    else:
        j = m
    # enumerate tuples in co-lex order
    while j < m - 1:
        # compute next tuple
        tj = t[j]
        if j:
            t[j] = 0
            t[0] = tj
        if tj > 1:
            t[j + 1] += 1
            j = 0
            start = 1
            v = 0
        else:
            j += 1
            start = j + 1
            v = r[tuple(t)]
            t[j] += 1
        # compute the value
        # NB: the initialization of v was done above
        for k in range(start, m):
            if t[k]:
                t[k] -= 1
                v += r[tuple(t)]
                t[k] += 1
        t[0] -= 1
        r[tuple(t)] = (v * tj) // (n - t[0])
    return r



def multinomial_coefficients_iterator(m, n, _tuple=tuple):
    """multinomial coefficient iterator

    This routine has been optimized for `m` large with respect to `n` by taking
    advantage of the fact that when the monomial tuples `t` are stripped of
    zeros, their coefficient is the same as that of the monomial tuples from
    ``multinomial_coefficients(n, n)``. Therefore, the latter coefficients are
    precomputed to save memory and time.

    >>> from sympy.ntheory.multinomial import multinomial_coefficients
    >>> m53, m33 = multinomial_coefficients(5,3), multinomial_coefficients(3,3)
    >>> m53[(0,0,0,1,2)] == m53[(0,0,1,0,2)] == m53[(1,0,2,0,0)] == m33[(0,1,2)]
    True

    Examples
    ========

    >>> from sympy.ntheory.multinomial import multinomial_coefficients_iterator
    >>> it = multinomial_coefficients_iterator(20,3)
    >>> next(it)
    ((3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 1)
    """
    if m < 2*n or n == 1:
        mc = multinomial_coefficients(m, n)
        for k, v in mc.items():
            yield(k, v)
    else:
        mc = multinomial_coefficients(n, n)
        mc1 = {}
        for k, v in mc.items():
            mc1[_tuple(filter(None, k))] = v
        mc = mc1

        t = [n] + [0] * (m - 1)
        t1 = _tuple(t)
        b = _tuple(filter(None, t1))
        yield (t1, mc[b])
        if n:
            j = 0  # j will be the leftmost nonzero position
        else:
            j = m
        # enumerate tuples in co-lex order
        while j < m - 1:
            # compute next tuple
            tj = t[j]
            if j:
                t[j] = 0
                t[0] = tj
            if tj > 1:
                t[j + 1] += 1
                j = 0
            else:
                j += 1
                t[j] += 1

            t[0] -= 1
            t1 = _tuple(t)
            b = _tuple(filter(None, t1))
            yield (t1, mc[b])
