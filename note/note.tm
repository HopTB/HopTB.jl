<TeXmacs|2.1.2>

<style|<tuple|generic|invisible-multiply|reduced-margins|old-lengths|number-long-article>>

<\body>
  <\hide-preamble>
    <assign|braket|<macro|a|b|<around*|\<langle\>|<arg|a><mid|\|><arg|b>|\<rangle\>>>>

    <assign|average|<macro|a|<around*|\<langle\>|<arg|a>|\<rangle\>>>>

    <assign|bracket|<macro|a|b|c|<around*|\<langle\>|<arg|a><mid|\|><arg|b><mid|\|><arg|c>|\<rangle\>>>>

    <assign|bra|<macro|a|<around*|\<langle\>|<arg|a>|\|>>>

    <assign|ket|<macro|a|<around*|\||<arg|a>|\<rangle\>>>>
  </hide-preamble>

  <doc-data|<doc-title|Note for HopTB.jl>>

  <section|Convention>

  Since HopTB.jl is doing real material calculations, this unit will always
  be standard units. Especially, it it important to keep in mind that

  <\itemize-dot>
    <item><math|e> is positive elementary charge. Electrons have charge
    <math|-e>.

    <item><math|\<b-v\>> is band velocity, whose unit is
    <math|\<up-m\>\<cdot\>\<up-s\><rsup|-1>>.
    <math|\<nabla\><rsub|\<b-k\>>\<varepsilon\>> is <math|\<hbar\>\<b-v\>>.

    <item><math|\<omega\>> is frequency, whose unit is
    <math|\<up-s\><rsup|-1>>. <math|\<hbar\>\<omega\>> is of energy unit.
  </itemize-dot>

  <section|Gauge choice for NoTB>

  Although the observables are gauge invariant, their ingredients are not.
  For example, to calculate the Berry connection <math|\<b-A\>>, a gauge must
  be fixed. Our gauge choice is done at each point in the Brillouin zone and
  is only a local gauge. The gauge choice starts with arbitrary states at
  <math|\<b-k\><rsub|0>> point (which is convenient since we won't need to
  modify the states obtained by direct diagonalization). For a specific band,
  we collect all the bands that are degenerate everywhere around
  <math|\<b-k\><rsub|0>> and call them a band set. Especially, if the band is
  nondegenerate, the band set contains only one band. For arbitrary band
  <math|n> and band <math|m> in a band set, we choose a gauge such that
  <math|\<b-A\><rsub|n\<nocomma\>m><around*|(|\<b-k\><rsub|0>|)>=\<b-0\>>. To
  show this is possible, we start with a random smooth gauge
  <math|<ket|<wide|u|~><rsub|n\<b-k\>>>>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|\<b-A\>|~><rsub|n\<nocomma\>m>>|<cell|=>|<cell|\<mathi\><braket|<wide|u|~><rsub|n\<b-k\>>|\<nabla\><rsub|\<b-k\>><wide|u|~><rsub|m\<b-k\>>>.<eq-number>>>>>
  </eqnarray*>

  We then do a transformation <math|<ket|u<rsub|n\<b-k\>>>=<big|sum><rsub|m>U<rsub|m\<nocomma\>n><around*|(|\<b-k\>|)><ket|<wide|u|~><rsub|m\<b-k\>>>>,
  where <math|U<around*|(|\<b-k\>|)>> is a smooth unitary matrix with
  <math|U<around*|(|\<b-k\><rsub|0>|)>> being identity matrix. We obtain

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-A\>>|<cell|=>|<cell|U<rsup|\<dagger\>><wide|\<b-A\>|~>*U+\<mathi\>U<rsup|\<dag\>>\<nabla\><rsub|\<b-k\>>U,<eq-number>>>>>
  </eqnarray*>

  and especially at <math|\<b-k\><rsub|0>>,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-A\><around*|(|\<b-k\><rsub|0>|)>>|<cell|=>|<cell|<wide|\<b-A\>|~><around*|(|\<b-k\><rsub|0>|)>*+<around*|\<nobracket\>|\<mathi\>\<nabla\><rsub|\<b-k\>>U|\|><rsub|\<b-k\>=\<b-k\><rsub|0>>.<eq-number>>>>>
  </eqnarray*>

  Since <math|<wide|\<b-A\>|~>> is Hermitian, we can choose
  <math|U<around*|(|\<b-k\>|)>=\<mathe\><rsup|\<mathi\>*\<b-k\>\<cdot\><wide|\<b-A\>|~><around*|(|\<b-k\><rsub|0>|)>>>
  to make <math|\<b-A\><around*|(|\<b-k\><rsub|0>|)>=0>, arriving at our
  desired gauge.<math|>

  <section|Derivative of generalized eigenvalues and generalized
  eigenvectors>

  The generalized eigenvalue equation is

  <\eqnarray*>
    <tformat|<table|<row|<cell|H<rsub|\<b-k\>>v<rsub|n\<b-k\>>>|<cell|=>|<cell|\<varepsilon\><rsub|n\<b-k\>>S<rsub|\<b-k\>>v<rsub|n\<b-k\>><eq-number>>>>>
  </eqnarray*>

  where <math|\<b-k\>> is a multidimensional parameter, <math|H> is a
  Hermitian matrix and <math|S> is a positive definite matrix. Written as
  matrices, the generalized eigenvalue equation becomes

  <\eqnarray*>
    <tformat|<table|<row|<cell|H<rsub|\<b-k\>>*V<rsub|\<b-k\>>>|<cell|=>|<cell|S<rsub|\<b-k\>>V<rsub|\<b-k\>>E<rsub|\<b-k\>>,<eq-number><label|gee>>>>>
  </eqnarray*>

  where the columns of <math|V<rsub|\<b-k\>>> are <math|v<rsub|n\<b-k\>>> and
  the diagonal elements of <math|E<rsub|\<b-k\>>> are
  <math|\<varepsilon\><rsub|n\<b-k\>>>. The eigenvectors are normalized
  according to <math|V<rsup|\<dag\>>S*V=I>, where <math|I> is identity
  matrix. With <math|\<partial\><rsub|\<alpha\>>\<assign\>\<partial\><rsub|\<b-k\><rsup|\<alpha\>>>>
  and omitting the <math|\<b-k\>> index, we take the derivative of Eq.
  (<reference|gee>) once and then multiply <math|V<rsup|\<dag\>>> to the
  left, obtaining

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|(|\<partial\><rsub|\<alpha\>>E|)>-<around*|[|E,*D<rsup|\<alpha\>>|]>>|<cell|=>|<cell|<around*|[|\<partial\><rsub|\<alpha\>>H|]>-<around*|[|\<partial\><rsub|\<alpha\>>S|]>*E,>>|<row|<cell|D<rsup|\<alpha\>>>|<cell|\<assign\>>|<cell|V<rsup|\<dag\>>S\<partial\><rsub|\<alpha\>>V,>>|<row|<cell|<around*|[|O|]>>|<cell|\<assign\>>|<cell|V<rsup|\<dag\>>O*V.<eq-number>>>>>
  </eqnarray*>

  The diagonal elements of the above equation gives the first order
  derivative of eigenvalues and the off-diagonal elements of the above
  equation gives the off-diagonal terms of <math|D<rsup|\<alpha\>>>. The
  diagonal terms of <math|D<rsup|\<alpha\>>> depends on the gauge choice of
  the eigenvectors. Taking the second order derivative of Eq.
  (<reference|gee>) and then multiply <math|V<rsup|\<dag\>>> to the left, we
  obtain

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|(|\<partial\><rsub|\<beta\>>\<partial\><rsub|\<alpha\>>E|)>-<around*|[|E,F<rsup|\<alpha\>\<beta\>>|]>>|<cell|=>|<cell|<around*|[|\<partial\><rsub|\<beta\>>\<partial\><rsub|\<alpha\>>H|]>-<around*|[|\<partial\><rsub|\<beta\>>\<partial\><rsub|\<alpha\>>S|]>**E+<around*|{|<around*|[|\<partial\><rsub|\<alpha\>>H|]>D<rsup|\<beta\>>-<around*|[|\<partial\><rsub|\<alpha\>>S|]>D<rsup|\<beta\>>**E-<around*|[|\<partial\><rsub|\<alpha\>>S|]>*<around*|(|\<partial\><rsub|\<beta\>>*E|)>-D<rsup|\<alpha\>><around*|(|\<partial\><rsub|\<beta\>>*E|)>+\<alpha\>\<leftrightarrow\>\<beta\>|}>,>>|<row|<cell|F<rsup|\<alpha\>\<beta\>>>|<cell|=>|<cell|V<rsup|\<dag\>>S\<partial\><rsub|\<beta\>>\<partial\><rsub|\<alpha\>>V.>>>>
  </eqnarray*>

  Again, the diagonal elements of the above equation are second order
  derivative of the eigenvalues and the off-diagonal elements gives off
  diagonal elements of <math|F<rsup|\<alpha\>\<beta\>>>. The diagonal terms
  of <math|F<rsup|\<alpha\>\<beta\>>> depends on the gauge choice of the
  eigenvectors. In the above equation, <math|<around*|{|\<alpha\>\<leftrightarrow\>\<beta\>|}>>
  means switching the <math|\<alpha\>> and <math|\<beta\>> of the previous
  term in the bracket.

  <section|Berry curvature>

  Berry curvature is defined as\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-Omega\>>|<cell|=>|<cell|\<nabla\><rsub|\<b-k\>>\<times\>\<b-A\>,<eq-number>>>>>
  </eqnarray*>

  where <math|\<b-A\>> is Berry connection. Explicitly,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<Omega\><rsup|\<gamma\>>>|<cell|=>|<cell|<frac|1|2>\<varepsilon\><rsup|\<alpha\>\<beta\>\<gamma\>>\<Omega\><rsup|\<alpha\>\<beta\>>,<eq-number>>>>>
  </eqnarray*>

  where <math|\<Omega\><rsup|\<alpha\>\<beta\>>=\<partial\><rsub|\<b-k\>><rsup|\<alpha\>>A<rsup|\<beta\>>-\<partial\><rsub|\<b-k\>><rsup|\<beta\>>A<rsup|\<alpha\>>>.
  HopTB.jl calculates <math|\<Omega\><rsup|\<alpha\>\<beta\>>>.

  <section|Berry curvature dipole>

  The Berry curvature dipole contribution to second order photocurrent
  (charge current, not particle current) is

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<sigma\><rsup|\<alpha\>\<beta\>\<gamma\>><rsub|BCD><around*|(|-\<omega\>,\<omega\>|)>>|<cell|=>|<cell|-<frac|\<mathi\>e<rsup|3>|2\<hbar\><around*|(|\<hbar\>\<omega\>+\<mathi\>\<eta\>|)>><big|sum><rsub|n><big|int><rsub|BZ><frac|\<mathd\>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>>\<Omega\><rsub|n\<nocomma\>n><rsup|\<alpha\>\<beta\>>\<partial\><rsub|\<b-k\>><rsup|\<gamma\>>f<rsub|n>+<around*|(|\<beta\>\<leftrightarrow\>\<gamma\>,\<omega\>\<leftrightarrow\>-\<omega\>|)>,<eq-number>>>>>
  </eqnarray*>

  where <math|\<Omega\><rsup|\<alpha\>\<beta\>>=\<partial\><rsub|\<b-k\>><rsup|\<alpha\>>A<rsup|\<beta\>>-\<partial\><rsub|\<b-k\>><rsup|\<beta\>>A<rsup|\<alpha\>>>,
  <math|\<b-A\>> is Berry connection, <math|\<eta\>> accounts for scattering.
  HopTB.jl calculates the following dimensionless tensor at zero temperature

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<Lambda\><rsup|\<alpha\>\<beta\>\<gamma\>>>|<cell|=>|<cell|-<big|sum><rsub|n><big|int><rsub|BZ><frac|\<mathd\>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>>\<Omega\><rsub|n\<nocomma\>n><rsup|\<alpha\>\<beta\>>\<partial\><rsub|\<b-k\>><rsup|\<gamma\>>f<rsub|n>>>|<row|<cell|>|<cell|=>|<cell|\<hbar\><big|sum><rsub|n><big|int><rsub|BZ><frac|\<mathd\>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>>\<Omega\><rsub|n\<nocomma\>n><rsup|\<alpha\>\<beta\>>v<rsub|n><rsup|\<gamma\>>\<delta\><around*|(|\<varepsilon\><rsub|n>-\<mu\>|)>,<eq-number>>>>>
  </eqnarray*>

  where <math|\<mu\>> is the Fermi energy.
  <math|\<Lambda\><rsup|\<alpha\>\<beta\>\<gamma\>>> can be expressed as a
  Fermi surface integral

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<Lambda\><rsup|\<alpha\>\<beta\>\<gamma\>>>|<cell|=>|<cell|<big|sum><rsub|n><big|int><rsub|FS<rsub|n>><frac|\<mathd\>\<sigma\>|<around*|(|2\<mathpi\>|)><rsup|3>>\<Omega\><rsub|n\<nocomma\>n><rsup|\<alpha\>\<beta\>><frac|v<rsub|n><rsup|\<gamma\>>|<around*|\||\<b-v\><rsub|n>|\|>>.<eq-number>>>>>
  </eqnarray*>

  <section|Drude weight>

  The optical Drude conductivity is defined as

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<sigma\><rsub|Drude><rsup|\<alpha\>\<beta\>><around*|(|\<omega\>|)>>|<cell|=>|<cell|-e<rsup|2>\<hbar\><big|sum><rsub|n><big|int><rsub|BZ><frac|\<mathd\>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>><frac|\<mathi\>v<rsup|\<alpha\>><rsub|n\<nocomma\>\<b-k\>>v<rsub|n\<b-k\>><rsup|\<beta\>>|\<hbar\>\<omega\>+\<mathi\>\<eta\>><frac|\<partial\>f<rsub|n\<b-k\>>|\<partial\>\<varepsilon\><rsub|n\<b-k\>>>,<eq-number>>>>>
  </eqnarray*>

  where <math|\<eta\>> accounts for scattering. The optical Drude
  conductivity is part of the contribution of the conductivity for metals
  defined as

  <\eqnarray*>
    <tformat|<table|<row|<cell|J<rsup|\<alpha\>><around*|(|\<omega\>|)>>|<cell|=>|<cell|\<sigma\><rsup|\<alpha\>\<beta\>><around*|(|\<omega\>|)>E<rsup|\<beta\>><around*|(|\<omega\>|)>,<eq-number>>>>>
  </eqnarray*>

  where <math|\<b-J\>> is charge current density (the electron has charge
  <math|-e> and this minus sign has been included) and <math|\<b-E\>> is the
  electric field. At zero frequency,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<sigma\><rsub|Drude><rsup|\<alpha\>\<beta\>><around*|(|0|)>>|<cell|=>|<cell|-<frac|e<rsup|2>\<hbar\>|\<eta\>><big|sum><rsub|n><big|int><rsub|BZ><frac|\<mathd\>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>><frac|\<partial\>f<rsub|n\<b-k\>>|\<partial\>\<varepsilon\><rsub|n\<b-k\>>>v<rsup|\<alpha\>><rsub|n\<nocomma\>\<b-k\>>v<rsub|n\<b-k\>><rsup|\<beta\>>.>>|<row|<cell|>|<cell|=>|<cell|<frac|D<rsup|\<alpha\>\<beta\>>|\<eta\>>,<eq-number>>>>>
  </eqnarray*>

  where\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|D<rsup|\<alpha\>\<beta\>>>|<cell|=>|<cell|-e<rsup|2>\<hbar\><big|sum><rsub|n><big|int><rsub|BZ><frac|\<mathd\>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>><frac|\<partial\>f<rsub|n\<b-k\>>|\<partial\>\<varepsilon\><rsub|n\<b-k\>>>v<rsup|\<alpha\>><rsub|n\<nocomma\>\<b-k\>>v<rsub|n\<b-k\>><rsup|\<beta\>><eq-number>>>>>
  </eqnarray*>

  is Drude weight.

  <section|OpenMX: the ordering of basis>

  The ordering of basis for OpenMX is

  <\itemize>
    <item><math|s> orbitals: <math|s>;

    <item><math|p> orbitals: <math|p<rsub|x>>, <math|p<rsub|y>>,
    <math|p<rsub|z>>;

    <item><math|d> orbitals: <math|d<rsub|z<rsup|2>>>,
    <math|d<rsub|x<rsup|2>-y<rsup|2>>>, <math|d<rsub|x\<nocomma\>y>>,
    <math|d<rsub|x\<nocomma\>z>>, <math|d<rsub|y\<nocomma\>z>>;

    <item><math|f> orbitals: <math|f<rsub|z<rsup|3>>>,
    <math|f<rsub|x\<nocomma\>z<rsup|2>>>,
    <math|f<rsub|y\<nocomma\>z<rsup|2>>>,
    <math|f<rsub|z<around*|(|x<rsup|2>-y<rsup|2>|)>\<nocomma\>>>,
    <math|f<rsub|x\<nocomma\>y\<nocomma\>z>>,
    <math|f<rsub|x<around*|(|x<rsup|2>-3y<rsup|2>|)>>>,
    <math|f<rsub|y<around*|(|3\<nocomma\>x<rsup|2>-y<rsup|2>|)>>>.
  </itemize>

  In <verbatim|group.jl>, <verbatim|Us_openmx> contains <math|U<rsub|l>>
  matrix, defined by

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|[|U<rsub|l>|]><rsub|n\<nocomma\>m>>|<cell|=>|<cell|<braket|Y<rsub|l><rsup|n>|X<rsub|l><rsup|m>>.<eq-number>>>>>
  </eqnarray*>

  Here, <math|Y<rsub|l><rsup|n>> is the complex spherical harmonics in the
  order of decreasing <math|n>; <math|X<rsub|l><rsup|m>> is the real
  spherical harmonics in the order of the OpenMX convention.

  <section|Wannier90: the ordering of basis>

  The ordering of basis for Wannier90 is

  <\itemize>
    <item><math|s> orbitals: <math|s>;

    <item><math|p> orbitals: <math|p<rsub|z>>, <math|p<rsub|x>>,
    <math|p<rsub|y>>;

    <item><math|d> orbitals: <math|d<rsub|z<rsup|2>>>,
    <math|d<rsub|x\<nocomma\>z>>, <math|d<rsub|y\<nocomma\>z>>,
    <math|d<rsub|x\<nocomma\><rsup|2>-y<rsup|2>>>,
    <math|d<rsub|x\<nocomma\>y>>;

    <item><math|f> orbitals: <math|f<rsub|z<rsup|3>>>,
    <math|f<rsub|x\<nocomma\>z<rsup|2>>>,
    <math|f<rsub|y\<nocomma\>z<rsup|2>>>,
    <math|f<rsub|z<around*|(|x<rsup|2>-y<rsup|2>|)>\<nocomma\>>>,
    <math|f<rsub|x\<nocomma\>y\<nocomma\>z>>,
    <math|f<rsub|x<around*|(|x<rsup|2>-3y<rsup|2>|)>>>,
    <math|f<rsub|y<around*|(|3\<nocomma\>x<rsup|2>-y<rsup|2>|)>>>.
  </itemize>

  In <verbatim|group.jl>, <verbatim|Us_wann> contains <math|U<rsub|l>>
  matrix, defined by

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|[|U<rsub|l>|]><rsub|n\<nocomma\>m>>|<cell|=>|<cell|<braket|Y<rsub|l><rsup|n>|X<rsub|l><rsup|m>>.<eq-number>>>>>
  </eqnarray*>

  Here, <math|Y<rsub|l><rsup|n>> is the complex spherical harmonics in the
  order of decreasing <math|n>; <math|X<rsub|l><rsup|m>> is the real
  spherical harmonics in the order of the Wannier90 convention.

  <section|Second order instrinsic nonlinear conductivity>

  The intrinsic nonlinear conductivity is defined as

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<sigma\><rsup|\<alpha\>\<beta\>\<gamma\>>>|<cell|=>|<cell|2e<rsup|3><big|sum><rsub|n,m><rsup|\<varepsilon\><rsub|m>\<neq\>\<varepsilon\><rsub|n>><big|int><frac|\<mathd\>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>>Re<around*|[|<frac|v<rsup|\<alpha\>><rsub|n>A<rsup|\<beta\>><rsub|n\<nocomma\>m>A<rsup|\<gamma\>><rsub|m\<nocomma\>n>-v<rsup|\<beta\>><rsub|n>A<rsup|\<alpha\>><rsub|n\<nocomma\>m>A<rsup|\<gamma\>><rsub|m\<nocomma\>n>|\<varepsilon\><rsub|n>-\<varepsilon\><rsub|m>>|]>f<rprime|'><around*|(|\<varepsilon\><rsub|n>|)>>>>>
  </eqnarray*>

  where

  <\eqnarray*>
    <tformat|<table|<row|<cell|f<around*|(|E|)>>|<cell|=>|<cell|<frac|1|\<mathe\><rsup|\<beta\><around*|(|E-\<mu\>|)>>+1><eq-number>>>>>
  </eqnarray*>

  \ is Fermi-Dirac distribution, <math|\<b-A\>> is Berry connection,
  <math|\<varepsilon\>> is band energy. At zero temperature,
  <math|f<rprime|'><around*|(|E|)>=-\<delta\><around*|(|E-\<mu\>|)>>, and
  therefore (<math|FS<rsub|n>> is the Fermi surface of band <math|n>)

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<sigma\><rsup|\<alpha\>\<beta\>\<gamma\>>>|<cell|=>|<cell|-2e<rsup|3><big|sum><rsub|n,m><rsup|\<varepsilon\><rsub|m>\<neq\>\<varepsilon\><rsub|n>><big|int><frac|\<mathd\>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>>Re<around*|[|<frac|v<rsup|\<alpha\>><rsub|n>A<rsup|\<beta\>><rsub|n\<nocomma\>m>A<rsup|\<gamma\>><rsub|m\<nocomma\>n>-v<rsup|\<beta\>><rsub|n>A<rsup|\<alpha\>><rsub|n\<nocomma\>m>A<rsup|\<gamma\>><rsub|m\<nocomma\>n>|\<varepsilon\><rsub|n>-\<varepsilon\><rsub|m>>|]>\<delta\><around*|(|\<varepsilon\><rsub|n>-\<mu\>|)>>>|<row|<cell|>|<cell|=>|<cell|-<frac|2e<rsup|3>|\<hbar\>><big|sum><rsub|n><big|int><rsub|FS<rsub|n>><frac|\<mathd\>\<sigma\>|<around*|(|2\<mathpi\>|)><rsup|3>><big|sum><rsub|m><rsup|\<varepsilon\><rsub|m>\<neq\>\<varepsilon\><rsub|n>>Re<around*|[|<frac|v<rsup|\<alpha\>><rsub|n>A<rsup|\<beta\>><rsub|n\<nocomma\>m>A<rsup|\<gamma\>><rsub|m\<nocomma\>n>-v<rsup|\<beta\>><rsub|n>A<rsup|\<alpha\>><rsub|n\<nocomma\>m>A<rsup|\<gamma\>><rsub|m\<nocomma\>n>|\<varepsilon\><rsub|n>-\<varepsilon\><rsub|m>>|]><frac|1|<around*|\||\<b-v\><rsub|n>|\|>>.<eq-number>>>>>
  </eqnarray*>

  <section|Second order Drude weight>

  Second order Drude conductivity is defined by

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<sigma\><rsup|\<alpha\>\<beta\>\<gamma\>>>|<cell|=>|<cell|-<frac|e<rsup|3>\<tau\><rsup|2>|\<hbar\><rsup|3>><big|sum><rsub|n><big|int><frac|\<mathd\><rsup|3>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>><around*|(|\<partial\><rsub|\<alpha\>>\<partial\><rsub|\<beta\>>\<partial\><rsub|\<gamma\>>\<varepsilon\><rsub|n>|)>f<around*|(|\<varepsilon\><rsub|n>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|e<rsup|3>\<tau\><rsup|2>|\<hbar\><rsup|3>><big|sum><rsub|n><big|int><frac|\<mathd\><rsup|3>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>><around*|(|\<partial\><rsub|\<alpha\>>\<varepsilon\><rsub|n>|)><around*|(|\<partial\><rsub|\<beta\>>\<partial\><rsub|\<gamma\>>\<varepsilon\><rsub|n>|)>f<rprime|'><around*|(|\<varepsilon\><rsub|n>|)><eq-number>>>>>
  </eqnarray*>

  where <math|\<varepsilon\><rsub|n>> is the band velocity, \ <math|f> is the
  Fermi-Dirac distribution, <math|\<tau\>> is the scattering lifetime. At
  zero temperature

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<sigma\><rsup|\<alpha\>\<beta\>\<gamma\>>>|<cell|=>|<cell|-<frac|e<rsup|3>\<tau\><rsup|2>|\<hbar\><rsup|3>><big|sum><rsub|n><big|int><frac|\<mathd\><rsup|3>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>><around*|(|\<partial\><rsub|\<alpha\>>\<varepsilon\><rsub|n>|)><around*|(|\<partial\><rsub|\<beta\>>\<partial\><rsub|\<gamma\>>\<varepsilon\><rsub|n>|)>\<delta\><around*|(|\<varepsilon\><rsub|n>-\<mu\>|)>,>>|<row|<cell|>|<cell|=>|<cell|-<frac|e<rsup|3>\<tau\><rsup|2>|\<hbar\><rsup|3>><big|sum><rsub|n><big|int><rsub|FS<rsub|n>><frac|\<mathd\>\<sigma\>|<around*|(|2\<mathpi\>|)><rsup|3>><around*|(|\<partial\><rsub|\<alpha\>>\<varepsilon\><rsub|n>|)><around*|(|\<partial\><rsub|\<beta\>>\<partial\><rsub|\<gamma\>>\<varepsilon\><rsub|n>|)><frac|1|<around*|\||\<nabla\>\<varepsilon\><rsub|n>|\|>>.<eq-number>>>>>
  </eqnarray*>

  HopTB.jl calculates the following quantity

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<Lambda\><rsup|\<alpha\>\<beta\>\<gamma\>>>|<cell|=>|<cell|-<frac|e<rsup|3>|\<hbar\>><big|sum><rsub|n><big|int><rsub|FS<rsub|n>><frac|\<mathd\>\<sigma\>|<around*|(|2\<mathpi\>|)><rsup|3>><around*|(|\<partial\><rsub|\<alpha\>>\<varepsilon\><rsub|n>|)><around*|(|\<partial\><rsub|\<beta\>>\<partial\><rsub|\<gamma\>>\<varepsilon\><rsub|n>|)><frac|1|<around*|\||\<nabla\>\<varepsilon\><rsub|n>|\|>>.<eq-number>>>>>
  </eqnarray*>

  <section|Injection current>

  The injection current conductivity is

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<sigma\><rsup|\<alpha\>\<beta\>\<gamma\>><around*|(|\<omega\>|)>>|<cell|=>|<cell|<frac|\<mathpi\>e<rsup|3>|\<eta\>><big|sum><rsub|n,m><big|int><frac|\<mathd\><rsup|3>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>><around*|(|v<rsup|\<alpha\>><rsub|n>-v<rsub|m><rsup|\<alpha\>>|)>A<rsup|\<beta\>><rsub|n\<nocomma\>m>A<rsup|\<gamma\>><rsub|m\<nocomma\>n>f<rsub|n\<nocomma\>m>\<delta\><around*|(|\<hbar\>\<omega\>-\<varepsilon\><rsub|m\<nocomma\>n>|)>.>>>>
  </eqnarray*>

  HopTB.jl calculates

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<Lambda\><rsup|\<alpha\>\<beta\>\<gamma\>><around*|(|\<omega\>|)>>|<cell|=>|<cell|\<mathpi\>e<rsup|3><big|sum><rsub|n,m><big|int><frac|\<mathd\><rsup|3>\<b-k\>|<around*|(|2\<mathpi\>|)><rsup|3>><around*|(|v<rsup|\<alpha\>><rsub|n>-v<rsub|m><rsup|\<alpha\>>|)>A<rsup|\<beta\>><rsub|n\<nocomma\>m>A<rsup|\<gamma\>><rsub|m\<nocomma\>n>f<rsub|n\<nocomma\>m>\<delta\><around*|(|\<hbar\>\<omega\>-\<varepsilon\><rsub|m\<nocomma\>n>|)>.>>>>
  </eqnarray*>

  For nonmagnetic materials, <math|\<Lambda\><around*|(|\<omega\>|)>> is
  purely imaginary. In addition, for nonmagnetic materials,
  <math|\<Lambda\><rsup|\<alpha\>\<beta\>\<gamma\>><around*|(|\<omega\>|)>=-\<Lambda\><rsup|\<alpha\>\<gamma\>\<beta\>><around*|(|\<omega\>|)>>.
</body>

<\initial>
  <\collection>
    <associate|font|stix>
    <associate|font-family|rm>
    <associate|math-font|math-stix>
    <associate|page-medium|automatic>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-10|<tuple|10|?>>
    <associate|auto-11|<tuple|11|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|3|?>>
    <associate|auto-4|<tuple|4|?>>
    <associate|auto-5|<tuple|5|?>>
    <associate|auto-6|<tuple|6|?>>
    <associate|auto-7|<tuple|7|?>>
    <associate|auto-8|<tuple|8|?>>
    <associate|auto-9|<tuple|9|?>>
    <associate|gee|<tuple|3.2|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Convention>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Gauge
      choice for NoTB> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Derivative
      of generalized eigenvalues and generalized eigenvectors>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Berry
      curvature> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Berry
      curvature dipole> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Drude
      weight> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|7<space|2spc>OpenMX:
      the ordering of basis> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|8<space|2spc>Wannier90:
      the ordering of basis> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|9<space|2spc>Second
      order instrinsic nonlinear conductivity>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|10<space|2spc>Second
      order Drude weight> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>