# compass
## A robust header-only curve fitting library.
<p align="center">
  <a href="#summary">Summary</a> •
  <a href="#installation">Installation</a> •
  <a href="#features">Features</a> •
  <a href="#documentation">Documentation</a> •
  <a href="#references">References</a> •
</p>

---

## Summary
**compass** is the first comprehensive open source circle-fitting library written in C++. 

Using an easily extensible interface based on the [curiously recurring template pattern](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern), **compass** implements (insert number) numerical algorithms for solving the canonical problem of fitting circles to given points in the plane. Of the algorithms included in the library, most contain multiple implementations so as to optimize for either numerical stability or computational efficiency.

**compass** is designed to be easily integrated into larger projects. To that end, its only dependecy is [Eigen 3.3](http://eigen.tuxfamily.org). Further, **compass** can be collapsed into a singular header file using the included Python script, though, admittedly, the benefits of this are minimal considering the small scale and incomplexity of the codebase.

**compass** is in the early stages of its development and has not been exhaustively tested and optimized. As of the current release, it is not recommended for use in production.

The taxonomy of **compass** comprises three broad categories of algorithms &mdash; algebraic, geometric, and specialized &mdash; as well as several subcatagories which are further discussed in the [features](#features) section.

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;F\left&space;(&space;a,&space;b,&space;R&space;\right&space;)&space;=&space;\sum_{i=1}^{n}\left&space;|&space;\sqrt{(x_{i}&space;-&space;a)^{2}&space;&plus;&space;(y_{i}&space;-&space;b)^{2}}&space;-&space;R&space;\right&space;|" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\large&space;F\left&space;(&space;a,&space;b,&space;R&space;\right&space;)&space;=&space;\sum_{i=1}^{n}\left&space;|&space;\sqrt{(x_{i}&space;-&space;a)^{2}&space;&plus;&space;(y_{i}&space;-&space;b)^{2}}&space;-&space;R&space;\right&space;|" title="\large F\left ( a, b, R \right ) = \sum_{i=1}^{n}\left | \sqrt{(x_{i} - a)^{2} + (y_{i} - b)^{2}} - R \right|" /></a>
</p>

- Algebraic algorithms are non-iterative procedures that reparameterize <a href="https://www.codecogs.com/eqnedit.php?latex=F" target="_blank"><img src="https://latex.codecogs.com/svg.latex?F" title="F" /></a> and minimize the corresponding objective function by way of orthogonal least squares.
- Geometric algorithms, supplied with an initial guess (usually produced by an algorithm of the Algebraic variety), iteratively converge to a minimum of <a href="https://www.codecogs.com/eqnedit.php?latex=F" target="_blank"><img src="https://latex.codecogs.com/svg.latex?F" title="F" /></a>.
- Specialized algorithms are a group of sophisticated mathematical routines that bear little similarity with their geometric and algebraic counterparts as well as with one another. They employ 

## Features
**compass** implements a module-centric organizational structure to reduce the amount of #include directives needed to write a minimally functional program. To avoid compilation overhead, files can be included individually.

<table>
  <tr>
  <th></th>
    <th scope="col">Algebraic.hpp</th>
    <th scope="col">Geometric.hpp</th>
    <th scope="col">Specialized.hpp</th>
    <th scope="col">Compass.hpp</th>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/GanderGolubStrebel.hpp">GanderGolubStrebel</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Kasa.hpp">Kasa</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
    <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Kasa.hpp">KasaConsistent</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Hyper.hpp">HyperSVD</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>

  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Hyper.hpp">HyperSimple</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
   
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/KukushMarkovskyHuffel.hpp">KukushMarkovskyHuffel</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Nievergelt.hpp">Nievergelt</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Pratt.hpp">PrattSVD</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Pratt.hpp">PrattNewton</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Pratt.hpp">PrattRobust</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Taubin.hpp">TaubinSVD</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Taubin.hpp">TaubinNystromSVD</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Algebraic/Taubin.hpp">TaubinNewton</a></th>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Geometric/Landau.hpp">Landau</a></th>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Geometric/LevenbergMarquardt.hpp">LevenbergMarquardt</a></th>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Geometric/LevenbergMarquardt.hpp">LevenbergMarquardtReduced</a></th>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Geometric/Spath.hpp">Spath</a></th>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
  <tr>
    <th scope="row", align="left"><a href="Compass/src/Geometric/Trust.hpp">Trust</a></th>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Specialized/Riemann.hpp">Riemann</a></th>
    <td align="center">:x:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Specialized/Riemann.hpp">RiemannAlgebraic</a></th>
    <td align="center">:x:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Specialized/Inversion.hpp">InversionIterative</a></th>
    <td align="center">:x:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>
   <tr>
    <th scope="row", align="left"><a href="Compass/src/Specialized/Inversion.hpp">InversionNonIterative</a></th>
    <td align="center">:x:</td>
    <td align="center">:x:</td>
    <td align="center">:heavy_check_mark:</td>
    <td align="center">:heavy_check_mark:</td>
  </tr>

</table>

## Installation

## Documentation
