# IPO - Investigating Polyhedra by Oracles # {#mainpage}

\section description Description

<span style="color: #007bff; font-weight: bold;">I</span><span style="font-weight: bold;">P</span><span style="color: #ff0000; font-weight: bold;">O</span></span>
is a software library for certain problems in polyhedral combinatorics
such as computing the affine hull, checking adjacency, or computing facets.
In contrast to other tools, it only requires the polyhedra to be given
implicitly by means of so-called <em>optimization oracles</em>.

\section motivation Motivation

In mathematical optimization, it is common to encode certain
combinatorial objects (e.g., routes, schedules, decisions) as points in the Euclidean space
and consider their convex hulls, which are usually polyhedra.
Using this approach, properties concerning the facial structures of these polyhedra become of central interest.
In general, classical enumeration-based algorithms investigating such properties turn out to be impractical already for small dimensions, e.g., <em>n=20</em>.

On the other hand, maximizing linear objective functions over these polyhedra (though most often NP-hard)
can be done very efficiently for moderate sizes (say <em>n=100</em>), e.g., by mixed-integer programming solvers.
<span style="color: #007bff; font-weight: bold;">I</span><span style="font-weight: bold;">P</span><span style="color: #ff0000; font-weight: bold;">O</span></span>
uses such <em>optimization oracles</em> and (partially) solves some of the problems mentioned above.
For this it utilizes the LP solver
<a href="http://soplex.zib.de/">SoPlex</a>,
in particular its ability to compute solutions in exact arithmetic.
        
\section functionality Functionality
        
Given an optimization oracle defining a polyhedron <em>P</em>,
<span style="color: #007bff; font-weight: bold;">I</span><span style="font-weight: bold;">P</span><span style="color: #ff0000; font-weight: bold;">O</span></span>
can...
\li compute the <b>affine hull</b> of <em>P</em>
    in the sense that it returns <em>(d+1)</em> many affinely independent points in <em>P</em> and
    <em>(n-d)</em> many linearly independent equations valid for <em>P</em>, where <em>d</em> is <em>P</em> 's dimension.
    See AffineHull for details.
\li compute a <b>facet-defining inequality</b> that is <b>violated</b> by a given point.<br/>
    See Separation for details.
\li compute <b>facets that are "helpful"</b> when maximizing a specified objective vector <em>c</em>.
    More precisely, it iteratively solves LPs whose inequalities correspond to some of <em>P</em> 's facets,
    until the current optimum is in <em>P</em>.
    As long as this is not the case, the procedure returns violated facet-defining inequalities and adds them to the LP.
\li <b>check</b> whether two given vertices of <em>P</em> are <b>adjacent</b>.<br/>
    See SmallestFace for details.
\li compute the <b>smallest face</b> of <em>P</em> <b>containing a given point</b>
    by means of an inequality defining this face.<br/>
    See SmallestFace for details.
\li <b>check</b> whether a given point is a <b>vertex</b> of <em>P</em></b>.<br/>
    See SmallestFace for details.


\section oracles Oracles

<span style="color: #007bff; font-weight: bold;">I</span><span style="font-weight: bold;">P</span><span style="color: #ff0000; font-weight: bold;">O</span></span>'s
requirements for an optimization oracle are very simple:

\li Essentially, for a given objective vector, it only has to return
    a point in <em>P</em> of maximum objective value.
\li If the maximum is not attained, it must return an unbounded direction proving this.
    Of course, if <em>P</em> is empty, it must claim "infeasible".
\li There exist several ways to speed up certain algorithms, e.g., by returning additional solutions, implementing
    another heuristic oracle that is faster but does not guarantee optimality, etc.

\section technicalDetails Technical Details

<span style="color: #007bff; font-weight: bold;">I</span><span style="font-weight: bold;">P</span><span style="color: #ff0000; font-weight: bold;">O</span></span>
is implemented in C++ and uses exact arithmetic. 
It depends on the LP solver
<a href="http://soplex.zib.de/">SoPlex</a>.
In order to implement an optimization oracle,
there are some base classes from which one must inherit:

\li 
  OracleBase: 
  The most general base class.
\li 
  FaceOracleBase: 
  Inherit from this class if your oracle is not capable of optimizing over a face of <em>P</em> explicitly.

For mixed-integer programs, several oracles are already implemented:
\li 
  SCIPOracle:
  This oracle is built from a 
  <a href="http://scip.zib.de/">SCIP</a>
  instance which is used for optimization.
  It inherits from MIPOracle which cares about the correct handling of the floating-point solutions returned by SCIP.
<!--
\li
  ExactSCIPOracle:
  This oracle is built from a MixedIntegerProgram instances which explicitly represents a MIP and optionally
  another oracle which may act as a heuristic.
  If non-optimality is permitted, it runs the heuristic oracle,
  and otherwise
  <a href="http://scip.zib.de/#exact">SCIP-ex</a>.
-->!
\li
  ExternalOracle:
  This oracle uses an external program using a defined interface.

\section tools Command Line Tool

<span style="color: #007bff; font-weight: bold;">I</span><span style="font-weight: bold;">P</span><span style="color: #ff0000; font-weight: bold;">O</span></span>
comes with a tool that use
<a href="http://scip.zib.de/">SCIP</a>
as an oracle.
It can handle all polyhedra that SCIP understands, e.g., mixed-integer hulls of ZIMPL models.
Here are some sample invocations:
\li
  \c ipo --dimension: Computes the affine hull.
\li
  \c ipo \c --facets: Uses the objective function present in the instance to computes facets.
\li
  \c ipo \c --facets \c --random \c 1000: Additionally usese 1000 randomly samples objective functions for facet generation..
\li
  \c ipo \c --smallest-face \c --point \c '(x1=1)': Computes the smallest face that contains the unit vector with x1=1.

\section about About

<span style="color: #007bff; font-weight: bold;">I</span><span style="font-weight: bold;">P</span><span style="color: #ff0000; font-weight: bold;">O</span></span>
is an open source project.
You are free to use it in your commercial or
non-commercial applications under very permissive
\ref license "license terms".

Any publication of results to which using 
<span style="color: #007bff; font-weight: bold;">I</span><span style="font-weight: bold;">P</span><span style="color: #ff0000; font-weight: bold;">O</span></span>
contributed must cite this project.

\code
@MISC{Walter16,
  author = {Walter, Matthias},
  title = {\emph{IPO -- Investigating Polyhedra by Oracles}, \\ \url{http://polyhedra-oracles.bitbucket.org/}},
  year = {2016},
}
\endcode

The project is maintained by Matthias Walter.

\section getting Getting it to work

Pull from the <a href="http://bitbucket.org/polyhedra-oracles/ipo/">git repository</a>.
Note that IPO is not available as in a packaged form (i.e., a tar.gz archive) since it
is still under heavy development. This will change once it reaches a more stable state.

After cloning, please follow the instructions in the INSTALL file.

*/

}
