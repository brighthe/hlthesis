\subsection{Management of Degrees of Freedom for BDM and N\'ed\'elec Finite Element Spaces}

One of the focal points in this section is the intricacies involved in the
management of DoFs, primarily revolving around ensuring the continuity of the
finite element space.

Both BDM and N\'ed\'elec are distinguished as vector finite element spaces.
Their unique approach involves designating DoFs of the vector type by
stipulating a vector frame for each interpolation point. Concurrently, the
formulation of their vector basis functions arises from a fusion of the Lagrange
basis function and the DoFs' dual frame.

To offer a different perspective, each DoF within the BDM or N\'ed\'elec space
correlates to a distinct interpolation point, denoted as \(p\), and an exclusive
vector, \(\boldsymbol e\). Furthermore, each basis function aligns with a
specific Lagrange basis function defined on \(p\) and a vector \(\boldsymbol
e'\), recognized as the dual counterpart of \(\boldsymbol e\).

Fundamentally, managing DoFs can be perceived as a systematic counting endeavor.
Initially, our task is to establish both global and local indexing protocols for
all existing DoFs.

Globally categorizing, DoFs can either be shared among simplices or remain
exclusive. Shared DoFs, depending on the dimensionality of their residing
sub-simplex, can be classified further as on-edge or on-face. It's worth noting
that in the realms of BDM and N\'ed\'elec spaces, there's a void of DoFs shared
on nodes. Moreover, in a 3D BDM space, shared DoFs on edges are absent. Hence,
the global numbering pattern echoes that of the Lagrange interpolation points:
initiate by counting shared DoFs on each edge based on edge order, progress to
count shared DoFs on each face by face order, and culminate with counting the
unshared DoFs within the cell. For each edge or face, the sequence of DoFs can
seamlessly align with the interpolation points' order.

Given the global indexing rules, one can derive an array named
\lstinline{dof2vector} of dimensions \lstinline{(gdof, GD)}, where
\lstinline{gdof} represents the global DoFs count, and \lstinline{GD} signifies
geometry dimensions. In this array, \lstinline{dof2vector[i, :]} encapsulates
the vector linked to the \(i\)-th DoF.

Our subsequent task involves defining local indexing norms and synthesizing an
array named \lstinline{cell2dof} with dimensions \lstinline{(NC,ldof)}, where
\lstinline{ldof} equates to the local DoFs count per cell. It's imperative to
realize that each DoF stems from a combination of an interpolation point and a
vector. With every interpolation point, there's an associated frame, housing
\lstinline{GD} vectors. For a DoF on the \(i\)-th cell, let the local index of
its interpolation point be \(p\) and its vector in the frame be \(q\). Such a
DoF can be represented by a unique local index \(j\), calculated as:
$$
j = n\cdot q + p
$$
where \(n\) is the tally of interpolation points in the \(i\)-th cell. Expanding
on this, the value of \lstinline{cell2dof[i,j]} can be determined using the
global index \lstinline{cell2ipoint[i, p]}, the sub-simplex location of the
interpolation point, and the global indexing standards.

\begin{remark}
It's pivotal to recognize that neither the local nor the global numbering
protocols elucidated above are the sole possible configurations. Additionally,
armed with the \lstinline{cell2dof} array, executing the higher-order finite
element methodologies mentioned in this paper differs negligibly from the
traditional finite element method, especially in terms of matrix-vector assembly
and handling boundary conditions.
\end{remark}

