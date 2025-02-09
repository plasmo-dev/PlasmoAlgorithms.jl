# Solution Algorithm for PlasmoSchwarz.jl

Schwarz Decomposition is an iterative algorithm takes advantage of Plasmo.jl's subgraph structuring capabilities by partitioning an optimization problem into separate subgraphs. The subgraphs can then be expanded to overlap with other subgraphs, and a single iteration of the algorithm solves each subgraph, which can be done in parallel. After optimization, solutions around the overlapping regions are shared between subgraphs and used within the algorithm to help the subgraphs converge to the solution of the overall problem. Schwarz decomposition does not support integer variables. Further, in practice, convergence of PlasmoSchwarz is dependent on problem formulation and the degree of overlap between problems. Convergence is not guaranteed, but primal and dual values can be used to detect whether a local solution has been reached. 

## Schwarz Decomposition

Schwarz Decomposition can be applied to graph-structured problems defined in Plasmo.jl The overall graph structure in Plasmo.jl has the form

```math
\begin{align*}
    \min&\;  \sum_{n \in \mathcal{N}(\mathcal{G})} f_n (x_n) \\
    \textrm{s.t.} &\; x_n \in \mathcal{X}_n, n \in \mathcal{N}(\mathcal{G}) \\
    &\; g_e(\{x_n \}_{n \in \mathcal{N}(e)}) = 0, e \in \mathcal{E}(\mathcal{G})
\end{align*}
```

Here, $x_n$ are the variables stored on node $n$, $\mathcal{N}(\mathcal{G})$ represents the set of nodes on graph $\mathcal{G}$, $\mathcal{E}(\mathcal{G})$ represents the set of edges on graph $\mathcal{G}$, and $\mathcal{N}(E)$ are the nodes connected by edge $e$. Thus, the first constraint set are the local constraints (stored on each node) and the second constraint set are the linking constraints for constraints on two or more nodes. 

Further details of the algorithm can be found in the original Plasmo.jl paper [here](https://arxiv.org/abs/2006.05378), but we will give an overview of the algorithm here. A graph $\mathcal{G}$ can be partitioned into $N$ subgraphs notated by $\{ \mathcal{SG}_i \}^N_{i=1}$. These subgraphs can be _expanded_ (meaning that they will now include nodes within a prescribed distance of nodes in the original subgraph) to obtain _expanded subgraphs_ notated by $\{ \mathcal{SG}_i' \}^N_{i=1}$. Under Schwarz Decomposition, the constraints on the incident edges of the expanded subgraphs can be enforced using the previous iteration's solutions from the adjacent subgraphs. The formulation for a given subgraph, $\mathcal{SG}_i'$ is given by:

```math
\begin{align*}
    \min&\;  \sum_{n \in \mathcal{N}(\mathcal{SG}_i')} f_n (x_n) - \sum_{e \in \mathcal{I}_1(\mathcal{SG}_i')} (\lambda_e^k)^\top g_e(\{x_n \}_{n \in \mathcal{N}(e) \cap\mathcal{N}(\mathcal{SG}_i')}, \{x_n^k \}_{n\in \mathcal{N}(e) \setminus \mathcal{N}(\mathcal{SG}_i')}) \\
    \textrm{s.t.} &\; x_n \in \mathcal{X}_n, n \in \mathcal{N}(\mathcal{SG}_i') \\
    &\; g_e(\{x_n \}_{n \in \mathcal{N}(e)}) = 0, e \in \mathcal{E}(\mathcal{SG}_i') \\
    &\; g_e(\{x_n \}_{n \in \mathcal{N}(e) \cap\mathcal{N}(\mathcal{SG}_i')}, \{x_n^k \}_{n\in \mathcal{N}(e) \setminus \mathcal{N}(\mathcal{SG}_i')}) = 0, e \in \mathcal{I}_2(\mathcal{SG}_i')
\end{align*}
```

Here, $(\cdot)^k$ denotes the iteration counter, $\lambda_e$ denotes the dual variable on their respective constraints, and $\mathcal{I}(\mathcal{SG}_i')$ denotes the set of incident edges (i.e., edges which contain nodes in the overlapped subgraph and outside of the overlapped subgraph). The sets $\mathcal{I}_1$ and $\mathcal{I}_2$ denote how the incident linking constraints are formulated in the subproblem (as either primal or dual coupled). 

In the [thesis](https://asset.library.wisc.edu/1711.dl/V2UHW7KSFIKBQ8Q/R/file-e04b2.pdf) of Prof. Sungho Shin, it was also shown to improve convergence to include an augmented penalty term in the objective, and this is included in PlasmoSchwarz. The objective therefore also gets updated with the augmented term 

```math
\begin{align*}
\sum_{e \in \mathcal{I}_1(\mathcal{SG}_i')} \mu \frac{1}{2} ||g_e(\{x_n \}_{n \in \mathcal{N}(e) \cap\mathcal{N}(\mathcal{SG}_i')}, \{x_n^k \}_{n\in \mathcal{N}(e) \setminus \mathcal{N}(\mathcal{SG}_i')})||
\end{align*}
```

where $\mu$ is a parameter that can be set by the user. 

Convergence of the algorithm is determined by primal and dual infeasibility which is computed at each iteration. The primal feasibility is determined by the incident linking constraints being feasible between the different subproblems, and the dual infeasibility is determined by the difference in the dual variables of the link constraints in the overlapping regions. 