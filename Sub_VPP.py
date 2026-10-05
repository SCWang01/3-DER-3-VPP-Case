# Sub_VPP.py — Gurobi best-response subproblems for the 3-VPPO, 3-DER case.
# Translated from Sub_VPP_1.m / Sub_VPP_2.m / Sub_VPP_3.m (YALMIP/Gurobi).
#
# Market parameters (¥ per contract period):
#   Switching costs   :  C1 = 0,   C2 = 30,  C3 = 60
#   Offer slopes      :  DER2→VPP1: 120,  DER1→VPP2: 120,  DER2→VPP3:  72
#                        DER3→VPP1: 240,  DER1→VPP3: 240,  DER3→VPP2:  72
#   Contributions     :  VPP1←DER2: 240,  VPP1←DER3: 480
#                        VPP2←DER1: 240,  VPP2←DER3: 144
#                        VPP3←DER1: 480,  VPP3←DER2: 144
#   Switching revenue :  VPP1: 0,  VPP2: 60*(1-y2[1]),  VPP3: 120*(1-y3[2])

import gurobipy as gp
from gurobipy import GRB

_M = 1e4   # big-M for complementarity encoding


def _require_solution(model):
    """Raise if Gurobi did not find a feasible solution."""
    if model.SolCount <= 0:
        raise RuntimeError(
            f"Gurobi found no feasible solution. Status: {model.Status}"
        )


def _der_kkt(model, y, mu, K, lam, epr, eps, tag):
    """Embed KKT stationarity and complementarity for one DER.

    Stationarity :  eps*y[j] - epr[j] + lam - mu[j] = 0   for j in {0,1,2}
    Compl. (big-M):  y[j] <= K[j]*M,   mu[j] <= (1-K[j])*M
    Feasibility  :  sum(y) = 1,  y >= 0,  mu >= 0  (bounds set at addVars)
    """
    for j in range(3):
        model.addConstr(y[j] <= K[j] * _M,        name=f"{tag}_yK{j}")
        model.addConstr(mu[j] <= (1 - K[j]) * _M, name=f"{tag}_muK{j}")
        model.addConstr(
            eps * y[j] - epr[j] + lam - mu[j] == 0,
            name=f"{tag}_stat{j}",
        )
    model.addConstr(y[0] + y[1] + y[2] == 1, name=f"{tag}_sum")


def _der_vars(model, tag):
    """Create (y, mu, K, lam) decision variables for one DER."""
    y   = model.addVars(3, lb=0.0, name=f"y_{tag}")
    mu  = model.addVars(3, lb=0.0, name=f"mu_{tag}")
    K   = model.addVars(3, vtype=GRB.BINARY, name=f"K_{tag}")
    lam = model.addVar(lb=-GRB.INFINITY, name=f"lam_{tag}")
    return y, mu, K, lam


# ---------------------------------------------------------------------------
# Subproblem 1: VPP 1 optimises x_1 given fixed x_2, x_3
# ---------------------------------------------------------------------------

def Sub_VPP_1(x_2, x_3, eps, verbose=False):
    """Best-response subproblem for VPP 1.

    Objective (from Sub_VPP_1.m):
        max  x_1 * (240*y_2[0] + 480*y_3[0])
             + 2*C1*(1 - y_1[0])          # C1 = 0, term vanishes
    where y_d[i] = P(DER d affiliates with VPP i+1).

    Offer slopes (¥):   DER1→VPP2: 120, DER1→VPP3: 240
                         DER2→VPP1: 120, DER2→VPP3:  72
                         DER3→VPP1: 240, DER3→VPP2:  72
    Switching costs:    C1=0, C2=30, C3=60
    """
    C1, C2, C3 = 0, 30, 60

    m = gp.Model("Sub_VPP_1")
    m.Params.NonConvex  = 2
    m.Params.OutputFlag = int(verbose)

    x1 = m.addVar(lb=0.0, ub=1.0, name="x_1")

    y1, mu1, K1, lam1 = _der_vars(m, "d1")   # DER 1
    y2, mu2, K2, lam2 = _der_vars(m, "d2")   # DER 2
    y3, mu3, K3, lam3 = _der_vars(m, "d3")   # DER 3

    # Evaluated profit vectors (EPR[j] = offer from VPP j+1 to this DER)
    #   DER 1 — incumbent VPP 1 (index 0 → EPR = 0)
    epr1 = [0,
            120*(1 - x_2) - C1,
            240*(1 - x_3) - C1]
    #   DER 2 — incumbent VPP 2 (index 1 → EPR = 0)
    epr2 = [120*(1 - x1) - C2,
            0,
            72*(1 - x_3) - C2]
    #   DER 3 — incumbent VPP 3 (index 2 → EPR = 0)
    epr3 = [240*(1 - x1) - C3,
            72*(1 - x_2) - C3,
            0]

    _der_kkt(m, y1, mu1, K1, lam1, epr1, eps, "d1")
    _der_kkt(m, y2, mu2, K2, lam2, epr2, eps, "d2")
    _der_kkt(m, y3, mu3, K3, lam3, epr3, eps, "d3")

    # Objective: maximise VPP 1 payoff (switching revenue = 0 since C1=0)
    m.setObjective(
        -(x1 * (240*y2[0] + 480*y3[0])),
        GRB.MINIMIZE,
    )
    m.optimize()
    _require_solution(m)

    return {
        "x_1": x1.X,
        "y1":  [y1[j].X for j in range(3)],
        "y2":  [y2[j].X for j in range(3)],
        "y3":  [y3[j].X for j in range(3)],
        # affiliation probability vectors indexed to VPP 1 (index 0)
        "y1_aff": y1[0].X,
        "y2_aff": y2[0].X,
        "y3_aff": y3[0].X,
        "obj": -m.ObjVal,
    }


# ---------------------------------------------------------------------------
# Subproblem 2: VPP 2 optimises x_2 given fixed x_1, x_3
# ---------------------------------------------------------------------------

def Sub_VPP_2(x_1, x_3, eps, verbose=False):
    """Best-response subproblem for VPP 2.

    Objective (from Sub_VPP_2.m):
        max  x_2 * (240*y_1[1] + 144*y_3[1])
             + 2*C2*(1 - y_2[1])          # C2 = 30, switching revenue
    """
    C1, C2, C3 = 0, 30, 60

    m = gp.Model("Sub_VPP_2")
    m.Params.NonConvex  = 2
    m.Params.OutputFlag = int(verbose)

    x2 = m.addVar(lb=0.0, ub=1.0, name="x_2")

    y1, mu1, K1, lam1 = _der_vars(m, "d1")
    y2, mu2, K2, lam2 = _der_vars(m, "d2")
    y3, mu3, K3, lam3 = _der_vars(m, "d3")

    epr1 = [0,
            120*(1 - x2) - C1,
            240*(1 - x_3) - C1]
    epr2 = [120*(1 - x_1) - C2,
            0,
            72*(1 - x_3) - C2]
    epr3 = [240*(1 - x_1) - C3,
            72*(1 - x2) - C3,
            0]

    _der_kkt(m, y1, mu1, K1, lam1, epr1, eps, "d1")
    _der_kkt(m, y2, mu2, K2, lam2, epr2, eps, "d2")
    _der_kkt(m, y3, mu3, K3, lam3, epr3, eps, "d3")

    # Switching revenue: VPP 2 earns C2 per period from DER 2 if it stays
    # The MATLAB objective uses 2*switching_2*(1-y_2[1]); here index 1 = VPP 2
    m.setObjective(
        -(x2 * (240*y1[1] + 144*y3[1]) + 2*C2*(1 - y2[1])),
        GRB.MINIMIZE,
    )
    m.optimize()
    _require_solution(m)

    return {
        "x_2": x2.X,
        "y1":  [y1[j].X for j in range(3)],
        "y2":  [y2[j].X for j in range(3)],
        "y3":  [y3[j].X for j in range(3)],
        "y1_aff": y1[1].X,
        "y2_aff": y2[1].X,
        "y3_aff": y3[1].X,
        "obj": -m.ObjVal,
    }


# ---------------------------------------------------------------------------
# Subproblem 3: VPP 3 optimises x_3 given fixed x_1, x_2
# ---------------------------------------------------------------------------

def Sub_VPP_3(x_1, x_2, eps, verbose=False):
    """Best-response subproblem for VPP 3.

    Objective (from Sub_VPP_3.m):
        max  x_3 * (480*y_1[2] + 144*y_2[2])
             + 2*C3*(1 - y_3[2])          # C3 = 60, switching revenue
    """
    C1, C2, C3 = 0, 30, 60

    m = gp.Model("Sub_VPP_3")
    m.Params.NonConvex  = 2
    m.Params.OutputFlag = int(verbose)

    x3 = m.addVar(lb=0.0, ub=1.0, name="x_3")

    y1, mu1, K1, lam1 = _der_vars(m, "d1")
    y2, mu2, K2, lam2 = _der_vars(m, "d2")
    y3, mu3, K3, lam3 = _der_vars(m, "d3")

    epr1 = [0,
            120*(1 - x_2) - C1,
            240*(1 - x3) - C1]
    epr2 = [120*(1 - x_1) - C2,
            0,
            72*(1 - x3) - C2]
    epr3 = [240*(1 - x_1) - C3,
            72*(1 - x_2) - C3,
            0]

    _der_kkt(m, y1, mu1, K1, lam1, epr1, eps, "d1")
    _der_kkt(m, y2, mu2, K2, lam2, epr2, eps, "d2")
    _der_kkt(m, y3, mu3, K3, lam3, epr3, eps, "d3")

    m.setObjective(
        -(x3 * (480*y1[2] + 144*y2[2]) + 2*C3*(1 - y3[2])),
        GRB.MINIMIZE,
    )
    m.optimize()
    _require_solution(m)

    return {
        "x_3": x3.X,
        "y1":  [y1[j].X for j in range(3)],
        "y2":  [y2[j].X for j in range(3)],
        "y3":  [y3[j].X for j in range(3)],
        "y1_aff": y1[2].X,
        "y2_aff": y2[2].X,
        "y3_aff": y3[2].X,
        "obj": -m.ObjVal,
    }


# ---------------------------------------------------------------------------
# Gauss–Seidel best-response iteration (mirrors Gauss_Iteration.m)
# ---------------------------------------------------------------------------

def gauss_iteration(
    x_init=(0.1, 0.1, 0.1),
    eps=0.1,
    max_iter=200,
    tol=1e-10,
    verbose=False,
):
    """Sequential best-response iteration over the three VPP subproblems.

    At each iteration:
        1. Solve Sub_VPP_1 with the current (x_2, x_3).
        2. Solve Sub_VPP_2 with the updated x_1 and current x_3.
        3. Solve Sub_VPP_3 with the updated (x_1, x_2).
    Stop when the L2-norm change in strategies falls below ``tol`` or
    ``max_iter`` iterations are reached.

    Returns
    -------
    history : list of (x_1, x_2, x_3) tuples, one per iteration.
    aff_history : list of dicts with affiliation probabilities per iteration.
    """
    x1, x2, x3 = x_init
    history     = [(x1, x2, x3)]
    aff_history = []

    for it in range(1, max_iter + 1):
        prev = (x1, x2, x3)

        sol1 = Sub_VPP_1(x2, x3, eps, verbose=verbose)
        x1   = sol1["x_1"]

        sol2 = Sub_VPP_2(x1, x3, eps, verbose=verbose)
        x2   = sol2["x_2"]

        sol3 = Sub_VPP_3(x1, x2, eps, verbose=verbose)
        x3   = sol3["x_3"]

        history.append((x1, x2, x3))
        aff_history.append({
            "iter": it,
            # affiliation probabilities with the winning VPP (own index)
            "VPP1": {"y1": sol1["y1_aff"], "y2": sol1["y2_aff"], "y3": sol1["y3_aff"]},
            "VPP2": {"y1": sol2["y1_aff"], "y2": sol2["y2_aff"], "y3": sol2["y3_aff"]},
            "VPP3": {"y1": sol3["y1_aff"], "y2": sol3["y2_aff"], "y3": sol3["y3_aff"]},
        })

        diff = sum((a - b)**2 for a, b in zip((x1, x2, x3), prev)) ** 0.5
        if diff <= tol:
            if verbose:
                print(f"Converged at iteration {it}  (Δ={diff:.2e})")
            break

    return history, aff_history


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import csv, sys

    eps      = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
    max_iter = int(sys.argv[2])   if len(sys.argv) > 2 else 200

    print(f"Running Gauss–Seidel  eps={eps}  max_iter={max_iter}")
    hist, aff = gauss_iteration(eps=eps, max_iter=max_iter, verbose=False)

    out = f"gs_trajectory_eps{eps}.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["iter", "x_1", "x_2", "x_3"])
        for i, (x1, x2, x3) in enumerate(hist):
            w.writerow([i, round(x1, 6), round(x2, 6), round(x3, 6)])
    print(f"Saved {len(hist)} iterations → {out}")
