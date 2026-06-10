import gurobipy as gp
from gurobipy import GRB


def _require_solution(model):
    if model.SolCount <= 0:
        raise RuntimeError(f"Gurobi did not find a feasible solution. Status: {model.Status}")


def Sub_VPP_1(x_2, x_3_or_epis=None, eps=None, verbose=False):
    """Gurobi version of the YALMIP Sub_VPP_1 model.

    Supports both call forms:
        Sub_VPP_1(x_2, epis)
        Sub_VPP_1(x_2, x_3, eps)

    The x_3 argument is accepted for compatibility with the existing skeleton,
    but it is not used by the MATLAB model shown for VPP 1.
    """
    if eps is None:
        epis = x_3_or_epis
    else:
        epis = eps

    if epis is None:
        raise ValueError("epis/eps must be provided")

    switching_1 = 40
    switching_2 = 35
    switching_3 = 45
    M = 1e4

    model = gp.Model("Sub_VPP_1")
    model.Params.NonConvex = 2
    model.Params.OutputFlag = 1 if verbose else 0

    x_1 = model.addVar(lb=0.0, ub=1.0, name="x_1")
    y_1 = model.addVars(3, lb=0.0, name="y_1")
    y_2 = model.addVars(3, lb=0.0, name="y_2")
    y_3 = model.addVars(3, lb=0.0, name="y_3")

    K = model.addVars(6, vtype=GRB.BINARY, name="K")
    mu_l = model.addVars(3, lb=0.0, name="mu_l")
    mu_h = model.addVars(3, lb=0.0, name="mu_h")

    # DER-1 Constraints:
    

    # Transferred KKT stationarity conditions.
    model.addConstr(
        epis / 2 * (4 * y_1 - 2)
        - 50 * (1 - x_2)
        + switching_1
        - mu_l[0]
        + mu_h[0]
        == 0,
        name="kkt_der1",
    )
    model.addConstr(
        epis / 2 * (4 * y_2 - 2)
        - 50 * (1 - x_1)
        + switching_2
        - mu_l[1]
        + mu_h[1]
        == 0,
        name="kkt_der2",
    )

    # Big-M complementarity constraints.
    model.addConstr(y_1 <= M * K[0], name="der1_y_lower_active")
    model.addConstr(mu_l[0] <= M * (1 - K[0]), name="der1_mu_l")
    model.addConstr(1 - y_1 <= M * K[1], name="der1_y_upper_active")
    model.addConstr(mu_h[0] <= M * (1 - K[1]), name="der1_mu_h")

    model.addConstr(y_2 <= M * K[2], name="der2_y_lower_active")
    model.addConstr(mu_l[1] <= M * (1 - K[2]), name="der2_mu_l")
    model.addConstr(1 - y_2 <= M * K[3], name="der2_y_upper_active")
    model.addConstr(mu_h[1] <= M * (1 - K[3]), name="der2_mu_h")

    model.setObjective(-100 * x_1 * y_2 - switching_1 * y_1, GRB.MINIMIZE)
    model.optimize()
    _require_solution(model)

    return {
        "x_1": x_1.X,
        "y_1": y_1.X,
        "y_2": y_2.X,
        "mu_l": [mu_l[i].X for i in range(2)],
        "mu_h": [mu_h[i].X for i in range(2)],
        "K": [round(K[i].X) for i in range(4)],
        "obj": model.ObjVal,
        "status": model.Status,
    }


def Sub_VPP_2(x_1, x_3, eps):
    raise NotImplementedError("Sub_VPP_2 has not been translated yet.")


def Sub_VPP_3(x_1, x_2):
    raise NotImplementedError("Sub_VPP_3 has not been translated yet.")
