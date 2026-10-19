from steadiercom.steadiercom import allocation_constraints, build_problem
from reframed.solvers.solution import Status
from reframed.community.solution import CommunitySolution
from reframed.solvers.solver import Parameter
from warnings import warn


def build_problem_extended(community, growth=0, abundance=None, constraints=None):
    """
    function for building steadiercom problem for steadiercomVA and binary_search. 
    
    Essentially: makes it possibiltiy to adjust growth rate.

    NB
    growth=0.0  just to avoid SteadierCom throwing an error (overidden by the objective)
    """
    solver = build_problem(community, growth=growth, abundance=abundance, constraints=constraints)
    solver.set_parameter(Parameter.FEASIBILITY_TOL,1e-7)

    def update_growth(value):
        grb_model = solver.problem

        for org_id in community.organisms:
            constr = grb_model.getConstrByName(f"g_{org_id}")
            var = grb_model.getVarByName(f"x_{org_id}")

            if constr is None:
                raise RuntimeError(f"Constraint g_{org_id} not found in Gurobi model.")

            if var is None:
                raise RuntimeError(f"Variable x_{org_id} not found in Gurobi model.")

            grb_model.chgCoeff(constr, var, float(value))

        grb_model.update()

    solver.update_growth = update_growth

    return solver


def SteadierComVA(community, obj_frac=1.0, constraints=None, w_e=0.002, w_r=0.2, max_growth=None,organisms=None):
    """
    SteadyComVA-inspired method for determining variability in abundance
    """
    solver = build_problem_extended(community, constraints=constraints)
    solver.problem.setParam("Method", 1)
    solver.problem.setParam("LPWarmStart", 2)

    allocation_constraints(community, solver, w_e=w_e, w_r=w_r)

    if max_growth is None:
        objective = {community.merged_model.biomass_reaction: 1} 
        sol = binary_search(solver, objective, constraints=constraints)
        growth = obj_frac * sol.values[community.merged_model.biomass_reaction]
    else:
        growth = obj_frac * max_growth

    solver.update_growth(growth)

    variability = {org_id: [None, None] for org_id in community.organisms}

    print("\t finding min abundance for all members...")

    if organisms is None:
        organisms = community.organisms.keys()


    for org_id in organisms:
        sol2 = solver.solve({f"x_{org_id}": 1}, minimize=True, get_values=False, constraints=constraints)
        variability[org_id][0] = sol2.fobj

    print("\t finding max abundance for all members...")
    for org_id in organisms:
        sol2 = solver.solve({f"x_{org_id}": 1}, minimize=False, get_values=False, constraints=constraints)
        variability[org_id][1] = sol2.fobj

    return variability



def SteadierCom_binary_search(community, constraints=None, w_e=0.002, w_r=0.2,max_iters=20, abs_tol=1e-3,test_zero_growth=False):
    """
    SteadyCom-inspired method for optimising growth and predicting abundance simultaneously.

    The SteadyCom implementation is in reframed
    """
    solver = build_problem_extended(community, constraints=constraints)
    solver.problem.setParam("TimeLimit", 180)
    allocation_constraints(community, solver, w_e=w_e, w_r=w_r)

    objective = {community.merged_model.biomass_reaction: 1} 
    
    sol = binary_search(solver, objective, constraints=constraints,max_iters=max_iters, abs_tol=abs_tol,test_zero_growth=test_zero_growth)

    if sol.status == Status.OPTIMAL:
        solution = CommunitySolution(community, sol.values)
        solution.solver = solver

        return solution,sol.status
    
    return None, sol.status



def binary_search(solver, objective, obj_frac=1, minimize=False, max_iters=20, abs_tol=1e-3, constraints=None,test_zero_growth=False):
    """
    Optimize growth iteratively and predict abundance.
    """

    # ------------------------------------------------------------
    # Test feasibility at zero growth
    # ------------------------------------------------------------

    if test_zero_growth:
        solver.update_growth(0)

        # Pure feasibility test.
        sol0 = solver.solve(
            {},
            get_values=False,
            minimize=True,
            constraints=constraints,
        )

        zero_feasible = sol0.status == Status.OPTIMAL

        if not zero_feasible:
            warn(
                "Problem is infeasible at growth = 0. "
                "No positive growth rate is expected to be feasible."
            )
            return sol0

    # ------------------------------------------------------------
    # Binary search
    # ------------------------------------------------------------

    previous_value = 0
    value = 1
    fold = 2
    feasible = False
    last_feasible = 0

    for i in range(max_iters):

        diff = value - previous_value

        if diff < abs_tol:
            break

        if feasible:
            last_feasible = value
            previous_value = value
            value = fold*diff + value
        else:
            if i > 0:
                fold = 0.5
            value = fold*diff + previous_value

        print(f"\t {i}: last_feasible={last_feasible}, testing value={value}")

        solver.update_growth(value)
        sol = solver.solve({}, get_values=False, minimize=minimize, constraints=constraints) # remove objective to determine if solution space is valid

        feasible = sol.status == Status.OPTIMAL

    if feasible:
        solver.update_growth(obj_frac * value)
    else:
        solver.update_growth(obj_frac * last_feasible)

    sol = solver.solve(objective, minimize=minimize, constraints=constraints)

    if i == max_iters - 1:
        warn("Max iterations exceeded.")

    return sol