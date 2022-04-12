import matlab.engine
import numpy as np
from pymoo.util.misc import stack
from pymoo.optimize import minimize
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.visualization.scatter import Scatter
from multiprocessing.pool import ThreadPool
from pymoo.factory import get_problem
    
def nsga(R_input, theta_input, D_input):
    # it takes the user input for radius, pitch angle and diameter and do an NSGA-II optimization and return the best solution index

    
    eng = matlab.engine.start_matlab()
    
    class MyProblem(Problem):
        def __init__(self):
            super().__init__(n_var=6, n_obj=3, n_constr=0, xl=1, xu=8, type_var=int, ElementwiseProblem = False)
        def _evaluate(self, x, out, *args, **kwargs):
            x = x.tolist()
            x = matlab.double(x)
            radius, angle, diameter = eng.vectorized_for_helix(x, nargout=3)
            radius = np.array(radius)
            angle = np.array(angle)
            diameter = np.array(diameter)
            radius = radius - R_input
            angle = angle - theta_input
            diameter = diameter - D_input
            out["F"] = np.column_stack([-radius, -angle, -diameter])
    problem = MyProblem()
    #problem = get_problem("bnh")
    method = NSGA2(
        pop_size=40,
        n_offsprings=30,
        sampling=get_sampling("int_random"),
        crossover=get_crossover("int_sbx", prob=0.9, eta=15),
        mutation=get_mutation("int_pm", eta=20),
        eliminate_duplicates=True
    )
    
    
    res = minimize(problem,
                   method,
                   termination=('n_gen', 100),
                   seed=1,
                   save_history=True,
                   verbose=True
                   )
    
    print("Best solution found: %s" % res.X)
    print("Function value: %s" % res.F)
    print("Constraint violation: %s" % res.CV)
    
    # Objective Space
    ps = problem.pareto_set(use_cache=False, flatten=False)
    pf = problem.pareto_front(use_cache=False, flatten=False)
    plot = Scatter(title = "Objective Space")
    plot.add(res.F)
    if pf is not None:
        plot.add(pf, plot_type="line", color="black", alpha=0.7)
    plot.show()
    
    from datetime import datetime
    
    now = datetime.now()
    
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    
    desired_result_index = int(abs((len(res.X)/2)+1))
    return res.X[desired_result_index]

def position(index):
    # it takes index of the best result and return the xyz coordinates of the cnc
    eng = matlab.engine.start_matlab()
    index = list(index)
    index = matlab.double(index)
    XYZ = eng.helix_pos(index, nargout=1)
    XYZ = np.asarray(XYZ)
    return XYZ