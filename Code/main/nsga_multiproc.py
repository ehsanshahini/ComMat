# -*- coding: utf-8 -*-
"""
Created on Thu May 27 22:38:19 2021

@author: shahini
"""


'''
eng = matlab.engine.start_matlab()
type_cnc = 2.0
index1 = 2.0
index2 = 2.0
index3 = 2.0
index4 = 2.0
shift_cnc = 7.0
radius,diameter,pitchangle = eng.helix(type_cnc, index1, index2, index3, index4, shift_cnc, nargout=3)

print(radius)
'''
import matlab.engine
import numpy as np
from pymoo.util.misc import stack
from pymoo.optimize import minimize
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.visualization.scatter import Scatter
from multiprocessing.pool import ThreadPool

eng = matlab.engine.start_matlab()

n_threads = 4
pool = ThreadPool(n_threads)

class MyProblem(Problem):
    def __init__(self):
        super().__init__(n_var=6, n_obj=2, n_constr=0, xl=1, xu=7, type_var=int, ElementwiseProblem=True, parallelization = ('starmap', pool.starmap))
    def _evaluate(self, x, out, *args, **kwargs):

        x = x.tolist()
        x = matlab.double(x)
        f1, f2 = eng.main(x, nargout=2)
        f1 = np.array(f1)
        f2 = np.array(f2)
        f1 = f1 - 5.0
        f2 = f2 - 45.0
        out["F"] = np.column_stack([-f1, -f2])
        
problem = MyProblem()

method = NSGA2(
    pop_size=40,
    n_offsprings=10,
    sampling=get_sampling("int_random"),
    crossover=get_crossover("int_sbx", prob=0.9, eta=15),
    mutation=get_mutation("int_pm", eta=20),
    eliminate_duplicates=True
)


res = minimize(problem,
               method,
               termination=('n_gen', 1),
               seed=1,
               save_history=True
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