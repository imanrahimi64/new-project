import autograd.numpy as anp
import numpy as np
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.factory import get_problem
from pymoo.optimize import minimize
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.problems.constr_as_penalty import ConstraintsAsPenalty
from pymoo.factory import get_problem, get_termination
import matplotlib.pyplot as plt
from pymoo.core.repair import Repair
from matplotlib import rc, rcParams
from pylab import *
from pymoo.factory import get_visualization

from pymoo.util.termination.default import MultiObjectiveDefaultTermination
from pymoo.core.population import Population
from pymoo.core.evaluator import Evaluator
class CantileveredBeam(Problem):

    def __init__(self):
        super().__init__(n_var=10, n_obj=2, n_constr=11, type_var=np.double)
        self.xl = np.array([1, 30, 2.4, 45, 2.4, 45, 1, 30, 1, 30])
        self.xu = np.array([5, 65, 3.1, 60, 3.1, 60, 5, 65, 5, 65])
        # self.h1 = np.array([0.1, 0.25, 0.35, 0.5, 0.65, 0.75, 0.9, 1.0])

    def _evaluate(self, x, out, *args, **kwargs):
        E, L, P, l, deltamax, sigmamax = 2e7, 500.0, 50000.0, 100.0, 2.7, 14000

        # b1, h1, b2, h2,b3,h3,b4,h4,b5,h5 = x[:, 0], x[:, 1], x[:, 2], x[:, 3],x[:, 4],x[:, 6],x[:, 7],x[:, 8],x[:, 9],x[:, 10]

        # I = 1 / 12 * b2 * (H - 2 * h1) ** 3 + 2 * (1 / 12 * b1 * h1 ** 3 + b1 * h1 * (H - h1) ** 2 / 4)
        # sigma = P * L * H / (2 * I)
        volume = (x[:, 0] * x[:, 1] + x[:, 2] * x[:, 3] + x[:, 4] * x[:, 5] + x[:, 6] * x[:, 7] + x[:, 8] * x[:, 9]) * l
        g6 = 10000.0 * ((244 / (x[:, 4] * x[:, 5] ** 3)) + (148 / (x[:, 6] * x[:, 7] ** 3)) + (
                    76 / (x[:, 8] * x[:, 9] ** 3)) + (28 / (x[:, 0] * x[:, 1] ** 3)) + (4 / (x[:, 2] * x[:, 3] ** 3)))
        # deflection=(sigma - 5000.0)/5000.0
        out["F"] = anp.column_stack([volume, g6])
        # volume

        # sigma = P * L * H / (2 * I)
        # delta = P * L ** 3 / (3 * E * I)

        g1 = 10.7141 - ((x[:, 4] * x[:, 5] ** 2) / 1000)
        g2 = 8.5714 - ((x[:, 6] * x[:, 7] ** 2) / 1000)
        g3 = 6.4286 - ((x[:, 8] * x[:, 9] ** 2) / 1000)
        g4 = 4.2857 - ((x[:, 0] * x[:, 1] ** 2) / 1000)
        g5 = 2.1429 - ((x[:, 2] * x[:, 3] ** 2) / 1000)
        # g6 = 10000.0*((244/(x[:, 4]*x[:, 5]**3))+(148/(x[:, 6]*x[:, 7]**3))+(76/(x[:, 8]*x[:, 9]**3))+(28/(x[:, 0]*x[:, 1]**3))+(4/(x[:, 2]*x[:, 3]**3)))-1086.0
        g7 = x[:, 5] - 20 * x[:, 4]
        g8 = x[:, 7] - 20 * x[:, 6]
        g9 = x[:, 9] - 20 * x[:, 8]
        g10 = x[:, 1] - 20 * x[:, 0]
        g11 = x[:, 3] - 20 * x[:, 2]
        out["G"] = anp.column_stack([g1, g2, g3, g4, g5, g7, g8, g9, g10, g11])


class MyRepair(Repair):
    def _do(self, Problem, pop, **kwargs):
        for k in range(len(pop)):
            x = pop[k].X

            xl = [1, 0, 2.4, 0, 2.4, 0, 1, 0, 1, 0]
            xu = [5, 1, 3.1, 1, 3.1, 1, 5, 1, 5, 1]
            xll = [1, 30, 2.4, 45, 2.4, 45, 1, 30, 1, 30]
            xuu = [5, 65, 3.1, 60, 3.1, 60, 5, 65, 5, 65]

            xll[1] = min(max([xll[1], 20 * x[0]]), xuu[1])
            xuu[1] = max(min([xuu[1], 20 * x[0]]), xll[1])
            xll[3] = min(max([xll[3], 20 * x[2]]), xuu[3])
            xuu[3] = max(min([xuu[3], 20 * x[2]]), xll[3])
            xll[5] = min(max([xll[5], 20 * x[4]]), xuu[5])
            xuu[5] = max(min([xuu[5], 20 * x[4]]), xll[5])
            xll[7] = min(max([xll[7], 20 * x[6]]), xuu[7])
            xuu[7] = max(min([xuu[7], 20 * x[6]]), xll[7])
            xll[9] = min(max([xll[9], 20 * x[8]]), xuu[9])
            xuu[9] = max(min([xuu[9], 20 * x[8]]), xll[9])
            # xuu[1]=max(min([xuu[1],20*x[0]]),xll[1])

            x[1] = xll[1] + x[1] * (xuu[1] - xll[1])
            x[3] = xll[3] + x[3] * (xuu[3] - xll[3])
            x[5] = xll[5] + x[5] * (xuu[5] - xll[5])
            x[7] = xll[7] + x[7] * (xuu[7] - xll[7])
            x[9] = xll[9] + x[9] * (xuu[9] - xll[9])

        return pop


problem = CantileveredBeam()
#problem = ConstraintsAsPenalty(CantileveredBeam(), penalty=1e6)
X = np.random.random((300, 10))
termination = MultiObjectiveDefaultTermination(

    cv_tol=1e-6,
    # f_tol=1e-6,
    n_max_gen=60

)

termination1 = get_termination("n_gen", 120)
termination2 = get_termination("n_gen", 60)

from pymoo.factory import get_problem, get_termination
from pymoo.util.termination.x_tol import DesignSpaceToleranceTermination
from matplotlib import rc, rcParams
from pylab import *
import math
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.factory import get_problem, get_reference_directions
from pymoo.optimize import minimize
from pymoo.algorithms.moo.moead import MOEAD
from pymoo.factory import get_visualization
from pymoo.algorithms.soo.nonconvex.brkga import BRKGA
# termination = DesignSpaceToleranceTermination(tol=0.0000000025, n_last=20)
# termination = get_termination("cv_tol", 1e-10)
algorithm2 = NSGA2(pop_size=200,  sampling=X)
res2 = minimize(problem, algorithm2, termination1, seed=1, verbose=True, save_history=True)
algorithm1 = NSGA2(pop_size=200,  repair=MyRepair(),sampling=X)
res1 = minimize(problem, algorithm1, termination,copy_algorithm=False, seed=1, verbose=True, save_history=True)
np.save("checkpoint", algorithm1)
checkpoint, = np.load("checkpoint.npy", allow_pickle=True).flatten()
print("Loaded Checkpoint:", checkpoint)
checkpoint.has_terminated = False
res3 = minimize(problem, checkpoint, termination1, seed=1, verbose=True, save_history=True)
#opt1 = res1.opt[0]
X = res1.X
pop = res1.pop
XX = pop.get("X")
pop = Population.new("X", XX)
Evaluator().eval(problem, pop)
algorithm = NSGA2(pop_size=200,sampling=pop)
res = minimize(problem, algorithm, termination2, seed=1, verbose=True, save_history=True)

from pymoo.visualization.scatter import Scatter
from pymoo.factory import get_visualization

# print(res.X,res.F)

# plt.figure(figsize=(10, 10))

plot = Scatter(legend=(True, {'loc': "upper right"}))
plot = get_visualization("scatter", legend=True, angle=(45, 30))
plot.add(res2.F, color="blue", marker="s", label="Without BU")
plot.add(res3.F, color="red", facecolor="none", marker="v", label="With BU")
plot.add(res.F, color="green", facecolor="none", marker="o", label="Hybrid")
plt.xlabel("$f1$", fontsize=12)
plt.ylabel("$f2$", fontsize=12)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
ax = gca()
ax.xaxis.set_tick_params(color='black', labelsize=12)
ax.yaxis.set_tick_params(color='black', labelsize=12)
#plot.legend(loc='upper right', prop={'size': 16})
# plt.rc('xtick', labelsize=16)
# plt.rc('ytick', labelsize=16)
plt.savefig('CantileveredBeam.jpg')
plot.show()
