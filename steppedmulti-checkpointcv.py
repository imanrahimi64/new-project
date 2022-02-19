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
from pymoo.visualization.scatter import Scatter
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
    n_max_gen=50

)

termination1 = get_termination("n_gen", 100)
termination2 = get_termination("n_gen", 50)

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
from pymoo.visualization.scatter import Scatter
from pymoo.factory import get_visualization
# termination = DesignSpaceToleranceTermination(tol=0.0000000025, n_last=20)
# termination = get_termination("cv_tol", 1e-10)
algorithm2 = NSGA2(pop_size=200,  sampling=X)
res2 = minimize(problem, algorithm2, termination1, seed=1, verbose=True, save_history=True)
algorithm1 = NSGA2(pop_size=200,  repair=MyRepair(),sampling=X)
res1 = minimize(problem, algorithm1, termination1,copy_algorithm=False, seed=1, verbose=True, save_history=True)




# print(res.X,res.F)

# plt.figure(figsize=(10, 10))

hist = res2.history
histt = res1.history

n_evals = []             # corresponding number of function evaluations\
hist_F = []              # the objective space values in each generation
hist_cv = []             # constraint violation in each generation
hist_cv_avg = []         # average constraint violation in the whole population

n_evalss = []             # corresponding number of function evaluations\
hist_FF = []              # the objective space values in each generation
hist_cvv = []             # constraint violation in each generation
hist_cv_avgg = []

for algo in hist:

    # store the number of function evaluations
    n_evals.append(algo.evaluator.n_eval)

    # retrieve the optimum from the algorithm
    opt = algo.opt

    # store the least contraint violation and the average in each population
    hist_cv.append(opt.get("CV").min())
    hist_cv_avg.append(algo.pop.get("CV").mean())

    # filter out only the feasible and append and objective space values
    feas = np.where(opt.get("feasible"))[0]
    hist_F.append(opt.get("F")[feas])
vals = hist_cv_avg


for algoo in histt:

    # store the number of function evaluations
    n_evalss.append(algoo.evaluator.n_eval)

    # retrieve the optimum from the algorithm
    optt = algoo.opt

    # store the least contraint violation and the average in each population
    hist_cvv.append(optt.get("CV").min())
    hist_cv_avgg.append(algoo.pop.get("CV").mean())

    # filter out only the feasible and append and objective space values
    feass = np.where(optt.get("feasible"))[0]
    hist_FF.append(optt.get("F")[feass])
valss = hist_cv_avgg




k = np.where(np.array(hist_cv) <= 0.0)[0].min()
print(f"At least one feasible solution in Generation {k} after {n_evals[k]} evaluations without BU.")

k1 = np.where(np.array(vals) <= 0.0)[0].min()
print(f"Whole population feasible in Generation {k1} after {n_evals[k1]} evaluations without BU.")


kk = np.where(np.array(hist_cvv) <= 0.0)[0].min()
print(f"At least one feasible solution in Generation {kk} after {n_evalss[kk]} evaluations with BU and hybrid.")

kk1 = np.where(np.array(valss) <= 0.0)[0].min()
print(f"Whole population feasible in Generation {kk1} after {n_evalss[kk1]} evaluations with BU and hybrid.")


plt.figure(figsize=(15, 10))

plt.plot(n_evals, vals,  color='black', linewidth=4,label="Avg. CV of Pop (without BU)")
plt.scatter(n_evals, vals,  facecolor="none", edgecolor='black', marker="p")
plt.plot(n_evalss, valss,  color='magenta', linewidth=4, label="Avg. CV of Pop (with BU and hybrid)")
plt.scatter(n_evalss, valss,  facecolor="none", edgecolor='black', marker="p")

plt.axvline(n_evals[k], color="red", label="First feasible solution found- without BU",linewidth=4, linestyle=":")
plt.axvline(n_evalss[kk], color="green", label="First  feasibile solution(s) found-with BU and hybrid",linewidth=4, linestyle="--")

#plt.plot(n_evals, vals,  color='black', lw=0.7, label="Avg. CV of Pop")
#plt.scatter(n_evals, vals,  facecolor="none", edgecolor='black', marker="p")
plt.axvline(n_evals[k1], color="blue", label="All Feasible solutions found-without BU",linewidth=4, linestyle="-.")
plt.axvline(n_evalss[kk1], color="cyan", label="All Feasible solutions found-with BU and hybrid", linewidth=4,linestyle="--")

#plt.title("Convergence")
plt.xlabel("Function Evaluations",fontsize=18)
plt.ylabel("Constraint Violation",fontsize=18)

plt.legend(loc='upper right',prop={'size': 18})
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
ax = gca()
ax.xaxis.set_tick_params(color='black',labelsize=18)
ax.yaxis.set_tick_params(color='black',labelsize=18)
plt.legend(fontsize=12)
plt.savefig('CV convergence.jpg')
plt.show()