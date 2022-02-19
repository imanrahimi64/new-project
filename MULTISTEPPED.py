import autograd.numpy as anp
import numpy as np
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.factory import get_problem
from pymoo.optimize import minimize
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
import matplotlib.pyplot as plt
from pymoo.core.repair import Repair


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

            xl = [0.01, 0.01, 0.0]
            xu = [0.45, 0.10, 1.0]
            xll = [0.01, 0.01, 0.01]
            xuu = [0.45, 0.10, 0.10]

            xuu[2] = max(
                min([xuu[2], (1.08 - 0.00139 / (x[0] * x[1])) / 4.94, (1.0986 - 0.000306 / (x[0] * x[1])) / 1.082,
                     (45948.98 - 12.307 / (x[0] * x[1])) / 49408.24, (16696.71 - 02.098 / (x[0] * x[1])) / 8046.33,
                     (10705.04 - 2.138 / (x[0] * x[1])) / 7883.39, (2136.54 - 0.417 * (x[0] * x[1])) / 1721.26,
                     (604.48 - 0.164 / (x[0] * x[1])) / 631.13]), xll[2])
            # xuu[0]=max(min([xuu[0],0.00139/(1.08-4.94*x[2])*x[1],0.000306/(1.0986-1.082*x[2])*x[1],12.307/(45948.98-49408.24*x[2])*x[1],2.098/(16696.71-8046.33*x[2])*x[1],2.138/(10705.04-7883.39*x[2])*x[1],0.417/(2136.54-1721.26*x[2])*x[1],0.164/(604.48-631.13*x[2])*x[1]]),xll[0])
            xll[2] = min(
                max([xll[2], (1.08 - 0.00139 / (x[0] * x[1])) / 4.94, (1.0986 - 0.000306 / (x[0] * x[1])) / 1.082,
                     (45948.98 - 12.307 / (x[0] * x[1])) / 49408.24, (16696.71 - 02.098 / (x[0] * x[1])) / 8046.33,
                     (10705.04 - 2.138 / (x[0] * x[1])) / 7883.39, (2136.54 - 0.417 * (x[0] * x[1])) / 1721.26,
                     (604.48 - 0.164 / (x[0] * x[1])) / 631.13]), xuu[2])

            x[2] = xll[2] + x[2] * (xuu[2] - xll[2])

        return pop


problem = CantileveredBeam()
algorithm1 = NSGA2(pop_size=100,eliminate_duplicates=True,repair=MyRepair())


algorithm = NSGA2(
    pop_size=100,
    eliminate_duplicates=True)
res = minimize(problem,
               algorithm,
               seed=1,
               verbose=False, save_history=True)

res1 = minimize(problem,algorithm1,seed=1,verbose=False,save_history=True)
from pymoo.visualization.scatter import Scatter
from pymoo.factory import get_visualization

# print(res.X,res.F)


plot = Scatter()
plot = get_visualization("scatter", legend=True, angle=(45, 30))
plot.add(res.F, color="blue", marker="s", label="CantileveredBeam coupled with BU")
plot.add(res1.F, color="red", facecolor="none",marker="v",label="CantileveredBeam without BU")
plot.show()
