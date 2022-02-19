import sys
from collections import defaultdict
import xlrd
import gurobipy as gp
from gurobipy import GRB


class Technician():
    def __init__(self, name, cap, depot):
        self.name = name
        self.cap = cap
        self.depot = depot

    def __str__(self):
        return f"Technician: {self.name}\n  Capacity: {self.cap}\n  Depot: {self.depot}"


class Job():
    def __init__(self, name, priority, duration, coveredBy):
        self.name = name
        self.priority = priority
        self.duration = duration
        self.coveredBy = coveredBy

    def __str__(self):
        about = f"Job: {self.name}\n  Priority: {self.priority}\n  Duration: {self.duration}\n  Covered by: "
        about += ", ".join([t.name for t in self.coveredBy])
        return about


class Customer():
    def __init__(self, name, loc, job, tStart, tEnd, tDue):
        self.name = name
        self.loc = loc
        self.job = job
        self.tStart = tStart
        self.tEnd = tEnd
        self.tDue = tDue

    def __str__(self):
        coveredBy = ", ".join([t.name for t in self.job.coveredBy])
        return f"Customer: {self.name}\n  Location: {self.loc}\n  Job: {self.job.name}\n  Priority: {self.job.priority}\n  Duration: {self.job.duration}\n  Covered by: {coveredBy}\n  Start time: {self.tStart}\n  End time: {self.tEnd}\n  Due time: {self.tDue}"


# Open Excel workbook
wb = xlrd.open_workbook('data-Sce0.xlsx')

# Read technician data
ws = wb.sheet_by_name('Technicians')
technicians = []
for i, t in enumerate(ws.col_values(0)[3:]):
    # Create Technician object
    thisTech = Technician(*ws.row_values(3 + i)[:3])
    technicians.append(thisTech)
# Read job data
jobs = []
for j, b in enumerate(ws.row_values(0)[3:]):
    coveredBy = [t for i, t in enumerate(technicians) if ws.cell_value(3 + i, 3 + j) == 1]
    # Create Job object
    thisJob = Job(*ws.col_values(3 + j)[:3], coveredBy)
    jobs.append(thisJob)
# Read location data
ws = wb.sheet_by_name('Locations')
locations = ws.col_values(0)[1:]
dist = {(l, l): 0 for l in locations}
for i, l1 in enumerate(locations):
    for j, l2 in enumerate(locations):
        if i < j:
            dist[l1, l2] = ws.cell_value(1 + i, 1 + j)
            dist[l2, l1] = dist[l1, l2]
# Read customer data
ws = wb.sheet_by_name('Customers')
customers = []
for i, c in enumerate(ws.col_values(0)[1:]):
    for b in jobs:
        if b.name == ws.cell_value(1 + i, 2):
            # Create Customer object using corresponding Job object
            rowVals = ws.row_values(1 + i)
            # print(rowVals)
            thisCustomer = Customer(*rowVals[:2], b, *rowVals[3:])
            customers.append(thisCustomer)
            break


def solve_trs0(technicians, customers, dist):
    # Build useful data structures
    K = [k.name for k in technicians]
    C = [j.name for j in customers]
    J = [j.loc for j in customers]
    L = list(set([l[0] for l in dist.keys()]))
    D = list(set([t.depot for t in technicians]))
    cap = {k.name: k.cap for k in technicians}
    loc = {j.name: j.loc for j in customers}
    depot = {k.name: k.depot for k in technicians}
    canCover = {j.name: [k.name for k in j.job.coveredBy] for j in customers}
    dur = {j.name: j.job.duration for j in customers}
    tStart = {j.name: j.tStart for j in customers}
    tEnd = {j.name: j.tEnd for j in customers}
    tDue = {j.name: j.tDue for j in customers}
    priority = {j.name: j.job.priority for j in customers}

    ### Create model
    m = gp.Model("trs0")

    ### Decision variables
    # Customer-technician assignment
    x = m.addVars(C, K, vtype=GRB.BINARY, name="x")

    # Technician assignment
    u = m.addVars(K, vtype=GRB.BINARY, name="u")

    # Edge-route assignment to technician
    y = m.addVars(L, L, K, vtype=GRB.BINARY, name="y")

    # Technician cannot leave or return to a depot that is not its base
    for k in technicians:
        for d in D:
            if k.depot != d:
                for i in L:
                    y[i, d, k.name].ub = 0
                    y[d, i, k.name].ub = 0

    # Start time of service
    t = m.addVars(L, ub=600, name="t")

    # Lateness of service
    z = m.addVars(C, name="z")

    # Artificial variables to correct time window upper and lower limits
    xa = m.addVars(C, name="xa")
    xb = m.addVars(C, name="xb")

    # Unfilled jobs
    g = m.addVars(C, vtype=GRB.BINARY, name="g")

    ### Constraints

    # A technician must be assigned to a job, or a gap is declared (1)
    m.addConstrs((gp.quicksum(x[j, k] for k in canCover[j]) + g[j] == 1 for j in C), name="assignToJob")

    # At most one technician can be assigned to a job (2)
    m.addConstrs((x.sum(j, '*') <= 1 for j in C), name="assignOne")

    # Technician capacity constraints (3)
    capLHS = {k: gp.quicksum(dur[j] * x[j, k] for j in C) + \
                 gp.quicksum(dist[i, j] * y[i, j, k] for i in L for j in L) for k in K}
    m.addConstrs((capLHS[k] <= cap[k] * u[k] for k in K), name="techCapacity")

    # Technician tour constraints (4 and 5)
    m.addConstrs((y.sum('*', loc[j], k) == x[j, k] for k in K for j in C), \
                 name="techTour1")
    m.addConstrs((y.sum(loc[j], '*', k) == x[j, k] for k in K for j in C), \
                 name="techTour2")

    # Same depot constraints (6 and 7)
    m.addConstrs((gp.quicksum(y[j, depot[k], k] for j in J) == u[k] for k in K), \
                 name="sameDepot1")
    m.addConstrs((gp.quicksum(y[depot[k], j, k] for j in J) == u[k] for k in K), \
                 name="sameDepot2")

    # Temporal constraints (8) for customer locations
    M = {(i, j): 600 + dur[i] + dist[loc[i], loc[j]] for i in C for j in C}
    m.addConstrs((t[loc[j]] >= t[loc[i]] + dur[i] + dist[loc[i], loc[j]] \
                  - M[i, j] * (1 - gp.quicksum(y[loc[i], loc[j], k] for k in K)) \
                  for i in C for j in C), name="tempoCustomer")

    # Temporal constraints (8) for depot locations
    M = {(i, j): 600 + dist[i, loc[j]] for i in D for j in C}
    m.addConstrs((t[loc[j]] >= t[i] + dist[i, loc[j]] \
                  - M[i, j] * (1 - y.sum(i, loc[j], '*')) for i in D for j in C), \
                 name="tempoDepot")

    # Time window constraints (9 and 10)
    m.addConstrs((t[loc[j]] + xa[j] >= tStart[j] for j in C), name="timeWinA")
    m.addConstrs((t[loc[j]] - xb[j] <= tEnd[j] for j in C), name="timeWinB")

    # Lateness constraint (11)
    m.addConstrs((z[j] >= t[loc[j]] + dur[j] - tDue[j] for j in C), \
                 name="lateness")

    ### Objective function
    M = 6100

    m.setObjective(z.prod(priority) + gp.quicksum(0.01 * M * priority[j] * (xa[j] + xb[j]) for j in C) +
                   gp.quicksum(M * priority[j] * g[j] for j in C), GRB.MINIMIZE)

    m.write("TRS0.lp")
    m.optimize()

    status = m.Status
    if status in [GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED]:
        print("Model is either infeasible or unbounded.")
        sys.exit(0)
    elif status != GRB.OPTIMAL:
        print("Optimization terminated with status {}".format(status))
        sys.exit(0)

    ### Print results
    # Assignments
    print("")
    for j in customers:
        if g[j.name].X > 0.5:
            jobStr = "Nobody assigned to {} ({}) in {}".format(j.name, j.job.name, j.loc)
        else:
            for k in K:
                if x[j.name, k].X > 0.5:
                    jobStr = "{} assigned to {} ({}) in {}. Start at t={:.2f}.".format(k, j.name, j.job.name, j.loc,
                                                                                       t[j.loc].X)
                    if z[j.name].X > 1e-6:
                        jobStr += " {:.2f} minutes late.".format(z[j.name].X)
                    if xa[j.name].X > 1e-6:
                        jobStr += " Start time corrected by {:.2f} minutes.".format(xa[j.name].X)
                    if xb[j.name].X > 1e-6:
                        jobStr += " End time corrected by {:.2f} minutes.".format(xb[j.name].X)
        print(jobStr)

    # Technicians
    print("")
    for k in technicians:
        if u[k.name].X > 0.5:
            cur = k.depot
            route = k.depot
            while True:
                for j in customers:
                    if y[cur, j.loc, k.name].X > 0.5:
                        route += " -> {} (dist={}, t={:.2f}, proc={})".format(j.loc, dist[cur, j.loc], t[j.loc].X,
                                                                              j.job.duration)
                        cur = j.loc
                for i in D:
                    if y[cur, i, k.name].X > 0.5:
                        route += " -> {} (dist={})".format(i, dist[cur, i])
                        cur = i
                        break
                if cur == k.depot:
                    break
            print("{}'s route: {}".format(k.name, route))
        else:
            print("{} is not used".format(k.name))

            # Utilization
    print("")
    for k in K:
        used = capLHS[k].getValue()
        total = cap[k]
        util = used / cap[k] if cap[k] > 0 else 0
        print("{}'s utilization is {:.2%} ({:.2f}/{:.2f})".format(k, util, \
                                                                  used, cap[k]))
    totUsed = sum(capLHS[k].getValue() for k in K)
    totCap = sum(cap[k] for k in K)
    totUtil = totUsed / totCap if totCap > 0 else 0
    print("Total technician utilization is {:.2%} ({:.2f}/{:.2f})".format(totUtil, totUsed, totCap))


def printScen(scenStr):
    sLen = len(scenStr)
    print("\n" + "*" * sLen + "\n" + scenStr + "\n" + "*" * sLen + "\n")


if __name__ == "__main__":
    # Base model
    printScen("Solving base scenario model")
    solve_trs0(technicians, customers, dist)