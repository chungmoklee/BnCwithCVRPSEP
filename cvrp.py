import cplex
import time
import numpy as np
from enum import Enum
from cvrpsep import lib, ffi
import networkx as nx
import sys
import traceback


class Prob:
    def __init__(self, datafile):

        prob_name, n, num_vehicles, capacity, pos_x, pos_y, demand, distance_mat, time_mat, round_mult = self.read_solomon(datafile)

        num_nodes = n + 2

        node_id = n + 1

        demand =  demand.astype(np.int64)

        N = [i for i in range(1, n+1)]
        N_s = [0] + N
        N_t = N + [n+1]
        N_st = [0] + N + [n+1]

        T = time_mat.astype(np.int64)
        D = distance_mat.astype(np.int64)

        A = [(i,j) for i in N_s for j in N_t if i != j and (i,j) != (0,n+1)]
        
        Out_i = [[a[1] for a in A if a[0] == i] for i in N_st]
        In_i = [[a[0] for a in A if a[1] == j] for j in N_st]

        self.datafile = datafile
        self.prob_name = prob_name

        self.n = n
        self.num_nodes = num_nodes
        self.q = demand
        
        self.N = N
        self.N_s = N_s
        self.N_t = N_t
        self.N_st = N_st

        self.A = A
        self.Out_i = Out_i
        self.In_i = In_i
        
        self.T = T
        self.D = D

        self.pos_x = pos_x
        self.pos_y = pos_y

                
        self.Q = capacity
        self.K = num_vehicles

        self.round_mult = round_mult

    def read_solomon(self, datafile):
        with open(datafile, 'r') as f:
            lines = f.readlines()

        lines = [l.strip() for l in lines]
        prob_name = lines[0]

        n = int([l.split()[0] for l in lines if len(l.split()) > 0][-1])
        num_nodes = n + 2

        num_vehicles, capacity = int(lines[4].split()[0]), int(lines[4].split()[1])
        
        pos_x = np.zeros(num_nodes)
        pos_y = np.zeros(num_nodes)
        demand = np.zeros(num_nodes)
        ready_time = np.zeros(num_nodes)
        due_time = np.zeros(num_nodes)
        service_time = np.zeros(num_nodes)


        round_mult = 100


        for l in lines[9:n+10]:
            toks = l.split()

            node_id = int(toks[0])
            pos_x[node_id] = int(toks[1])
            pos_y[node_id] = int(toks[2])
            demand[node_id] = int(toks[3])
            ready_time[node_id] = int(toks[4]) * round_mult
            due_time[node_id] = int(toks[5]) * round_mult
            service_time[node_id] = int(toks[6]) * round_mult

        t = n + 1
        pos_x[t] = pos_x[0]
        pos_y[t] = pos_y[0]
        demand[t] = demand[0]
        ready_time[t] = ready_time[0]
        due_time[t] = due_time[0]
        service_time[t] = service_time[0]


        distance_mat = np.zeros((num_nodes, num_nodes))
        for i in range(num_nodes):
            for j in  range(num_nodes):
                distance_mat[i,j] = round(round_mult * np.hypot(pos_x[i] - pos_x[j], pos_y[i] - pos_y[j]))

        time_mat = distance_mat.copy()


        return prob_name, n, num_vehicles, capacity, pos_x, pos_y, demand, distance_mat, time_mat, round_mult
    
    def __repr__(self):
        return f'Prob("{self.datafile}")'

    def prob_info(self):
        return {
                'name': self.prob_name,
                'datafile': self.datafile,
                'num_cust_nodes': self.n,
                'num_tot_nodes': self.num_nodes
            }

    def get_file_name(self):
        return f'{self.prob_name}_N[{self.n}]'
    
    def draw_solution(self, sol):
        # sol = [(i,j) for (i,j) in A if cpx.solution.get_values(f'x_{i}_{j}') > 0.9]

        G = nx.Graph()

        for i in self.N_st:
            G.add_node(i)

        pos = {i: np.array([self.pos_x[i], self.pos_y[i]]) for i in range(self.num_nodes)}

        nx.draw_networkx_nodes(G, pos=pos, nodelist=[0], node_size=100, node_shape='s', node_color='gray')
        nx.draw_networkx_nodes(G, pos=pos, nodelist=self.N, node_size=100, node_color='red')
        nx.draw_networkx_labels(G, pos=pos, labels={i:f'{i}' for i in self.N_s}, font_size=7, font_weight='bold', font_color='white')

        G.add_edges_from(sol)

        nx.draw_networkx_edges(G, pos=pos, edge_color='g')


    


class MIPSolver:
    def __init__(self, prob):
        self.prob = prob
        self.formulate()

    def formulate(self):
        prob = self.prob

        cpx = cplex.Cplex()

        N = prob.N
        N_st = prob.N_st
        A = prob.A
        T = prob.T
        Out_i = prob.Out_i
        In_i = prob.In_i
        K = prob.K
        Q = prob.Q
        s = 0
        t = len(N)+1
        q = prob.q

        cpx.variables.add(
            names = [f'x_{i}_{j}' for i,j in A],
            types = ['B'] * len(A),
            obj = [float(T[i,j]) for i,j in A]
        )

        cpx.variables.add(
            names = [f'u_{i}' for i in N_st],
            types = ['C'] * len(N_st),
            lb = [float(q[i]) for i in N_st],
            ub = [Q for i in N_st]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in Out_i[i]], ([1] * (len(Out_i[i])))
                ) for i in N
            ],
            senses=['E'] * len(N),
            rhs=[1] * len(N),
            names=[f'out_degree_{i}' for i in N])


        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for i in In_i[j]], ([1] * (len(In_i[j])))
                ) for j in N
            ],
            senses=['E'] * len(N),
            rhs=[1] * len(N),
            names=[f'in_degree_{j}' for j in N])


        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{s}_{j}' for j in Out_i[s]], ([1] * (len(Out_i[s])))
                ) 
            ],
            senses=['L'],
            rhs=[K],
            names=['num_vehicles'])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'u_{i}', f'u_{j}', f'x_{i}_{j}'], ([1, -1, Q])
                ) for (i,j) in A
            ],
            senses=['L'] * len(A),
            rhs=[float(Q - q[j]) for (i,j) in A],
            names=[f'MTZ_{i}_{j}' for (i,j) in A])
        

        self.cpx = cpx



    def solve(self, timelimit=3600, cpu=4):

        cpx = self.cpx

        start_time = time.time()
        cpx.parameters.threads.set(cpu)
        cpx.parameters.timelimit.set(timelimit)

        # cpx.set_log_stream(None)
        # cpx.set_results_stream(None)

        cpx.solve() 

        num_bnb_nodes = cplex._internal._procedural.getnodecnt(cpx._env._e, cpx._lp)
        num_gap = cplex._internal._procedural.getmiprelgap(cpx._env._e, cpx._lp)
        solve_time = time.time() - start_time

        results = {
            'status': cpx.solution.get_status_string(),
            'obj': cpx.solution.get_objective_value() / self.prob.round_mult,
            'time': solve_time,
            'num_bnb_nodes': num_bnb_nodes,
            'num_gap': num_gap,
            'cpu': cpu,
            'timelimit': timelimit,
            'solution': [(i,j) for (i,j) in self.prob.A if cpx.solution.get_values(f'x_{i}_{j}') > 0.9]

        }

        self.results = results

        return results




class BnCSolver:
    def __init__(self, prob):
        self.prob = prob
        self.formulate()

    def formulate(self):
        prob = self.prob

        cpx = cplex.Cplex()

        N = prob.N
        N_s = prob.N_s
        N_t = prob.N_t
        N_st = prob.N_st
        A = prob.A
        T = prob.T
        Out_i = prob.Out_i
        In_i = prob.In_i
        K = prob.K
        Q = prob.Q
        s = 0
        t = len(N)+1
        q = prob.q

        cpx.variables.add(
            names = [f'x_{i}_{j}' for i,j in A],
            types = ['B'] * len(A),
            obj = [float(T[i,j]) for i,j in A]
        )


        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in Out_i[i]], ([1] * (len(Out_i[i])))
                ) for i in N
            ],
            senses=['E'] * len(N),
            rhs=[1] * len(N),
            names=[f'out_degree_{i}' for i in N])


        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for i in In_i[j]], ([1] * (len(In_i[j])))
                ) for j in N
            ],
            senses=['E'] * len(N),
            rhs=[1] * len(N),
            names=[f'in_degree_{j}' for j in N])


        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{s}_{j}' for j in Out_i[s]], ([1] * (len(Out_i[s])))
                ) 
            ],
            senses=['L'],
            rhs=[K],
            names=['num_vehicles'])
        
        self.cpx = cpx


    def solve(self, timelimit=3600, cpu=4, sep_frac_sols=False, max_num_cuts=10):

        cpx = self.cpx
        prob = self.prob

        start_time = time.time()
        cpx.parameters.timelimit.set(timelimit)

        num_threads = cpu
        cpx.parameters.threads.set(num_threads)


        cb = CVRPCallback(num_threads, prob, max_num_cuts=max_num_cuts)
        contextmask = cplex.callbacks.Context.id.thread_up
        contextmask |= cplex.callbacks.Context.id.thread_down
        contextmask |= cplex.callbacks.Context.id.candidate
        if sep_frac_sols:
            contextmask |= cplex.callbacks.Context.id.relaxation
        cpx.set_callback(cb, contextmask)


        cpx.solve()


        num_bnb_nodes = cplex._internal._procedural.getnodecnt(cpx._env._e, cpx._lp)
        num_gap = cplex._internal._procedural.getmiprelgap(cpx._env._e, cpx._lp)
        solve_time = time.time() - start_time

        results = {
            'status': cpx.solution.get_status_string(),
            'obj': cpx.solution.get_objective_value() / prob.round_mult,
            'time': solve_time,
            'num_bnb_nodes': num_bnb_nodes,
            'num_gap': num_gap,
            'cpu': cpu,
            'timelimit': timelimit,
            'num_sep_called': cb.num_separation_called,
            'num_added_cuts': cb.num_added_cuts,
            'total_sep_time': cb.total_separation_time,
            'solution': [(i,j) for (i,j) in self.prob.A if cpx.solution.get_values(f'x_{i}_{j}') > 0.9]

        }

        return results


class SepWorker:

    def __init__(self, prob, max_num_cuts=10):
        self.prob = prob
        self.s = 0
        self.t = len(prob.N)+1

        # CVEPSEP-related stuffs
        p_OldCutsCMP = ffi.new('CnstrMgrRecord **')
        p_NewCutsCMP = ffi.new('CnstrMgrRecord **')

        lib.CMGR_CreateCMgr(p_OldCutsCMP, 100)
        lib.CMGR_CreateCMgr(p_NewCutsCMP, 100)

        OldCutsCMP = p_OldCutsCMP[0]
        NewCutsCMP = p_NewCutsCMP[0]

        self.p_OldCutsCMP = p_OldCutsCMP
        self.p_NewCutsCMP = p_NewCutsCMP

        self.OldCutsCMP = OldCutsCMP
        self.NewCutsCMP = NewCutsCMP

        self.p_IntegerAndFeasible = ffi.new('char *')
        self.NoOfCustomers = len(prob.N)
        self.CAP = prob.Q
        self.p_Demand = ffi.new('int[]', list(prob.q[:-1]))
        self.MaxNoOfCuts = max_num_cuts
        self.EpsForIntegrality = 0.01
        self.p_MaxViolation = ffi.new('double *')


    def separate_capcuts(self, positive_xstar):

        # Call CVRPSEP's CAPSEP_SeparateCapCuts()

        NoOfEdges = len(positive_xstar)
        p_EdgeTail = ffi.new('int[]', NoOfEdges+1)
        p_EdgeHead = ffi.new('int[]', NoOfEdges+1)
        p_EdgeX = ffi.new('double[]', NoOfEdges+1)

        p_EdgeTail[0] = 0
        p_EdgeHead[0] = 0
        p_EdgeX[0] = 0

        for idx, ((i,j), val) in enumerate(positive_xstar):
            p_EdgeTail[idx+1] = i if i < self.t else 0
            p_EdgeHead[idx+1] = j if j < self.t else 0
            p_EdgeX[idx+1] = val


        lib.CAPSEP_SeparateCapCuts(
            self.NoOfCustomers,
            self.p_Demand,
            self.CAP,
            NoOfEdges,
            p_EdgeTail,
            p_EdgeHead,
            p_EdgeX,
            self.OldCutsCMP,
            self.MaxNoOfCuts,
            self.EpsForIntegrality,
            self.p_IntegerAndFeasible,
            self.p_MaxViolation,
            self.NewCutsCMP
        )

        cuts = []

        if self.NewCutsCMP.Size == 0:
            return cuts

        for idx in range(self.NewCutsCMP.Size):
            S = set()
            cut_info = self.NewCutsCMP[0].CPL[idx]
            if cut_info.CType == lib.CMGR_CT_CAP:
                for j in range(1, cut_info.IntListSize+1):
                    S.add(cut_info.IntList[j])

                LHS = [(i,j) for i in S for j in S if i!=j]
                RHS = round(cut_info.RHS)

                if len(S) > 0: # Sometimes CVRPSEP returns empry sets... is this a bug??????
                    cuts.append((LHS, RHS, S))

        for idx in range(self.NewCutsCMP.Size):
            lib.CMGR_MoveCnstr(self.NewCutsCMP, self.OldCutsCMP, idx, 0)

        self.NewCutsCMP.Size = 0

        return cuts





class CVRPCallback():

    def __init__(self, num_threads, prob, max_num_cuts=10):
        self.num_threads = num_threads
        self.prob = prob
        self.max_num_cuts = max_num_cuts
        self.workers = [None] * num_threads
        self.addedcuts = []

        self.num_separation_called = 0
        self.num_added_cuts = 0
        self.total_separation_time = 0
        self.generated_cuts = []

    def get_positive_xstar(self, get_sol_value_fn):
        # Convert the solution to symmetry form...
        return [
            ((i,j), get_sol_value_fn(f'x_{i}_{j}') + get_sol_value_fn(f'x_{j}_{i}')) for i in self.prob.N for j in self.prob.N if i<j and (get_sol_value_fn(f'x_{i}_{j}') + get_sol_value_fn(f'x_{j}_{i}') > 0.001)
        ]

    def separate_user_cuts(self, context, worker):

        get_sol_value_fn = context.get_relaxation_point
        positive_xstar = self.get_positive_xstar(get_sol_value_fn)

        cuts = worker.separate_capcuts(positive_xstar)

        if len(cuts) > 0:
            cutmanagement = cplex.callbacks.UserCutCallback.use_cut.purge
            context.add_user_cuts(
                cuts=[
                        cplex.SparsePair(
                            [f'x_{i}_{j}' for (i,j) in LHS], [1] * len(LHS)
                        )
                        for LHS, RHS, S in cuts
                ],
                senses=['L'] * len(cuts), 
                rhs=[RHS for LHS, RHS, S in cuts],
                cutmanagement=[cutmanagement] * len(cuts), 
                local=[False] * len(cuts)
            )

            self.num_added_cuts += len(cuts)

    def separate_lazy_constraints(self, context, worker):

        get_sol_value_fn = context.get_candidate_point
        positive_xstar = self.get_positive_xstar(get_sol_value_fn)

        cuts = worker.separate_capcuts(positive_xstar)
        
        if len(cuts) > 0:
            context.reject_candidate(
                constraints=[
                    cplex.SparsePair(
                        [f'x_{i}_{j}' for (i,j) in LHS], [1] * len(LHS)
                    )
                    for LHS, RHS, S in cuts
                ],
                senses=['L'] * len(cuts), 
                rhs=[RHS for LHS, RHS, S in cuts]
            )

            self.num_added_cuts += len(cuts)

    def invoke(self, context):
        """Whenever CPLEX needs to invoke the callback it calls this
        method with exactly one argument: an instance of
        cplex.callbacks.Context.
        """
        try:
            thread_id = context.get_int_info(
                cplex.callbacks.Context.info.thread_id)
            if context.get_id() == cplex.callbacks.Context.id.thread_up:
                self.workers[thread_id] = SepWorker(self.prob, self.max_num_cuts)
            elif context.get_id() == cplex.callbacks.Context.id.thread_down:
                self.workers[thread_id] = None
            elif context.get_id() == cplex.callbacks.Context.id.relaxation:
                self.num_separation_called += 1
                start_time = time.time()
                self.separate_user_cuts(context, self.workers[thread_id])
                self.total_separation_time += time.time() - start_time
            elif context.get_id() == cplex.callbacks.Context.id.candidate:
                self.num_separation_called += 1
                start_time = time.time()
                self.separate_lazy_constraints(context, self.workers[thread_id])
                self.total_separation_time += time.time() - start_time
            else:
                print("Callback called in an unexpected context {}".format(
                    context.get_id()))
        except:
            info = sys.exc_info()
            print('#### Exception in callback: ', info[0])
            print('####                        ', info[1])
            print('####                        ', info[2])
            traceback.print_tb(info[2], file=sys.stdout)
            raise



