from random import uniform
from math import log, exp, sqrt

state = (0, 1, 2, 3)

rates0 = [[1 for x in state] for y in state]
for x in state:
    rates0[x][x] = 0
   
        
def diag(array):
    
    diag = [[0 for x in array] for y in array]
    
    for x in range(len(array)):
        
        diag[x][x] = array[x]
        
    return diag
    
    
def norm(data):
 
    norm = sum(data)
    return [x / norm for x in data]
    

def steady(rates):

    import numpy as np

    exits = [sum(rates[x]) for x in state]
    steady = []

    for x in state:
        R = np.matrix(rates)-np.diag(exits)
        R = np.delete(R, x, 0)
        R = np.delete(R, x, 1)
        steady.append(-np.linalg.det(R))
        
    steady = [x / sum(steady) for x in steady]

    return steady
    
def chop(expr, delta = 10**-10):
    if isinstance(expr, (int, float, complex)):
        return 0 if -delta <= expr <= delta else expr
    else:
        return [chop(x) for x in expr]
        
        
def x_choice(prob):
    
    l = len(prob)
    s = [sum(prob[0:x]) for x in range(0,l+1)]
    y = uniform(0, 1)
    for x in range(0,l):
        if y >= s[x] and y < s[x+1]:
            break
    return x

def t_choice(x, exits):

    return log(1 / uniform(0, 1)) / exits[x]


class Process:
        
    def __init__(self, rates = rates0, start = None):

        self.rates = [x.copy() for x in rates]
        self.type = type(start)
        
        if start is not None:
           
            if self.type is list:
                norm = sum(start)
                self.start = [x / norm for x in start]
            elif self.type is int:
                if start in state:
                    self.init = start
                    self.start = [0 for x in state]
                    self.start[start] = 1
                else:
                    raise TypeError("Only states in " + str(state) + " allowed")
                                    
        else:
            
            self.start = steady(self.rates)
            self.type = list
            
        self.exits = [sum(self.rates[x]) for x in state]
        self.chain = [[self.rates[x][y] / self.exits[x] for y in state] for x in state]   
        self.generator = [[self.rates[x][y] for y in state] for x in state]  
        for x in state:
            self.generator[x][x] = - self.exits[x]
           
            
    def steady(self):
        
        return steady(self.rates)
            
    def set_lvec(self):
        
        self.lvec = [1 for x in state]
        
    def set_taboo(self, transitions):
        
        self.taboo = [x.copy() for x in self.generator]
        for k in transitions:
            self.taboo[k[0]][k[1]] = 0
        
    def tilted(self, q, transitions = [(0, 1)] , weights = [1]):
           
        tilted = [x.copy() for x in p.generator]
                       
        for i in range(len(transitions)):
            
            da = transitions[i][0]
            ad = transitions[i][1]
            tilted[ad][da] *= exp(q) * weights[i]
        
        return tilted

    def scgf(self, q, transitions = [(0, 1)] , weights = [1]):
                
        import sympy
        
        qq = sympy.Symbol('qq')        
        tilted = [x.copy() for x in self.generator]
                       
        for i in range(len(transitions)):
            
            da = transitions[i][0]
            ad = transitions[i][1]
            tilted[ad][da] *= sympy.exp(qq) * weights[i]
        
        m = sympy.Matrix(tilted)  

        for x in state:
            dominant = m.eigenvals(multiple = list)[x] 
            n = float(sympy.re(dominant.subs({qq:0}).evalf())) 
            if chop(n) == 0:
                break

        return sympy.re(dominant.subs({qq:q}).evalf())
    
       
    def normcomm(self, q, transitions = [(0, 1)], weights = [1]):
        
        import numpy as np
        
        a1 = np.matrix(self.tilted(q, transitions, weights))
        b1 = np.matrix(diag(self.exits))
        c1 = np.matmul(a1, b1) - np.matmul(b1, a1)
        
        a2 = np.matrix(p.generator)
        b2 = np.matrix(diag(self.exits))
        c2 = np.matmul(a2, b2) - np.matmul(b2, a2)
        
        return np.linalg.norm(c1) / np.linalg.norm(c2)

    
    
def path(n_or_t, p = Process()):
    
    '''Accepts both an integer number of steps or a real-valued final time for stopping.'''
    
    states = []
    times = []
    x = x_choice(p.start) if p.type is list else p.init

    if type(n_or_t) is int:

        for k in range(n_or_t):
            
            t = t_choice(x, p.exits)
            states.append(x)
            times.append(t)          
            x = x_choice(p.chain[x])  

    elif type(n_or_t) is float:
        
        t_tot = 0
        
        while t_tot < n_or_t:
            t = t_choice(x, p.exits)         
            states.append(x)
            times.append(t)  
            x = x_choice(p.chain[x])
            t_tot += t
            
        times[-1] += n_or_t - t_tot

    return states, times
    

class Observe:
    
    def __init__(self, n_or_t, p = Process()):
        
        self.path = path(n_or_t, p)
        self.states = self.path[0]
        self.times = self.path[1]
        self.steps = len(self.states)
        self.duration = sum(self.times)

    def souj(self):
        
        souj = [0 for x in state]
        
        for k in range(self.steps):
            souj[self.states[k]] += self.times[k] 
                        
        return [x / self.duration for x in souj]
    
    def fluxes(self, transitions = None):
        
        if transitions == None:
            transitions = [(x, y) for x in state for y in state if x != y]
            
        fluxes = {}
        for k in range(self.steps - 1):
            trans = (self.states[k], self.states[k + 1])
            if trans in transitions:
                fluxes[trans] = fluxes.get(trans, 0) + 1

        return fluxes      
        
    def currents(self, transitions = None):
        
        if transitions == None:
            transitions = [(x, y) for x in state for y in state if x < y]

            
        sym_trans = transitions.copy()
        for k in transitions:
            sym_trans.append((k[1], k[0]))
            
        f = self.fluxes(sym_trans)
        
        c = {}
        for k in transitions:
            c[k] = f[(k[0], k[1])] - f[(k[1], k[0])]
        
        return c

    
### FROM HERE ON TO BE FIXED
    
    def emp_ep(self):
        
        emp_ep = 0
        f = self.fluxes()
        for x in state:
            for y in state:
                if y != x:
                    emp_ep += f[x][y] * log(f[x][y] / f[y][x])
        return emp_ep
    
    
        
    def emp_gen(self):
        
        emp_gen = [[0 for x in state] for y in state]
        f = self.fluxes()
        
        for k in f:
            emp_gen[k[0]][k[1]] = f[k] / self.steps # Notice transposition
            
        return emp_gen
    
    def cycle(self, name): ### THERE ARE STILL PROBLEMS WITH THIS FUNCTION

        # This is to account for all cyclic permutations of the states:
        
        letters = list(name)
        l = len(letters)
        name_perm = [[letters[i - j] for i in range(l)] for j in range(l)]
        name_list = ["".join(i) for i in name_perm]

        x_list = []
        t_list = []
        times = []
        
        for xt in self.path:
            
            x, t = xt[0], xt[1]

            if x in x_list:
                i = x_list.index(x)
                cname = ''
                ctime = 0
                for k in range(i,len(x_list)):
                    cname += str(x_list.pop(k))
                    ctime += t_list.pop(k)
                if cname in name_list:
                    times.append(ctime)

            x_list.append(x)
            t_list.append(t)
                
        return times

