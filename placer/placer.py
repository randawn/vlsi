#! /usr/bin/python

import copy
import re
import numpy
from numpy import array
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from itertools import combinations

class NetList(object):
    def __init__(self):
        #self.gate = {}#{'nets':[(wire, weight)], 'coord':(None, None)}
        #self.pad = {}#{'wire':None, 'coord':(None, None)}
        self.gates = {}#{'ID':gate, }
        self.pads = {}#{'ID':wire, }
        self.gatesk = {}
        self.x = []
        self.y = []

    def AddGate(self, gate_ID, gate):
        self.gates[gate_ID] = gate

    def AddPad(self, pad_ID, pad):
        self.pads[pad_ID] = pad

    def DisplayGate(self):
        for gate in self.gate_lst:
            print gate

    def DisplayPad(self):
        for pad in self.pad_lst:
            print pad

    def InputParser(self):
        gate = {}
        pad = {}
        fh = open('struct')   # no comment for parsing!!!
        pattern0 = re.compile("(\d+) (\d+)")
        pattern1 = re.compile("(\d+)")  #TODO: we do know how many of \d+
        pattern2 = re.compile("(\d+)")
        pattern3 = re.compile("(\d+) (\d+) (\d+) (\d+)")
        first_line = fh.readline()
        match = pattern0.search(first_line)
        (gate_sum, net_sum) = match.groups()
        # first line: two int G: num of Gates N:num of nets
        for line_cnt in range(int(gate_sum)):
            line = fh.readline()
            all_lst = pattern1.findall(line)
            gate_ID = all_lst[0]
            for i in range(int(all_lst[1])):
                try:
                    gate['nets'].append((all_lst[2+i], 1))
                except KeyError as e:
                    gate['nets'] = [(all_lst[2+i], 1)]
            self.AddGate(gate_ID, gate)
            gate = {}
        #next line: one int P: num of pad P
        next_line = fh.readline()
        match = pattern2.search(next_line)
        (pad_sum,) = match.groups()
        for line_cnt in range(int(pad_sum)):
            line = fh.readline()
            match = pattern3.search(line)
            (pin_ID, net_num, X, Y) = match.groups()
            pad_ID = pin_ID
            pad['wire'] = net_num
            pad['coord'] = (int(X), int(Y))
            self.AddPad(pad_ID, pad)
            pad = {}
        fh.close()
        self.gatesk = copy.deepcopy(self.gates)
        self.KTo2(int(net_sum))

    def KTo2(self, net_sum):
        '''trans k-point net to 2-point'''
        gates = self.gates
        # decompose
        for i in range(1,net_sum+1):
            gate_rcd = []
            net_rcd = []
            for ID, gate in gates.items():
                for net_name, weight in gate['nets']:
                    if (str(i) == net_name):
                        gate_rcd.append(ID)
                        net_rcd.append(net_name)
            if (len(gate_rcd) > 2):
                weight = 1.0/(len(gate_rcd)-1)
                # add new
                for (id0, id1) in list(combinations(gate_rcd, 2)):
                    net_name = str(id0)+'_'+str(id1)+'tmp'
                    gates[id0]['nets'].append((net_name, weight))
                    gates[id1]['nets'].append((net_name, weight))
                # remove old
                for i in range(len(gate_rcd)):
                    ID = gate_rcd[i]
                    net_name = net_rcd[i]
                    gates[ID]['nets'].remove((net_name, 1))
        #merge
        for ID, ID_ex in list(combinations(gates.keys(), 2)):
            net_rcd = []
            weight_rcd = []
            for net, weight in gates[ID]['nets']:
                for net_ex, weight_ex in gates[ID_ex]['nets']:
                    if (net == net_ex):
                        net_rcd.append(net)
                        weight_rcd.append(weight)
            if (len(net_rcd) > 1):
                # add new
                net_name = str(ID)+'_'+str(ID_ex)
                weight_sum = sum(weight_rcd)
                gates[ID]['nets'].append((net_name, weight_sum))
                gates[ID_ex]['nets'].append((net_name, weight_sum))
                # remove old
                for i in range(len(net_rcd)):
                    try:
                        gates[ID]['nets'].remove((net_rcd[i], weight_rcd[i]))
                        gates[ID_ex]['nets'].remove((net_rcd[i], weight_rcd[i]))
                    except ValueError as e:
                        pass

    def OutputGen(self):
        x = self.x
        y = self.y
        fh = open('toy1.result','w')
        for i in range(len(x)):
            fh.write(str(i+1)+' '+str(x[i])+' '+str(y[i])+'\n')
        fh.close()
        

    def GetMatC(self):
        gates = self.gates
        R_lst = []
        C_lst = []
        V_lst = []
        for ID, gate in gates.items():
            for net, weight in gate['nets']:
                for ID_ex, gate_ex in gates.items():
                    if (ID == ID_ex):
                        continue
                    for net_ex, weight_ex in gate_ex['nets']:
                        if (net == net_ex):
                            R_lst.append(int(ID)-1)
                            C_lst.append(int(ID_ex)-1)
                            V_lst.append(weight)
        return (R_lst, C_lst, V_lst)

    def GetMatA(self):
        pads = self.pads
        gatesk = self.gatesk
        (R_lst, C_lst, V_lst) = self.GetMatC()
        V_lst_n = [-1*value for value in V_lst]
        rows = list(set(R_lst))
        for row in rows:
            r_sum = 0
            for i in range(len(R_lst)):
                if (row == R_lst[i]):
                    r_sum += V_lst[i]
            for padID, pad in pads.items():
                for net, weight in gatesk[str(row+1)]['nets']:
                    if (net == pad['wire']):
                        r_sum += 1.0
            R_lst.append(row)
            C_lst.append(row)
            V_lst_n.append(r_sum)
        #for i in range(len(R_lst)):
        #    try:
        #        Diag_dict[R_lst[i]] += V_lst[i]
        #    except KeyError as e:
        #        Diag_dict[R_lst[i]] = V_lst[i]
        #for diag in Diag_dict.keys():
        #    for pad_ID, pad in pads.items():
        #        if (diag == pad['wire']):
        #            Diag_dict[diag] += 1.0
        #for diag in Diag_dict.keys():
        #    R_lst.append(diag)
        #    C_lst.append(diag)
        #    V_lst_n.append(Diag_dict[diag])
        R = array(R_lst)
        C = array(C_lst)
        V = array(V_lst_n)
        return (R, C, V)

    def GetVecb(self):
        gates = self.gatesk
        pads = self.pads
        bx = [0]*len(gates)
        by = [0]*len(gates)
        for pad_ID, pad in pads.items():
            for gate_ID, gate in gates.items():
                for net, weight in gate['nets']:
                    if (pad['wire'] == net):
                        bx[int(gate_ID)-1] += pad['coord'][0]
                        by[int(gate_ID)-1] += pad['coord'][1]
        return (array(bx), array(by))

    def SolveMat(self):
        (R, C, V) = self.GetMatA()
        (bx, by) = self.GetVecb()
        A = coo_matrix((V, (R, C)), shape=(len(bx), len(bx)))
        # convert to csr format for efficiency
        x = spsolve(A.tocsr(), bx)
        y = spsolve(A.tocsr(), by)
        self.x = list(x)
        self.y = list(y)
        #print "x = ", x
        #print "y = ", y

if __name__ == '__main__':
    my_netlist = NetList()
    my_netlist.InputParser()
    #my_netlist.DisplayGate()
    #my_netlist.DisplayPad()
    my_netlist.SolveMat()
    my_netlist.OutputGen()

