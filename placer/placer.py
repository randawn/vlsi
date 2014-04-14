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
        self.gates = {}
        self.pads = {}
        #self.x = []
        #self.y = []
    def PrintGates(self, gates):
        for ID, gate in gates.items():
            print ID
            print gate
    def PrintPads(self, pads):
        for ID, pad in pads.items():
            print ID
            print pad
    def Place(self):
        gates = self.gates
        pads = self.pads
        (x, y) = self.SolveMat(gates, pads)
        print x
        print y
        # add coord TODO make it nicer
        for ID, gate in gates.items():
            gate['coord'] = (x[int(ID)-1], y[int(ID)-1])
        (gatesL, gatesR, padsL, padsR) = self.Assign(x, y, gates, pads)
        pads = copy.deepcopy(padsL)
        padsC = copy.deepcopy(padsR)
        pads_L = self.Contain(gatesL, gatesR, pads, padsC)
        pads = copy.deepcopy(padsR)
        padsC = copy.deepcopy(padsL)
        pads_R = self.Contain(gatesR, gatesL, pads, padsC)
        (xL, yL) = self.SolveMat(gatesL, pads_L)
        (xR, yR) = self.SolveMat(gatesR, pads_R)
        print xL
        print yL
        print xR
        print yR

    def Assign(self, x, y, gates, pads):
        tuple_for_sort = [ (x[i], y[i], i) for i in range(len(x))]
        tuple_for_sort.sort(key=lambda item: (item[0], item[1]))
        xsorted = [ str(i+1) for xi, yi, i in tuple_for_sort]
        left_lst = xsorted[:len(xsorted)/2]
        gatesL = {}
        gatesR = {}
        for gate_ID in gates.keys():
            if (gate_ID in left_lst):
                gatesL[gate_ID] = gates[gate_ID]
            else:
                gatesR[gate_ID] = gates[gate_ID]
        padsL = {}
        padsR = {}
        for pad_ID, pad in pads.items():
            if (pad['coord'][0]<50):
                padsL[pad_ID] = pads[pad_ID]
            else:
                padsR[pad_ID] = pads[pad_ID]
        return (gatesL, gatesR, padsL, padsR)

    def Contain(self, gates, gatesC, pads, padsC):
        '''propagate the \w+C to \w+ part '''
        # trans right side gate to pad
        for gate_IDC, gateC in gatesC.items():
            pad = {}
            gate_rcd = []
            for wire, weight in gateC['nets']:
                for gate_ID, gate in gates.items():
                    if (wire in [ net[0] for net in gate['nets']]):
                        gate_rcd.append(gate_ID)
            if (len(gate_rcd)>0):
                pad['gates'] = gate_rcd
                pad['coord'] = gateC['coord']
                padsC[gate_IDC+'gate'] = pad
        # do the propagate
        for pad_ID, pad in padsC.items():
            # we do not kick out the one have no connection on the other side
            pad['coord'] = (50, pad['coord'][1])
            pads['P'+pad_ID] = pad
        return pads

    def DisplayGate(self):
        for gate in self.gate_lst:
            print gate

    def DisplayPad(self):
        for pad in self.pad_lst:
            print pad

    def InputParser(self):
        '''get the pads and gates from file'''
        # gates            = {'ID':{'num':num, 'nets':[(wire, weight)], 'coord':(X, Y)}}
        # pads(in file)    = {'ID':{'wire':wire,      'coord':(X, Y)}}
        # pads(after Kto2) = {'ID':{'gates':[gate_ID], 'coord':(X, Y)}}
        gates = {}
        pads = {}
        fh = open('toy1')   # no comment for parsing!!!
        pattern0 = re.compile("(\d+) (\d+)")
        pattern1 = re.compile("(\d+)")  #TODO: we do know how many of \d+
        pattern2 = re.compile("(\d+)")
        pattern3 = re.compile("(\d+) (\d+) (\d+) (\d+)")
        first_line = fh.readline()
        match = pattern0.search(first_line)
        (gate_sum, net_sum) = match.groups()
        # first line: two int G: num of Gates N:num of nets
        for line_cnt in range(int(gate_sum)):
            gate = {}
            line = fh.readline()
            all_lst = pattern1.findall(line)
            gate_ID = all_lst[0]
            for i in range(int(all_lst[1])):
                try:
                    gate['nets'].append((all_lst[2+i], 1))
                except KeyError as e:
                    gate['nets'] = [(all_lst[2+i], 1)]
            gates[gate_ID] = gate
        #next line: one int P: num of pad P
        next_line = fh.readline()
        match = pattern2.search(next_line)
        (pad_sum,) = match.groups()
        for line_cnt in range(int(pad_sum)):
            pad = {}
            line = fh.readline()
            match = pattern3.search(line)
            (pin_ID, net_num, X, Y) = match.groups()
            pad_ID = pin_ID
            pad['wire'] = net_num
            pad['coord'] = (int(X), int(Y))
            pads[pad_ID] = pad
        fh.close()
        self.KTo2(int(net_sum), gates, pads)

    def KTo2(self, net_sum, gates, pads):
        '''trans k-point netlist to 2-point'''
        #trans pads first
        for pad_ID, pad in pads.items():
            pad['gates'] = []
            for gate_ID, gate in gates.items():
                for wire, weight in gate['nets']:
                    if (pad['wire'] == wire):# this gate and pad is connected
                        pad['gates'].append(gate_ID)
            del pad['wire']
        # now trans k-point netlist to 2-point
        # decompose the wire connected with more the 2 gates
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
                    net_name = str(id0)+'_'+str(id1)+'decompose'
                    gates[id0]['nets'].append((net_name, weight))
                    gates[id1]['nets'].append((net_name, weight))
                # remove old
                for i in range(len(gate_rcd)):
                    ID = gate_rcd[i]
                    net_name = net_rcd[i]
                    gates[ID]['nets'].remove((net_name, 1))
        # merge if there are more than 1 wire connect two gates 
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
                net_name = str(ID)+'_'+str(ID_ex)+'merge'
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
        self.gates = gates
        self.pads = pads

    def OutputGen(self, x, y):
        fh = open('toy1.result','w')
        for i in range(len(x)):
            fh.write(str(i+1)+' '+str(x[i])+' '+str(y[i])+'\n')
        fh.close()

    def GetMatC(self, gates):
        R_lst = []
        C_lst = []
        V_lst = []
        # assign a num
        IDseq = sorted(gates.keys(), key=lambda id: int(id));
        for i in range(len(IDseq)):
            gates[IDseq[i]]['num'] = i
        for ID, gate in gates.items():
            for net, weight in gate['nets']:
                for ID_ex, gate_ex in gates.items():
                    if (ID == ID_ex):
                        continue
                    for net_ex, weight_ex in gate_ex['nets']:
                        if (net == net_ex):
                            R_lst.append(gate['num'])
                            C_lst.append(gate_ex['num'])
                            V_lst.append(weight)
        return (R_lst, C_lst, V_lst)

    def GetMatA(self, gates, pads):
        (R_lst, C_lst, V_lst) = self.GetMatC(gates)
        V_lst_n = [-1*value for value in V_lst]
        rows = list(set(R_lst))         # unique the row
        for row in rows:                # gate_ID = row + 1
            r_sum = 0
            for i in range(len(R_lst)): # add up every row
                if (row == R_lst[i]):
                    r_sum += V_lst[i]
            for padID, pad in pads.items():# add pad if there is 
                if (str(row+1) in pad['gates']):
                    r_sum += 1.0
            R_lst.append(row)
            C_lst.append(row)
            V_lst_n.append(r_sum)
        R = array(R_lst)
        C = array(C_lst)
        V = array(V_lst_n)
        return (R, C, V)

    def GetVecb(self, gates, pads):
        bx = [0]*len(gates)
        by = [0]*len(gates)
        print '==============='
        #self.PrintGates(gates)
        #self.PrintPads(pads)
        for pad_ID, pad in pads.items():
            for gate_ID, gate in gates.items():
                if (gate_ID in pad['gates']):
                    bx[gate['num']] += pad['coord'][0]
                    by[gate['num']] += pad['coord'][1]
                    #by[int(gate_ID)-1] += pad['coord'][1]
        return (array(bx), array(by))

    def SolveMat(self, gates, pads):
        (R, C, V) = self.GetMatA(gates, pads)
        (bx, by) = self.GetVecb(gates, pads)
        A = coo_matrix((V, (R, C)), shape=(len(bx), len(bx)))
        # convert to csr format for efficiency
        x = spsolve(A.tocsr(), bx)
        y = spsolve(A.tocsr(), by)
        return (list(x), list(y))
        #print "x = ", x
        #print "y = ", y

if __name__ == '__main__':
    my_netlist = NetList()
    my_netlist.InputParser()
    #my_netlist.DisplayGate()
    #my_netlist.DisplayPad()
    my_netlist.Place()

