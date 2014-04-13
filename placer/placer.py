#! /usr/bin/python

import re
import numpy
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spslove

#class NetList(object):
class NetList(object):
    def __init__(self):
        #self.gate = {'ID':None, 'nets':[], 'coord':(None, None)}
        #self.pad = {'ID':None, 'wire':None, 'coord':(None, None)}
        self.gate = {}
        self.pad = {}
        self.gate_lst = []
        self.pad_lst = []

    def AddGate(self, gate):
        self.gate_lst.append(gate)

    def AddPad(self, pad):
        self.pad_lst.append(pad)

    def DisplayGate(self):
        for gate in self.gate_lst:
            print gate

    def DisplayPad(self):
        for pad in self.pad_lst:
            print pad

    def InputParser(self):
        gate = self.gate
        pad = self.pad
        fh = open('toy1')   # no comment for parsing!!!
        # first line: two int G: num of Gates N:num of nets
        pattern0 = re.compile("(\d+) (\d+)")
        pattern1 = re.compile("(\d+)")  #TODO: we do know how many of \d+
        pattern2 = re.compile("(\d+)")
        pattern3 = re.compile("(\d+) (\d+) (\d+) (\d+)")
        first_line = fh.readline()
        match = pattern0.search(first_line)
        (gate_sum, net_sum) = match.groups()
        for line_cnt in range(int(gate_sum)):
            line = fh.readline()
            all_lst = pattern1.findall(line)
            gate['ID'] = int(all_lst[0])
            for i in range(2+int(all_lst[1])):
                try:
                    gate['nets'].append(all_lst[i])
                except KeyError as e:
                    gate['nets'] = [all_lst[i]]
            self.AddGate(gate)
            gate = {}
        #next line: one int P: num of pad P
        next_line = fh.readline()
        match = pattern2.search(next_line)
        (pad_sum,) = match.groups()
        for line_cnt in range(int(pad_sum)):
            line = fh.readline()
            match = pattern3.search(line)
            (pin_ID, net_num, X, Y) = match.groups()
            pad['ID'] = int(pin_ID)
            pad['wire'] = int(net_num)
            pad['coord'] = (int(X), int(Y))
            self.AddPad(pad)
            pad = {}
    def OutputGen():
        pass

    def GetMatC(self):

    def GetMatA(self):





if __name__ == '__main__':
    my_netlist = NetList()
    my_netlist.InputParser()
