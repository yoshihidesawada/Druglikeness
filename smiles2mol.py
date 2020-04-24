import os
import sys
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import QED, rdMolDescriptors
import sascorer
import networkx as nx
from rdkit.Chem import rdmolops

def rule_of_five(m):
    #分子量
    mw = Descriptors.MolWt(m)
    #値が高いほど疎水性が高い
    logp = Descriptors.MolLogP(m)
    #水素結合ドナー（OH，NHの数）
    hbd = rdMolDescriptors.CalcNumLipinskiHBD(m)
    #水素結合アクセプター（O，N原子の数）
    hba = rdMolDescriptors.CalcNumLipinskiHBA(m)
    if (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10):
        return 1
    else:
        return 0

def ComputeRewardConvert(data, num):

    print(num, data)
    m = Chem.MolFromSmiles(data)

    if m is not None:
        #値が高いほど疎水性が高い
        logp = Descriptors.MolLogP(m)
        print("logP:%f"%(logp))

        # 分子構造の複雑さを基準とした合成難易度
        SA_score = -sascorer.calculateScore(m)
        print("SA:%f"%(SA_score))

        # ring penality. 6 is hyperparameter
        cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(m)))
        if len(cycle_list) == 0:
            cycle_length = 0
        else:
            cycle_length = max([ len(j) for j in cycle_list ])
        if cycle_length <= 6:
            cycle_length = 0
        else:
            cycle_length = cycle_length - 6
        cycle_score = -cycle_length
        print("cycle:%d"%(cycle_score))

        qed = QED.qed(m)
        print("QED:%f"%(qed))

        flag = rule_of_five(m)
        print("Rule of Five:%d"%(flag))

        savefile = 'mol%d_logP%f_SA%f_cycle%d_QED%f_ruleoffive%d.png'%(num,logp,SA_score,cycle_score,qed,flag)
        Draw.MolToFile(m,savefile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='convert smiles to png')
    parser.add_argument('--file', type=str, default=None, metavar='type',help='filename (default: None)')
    args = parser.parse_args()

    if args.file:
        data = pd.read_csv(args.file,delimiter=',',engine="python",header=None)
        data = np.array(data)
        for i in range(len(data)):
            data[i][0] = data[i][0].lstrip('\'')
            data[i][0] = data[i][0].rstrip('\'')
            ComputeRewardConvert(data[i][0], i)
    else:
        data = 'OCc1ccccc1'
        ComputeRewardConvert(data, 0)
