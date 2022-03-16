#!/usr/bin/env python
import sys
import argparse
import numpy as np
from ubergauss import tools
from scipy.sparse import csr_matrix, hstack, load_npz, save_npz, vstack
partner =  {a:b for a,b in zip("({[<",")}]>")  }
import pandas as pd
import eden.graph as eg
import networkx as nx
from collections import defaultdict
from lmz import *
import subprocess
import os
import rrieval.lib as rl



def mkgraph(sequence, structure):
    graph = nx.Graph()
    lifo = defaultdict(list)
    cut = structure.index("&")

    structure = structure.replace("&","")
    sequence = sequence.replace("&","")
    for i,(s,n) in enumerate(zip(structure, sequence)):
        graph.add_node(i, label=n)
        if i > 0 and  i != cut:
            graph.add_edge(i, i-1, label='-')

        # ADD PAIRED BASES
        if s in ['(','[','<']:
            lifo[partner[s]].append(i)
        if s in [')',']','>']:
            j = lifo[s].pop()
            graph.add_edge(i, j, label='=')
    return graph
    #return eg.vectorize([graph], discrete = False) # keep here in case i want nested edges ...

def mkgr(x):
    return mkgraph(*x)


def convert(negname, posname, outname, graphfeatures=True):
    # d1 = cri.loadDF(negname)
    # d2 = cri.loadDF(posname)
    d1 = pd.read_csv(negname)
    d2 = pd.read_csv(posname)

    X, y = rl.read_pos_neg_data(posname, negname)
    #print(X.info())

    if graphfeatures:
        # makes list of subseqDP and hybridDP tupels
        data = [a for a in zip(X['subseqDP'],X['hybridDP'])]
        print(len(data))
        graphs = tools.xmap(mkgr, data,32)
        #print(graphs)
        #X2 = eg.vectorize(graphs)
        X2 = csr_matrix(vstack(tools.xmap(eg.vectorize,[[g] for g in graphs])))
        #print(len(X2))
        df_X = pd.DataFrame(X2.todense())
        print(len(df_X))

        #X = pd.concat([X, df_X], axis=1)
        X.join(df_X)
    X = X.drop(columns="subseqDP")
    X = X.drop(columns="hybridDP")

    # Whath dose Range do?
    #list = [ X,y, col_namez + Range(X.shape[1] - len(col_namez))]
    # saves object as 'np.savez_compressed'
    #tools.ndumpfile(list , outname)
    return X,y



def main():
    #parser = argparse.ArgumentParser(description='')
    #parser.add_argument("-i", "--input_file",
    #                    help= "path to input object",
    #                    default="/vol/scratch/output/modle_test/Stefan/tesmodel1C09.model")


    #args = parser.parse_args()
    #input_file = args.input_file
    #feature_file = '/vol/scratch/data/test/20211213_Cherri_evaluating_RRIs/feature_files/feature_filtered_test_data_context_150_pos.csv'
    what = 'selectft'
    ## DATA ################################################################
    out_dir = '/vol/scratch/20220126_Cherri_model_build/'
    model_stefan = f'{out_dir}model/optimized/model_test.model'
    model_neu = f'{out_dir}model.obj'

    dataset= 'model_test'
    d1 = f'{out_dir}/feature_files/feature_filtered_test_train_context_50_pos_occ_neg.csv'
    d2 = f'{out_dir}/feature_files/feature_filtered_test_train_context_50_pos_occ_pos.csv'

## Structure ###############################################################


    X,y = convert(d1,d2,f'{out_dir}/model/data/{dataset}', graphfeatures=True)

    model_neu, y_pred, y_score = rl.classify2(X,model_neu,'proba')
    accuracy_neu = rl.evaluate(model_neu, X, y)

    #model_stefan, y_pred, y_score = rl.classify2(X,model_stefan,'proba')
    #accuracy_stefan = rl.evaluate(model_neu, X, y)

    print(f'New accuray:{accuracy_neu}')
    #print(f'Stefan accuray:{accuracy_stefan}')






if __name__ == '__main__':
    main()
