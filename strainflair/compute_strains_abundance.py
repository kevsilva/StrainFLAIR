#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys # manage arguments
import getopt # manage arguments
import pandas as pd # read csv and manipulate dataframes

def usage():
    print(f"Usage: python {sys.argv[0]} -i input_file (csv) -o out_dir -t thr")

#if __name__ == "__main__":
def compute_strains_abundance_main():

    # check arguments

    input_file = "/home/kdasilva/master_project/project_mock_v3/counting_tables/table_mock1a_complete_allgraphs_20201119.csv" #None 
    out_dir = "/home/kdasilva/strain_profiling_vg/" #None 
    thr = 0.5
    
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hi:o:t:")
    
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
        elif o in ("-i"):
            input_file = a
        elif o in ("-o"):
            out_dir = a
        elif o in ("-t"):
            thr = a
        
        else:
            assert False, "unhandled option"
    if not input_file or not out_dir: 
        usage()
        exit()

    # start

    # read csv
    input_df = pd.read_csv(input_file, sep=";")
    # replace NaN with zeros
    input_df.fillna(0,inplace=True)
    # working only on specific genes
    input_df = input_df.drop( input_df[(input_df.iloc[:,0:(len(input_df.columns)-12)] > 0).sum(axis=1) != 1].index )

    # strain-level computation

    # empty dataframe
    strain_profile = pd.DataFrame(index=list((input_df.iloc[:,0:(len(input_df.columns)-12)]).columns) , columns=["detected_genes", "mean_abund", "mean_abund_nz", "median_abund", "median_abund_nz"])
    strain_profile.fillna(0,inplace=True)
    # detected genes
    strain_profile["detected_genes"] = [ input_df.loc[input_df["ratio_covered_nodes"] > 0,columnName].sum()/input_df.loc[:,columnName].sum() for (columnName, columnData) in (input_df.iloc[:,0:(len(input_df.columns)-12)]).iteritems() ]