#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys # manage arguments
import getopt # manage arguments
import pandas as pd # read csv and manipulate dataframes
import numpy as np # basic operations

def usage():
    print(f"Usage: python {sys.argv[0]} -i input_file (csv) -o out_dir -t thr")

#if __name__ == "__main__":
def compute_strains_abundance_main():

    # check arguments

    input_file = None
    out_dir = None
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
            thr = float(a)
        
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
    strains_profile = pd.DataFrame(index=list((input_df.iloc[:,0:(len(input_df.columns)-12)]).columns) , columns=["detected_genes", "mean_abund", "mean_abund_nz", "median_abund", "median_abund_nz"])
    strains_profile.fillna(0,inplace=True)
    # detected genes
    strains_profile["detected_genes"] = [ input_df.loc[input_df["ratio_covered_nodes"] > 0,columnName].sum()/input_df.loc[:,columnName].sum() for (columnName, columnData) in (input_df.iloc[:,0:(len(input_df.columns)-12)]).iteritems() ]
    # abundances
    strains_profile["mean_abund"] = [ np.mean(input_df.loc[input_df[columnName] > 0,"mean_abund_multiple"]/input_df.loc[input_df[columnName] > 0,columnName]) if strains_profile.loc[columnName,"detected_genes"] > thr else 0 for (columnName, columnData) in (input_df.iloc[:,0:(len(input_df.columns)-12)]).iteritems() ]
    if strains_profile["mean_abund"].sum() != 0: strains_profile["mean_abund"] /= strains_profile["mean_abund"].sum()/100
    strains_profile["mean_abund_nz"] = [ np.mean(input_df.loc[input_df[columnName] > 0,"mean_abund_multiple_nz"]/input_df.loc[input_df[columnName] > 0,columnName]) if strains_profile.loc[columnName,"detected_genes"] > thr else 0 for (columnName, columnData) in (input_df.iloc[:,0:(len(input_df.columns)-12)]).iteritems() ]
    if strains_profile["mean_abund_nz"].sum() != 0: strains_profile["mean_abund_nz"] /= strains_profile["mean_abund_nz"].sum()/100
    strains_profile["median_abund"] = [ np.median(input_df.loc[input_df[columnName] > 0,"mean_abund_multiple"]/input_df.loc[input_df[columnName] > 0,columnName]) if strains_profile.loc[columnName,"detected_genes"] > thr else 0 for (columnName, columnData) in (input_df.iloc[:,0:(len(input_df.columns)-12)]).iteritems() ]
    if strains_profile["median_abund"].sum() != 0: strains_profile["median_abund"] /= strains_profile["median_abund"].sum()/100
    strains_profile["median_abund_nz"] = [ np.median(input_df.loc[input_df[columnName] > 0,"mean_abund_multiple_nz"]/input_df.loc[input_df[columnName] > 0,columnName]) if strains_profile.loc[columnName,"detected_genes"] > thr else 0 for (columnName, columnData) in (input_df.iloc[:,0:(len(input_df.columns)-12)]).iteritems() ]
    if strains_profile["median_abund_nz"].sum() != 0: strains_profile["median_abund_nz"] /= strains_profile["median_abund_nz"].sum()/100

    # output
    strains_profile.to_csv(f"{out_dir}/strains_profile.csv")