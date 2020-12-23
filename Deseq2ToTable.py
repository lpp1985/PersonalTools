#!/usr/bin/env python2
# coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2020/12/20
"""

import pandas as pd
from io import BytesIO,StringIO
from itertools import combinations
from lpp import *
import numpy
if __name__ == '__main__':
    usage = '''usage: python2.7 %prog [options] 
'''
    parser = OptionParser(usage=usage)

    parser.add_option("-i", "--INPUT", action="store",
                      dest="INPUT",
                      help="input read count number File ")

    parser.add_option("-o", "--Output", action="store",
                      dest="Output",

                      help="Output Path")

    parser.add_option("-x", "--XLS", action="store",
                      dest="xls",

                      help="Excel File")

    (options, args) = parser.parse_args()
    input_path = options.INPUT
    output = options.Output
    xls = options.xls
    all_need = {}
    for a, b, c in os.walk(input_path):
        for f in c:
            if f == "AllDifferentialGene.tsv":
                data = pd.read_table(a + '/' + f)
                name = a.split("/")[-1].replace("__", " vs ")
                data[name + " Status"] = data["Status"]
                data[name + " padj"] = data["padj"]

                data[name + " padj"] = data["padj"]
                data = data[["Name", name + " padj", name + " Status", "SYMBOL"]]
                all_need[name] = data

    k = all_need.keys()[0]
    combined_data = all_need[k]
    # data_all = {}
    string_need = ""
    f = open(  output + ".venn",'w')
    for k1 in all_need:
        # data_all[k1] =all_need[ k1 ][ "Name" ]
        namelist = list(all_need[k1]["Name"])
        f.write(k1 + "\t" + "\t".join(namelist) + "\n" )

        if k1 != k:
            combined_data = pd.merge(all_need[k1], combined_data ,how='outer', on=["Name"])
    combined_data.to_csv(output + ".see", index=False, sep="\t")
    f.close()
    RAW= open( output + ".see",'rU' )
    has = {}
    END = open(  output + ".xls" ,'w' )
    for line in RAW:
	l_l = line.split("\t")
	if l_l[0] not in has:
		END.write(line)
		has[ l_l[0]  ] = ""

	
    #data = numpy.genfromtxt('cache',delimiter='\t',dtype=str)
    #data = numpy.transpose(  data )
    #df = pd.read_csv( "cache" , header=None,error_bad_lines=False, delimiter="\t")
    #df = df.stack()

    #pd.to_csv( data , output + ".venn", index=False, sep="\t")
