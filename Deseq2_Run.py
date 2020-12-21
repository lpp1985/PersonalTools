#!/usr/bin/env python2
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2020/12/20
"""

from lpp import *
import pandas as pd


if __name__ == '__main__':
    usage = '''usage: python2.7 %prog [options] 
'''
    parser = OptionParser(usage =usage )



    parser.add_option("-i", "--INPUT", action="store",
                      dest="INPUT",
                      help="input read count number File ")  
    
    parser.add_option("-s", "--Sample", action="store",
                      dest="Sample",
                      help="Sample file ,Sample\tGroup")
    
    parser.add_option("-t", "--Threshold", action="store",
                      dest="Threshold",
                      type='int',
                      default = 1,
                      help="Threshold ,LogFC ,Default 1")    
    
    parser.add_option("-r", "--Run", action="store",
                      dest="Run",
                      help="Run file ,Control\tCase\n")

    parser.add_option("-o", "--Output", action="store",
                      dest="Output",

                      help="Output Path")

    parser.add_option("-x", "--XLS", action="store",
                      dest="XLS",
                      default="",
                      help="xls Annotation file")
    (options, args) = parser.parse_args()
    reads_file = options.INPUT
    sample_file = options.Sample
    threshold = options.Threshold
    output_path = options.Output
    xls = options.XLS
    RUN_FILE = open(  options.Run,'rU' )
    sample_table = pd.read_table( sample_file )
    group_sample = Ddict()
    for _,i in sample_table.iterrows():
        group_sample[   i["Condition"]  ][  i[ "Sample"]  ] = ""     
    
    RUN_FILE.next()
    for line in RUN_FILE:
        line = line.strip()
        l_l = line.split("\t")
        case_group = l_l[1]
        control_group = l_l[ 0 ]
        case_sample_list =  ",".join(  group_sample[  case_group ]    )
        control_sample_list = ",".join(  group_sample[  control_group ]    )
        print(  """deseq2.r  -i {count}  -o {output_path}/Differential/{control_group}__{case_group} --case_group_name {case_group} --control_group_name {control_group}  --case_sample_list {case_sample_list} --control_sample_list {control_sample_list}  --log2fc_cutoff {threshold} -x {xls} """.format(xls=xls, count =reads_file,  output_path = output_path,control_group = control_group, case_group = case_group, case_sample_list= case_sample_list,control_sample_list = control_sample_list,threshold = threshold  )  )
        os.system(  """deseq2.r  -i {count}  -o {output_path}/Differential/{control_group}__{case_group} --case_group_name {case_group} --control_group_name {control_group}  --case_sample_list {case_sample_list} --control_sample_list {control_sample_list}  --log2fc_cutoff {threshold} -x {xls} """.format(xls=xls, count =reads_file,  output_path = output_path,control_group = control_group, case_group = case_group, case_sample_list= case_sample_list,control_sample_list = control_sample_list,threshold = threshold  )  )
        
        
    os.system(
        "Deseq2Cluster.R  -i {count}  -s {sample_file}  -o {output_path}/Cluster".format( 
        count= reads_file,
        sample_file = sample_file,
        output_path = output_path
        ) 
    )
    
