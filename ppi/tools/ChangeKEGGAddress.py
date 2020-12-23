#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2020/12/20
"""
import re,os,sys

for a,b,c in os.walk(sys.argv[1]):
	for f in c:
		if ".htm" in f:
			data = open(a+'/'+f).read()
			data = re.sub("\/tmp\/[^/]+\/","../png/",data)
			data = re.sub("\"/","\"http://www.genome.jp/",data)
			END = open(a+'/'+f,'w')
			END.write(data)
