# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import re

def chr_mapping(file, chr_column=0, ncnumber_column=1):
    """
        Creates a dictionary which can be used to convert NC number to a simple
        identifier, the simplier version is commonly used in for example VCF
        files.

        Parameters:
        file: (string) path to a tab seperate file containing mapping
              information, ex
              #Chr name    NC ID           Length
              chr1         NC_000001.9     Chr1#NC_000001.9#1#247249719#-1 24724971
        chr_column: (int) index (zero based) of the column containing chr
              information
        ncnumber_column: (int) index (zero based) of the column containing nc
              number information

        Returns:
        dict: ex {'NC_000001.9': '1', 'NC_000024.8': 'Y'}
    """
    mapping_dictonary = dict()
    with open(file, 'r') as conversion_file:
        for line in conversion_file.readlines():
            if not line.startswith("#"):
                columns = line.split("\t")
                mapping_dictonary[columns[ncnumber_column]] = \
                    re.search("chr(\w+)", columns[chr_column]).group(1) #Extract only the identifier
    return mapping_dictonary
