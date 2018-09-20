#!/bin/bash
#
# Script for running local GATK version
#

# Make script crontab friendly:
PATH=$(dirname $0)

/usr/bin/java -jar $PATH/GenomeAnalysisTKLite_molecules.jar $@
