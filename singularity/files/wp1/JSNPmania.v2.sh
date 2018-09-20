#!/bin/bash
#
# Script for running jSNPmania
#

# Make script crontab friendly:
JSNPMANIA_LIB_PATH=$(dirname $0);

java -Xmx48G -classpath $JSNPMANIA_LIB_PATH/config:$JSNPMANIA_LIB_PATH/lib/jSNPmania.jar org.jsnpmania.JSNPmania $@

