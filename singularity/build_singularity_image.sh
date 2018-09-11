#!/usr/bin/bash

GREEN='\033[1;32m';
YELLOW='\033[1;33m';
RED='\033[0;31m';
COLOR_OFF='\033[0m';
STATUS=0;

if [ -f "caps.simg" ]; then
  echo -e $YELLOW "Caps image already exists!";
  echo -e $COLOR_OFF;
  exit 1;
fi

echo -e $GREEN "Checking for dependent files.";

if [ ! -f "../local_dependencies/pindel" ]; then
    echo -e $RED "\t- Missing pindel executable!";
    STATUS=1;
fi

if [ ! -f "../local_dependencies/pindel" ]; then
    echo -e $RED "\t- Missing pindel executable!";
    STATUS=1;
fi

if [ ! -f "../local_dependencies/pindel2vcf" ]; then
    echo -e $RED "\t- Missing pindel2vcf executable!";
    STATUS=1;
fi

if [ ! -f "../local_dependencies/sam2pindel" ]; then
    echo -e $RED "\t- Missing sam2pindel executable!";
    STATUS=1;
fi

if [ ! -d "../local_dependencies/jSNPmania-0.0.7-SNAPSHOT-amplicons_v4" ]; then
    echo -e $RED "\t- Missing jSNPmania-0.0.7-SNAPSHOT-amplicons_v4 folder!";
    STATUS=1;
else
    if [ ! -f "../local_dependencies/jSNPmania-0.0.7-SNAPSHOT-amplicons_v4/JSNPmania.v2.sh" ]; then
        echo -e $RED "\t- Missing JSNPmania.v2.sh executable!";
        STATUS=1;
    fi
    if [ ! -f "../local_dependencies/jSNPmania-0.0.7-SNAPSHOT-amplicons_v4/lib/jSNPmania-0.0.7-SNAPSHOT.jar" ]; then
        echo -e $RED "\t- Missing lib/jSNPmania-0.0.7-SNAPSHOT.jar!";
        STATUS=1;
    fi
fi

if [ ! -d "../local_dependencies/AmpliconAlignment" ]; then
    echo -e $RED "\t- Missing AmpliconAlignment folder!";
    STATUS=1;
else
  if [ ! -f "../local_dependencies/AmpliconAlignment/AmpliconMapping.sh" ]; then
      echo -e $RED "\t- Missing AmpliconMapping.sh executable!";
      STATUS=1;
  fi
  if [ ! -f "../local_dependencies/AmpliconAlignment/GenomeAnalysisTKLite_molecules.jar" ]; then
      echo -e $RED "\t- Missing GenomeAnalysisTKLite_molecules.jar!";
      STATUS=1;
  fi
fi

if [ ! -f "../local_dependencies/table_annovar.pl" ]; then
    echo -e $RED "\t- Missing table_annovar.pl scripts!";
    STATUS=1;
fi

if [ ! -f "../local_dependencies/jdk-7u80-linux-x64.rpm" ]; then
    echo -e $RED "\t- Missing oracle jdk rpm!";
    STATUS=1;
fi

if [ $STATUS -ne 0 ]; then
    echo -e $RED "One ore more local files are missing!!!";
else
  echo -e $GREEN "Building singularity image.";
  singularity build caps.simg caps.def;
fi

echo -e "${COLOR_OFF}";
