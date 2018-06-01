#!/bin/sh

LIBS=lib/colt.jar:\
lib/jackson-core-asl-1.9.2.jar:\
lib/jackson-mapper-asl-1.9.2.jar:\
lib/jchem.jar:\
lib/jnlp.jar:\
lib/jcommon-1.0.17.jar:\
lib/jfreechart-1.0.14.jar:\
lib/jgoodies-common-1.7.0.jar:\
lib/jgoodies-looks-2.5.3.jar:\
lib/jide-components.jar:\
lib/jide-oss-3.5.9.jar:\
lib/jxlayer.jar:\
lib/prefuse.jar:\
lib/swingworker-0.8.0.jar:\
lib/swingx.jar

# use these as arguments: http://localhost:8080 chembl-v24
java -cp ${LIBS}:rgroup.jar gov.nih.ncgc.rgroup.RGroupTool $*
