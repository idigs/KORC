#!/bin/bash

MACHINE_ID=`uname -n`

if [ $MACHINE_ID == "MBP115573"     ]
then

   git -C /Users/21b/Desktop/KORC rev-parse HEAD > git_hash.txt

   git -C /Users/21b/Desktop/KORC diff HEAD > git_diff.txt
   
elif [ $MACHINE_ID == "cori01"   ] || \
	 [ $MACHINE_ID == "cori02"   ] || \
	 [ $MACHINE_ID == "cori03"   ] || \
	 [ $MACHINE_ID == "cori04"   ] || \
	 [ $MACHINE_ID == "cori05"   ] || \
	 [ $MACHINE_ID == "cori06"   ] || \
	 [ $MACHINE_ID == "cori07"   ] || \
	 [ $MACHINE_ID == "cori08"   ] || \
	 [ $MACHINE_ID == "cori09"   ] || \
	 [ $MACHINE_ID == "cori10"   ] || \
	 [ $MACHINE_ID == "cori11"   ] || \
	 [ $MACHINE_ID == "cori12"   ]
then

   git -C /global/u1/m/mbeidler/KORC rev-parse HEAD > git_hash.txt

   git -C /global/u1/m/mbeidler/KORC diff HEAD > git_diff.txt
    
else
# MACHINE_ID is new and unknown. Inform the user how to add support for this new machine.
    echo $MACHINE_ID not suported by this script.
    echo To support this machine, add a new elif statement of the form
    echo
    echo elif [ \$MACHINE_ID == \"$MACHINE_ID\" ]
fi

exit $?
