#/bin/bash

DIR=$1
for FILE in $DIR/*
do
    FILENAME=$(basename $FILE)
    BASE=$DIR/results/${FILENAME%.tif}
    OUT=$BASE"_ensured.tif"
    OUTSDP=$BASE"_ensured.sdp"
    OUTSVG=$BASE"_curvature.svg"
    # echo $FILE
    ./build/examples/EnsureConnectivity2D -i $FILE -o $OUT -m 1 2>$OUTSDP
    ./build/examples/vcm-curvature -i $OUTSDP -o $OUTSVG -r 4 -R 5 -t 0.25 2>/dev/null
done
