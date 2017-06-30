
if [ -z $CTAT_GENOME_LIB ]; then
    echo Error, set env var CTAT_GENOME_LIB to the CTAT genome lib directory before running test.
fi

../FusionAnnotator --genome_lib_dir $CTAT_GENOME_LIB --annotate test.fusions
