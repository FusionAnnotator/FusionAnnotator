
if [ -z $CTAT_GENOME_LIB ]; then
    echo Error, set env var CTAT_GENOME_LIB to the CTAT genome lib directory before running test.
fi

echo "######################"
echo "## Simple annotations"
echo "#####################"

../FusionAnnotator --genome_lib_dir $CTAT_GENOME_LIB --annotate test.fusions

echo "######################"
echo "## Complex annotations"
echo "#####################"

../FusionAnnotator --genome_lib_dir $CTAT_GENOME_LIB --annotate test.fusions --full


