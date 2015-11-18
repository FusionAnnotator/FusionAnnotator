all:
	tar xvf RESOURCES/Hg19_CTAT_fusion_annotator_lib.tar.gz
	./FusionAnnotator --fusion_annot_lib Hg19_CTAT_fusion_annotator_lib --build

