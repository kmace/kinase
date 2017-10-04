intermediate/images/normalized_data.RData: $(shell find input/fastq -type f)
	cd src/utils
	Rscript create_normalized_data_object.R
	cd ../../

