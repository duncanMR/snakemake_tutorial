def get_ENA_from_ID ( ID ):
    return samples.loc[samples["Sample"] == ID, "ENA"].item()

def get_input_filename_from_ID ( ID, read ):
    return config["fastq_filename_template"].format(
        ENA = get_ENA_from_ID( ID ),
        read = str(read)
    )

def get_read_group_line( ID ):
    ena = get_ENA_from_ID( ID )
    return "@RG'\\'tID:{ID}'\\'tSM:{ena}'\\'tPL:ILLUMINA".format(
        ID = ID,
        ena = ena
    )
