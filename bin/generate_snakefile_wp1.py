# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import pandas as pd

samples = pd.read_table("samples.tsv", index_col="sample")

result_files = {'common': {"fastqc/{0}.html","fastqc/{0}.zip", "reports/{0}.filteredMutations.tsv"},
    'haloplex': {
        "colon": {"reports/extracted_sample_info/{0}.msiMarkers.tsv"},
        "lung": {"reports/EGFR_T790M.allSamples.ampliconmapped.txt", "reports/EGFR_G719.allSamples.ampliconmapped.txt"},
        "gist": None,
        "melanom": None},
    'swift': {
        "lung": {"reports/EGFR_T790M.allSamples.ampliconmapped.txt", "reports/EGFR_G719.allSamples.ampliconmapped.txt"},
        "ovarial": None
    }
}

expected_result_files = []

for row in samples.itertuples():
    for expected_file in result_files['common']:
        expected_result_files.append(expected_file.format(row.Index,))
    if row['panel_type'][row.Index] in result_files:
        if result_files[row['panel_type']] is not None:
            for expected_file in result_files[row['panel_type']]:
                expected_result_files.append((expected_file.format((row.Index,))))

expected_result_files = set(expected_result_files)

with open("Snakefile", 'w') as run_file:
    run_file.write("import pandas as pd")
    run_file.write("\n\nconfigfile: \"config.yaml\"")
    run_file.write("\nsamples = pd.read_table(config[\"samples\"], index_col=\"sample\")")
    run_file.write("\nunits = pd.read_table(config[\"units\"], index_col=[\"sample\", \"unit\"], dtype=str)")
    run_file.write("\nunits.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index")
    run_file.write("\n\nrule all:")
    run_file.write("\n\tinput:")
    run_file.write("\n\t\t\"" + "\", \"".join(expected_result_files) + "\"")
    run_file.write("\n\ninclude: \"workflows/wp1.snakemake\"")
