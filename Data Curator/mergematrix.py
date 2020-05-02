import pandas
from vpolo.alevin import parser

run1_df = parser.read_quants_bin("./SRR387604_salmon")
run2_df = parser.read_quants_bin("./SRR387605_salmon")
run3_df = parser.read_quants_bin("./SRR387606_salmon")

run1_df = run1_df.T
run2_df = run2_df.T
run3_df = run3_df.T

umi_count = pandas.concat([run1_df, run2_df, run3_df], axis=1)

umi_count.to_csv("umi_matrix.csv")

