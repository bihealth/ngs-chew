#%%
from glob import glob
from fingerprint import load_fingerprints
from fingerprint import load_fingerprint
import pandas as pd
import vcfpy
import numpy as np
import matplotlib.pyplot as plt
# %%
p = "/home/memsonmi/DEV/ngs-chew/tests/ddata/fingerprints"
paths = glob("/home/memsonmi/DEV/ngs-chew/tests/ddata/fingerprints/*.npz")
vcf_paths = glob("/home/memsonmi/DEV/ngs-chew/tests/ddata/vcfs/*.vcf.gz")

fps = load_fingerprints(paths)
# %%

fps.keys()
# %%

het_aabs = {}
# if args.fps:
#     fps = load_fingerprint(args.fps)
#if args.vcf:
for path_vcf in vcf_paths:
    with vcfpy.Reader.from_path(path_vcf) as vcf_reader:
        for sample in vcf_reader.header.samples.names:
            het_aabs.setdefault(sample, {})
        for record in vcf_reader:
            key = "%s:%s" % (record.CHROM, record.POS)
            for call in record.calls:
                if call.data["GT"].replace("|", "/") in ("0/1", "1/0"):
                    if sum(call.data.get("AD", [0])):
                        het_aabs[call.sample][key] = call.data["AD"][0] / sum(call.data["AD"])
                    else:
                        het_aabs[call.sample][key] = 0
#%%

aabs = {}
for path_vcf in vcf_paths:
    with vcfpy.Reader.from_path(path_vcf) as vcf_reader:
        for sample in vcf_reader.header.samples.names:
            aabs.setdefault(sample, {})
        for record in vcf_reader:
            key = "%s:%s" % (record.CHROM, record.POS)
            for call in record.calls:
                if sum(call.data.get("AD", [0])):
                    aabs[call.sample][key] = call.data["AD"][0] / sum(call.data["AD"])
                else:
                    aabs[call.sample][key] = 0

df = pd.DataFrame.from_dict(aabs)

plot_data = {"sample": [], "bin": [], "count": []}
for c in df.columns:
    hist = np.histogram(df[c], 50, (0.0, 1.0))
    plot_data["count"] += hist[0].tolist()
    plot_data["sample"] += [c] * hist[0].shape[0]
    plot_data["bin"] += hist[1].tolist()[:-1]

plot_df = pd.DataFrame.from_dict(plot_data)
#fig = px.line(plot_df, x="bin", y="count", color="sample")
#fig.update_layout(title=args.title)
#pio.write_html(fig, file=args.out_html)


# %%
fps = load_fingerprints(paths)
plot_data = {"sample": [], "metric": [], "value": []}
for sample in fps.keys():
    is_alt = fps[sample][0][1]
    is_alt_hom = fps[sample][0][2]
    het_vars = fps[sample][1][is_alt & ~is_alt_hom]
    plot_data["sample"].append(sample)
    plot_data["metric"].append("var(het)")
    plot_data["value"].append(np.var(het_vars))
# %%
fps = load_fingerprints(paths)
D = pd.DataFrame.from_dict({key: fps[key][1] for key in fps.keys()})
D.drop(D[D.sum(axis=1) == 0].index,inplace=True)
#%%
for i,k in enumerate(fps.keys()):
    fps[k][1][fps[k][0][1]]
# %%
