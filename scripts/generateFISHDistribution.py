#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import os

#import pyasn1_modules.rfc2634
import tifffile
import plotly.graph_objs as go
import plotly.express as px
import copy
import multiprocessing as mp
from jinja2.compiler import find_undeclared


# In[9]:


#viewer = napari.Viewer()
#napari.run()


# In[2]:


import plotly.io as pio
pio.renderers.default = "notebook"


# In[4]:


def yield_image(images):
    for figure, probe_dict in images.items():
        for probe, treat_dict in probe_dict.items():
            for treatment, image_dict in treat_dict.items():
                for image_dict in images[figure][probe][treatment]:
                    yield figure, probe, treatment, image_dict


# In[5]:


import skimage

from plotly.subplots import make_subplots


# In[6]:


def show_images(images, titles=None, color_continuous_scale='gray', max_height=300, spacing=0.01, norm=False):
    """
    Displays multiple NumPy arrays (grayscale or RGB) using px.imshow inside subplots.

    Parameters:
    - images: list of 2D or 3D np.ndarray (grayscale or RGB)
    - titles: optional list of titles
    - color_continuous_scale: for grayscale images (e.g., 'gray', 'viridis', etc.)
    - max_height: height of figure in pixels
    - spacing: horizontal spacing between images
    """
    from plotly.subplots import make_subplots

    n = len(images)
    if titles is None:
        titles = [f"Image {i+1}" for i in range(n)]

    fig = make_subplots(rows=1, cols=n, subplot_titles=titles,
                        horizontal_spacing=spacing)

    for i, (key, img) in enumerate(images.items()):
        # Generate a single-image figure using px.imshow
        img = tifffile.imread(img)
        if img.ndim == 3:
            img = img[0, :]
        if norm:
            img = img / img.max()
        temp_fig = px.imshow(img,
                             color_continuous_scale=color_continuous_scale if img.ndim == 2 else None,
                             aspect='equal')
        # Extract the trace and layout from the temp fig
        for trace in temp_fig.data:
            fig.add_trace(trace, row=1, col=i+1)

        # Hide colorbar for all but the last subplot (optional)
        if i < n - 1:
            fig.update_traces(colorbar=dict(len=0), selector=dict(type='heatmap'))

    fig.update_layout(
        height=max_height,
        width=max_height * n,
        showlegend=False,
        margin=dict(t=40, b=10, l=10, r=10),
        title_text="Image Gallery"
    )

    # Turn off axis ticks
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)

    return fig


# In[7]:


def match_to_to_total_cell(tc_image, subimage):

    tc_labels = {l: None for l in np.unique(tc_image)}
    for tc_label in tc_labels:
        mask = (tc_image == tc_label)
        ulabels = np.unique(subimage[mask])
        ulabels = np.delete(ulabels, 0)
        tc_labels[tc_label] = ulabels
    del tc_labels[0]
    return tc_labels

def refactor_submask_indices(tc_image, submask, tc_labels):
    new_mask = np.zeros_like(submask)
    for key, val_list in tc_labels.items():
        for val in val_list:
            tc_mask = tc_image == key
            m = (submask == val) & tc_mask
            new_mask[m] = key
    return new_mask


# In[9]:


#import os
#os.chdir("../")
# os.getcwd()


# In[10]:


probes =  ["psaA", "cupA", "rnpB", "psbA2", "rpl10_12",  "16S rRNA", "fusB"]
treatments = ["untreated", "Kasugamycin", "Puromycin"]
images = {"Replicate 1": {}, "Replicate 2": {}, "Replicate 3": {}}
for probe in probes:
    subdirs = f"initialData/AllFISHData/{probe}"
    replicates = [r for r in os.listdir(subdirs) if "Replicate" in r and "old" not in r]
    for replicate in replicates:
        images[replicate][probe] = {}
        for treatment in treatments:
            images[replicate][probe][treatment] = {}
            path = os.path.join(subdirs, f"{replicate}/Segmentation_{treatment}")

            if os.path.exists(path):
                files = sorted([os.path.join(path, file) for file in os.listdir(path)])
                for filename in files:
                    if filename.endswith(".tiff") or filename.endswith(".tif"):
                        im_name = os.path.basename(filename).split("_")[1]

                        if "{Membrane}" in filename:
                            kind = "Membrane"
                        elif "_corrected_total_cell" in filename:
                            kind = "TotalCell"
                        elif "_label_image_corrected" in filename or "_segmentation_corrected" in filename:
                            kind = "Probe"
                        elif "{Cytosol}" in filename:
                            kind = "Cytosol"
                        elif "{FakeCytosol}" in filename:
                            continue
                        elif "{FixedMembrane}" in filename:
                            continue
                        elif "{FixedCytosol}" in filename:
                            continue
                        elif "Normalized" in filename:
                            continue
                        elif "Standardized" in filename:
                            continue
                        else:
                            kind = "Intensity"
                        if im_name not in images[replicate][probe][treatment]:
                            images[replicate][probe][treatment][im_name] = {}

                        images[replicate][probe][treatment][im_name][kind] = filename



# In[8]:


to_test = ["Intensity", "TotalCell", "Probe"]
to_del = []
for replicate, probe, treatment, image in yield_image(images):
    if not all([i in images[replicate][probe][treatment][image] for i in to_test]):
        l = (replicate, probe, treatment, image)
        if l not in to_del:
            to_del.append(l)
        else:
            b = 0
            h = 0
            print(l)
#to_del = set(to_del)
for replicate, probe, treatment, image in to_del:
    del images[replicate][probe][treatment][image]


# In[19]:


imgs = images["Replicate 1"]["psaA"]["untreated"]["Image 1"]
show_images(imgs, titles=list(imgs.keys()), norm=True)


# In[20]:


px.imshow(tifffile.imread(imgs["Intensity"])[0])



# In[21]:


px.imshow(tifffile.imread(imgs["Intensity"])[0])





# In[19]:


from scripts.foo import median_intensity

def _find_spots_center(images):
    dfs = []
    for figure, probe, treatment, image in yield_image(images):
        labeled_spots = tifffile.imread(images[figure][probe][treatment][image]["Probe"])
        tc_image = tifffile.imread(images[figure][probe][treatment][image]["TotalCell"])
        probe_probs = skimage.measure.regionprops_table(labeled_spots, intensity_image=tc_image, properties=("label", "centroid", "area"), extra_properties= ( median_intensity, ))
        intensity_image = tifffile.imread(images[figure][probe][treatment][image]["Intensity"])[1]
        probe_intensity = skimage.measure.regionprops_table(labeled_spots, intensity_image=intensity_image, properties=("label",), extra_properties= ( median_intensity, ))
        bg_intensity = skimage.measure.regionprops_table((labeled_spots == 0).astype(int), intensity_image=intensity_image, properties=("label",), extra_properties= ( median_intensity, ))

        bg_intensity = bg_intensity["median_intensity"]
        probe_intensity = pd.DataFrame(probe_intensity)
        probe_intensity = probe_intensity.rename({"median_intensity": "ProbeIntensity"}, axis=1)


        df = pd.DataFrame(probe_probs)
        df = df.rename({"median_intensity": "CellIdx", "centroid-0": "x", "centroid-1": "y"}, axis=1)
        df = df.merge(probe_intensity, on="label", how="inner")
        df["background"] = float(bg_intensity)
        df["CellIdx"] = df["CellIdx"].astype(int)
        if probe == "16S rRNA" and figure == "Replicate 1":
            print("lala")
            p = 0
        df["Image"] = image
        df["Probe"] = probe
        df["Treatment"] = treatment
        df["Figure"] = figure
        df = df[df["CellIdx"] != 0]
        dfs.append(df)
    df = pd.concat(dfs)
    return df

spots_df = _find_spots_center(images)
#spots_df = spots_df[(spots_df["area"] >= 20) & (spots_df["area"] <= 100)]


# In[10]:


spots_df.to_csv(snakemake.output.spots_df, sep="\t", index=False)


# In[24]:


import numpy as np


imgs = []
measured = []
probes = []
treatments = []
figures = []
imgs_to_norm = []
normalized = []
normalized_images = []


for figure, probe, treatment, image in yield_image(images):
    intensities = images[figure][probe][treatment][image]["Intensity"]
    tc_file = images[figure][probe][treatment][image]["TotalCell"]
    tc = tifffile.imread(tc_file)
    tc_mask = tc.astype(bool)
    intensities = tifffile.imread(intensities)
    if intensities.ndim > 2:
        intensities = intensities[0]
    mean = np.nanmean(intensities * tc_mask)
    intensities = intensities.astype(float) * tc_mask
    intensities[~tc_mask ] = np.nan
    imgs.append(image)
    probes.append(probe)
    treatments.append(treatment)
    measured.append(mean)
    figures.append(figure)
    normalized.append(False)
    nim = np.zeros_like(intensities) * np.nan
    nim2 = np.zeros_like(intensities) * np.nan

    u_labels = np.delete(np.unique(tc), 0)
    for u_label in u_labels:
        mask = tc == u_label
        vals = intensities[mask]
        nim[mask] = (vals - vals.min()) / (vals.max() - vals.min() )
        nim2[mask] = (vals - vals.mean()) / vals.std()

    full_normalized = np.zeros_like(intensities)
    vals = intensities[tc_mask]
    full_normalized[tc_mask] = (vals - vals.min()) / (vals.max() - vals.min() )
    dir = os.path.dirname(tc_file)
    name = f"{image}_IntensityNormalized.tif"
    name2 = f"{image}_IntensityStandardized.tif"
    filename = os.path.join(dir, name)
    filename2 = os.path.join(dir, name2)
    images[figure][probe][treatment][image]["NormalizedIntensity"] = filename
    images[figure][probe][treatment][image]["StandardizedIntensity"] = filename2
    tifffile.imwrite(filename, full_normalized)
    tifffile.imwrite(filename2, nim2)
    mean = np.nanmean(nim * tc_mask)
    imgs.append(image)
    probes.append(probe)
    treatments.append(treatment)
    measured.append(mean)
    figures.append(figure)
    normalized.append(True)






intensity_df = {
    "Image": imgs,
    "Probe": probes,
    "Treatment": treatments,
    "MeanIntensity": measured,
    "Figure": figures,
    "Normalized": normalized
}
intensity_df = pd.DataFrame(intensity_df)





# In[25]:


px.imshow(nim)


# In[26]:


px.imshow(full_normalized)


# In[27]:


intensity_df


# In[28]:


intensity_df["group"] = intensity_df["Probe"] + intensity_df["Figure"] + intensity_df["Treatment"]
pdf = intensity_df[intensity_df["Normalized"] == False]
px.box(pdf, x="group", y="MeanIntensity", color="Normalized", hover_data="Image")


# In[29]:


# = spots_df[spots_df["Treatment"] == "untreated"]
px.box(spots_df, x="Probe", y="area", color="Treatment", hover_data="Image")


# In[35]:


def _calc_image_data_only_probe(images, figure, probe, treatment, image, sub_spots_df, rad_step, max_rad, border_distance,):
    sub_dfs = []
    colnames = []
    n_steps = max_rad // rad_step

    for type in ("mean", "sum", "maxima"):
        colnames += [f"ring_{i}_{type}" for i in range(n_steps)]
    parts = ("Measured", )

    tc = tifffile.imread(images[figure][probe][treatment][image]["TotalCell"])

    intensity = tifffile.imread(images[figure][probe][treatment][image]["Intensity"])
    if intensity.ndim > 2:
        intensity = intensity[0, :]
    tc_labels = np.delete(np.unique(tc), 0)
    for cell_label in tc_labels:
        cell_df_mask = sub_spots_df["CellIdx"] == cell_label
        current_df = sub_spots_df[cell_df_mask]

        n_spots_per_cell = len(current_df)
        tc_masked = tc == cell_label

        intensity_masked = intensity * tc_masked

        mean_background = skimage.measure.regionprops_table((~(tc.astype(bool))).astype(int), intensity_image=intensity, properties=("mean_intensity", ), extra_properties=(median_intensity, foo.quantile_intensity ))
        bg_mean = mean_background["mean_intensity"][0]


        for part in parts:
            if part == "TotalCell":
                labels, points = generate_random_spots(tc_masked.astype(int), n_points=n_spots_per_cell, min_distance=border_distance)
            elif part == "Measured":
                points = current_df[["x", "y"]].values.round().astype(int)
                height, width = tc_masked.shape
                points = [
                    point for point in points
                    if (
                            border_distance <= point[0] < (height - border_distance) and
                            border_distance <= point[1] < (width - border_distance) and
                            tc_masked[point[0], point[1]]
                    )
                ]

            else:
                raise ValueError("Wrong cell part")





            means, sums, maxima = foo.calc_distribution(intensity_masked, points, rad_step, max_rad)

            values = np.concatenate(( means, sums, maxima), axis=1)
            df = pd.DataFrame(
                values,
                columns=colnames
            )
            df["Mask"] = part
            df["Image"] = image
            df["Probe"] = probe
            df["Treatment"] = treatment
            df["Cell"] = image + "_" + str(cell_label)
            df["Figure"] = figure
            df["bg_mean"] = bg_mean
            sub_dfs.append(df)
    return sub_dfs



# In[36]:


import foo

max_rad = 30
border_distance = max_rad
rad_step = 1
n_steps = max_rad // rad_step
dfs = []
ct = 0



calls = []
threads = snakemake.threads if "snakemake" in globals() else 1
for figure, probe, treatment, image in yield_image(images):
    image_mask = spots_df["Image"] == image
    treat_mask = spots_df["Treatment"] == treatment
    probe_mask = spots_df["Probe"] == probe
    sub_spots_df = spots_df[image_mask & treat_mask & probe_mask]
    calls.append((images, figure, probe, treatment, image, sub_spots_df, rad_step, max_rad, border_distance ))

if threads == 1:
    dfs = [_calc_image_data_only_probe(*call) for call in calls]
else:
    with mp.Pool(threads) as pool:
        dfs = pool.starmap(_calc_image_data_only_probe, calls)

dfs = [j for sub in dfs for j in sub]
print(dfs)
df_all = pd.concat(dfs)
print(ct)


# In[36]:





# In[37]:


def _color_to_calpha(color: str, alpha: float = 0.2):
    color = color.split("(")[-1].split(")")[0]
    return f"rgba({color}, {alpha})"


# In[38]:


# p = 0
# lines = [None, "dash", "dot"]
# probes = df_all["Probe"].unique().tolist()
# fig = make_subplots(cols=len(probes), column_titles=probes)
# pio.renderers.default = "browser"
# line_map = {l: lines[i] for i, l in enumerate(df_all["Treatment"].unique())}
# df_i = df_all
# for idx, image in enumerate(probes):
#     for treatment in df_i["Treatment"].unique():
#         sdf = df_i[df_i["Probe"] == image]
#         sdf = sdf[sdf["Treatment"] == treatment]
#
#         argmax_df = sdf["thresh_idx"]
#         pdf = sdf.loc[:, ["Probe", "Mask", "Cell"]]
#         pdf.loc[:, "argmax"] = argmax_df
#         pdf = pdf[~pdf["argmax"].isna()]
#         # Compute 1 - CDF
#         pdf = pdf.sort_values("argmax")
#         counts = (
#             pdf.groupby("Mask")["argmax"]
#             .apply(lambda x: pd.cut(x, bins=range(max_rad)).value_counts(sort=False))
#             .rename_axis(index=["Mask", "bin"])
#             .reset_index(name="count")
#         )
#         counts["bin"] = counts["bin"].apply(lambda cell: int(cell.left)).astype(int)
#         counts["cdf"] = counts.groupby("Mask")["count"].cumsum() / counts.groupby("Mask")["count"].transform("sum")
#
#         # Compute 1 - CDF
#         counts["1_minus_cdf"] = 1 - counts["cdf"]
#
#         for mask, color in [("Measured", "rgb(125, 74, 154)"), ("Cytosol", "rgb(52, 74, 154)"), ("Membrane", "rgb( 0, 160, 130)"), ]:
#             x = counts.reset_index()
#
#             sdf = counts[counts["Mask"] == mask]
#
#             fig.add_trace(
#                 go.Scatter(
#                     x=[0] + (sdf["bin"] + 1).tolist(),
#                     y=[1] + (sdf["1_minus_cdf"]).tolist(),
#                     marker=dict(color=color),
#                     line=dict(color=color, dash=line_map[treatment]),
#                     name=image + mask + treatment,
#                     legendgroup="Survival function", legendgrouptitle=dict(text="Survival function")
#
#                 ), col=idx+1, row=1
#
#
#             )
#             fig.add_trace(
#                 go.Histogram(
#                     x=pdf[pdf["Mask"] == mask]["argmax"], histnorm="probability",
#                     marker=dict(color=_color_to_calpha(color, 0.4)),
#                     name=image + mask,
#                     legendgroup="Probability", legendgrouptitle=dict(text="Probability")
#                 ), col=idx+1, row=1
#             )
#
#
#     fig.update_layout( xaxis_title="argmax", yaxis_title="Survival | Probability",)
# fig.show()


# In[39]:


df_all.to_csv(snakemake.output.tsv, sep="\t", index=False)


# In[ ]:


# def show_spots(imgs):
#     probe = tifffile.imread(imgs["Probe"])
#     intensity = tifffile.imread(imgs["Intensity"])
#     if intensity.ndim > 2:
#         intensity = intensity[0]
#     mem = tifffile.imread(imgs["FixedMembrane"])
#     b_mem = mem.astype(bool)
#     bprobe = probe.astype(bool)
#     m_int = intensity.max()
#     intensity[b_mem] = m_int + 10
#     intensity[bprobe] = m_int * 1.5
#
#
#     return px.imshow(intensity)


# 

# In[ ]:


# figure = "Figure5"
# probe = "rpl10_12"
# n = len(images[figure][probe]["untreated"])
# fig = make_subplots(cols=n)
# for i, img in enumerate(images[figure][probe]["untreated"]):
#     imgs = images[figure][probe]["untreated"][img]
#     sub_fig = show_spots(imgs)
#     fig.add_traces(sub_fig.data, cols=i+1, rows=1)
#
# fig.update_layout(width=n*1000, height=1000)
# fig.show()


# In[ ]:





# In[40]:


from scipy import stats

# id_vars = ["Mask", "Image", "Probe", "Figure", "Treatment"]
# mode = "mean"
# cols = [col for col in df_all.columns if mode in col and "bg_mean" not in col]
# df = df_all
# df = df[cols + id_vars]
# probes = df["Probe"].unique().tolist()
#
# # df[cols] = df[cols].values  / df[cols].values.max(axis=1)[:, None]
# df = df.rename({col: col.split("_")[1] for col in cols}, axis=1)
# df = df.melt(id_vars=id_vars)
# summary = df.groupby(['variable', 'Mask', "Probe", "Figure", "Treatment"])['value'].agg(
#  mean='mean',
#  q25=lambda x: np.nanquantile(x, 0.25),
#  q75=lambda x: np.nanquantile(x, 0.75)
# ).reset_index()
# print(summary)
# row_titles = [f"{probe}-{figure}" for (probe, figure), _ in summary.groupby(["Probe", "Figure"])]
# col_titles = ["untreated", "Puromycin", "Kasugamycin"]
#
# fig = make_subplots(rows=len(row_titles) + 1, row_titles=row_titles + ["All"], cols=3, column_titles=col_titles, vertical_spacing=0.005)
# # pio.renderers.default = "browser"
# # Add traces for each Mask group
# colors = px.colors.qualitative.Plotly + px.colors.qualitative.Alphabet
# print(row_titles)
# cmap = {row_tile: colors[i] for i, row_tile in enumerate(row_titles)}
# for k, treat in enumerate(col_titles):
#  for i, ((probe, figure), group_df) in enumerate(summary.groupby(["Probe", "Figure"])):
#      traces = [
#
#      ]
#      errors = []
#      sdf = summary[summary["Treatment"] == treat]
#
#      #for mask, color in [("Measured", "rgb(125, 74, 154)"), ("Cytosol", "rgb(52, 74, 154)"), ("Membrane", "rgb( 0, 160, 130)"), ("TotalCell", "rgb( 111, 160, 130)")]:
#      for mask, color in [("Measured", "rgb(125, 74, 154)")]:
#
#          group = sdf[sdf["Mask"] == mask]
#          group = group[group["Probe"] == probe]
#          group = group[group["Figure"] == figure]
#          group["variable"] = group.variable.astype(int)
#          group = group.sort_values("variable")
#          if not len(group):
#              continue
#          p = go.Scatter(
#              x=group['variable'],
#              y=group['mean'],
#              mode='lines',
#              name=f'{probe} {mask} Mean',
#              line=dict(width=2, color=color)
#          )
#          traces.append(p)
#          if mask == "Measured":
#              fig.add_trace(
#                  go.Scatter(
#                      x=group['variable'],
#                      y=group['mean'],
#                      mode='lines',
#                      name=f'{probe} {mask} Mean',
#                      line=dict(width=2, color=cmap[row_titles[i]])
#                  )
#                  , row=len(row_titles) + 1, col=k+1
#              )
#
#          print(probe, mask, idx_m)
#
#          idx_m = group["mean"].argmax()
#
#
#
#          # Add shaded area between q25 and q75
#          errors.append(go.Scatter(
#              x=pd.concat([group['variable'], group['variable'][::-1]]),
#              y=pd.concat([group['q75'], group['q25'][::-1]]),
#              fill='toself',
#              fillcolor=_color_to_calpha(color, 0.2),  # You can make this color dynamic per mask
#              line=dict(color='rgba(255,255,255,0)'),
#              hoverinfo="skip",
#              name=f'{mask} IQR'
#          ))
#      fig.add_traces(errors, rows=i+1, cols=k+1)
#      fig.add_traces(traces, rows=i+1, cols=k+1)
# fig.update_layout(height=6000)
# fig.show()


# In[ ]:


# from scipy import stats
# id_vars = ["Mask", "Image", "Probe", "Figure", "Treatment"]
# mode = "mean"
# cols = [col for col in df_all.columns if mode in col and "bg_mean" not in col]
# df = df_all[cols + id_vars]
# probes = df["Probe"].unique().tolist()
#
# # df[cols] = df[cols].values  / df[cols].values.max(axis=1)[:, None]
# df = df.rename({col: col.split("_")[1] for col in cols}, axis=1)
# df = df.melt(id_vars=id_vars)
# summary = df.groupby(['variable', 'Mask', "Probe", "Figure", "Treatment"])['value'].agg(
#     mean='mean',
#     q25=lambda x: np.nanquantile(x, 0.25),
#     q75=lambda x: np.nanquantile(x, 0.75)
# ).reset_index()
#
# fig = make_subplots(rows=1, cols=3)
# for i, treatment in enumerate(["untreated", "Puromycin", "Kasugamycin"]):
#     for _, ((probe, figure), group_df) in enumerate(summary.groupby(["Probe", "Figure"])):
#         traces = [
#
#         ]
#         errors = []
#         row_title = f"{probe}-{figure}"
#
#         group = summary[summary["Mask"] == "Measured"]
#         group = group[group["Treatment"] == treatment]
#         group = group[group["Probe"] == probe]
#         group = group[group["Figure"] == figure]
#         group["variable"] = group.variable.astype(int)
#         group = group.sort_values("variable")
#
#         fig.add_trace(
#             go.Scatter(
#                 x=group['variable'],
#                 y=group['mean'],
#                 mode='lines',
#                 name=f'{probe} {mask} Mean',
#                 line=dict(width=2, color=cmap[row_title])
#             )
#             , row=1, col=i+1
#         )
# fig.show()


# In[ ]:


# fig2 = go.Figure()
# fig2.add_traces([f for f in fig.data if f.xaxis == "x10"])
# fig2.show()


# In[ ]:


# mask = "Membrane"
# group = summary[summary["Mask"] == mask]
# group["variable"] = group.variable.astype(int)
# group = group.sort_values("variable")


# 

# In[ ]:


def _count_membrane_overlap(images):
    counts = []
    total_counts = []
    probez = []
    figurez = []
    treatmentz = []
    imagez = []
    for figure, probes in figures.items():
        for probe in probes:
            probe_mask = spots_df["Probe"] == probe
            for treatment in treatments:
                treat_mask = spots_df["Treatment"] == treatment
                for image in images[figure][probe][treatment]:
                    image_mask = spots_df["Image"] == image

                    membrane = tifffile.imread(images[figure][probe][treatment][image]["FixedMembrane"])
                    probe_image = tifffile.imread(images[figure][probe][treatment][image]["Probe"])
                    cytosol = tifffile.imread(images[figure][probe][treatment][image]["FixedCytosol"])

                    current_df = spots_df[probe_mask & treat_mask & image_mask]
                    spots = set(current_df["label"].unique())
                    f = np.unique(probe_image)
                    filter = set(f) - spots

                    probe_image[:border_distance, -border_distance:] = 0
                    probe_image[-border_distance:, :border_distance] = 0
                    before = set(np.unique(probe_image))
                    before = before - filter
                    total_count = len(before)
                    probe_image[~membrane.astype(bool)] = 0
                    after = set(np.unique(probe_image))
                    after = after - filter
                    count = len(after)
                    counts.append(count)
                    total_counts.append(total_count)
                    probez.append(probe)
                    figurez.append(figure)
                    treatmentz.append(treatment)
                    imagez.append(image)

    df = {
        "Treatment": treatmentz,
        "Probe": probez,
        "Figures": figurez,
        "Image": imagez,
        "counts": counts,
        "total_counts": total_counts,

    }
    return pd.DataFrame(df)

# count_df = _count_membrane_overlap(images)
# count_df["at Membrane [%]"] = count_df["counts"] / count_df["total_counts"] * 100


