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


if __name__ == '__main__':

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


    def median_intensity(mask, image):
        return np.median((image[mask]))

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

    def quantile_intensity(mask, image):
        vals = image[mask]
        return np.quantile(vals, 0.25), np.quantile(vals, 0.75)



    def calc_distribution(image, points, rad_step, max_rad):
        means = np.zeros(
            (len(points), max_rad // rad_step)
        )
        sums = np.zeros(
            (len(points), max_rad // rad_step)
        )
        maxima = np.zeros(
            (len(points), max_rad // rad_step)
        )
        mask = _generate_ring_mask(max_rad, rad_step)

        for idx, center in enumerate(points):
            x, y = center
            current_mask = np.zeros_like(image)
            try:
                current_mask[x-max_rad:x+max_rad+1, y-max_rad:y+max_rad+1] = mask
            except ValueError:
                p= 0
            ringprops = skimage.measure.regionprops_table(current_mask.astype(int), intensity_image=image,
                                                          properties= ('label', 'mean_intensity', 'area', "intensity_max"),  #Fragen an Thorsten: WEeelche Mengenangabe hat Area?
                                                          extra_properties = (intensity_sum,))
            means[idx] = ringprops["mean_intensity"]
            sums[idx] = ringprops["intensity_sum"]
            maxima[idx] = ringprops["intensity_max"]

        return means, sums, maxima


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

            mean_background = skimage.measure.regionprops_table((~(tc.astype(bool))).astype(int), intensity_image=intensity, properties=("mean_intensity", ), extra_properties=(median_intensity, quantile_intensity ))
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





                means, sums, maxima = calc_distribution(intensity_masked, points, rad_step, max_rad)

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

    df_all = pd.concat(dfs)



    # In[36]:





    # In[37]:


    def _color_to_calpha(color: str, alpha: float = 0.2):
        color = color.split("(")[-1].split(")")[0]
        return f"rgba({color}, {alpha})"


    df_all.to_csv(snakemake.output.tsv, sep="\t", index=False)


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


