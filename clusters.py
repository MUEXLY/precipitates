import re
import warnings

import ovito
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import linregress


# input stream of files, can be a single file or a pattern defining a sequence of files
INPUT_STREAM = "example/conf_*.txt"

# rmsd cutoff for the polyhedral template matching, same as in OVITO GUI
RMSD_CUTOFF = 0.15

# structures enabled for polyhedral template matching, same as in OVITO GUI
STRUCTURES_ENABLED = {
    "FCC": False,
    "HCP": False,
    "BCC": True,
    "ICO": False,
    "SC": False,
    "CUBIC_DIAMOND": False,
    "HEX_DIAMOND": False,
    "GRAPHENE": False
}

# expression to select atoms to delete after polyhedral template matching, same as OVITO GUI
EXPRESSION = "StructureType == 0"

# header pattern grabbing time and dose
HEADER_PATTERN = re.compile(r"t=\s*([\d.]+).*dose\s*:\s*([\d.]+)")

# cluster analysis parameters
CLUSTER_CUTOFF = 10.0
MIN_PRECIPITATE_SIZE = 10

# plot settings
NUM_ERROR_BARS = 10
OUTPUT_PLOT_FILE = "clusters.pdf"


def main():

    # load in files from input stream
    pipeline = ovito.io.import_file(INPUT_STREAM)

    # perform polyhedral template matching detecting patterns we specified above
    structure_modifier = ovito.modifiers.PolyhedralTemplateMatchingModifier(rmsd_cutoff=RMSD_CUTOFF)
    for structure_type, enabled in STRUCTURES_ENABLED.items():
        attr = getattr(ovito.modifiers.PolyhedralTemplateMatchingModifier.Type, structure_type)
        structure_modifier.structures[attr].enabled = enabled
    pipeline.modifiers.append(structure_modifier)

    # delete atoms matching expression defined above
    pipeline.modifiers.append(ovito.modifiers.ExpressionSelectionModifier(expression=EXPRESSION))
    pipeline.modifiers.append(ovito.modifiers.DeleteSelectedModifier())

    # perform cluster analysis
    pipeline.modifiers.append(ovito.modifiers.ClusterAnalysisModifier(cutoff=CLUSTER_CUTOFF, sort_by_size=True, compute_com=True))

    # initialize array of times
    num_frames = pipeline.source.num_frames
    times = np.zeros(num_frames)

    # initialize time series arrays
    doses = np.zeros_like(times)
    mean_cluster_sizes = np.zeros_like(times)
    std_cluster_sizes = np.zeros_like(times)
    num_clusters = np.zeros_like(times)

    # loop through frames in pipeline
    for frame in range(num_frames):

        # try to get data from frame
        # if OVITO throws a RuntimeError there was some issue with the file's formatting
        # set to NaN to easily get rid of values at those frames later
        try:
            data = pipeline.compute(frame)
        except RuntimeError as e:

            warnings.warn(f"Error raised for frame {frame} in source file {data.attributes['SourceFile']} with message:\n{e}\nSetting values to nan")

            mean_cluster_sizes[frame] = np.nan
            std_cluster_sizes[frame] = np.nan
            num_clusters[frame] = np.nan
            times[frame] = np.nan
            doses[frame] = np.nan

            continue

        # grab the cluster table and only consider clusters above a certain size
        cluster_table = data.tables["clusters"]
        sizes = cluster_table["Cluster Size"][...]
        sizes = sizes[sizes >= MIN_PRECIPITATE_SIZE]

        # populate time series from size table
        mean_cluster_sizes[frame] = np.mean(sizes)
        std_cluster_sizes[frame] = np.std(sizes)
        num_clusters[frame] = len(sizes)
        
        # get time and dose from file
        with open(data.attributes["SourceFile"], "r") as source_file:
            for line in source_file:
                if not (match := HEADER_PATTERN.search(line)):
                    continue
                t, dose = match.groups()
                t, dose = float(t), float(dose)
                times[frame] = t
                doses[frame] = dose
                break

        print(f"{(frame + 1) / num_frames:.2%} complete", end="\r")

    # initialize figure
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)

    # convert seconds to days
    times *= 1.157e-5

    # plot non-nan values for time series
    axs[0].errorbar(
        times[~np.isnan(times)],
        mean_cluster_sizes[~np.isnan(mean_cluster_sizes)],
        yerr=std_cluster_sizes[~np.isnan(std_cluster_sizes)],
        color="steelblue",
        zorder=6,
        capsize=5,
        errorevery=num_frames // NUM_ERROR_BARS
    )
    axs[0].set_ylabel(r"mean cluster size $\pm 1\sigma$")
    axs[0].grid()

    axs[1].plot(times[~np.isnan(times)], num_clusters[~np.isnan(num_clusters)], zorder=6, color="orangered")
    axs[1].set_xlabel("time (days)")
    axs[1].set_ylabel("number of clusters")
    axs[1].grid()

    time_dose_regression = linregress(times[~np.isnan(times)], doses[~np.isnan(doses)])
    top_ax = axs[0].secondary_xaxis("top", functions=(
        lambda t: time_dose_regression.intercept + time_dose_regression.slope * t,
        lambda dose: (dose - time_dose_regression.intercept) / time_dose_regression.slope
    ))
    top_ax.set_xlabel("dose (dpa)")

    fig.tight_layout()
    fig.savefig(OUTPUT_PLOT_FILE)


if __name__ == "__main__":

    mpl.use("Agg")
    main()
