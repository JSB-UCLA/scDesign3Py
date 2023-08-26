from typing import Optional, Union

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2.robjects as ro
import seaborn as sns
from matplotlib.patches import Patch

from .._utils._format import _anndata2sce, _other2list, _strvec_none2ri, convert


def plot_reduceddim(
    ref_anndata: ad.AnnData,
    anndata_dict: dict,
    name_list: list,
    color_by: str,
    assay_use: Optional[str] = None,
    default_assay_name: Optional[str] = None,
    n_pc: int = 50,
    pc_umap: bool = True,
    center: bool = True,
    scale: bool = True,
    if_plot: bool = True,
    shape_by: Optional[str] = None,
    point_size: Union[int, float] = 1,
):
    """Dimensionality reduction and visualization

    This function takes a reference sce and a list of new sces, performs the dimensionality reduction on the reference data, projects the synthetic datasets on the same low dimensional space, then visualize the results.

    Arguments:
    ----------
    ref_anndata: `anndata.AnnData`
        The reference `anndata.AnnData`

    anndata_dict: `dict`
        A dict of `anndata.AnnData` which is synthetic

    name_list: `list`
        A list of the names of each dataset. The length should be len(anndata_dict)+1, where the first name is for ref_sce.

    color_by: `str`
        The name in the `anndata.AnnData.obs` used for color.

    assay_use: `str` (default: None)
        Indicates the assay you will use. If None, please specify a name for the assay stored in `anndata.AnnData.X` in @default_assay_name. The @assay_use is both used in @ref_anndata and @anndata_dict.

    default_assay_name: `str` (default: None)
        Specified only when @assay_use is None. Asign a name to your default single cell experiment.

    n_pc: `int` (default: 50)
        The number of PCs.

    pc_umap: `bool` (default: True)
        Whether using PCs as the input of UMAP.

    center: `bool` (default: True)
        Whether centering the data before PCA.

    scale: `bool` (default: True)
        Whether scaling the data before PCA.

    if_plot: `bool` (default: True)
        Whether returning the plot. If False, return the reduced dimensions of each dataset.

    shape_by: `str` (default: None)
        The name in the `anndata.AnnData.obs` used for shape.

    point_size: `int` or `float` (default: 1)
        The point size in the final plot.
    """

    ref_sce, final_assay_use = _anndata2sce(
        data=ref_anndata,
        assay_use=assay_use,
        default_assay_name=default_assay_name,
    )

    sce_list = ro.ListVector(
        {
            key: _anndata2sce(
                data=value,
                assay_use=assay_use,
                default_assay_name=default_assay_name,
            )[0]
            for key, value in anndata_dict.items()
        }
    )

    name_list_, shape_by_, color_by_ = _other2list(name_list, shape_by, color_by)
    name_list_, shape_by_, color_by_ = _strvec_none2ri(name_list_, shape_by_, color_by_)

    with convert.context():
        plot_reduceddim_func = ro.r("scDesign3::plot_reduceddim")
        res = plot_reduceddim_func(
            ref_sce=ref_sce,
            sce_list=sce_list,
            name_vec=name_list_,
            assay_use=final_assay_use,
            n_pc=n_pc,
            pc_umap=pc_umap,
            center=center,
            scale_=scale,
            if_plot=False,
            shape_by=shape_by_,
            color_by=color_by_,
            point_size=point_size,
        )

    if if_plot:
        sns.set(style="whitegrid")
        methods = res.groupby("Method")
        colors = plt.get_cmap("viridis")

        fig_umap, axes_umap = plt.subplots(
            1, len(methods), figsize=(len(methods) * 5, 1 * 5), sharey=True, sharex=True
        )
        fig_umap.tight_layout()

        fig_pca, axes_pca = plt.subplots(
            1, len(methods), figsize=(len(methods) * 5, 1 * 5), sharey=True, sharex=True
        )
        fig_pca.tight_layout()

        if not (shape_by is None):
            shape_list = ["o", "s", "p", "P", "*", "h", "H", "X", "D", "d"] + [i for i in range(12)]
            shape_dict = {name: shape_list[i] for i, name in enumerate(res[shape_by].unique())}
        else:
            shape_dict = {}

        # discrete
        if res[color_by].dtype == "category":
            # plot for umap
            leg_dict = {}
            for i, (method, data) in enumerate(methods):
                ax = axes_umap[i]
                categories = data.groupby(color_by)
                color_list = np.linspace(0, 1, len(categories))

                if not shape_dict:
                    for j, (category, value) in enumerate(categories):
                        scatter = ax.scatter(
                            value["UMAP1"],
                            value["UMAP2"],
                            color=colors(color_list[j]),
                            alpha=0.5,
                            s=point_size,
                        )
                        leg_dict[category] = colors(color_list[j])
                else:
                    for j, (category, value) in enumerate(categories):
                        for name in value[shape_by].unique():
                            tmp = value[value[shape_by] == name]
                            scatter = ax.scatter(
                                tmp["UMAP1"],
                                tmp["UMAP2"],
                                color=colors(color_list[j]),
                                alpha=0.5,
                                s=point_size,
                                marker=shape_dict[name],
                            )
                        leg_dict[category] = colors(color_list[j])

                ax.set_title(method)

            # legend
            leg_patch = [Patch(facecolor=value, edgecolor=value, label=key) for key, value in leg_dict.items()]
            fig_umap.legend(
                handles=leg_patch,
                loc=[0.2, 0],
                framealpha=0,
                ncol=6,
            )
            if shape_dict:
                fig_umap.legend(
                    handles=[
                        plt.scatter([], [], marker=v, label=k, color="#74add1") for k, v in shape_dict.items()
                    ],
                    loc=2,
                    bbox_to_anchor=(0.97, 0.5),
                    framealpha=0,
                    ncol=2,
                    title=shape_by,
                )

            # text
            fig_umap.text(0.5, 0, "UMAP1", ha="center")
            fig_umap.text(0, 0.5, "UMAP2", va="center", rotation="vertical")
            fig_umap.text(0.1, -0.07, color_by)

            # plot for PCA
            for i, (method, data) in enumerate(methods):
                ax = axes_pca[i]
                categories = data.groupby(color_by)
                if not shape_dict:
                    for j, (category, value) in enumerate(categories):
                        scatter = ax.scatter(
                            value["PC1"],
                            value["PC2"],
                            color=leg_dict[category],
                            alpha=0.5,
                            s=point_size,
                        )
                else:
                    for j, (category, value) in enumerate(categories):
                        for name in value[shape_by].unique():
                            tmp = value[value[shape_by] == name]
                            scatter = ax.scatter(
                                tmp["PC1"],
                                tmp["PC2"],
                                color=leg_dict[category],
                                alpha=0.5,
                                s=point_size,
                                marker=shape_dict[name],
                            )

                ax.set_title(method)

            # legend
            fig_pca.legend(
                handles=leg_patch,
                loc=[0.2, 0],
                framealpha=0,
                ncol=6,
            )
            if shape_dict:
                fig_pca.legend(
                    handles=[
                        plt.scatter([], [], marker=v, label=k, color="#74add1") for k, v in shape_dict.items()
                    ],
                    loc=2,
                    bbox_to_anchor=(0.97, 0.5),
                    framealpha=0,
                    ncol=2,
                    title=shape_by,
                )

            # text
            fig_pca.text(0.5, 0, "PCA1", ha="center")
            fig_pca.text(0, 0.5, "PCA2", va="center", rotation="vertical")
            fig_pca.text(0.1, -0.07, color_by)

        # continuous
        elif np.issubdtype(res[color_by].dtype, np.number):
            # plot for UMAP
            for i, (method, data) in enumerate(methods):
                ax = axes_umap[i]
                # get color range
                norm = plt.Normalize(vmax=data[color_by].max(), vmin=data[color_by].min())

                if not shape_dict:
                    scatter = ax.scatter(
                        data["UMAP1"],
                        data["UMAP2"],
                        c=colors(norm(data[color_by])),
                        alpha=0.5,
                        s=point_size,
                    )
                else:
                    for name in data[shape_by].unique():
                        tmp = data[data[shape_by] == name]
                        scatter = ax.scatter(
                            tmp["UMAP1"],
                            tmp["UMAP2"],
                            c=colors(norm(tmp[color_by])),
                            alpha=0.5,
                            s=point_size,
                            marker=shape_dict[name],
                        )

                ax.set_title(method)

            fig_umap.text(0.5, 0, "UMAP1", ha="center")
            fig_umap.text(0, 0.5, "UMAP2", va="center", rotation="vertical")
            position = fig_umap.add_axes([0.2, -0.07, 0.60, 0.025])
            fig_umap.colorbar(
                mpl.cm.ScalarMappable(norm=norm, cmap=colors),
                cax=position,
                orientation="horizontal",
                label=color_by,
            )
            if shape_dict:
                fig_umap.legend(
                    handles=[
                        plt.scatter([], [], marker=v, label=k, color="#74add1") for k, v in shape_dict.items()
                    ],
                    loc=2,
                    bbox_to_anchor=(0.97, 0.5),
                    framealpha=0,
                    ncol=2,
                    title=shape_by,
                )
            # fig_umap.text(0.1, -0.07, color_by)

            # plot for PCA
            for i, (method, data) in enumerate(methods):
                ax = axes_pca[i]
                # get color range
                norm = plt.Normalize(vmax=data[color_by].max(), vmin=data[color_by].min())

                if not shape_dict:
                    scatter = ax.scatter(
                        data["PC1"],
                        data["PC2"],
                        c=colors(norm(data[color_by])),
                        alpha=0.5,
                        s=point_size,
                    )
                else:
                    for name in data[shape_by].unique():
                        tmp = data[data[shape_by] == name]
                        scatter = ax.scatter(
                            tmp["PC1"],
                            tmp["PC2"],
                            c=colors(norm(tmp[color_by])),
                            alpha=0.5,
                            s=point_size,
                            marker=shape_dict[name],
                        )
                ax.set_title(method)

            fig_pca.text(0.5, 0, "PCA1", ha="center")
            fig_pca.text(0, 0.5, "PCA2", va="center", rotation="vertical")
            position = fig_pca.add_axes([0.2, -0.07, 0.60, 0.025])
            fig_pca.colorbar(
                mpl.cm.ScalarMappable(norm=norm, cmap=colors),
                cax=position,
                orientation="horizontal",
                label=color_by,
            )
            if shape_dict:
                fig_pca.legend(
                    handles=[
                        plt.scatter([], [], marker=v, label=k, color="#74add1") for k, v in shape_dict.items()
                    ],
                    loc=2,
                    bbox_to_anchor=(0.97, 0.5),
                    framealpha=0,
                    ncol=2,
                    title=shape_by,
                )
            # fig_pca.text(0.1, -0.07, color_by)

        # if not numeric or factor, then try to handle as the discrete one
        else:
            # plot for umap
            leg_dict = {}
            for i, (method, data) in enumerate(methods):
                ax = axes_umap[i]
                categories = data.groupby(color_by)
                color_list = np.linspace(0, 1, len(categories))

                if not shape_dict:
                    for j, (category, value) in enumerate(categories):
                        scatter = ax.scatter(
                            value["UMAP1"],
                            value["UMAP2"],
                            color=colors(color_list[j]),
                            alpha=0.5,
                            s=point_size,
                        )
                        leg_dict[category] = colors(color_list[j])
                else:
                    for j, (category, value) in enumerate(categories):
                        for name in value[shape_by].unique():
                            tmp = value[value[shape_by] == name]
                            scatter = ax.scatter(
                                tmp["UMAP1"],
                                tmp["UMAP2"],
                                color=colors(color_list[j]),
                                alpha=0.5,
                                s=point_size,
                                marker=shape_dict[name],
                            )
                        leg_dict[category] = colors(color_list[j])

                ax.set_title(method)

            # legend
            leg_patch = [Patch(facecolor=value, edgecolor=value, label=key) for key, value in leg_dict.items()]
            fig_umap.legend(
                handles=leg_patch,
                loc=[0.2, 0],
                framealpha=0,
                ncol=6,
            )
            if shape_dict:
                fig_umap.legend(
                    handles=[
                        plt.scatter([], [], marker=v, label=k, color="#74add1") for k, v in shape_dict.items()
                    ],
                    loc=2,
                    bbox_to_anchor=(0.97, 0.5),
                    framealpha=0,
                    ncol=2,
                    title=shape_by,
                )

            # text
            fig_umap.text(0.5, 0, "UMAP1", ha="center")
            fig_umap.text(0, 0.5, "UMAP2", va="center", rotation="vertical")
            fig_umap.text(0.1, -0.07, color_by)

            # plot for PCA
            for i, (method, data) in enumerate(methods):
                ax = axes_pca[i]
                categories = data.groupby(color_by)
                if not shape_dict:
                    for j, (category, value) in enumerate(categories):
                        scatter = ax.scatter(
                            value["PC1"],
                            value["PC2"],
                            color=leg_dict[category],
                            alpha=0.5,
                            s=point_size,
                        )
                else:
                    for j, (category, value) in enumerate(categories):
                        for name in value[shape_by].unique():
                            tmp = value[value[shape_by] == name]
                            scatter = ax.scatter(
                                tmp["PC1"],
                                tmp["PC2"],
                                color=leg_dict[category],
                                alpha=0.5,
                                s=point_size,
                                marker=shape_dict[name],
                            )

                ax.set_title(method)

            # legend
            fig_pca.legend(
                handles=leg_patch,
                loc=[0.2, 0],
                framealpha=0,
                ncol=6,
            )
            if shape_dict:
                fig_pca.legend(
                    handles=[
                        plt.scatter([], [], marker=v, label=k, color="#74add1") for k, v in shape_dict.items()
                    ],
                    loc=2,
                    bbox_to_anchor=(0.97, 0.5),
                    framealpha=0,
                    ncol=2,
                    title=shape_by,
                )

            # text
            fig_pca.text(0.5, 0, "PCA1", ha="center")
            fig_pca.text(0, 0.5, "PCA2", va="center", rotation="vertical")
            fig_pca.text(0.1, -0.07, color_by)

        return {"p_pca": fig_pca, "p_umap": fig_umap}

    else:
        return res
