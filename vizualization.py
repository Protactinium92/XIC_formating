# -*- coding: utf-8 -*-

"""
    Part for the data vizualization
"""

from import_package import *


def fullheatmap(quanti, save_dir):
    # Clustered heat map from all normalized abundance

    heatmap_all = sns.clustermap(quanti, method='average', metric='euclidean',
                                 cbar_kws={"label": "Normalized Abundance", "orientation": "horizontal"},
                                 cbar_pos=(0.75, 0.90, 0.1, 0.05))

    # Add a title to the graph
    heatmap_all.fig.suptitle("Heatmap of all proteins", fontsize=20, y=1.05)
    heatmap_all.savefig(f"{save_dir}/heatmap_all.svg", bbox_inches='tight')
    plt.clf()  # Clear the current figure to start a new one


def volcano_plot(dataframe, logfc, adjpva, color_group, threshold_pval, threshold_fc, directory):
    """
        Generate a volcano plot of a 2 by 2 comparison
    :param dataframe: the dataframe for the vulcanoplot with the accession for the name, description,
                        and the logFC adjusted p-val and significant for each comparison
    :param logfc: column of the LogFC from a dataframe
    :param adjpva: column of the adjusted p-value from a dataframe
    :param color_group: Color groups if they are significant
    :param threshold_pval: threshold of the p-value
    :param threshold_fc: threshold of the LogFC
    :param directory: Path of the directory to save the graph
    :return: Save the graph
    """

    colors = dataframe[color_group].unique()
    colors = sorted(colors, reverse=True)

    # Create volcano plot with Plotly Express
    fig = px.scatter(dataframe, x=logfc, y=adjpva,
                     color=color_group, color_discrete_sequence=['orange', 'lightgrey'],
                     category_orders={color_group: colors},
                     labels={'color': 'Significant'},
                     hover_data={'description': True},
                     hover_name='accession')

    # Add threshold lines
    fig.add_hline(y=threshold_pval, line_width=2, line_dash="dash", line_color="black", layer="below")
    fig.add_vline(x=threshold_fc, line_width=2, line_dash="dash", line_color="black", layer="below")
    fig.add_vline(x=-threshold_fc, line_width=2, line_dash="dash", line_color="black", layer="below")

    # Styling & legend part
    fig.update_xaxes(title_text='<b>Log<sub>2</sub> Fold Change</b>', showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(title_text='<b>Log<sub>10</sub> Adjusted <i>p</i>-value</b>', showline=True, linewidth=2,
                     linecolor='black')
    fig.update_layout(
        legend_title_text='<b>Significant</b>',
        title=f'<b>{logfc.replace("(", "").replace(")", "").replace("_", " ").replace(" logFC", "")}</b>',
        plot_bgcolor='white'
        )

    # Save as HTML (enable interactive)
    html_file = f"{directory}/volcano_{'_'.join(logfc.split('_')[:-1])}.html"
    fig.write_html(html_file)


def volcanos(looking_FC, looking_2x2compare, threshold_pvalue, threshold_FC, save_dir):
    # Volcano plot for all comparisons
    for col in looking_FC:
        corresponding_pval = looking_2x2compare.columns.get_loc(col) + 1
        corresponding_pval = looking_2x2compare.columns[corresponding_pval]
        group = looking_2x2compare.columns.get_loc(col) + 2
        group = looking_2x2compare.columns[group]
        volcano_plot(dataframe=looking_2x2compare,logfc=col, adjpva=corresponding_pval, color_group=group,
                     threshold_pval=threshold_pvalue, threshold_fc=threshold_FC,
                     directory=save_dir)
        print(f"Volcano plot {looking_FC.columns.get_loc(col) + 1}/{len(looking_FC.columns)}")


def scree_plot(scree, save_dir):

    # 1 ) Scree plots
    scree.plot.bar(x="Dimension", y="% Explaine Variance")
    plt.text(5, 18, "5%")  #
    plt.axhline(y=5, linewidth=0.5, color="dimgray", linestyle="--")
    plt.savefig(f"{save_dir}/scree_polt_ExplainVar.svg")
    plt.clf()  # Clear the current figure to start a new one
    print("Plot Explain variance done")

    scree.plot.bar(x="Dimension", y="% Cumulative")
    plt.text(0.5, 90, "95%")
    plt.axhline(y=95, linewidth=0.5, color="dimgray", linestyle="--")
    plt.savefig(f"{save_dir}/scree_polt_CumVar.svg")
    plt.clf()  # Clear the current figure to start a new one
    print("Plot Cumumul Explain variance done")


def plot_PCA(table_pca, scree, save_dir, pca, table_stats):
    unique_conditions = table_pca['Condition'].unique()

    # Graph
    fig, ax = plt.subplots()

    scatter_handles = []
    for condition in unique_conditions:
        subset_df = table_pca[table_pca['Condition'] == condition]
        handle = plt.scatter(subset_df['PCA1'],
                             subset_df['PCA2'],
                             label=f'{condition}',
                             alpha=0.7, edgecolor='k')
        scatter_handles.append(handle)

    for index, sample in table_pca.iterrows():
        ax.annotate(" ".join(index.split("_")[1:]), (sample['PCA1'], sample['PCA2']), fontsize=4)

    plt.legend(handles=scatter_handles, loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.title('PCA Analysis')
    plt.axhline(y=0, linewidth=0.5, color="dimgray", linestyle="--")
    plt.axvline(x=0, linewidth=0.5, color="dimgray", linestyle="--")
    plt.xlabel(f'1st dim ({scree.at[0, "% Explaine Variance"]}%)')
    plt.ylabel(f'2nd dim ({scree.at[1, "% Explaine Variance"]}%)')
    plt.tight_layout()

    plt.savefig(f"{save_dir}/PCA.svg", dpi=300)
    plt.clf()  # Clear the current figure to start a new one
    print("Plot PCA done")

    # Multiple PCA

    total_5 = np.sum(pca.explained_variance_ratio_[:5]*100).round(1)
    label = {
        "PCA"+str(i+1): f"PC {i+1} ({var:.1f}%)"
        for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }
    fig = px.scatter_matrix(
        table_pca,
        dimensions=['PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5'],
        color=table_pca["Condition"],
        title=f'Total Explained Variance: {total_5}%',
        labels=label,
        hover_name=table_pca.index.str.replace("_", " ").str.replace("abundance", ""))

    fig.write_html(f"{save_dir}/PCA_multiple.html")
    print("Plot Multi PCA done")

    # 3D

    fig = px.scatter_3d(
        table_pca, x='PCA1', y='PCA2', z='PCA3',
        color=table_pca['Condition'],
        title=f'Total Explained Variance: {total_5}%',
        labels=label,
        hover_name=table_pca.index.str.replace("_", " ").str.replace("abundance", ""))
    fig.write_html(f"{save_dir}/PCA_3D.html")
    print("Plot 3D PCA done")

    # 3 ) Loading

    table_stats.drop('Condition', axis=1, inplace=True)

    table_stats = table_stats.T
    loadings_matrix = pca.components_.T * np.sqrt(pca.explained_variance_)

    table_stats = table_stats.T
    loadings_df = pd.DataFrame(loadings_matrix, columns=[f'PC{j+1}' for j in range(pca.n_components_)],
                               index=table_stats.columns)

    loadings_df.reset_index(inplace=True)

    for i in range(5):
        # plotly graph object
        trace = go.Bar(
            x=loadings_df['index'],
            y=loadings_df[f'PC{i+1}'],
        )

        # disposition
        layout = go.Layout(
            title=f'Loading PC {i+1}',
            xaxis=dict(title='Proteins'),
            yaxis=dict(title='Loading'),
            plot_bgcolor='rgba(255, 255, 255, 0)')

        # plot the figure
        fig = go.Figure(data=[trace], layout=layout)
        fig.write_html(f"{save_dir}/PCA_loading{i+1}.html")
        print(f"Plot loading PC{i + 1}/5 done")

