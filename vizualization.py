# -*- coding: utf-8 -*-

"""
    Part for the data vizualization
"""

from import_package import *

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


def volcano_plt(looking_FC, looking_2x2compare, threshold_pvalue, threshold_FC, save_dir):
    # Volcano plot
    for col in looking_FC:
        corresponding_pval = looking_2x2compare.columns.get_loc(col) + 1
        corresponding_pval = looking_2x2compare.columns[corresponding_pval]
        group = looking_2x2compare.columns.get_loc(col) + 2
        group = looking_2x2compare.columns[group]
        volcano_plot(dataframe=looking_2x2compare,logfc=col, adjpva=corresponding_pval, color_group=group,
                     threshold_pval=threshold_pvalue, threshold_fc=threshold_FC,
                     directory=save_dir)
        print(f"Volcano plot {looking_FC.columns.get_loc(col) + 1}/{len(looking_FC.columns)}")
