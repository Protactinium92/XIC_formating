# -*- coding: utf-8 -*-
"""
    Part dedicated for the transformation and the manipulations of dataframe

"""

from import_package import *

def keep_only_col(name_col, dataframe, compare=False):
    """
        Function to format the table to keep only desired column which start by a specific name
        Then replace empty cells by 0

    :param name_col: All the column which are started with this name, None: if it's not needed to have any column but
                    only the first usual
    :param dataframe: Table to format
    :param compare If True take the last columns which have no recurrent names, by default is false to skip it
    :return: New dataframe with the selected column and empty cell replace by 0
    """
    # test is the name is not empty
    if name_col:
        keep_col = [col for col in dataframe.columns if col.startswith(name_col)]
    else:
        keep_col = []

    # Addition of the conditions columns (all last column) if we want it
    if compare is not False:
        # Conditions columns are exactly after the last column of metacell (the last value of the keep_col list)
        pos_first_cond_col = dataframe.columns.get_loc(keep_col[-1])
        keep_col = keep_col + dataframe.columns.tolist()[pos_first_cond_col + 1:]

    keep_col = ['accession', 'description', 'samesets_accessions', 'subsets_accessions'] + keep_col
    return dataframe[keep_col].fillna(0)


def scientific(value):

    """
        Transform the vailue into the scientific format for all value under 0.05
        The value will be transform in str
    :param value: Value to be transform
    :return: return the value itself
    """

    if isinstance(value, float):
        if abs(value) < 0.05:
            value = np.format_float_scientific(value, precision=2, exp_digits=1)
        else:
            value = round(value, 2)  # still a float
    return value


def manipulation(feature, prot_set):
    # Keep only interested columns
    metacell_comp = keep_only_col('metacell', feature, True)
    psm = keep_only_col('psm', dataframe=feature)
    begining = keep_only_col(name_col=None, dataframe=feature)  # To order the table and have only the description...
    if prot_set is not None:
        spc = keep_only_col('psm_count', prot_set)
    else:
        spc = None

    # Simplification of the value on the df metacell
    for col in metacell_comp.columns[4:]:
        metacell_comp[col] = metacell_comp[col].replace(
            {'Quant. by direct id': 'D',
             'Quant. by recovery': 'R',
             'Imputed MEC': 'M',
             'Imputed POV': 'P'})
        metacell_comp[col] = metacell_comp[col].replace(".", ",")

    # Addition of the adjusted p-value
    looking_pval = metacell_comp.filter(regex=f"pval$", axis=1)  # keep only the p-value col (end of the sentence)
    for col in looking_pval.columns:
        adjusted_pval = multipletests(metacell_comp[col], method='fdr_bh')[1]
        metacell_comp.insert(metacell_comp.columns.get_loc(col) + 1, f"{col}adjusted", adjusted_pval)
    print(f"adjusted p-value: {len(looking_pval.columns)}/{len(looking_pval.columns)}")

    looking_adjpv = metacell_comp.filter(regex=f"adjusted$",
                                         axis=1)  # Keep non-modified adjusted p-value for stats part
    metacell_comp = metacell_comp.map(scientific)  # transform into scientific format (str)

    return metacell_comp, psm, begining, spc, looking_adjpv

def quanti_meta_transform (quanti, meta):

    """
        Transpose the quanti dataframe, et del the replicat column of the meta dataframe
    :param quanti: dataframe quanti
    :param meta: dataframe meta
    :return: transposed of quanti and meta simplified
    """

    # Transpose quanti
    trans_quanti = quanti.T
    trans_quanti.columns = trans_quanti.iloc[0]
    trans_quanti = trans_quanti[1:]
    trans_quanti.reset_index(inplace=True)
    if trans_quanti.columns[0] != 'Sample':
        trans_quanti.rename(columns={trans_quanti.columns[0]: 'Sample'}, inplace=True)

    # Preparation Meta
    del meta['Bio.Rep']
    if meta.columns[0] != 'Sample':
        meta.rename(columns={meta.columns[0]: 'Sample'}, inplace=True)

    return trans_quanti, meta

def merging(begining, quanti, psm, go_set, metacell_comp, meta, trans_quanti):
    """
        Merge all the dataframe in one to add to the excel file
        Merge meta + trans_quanti dataframe to the stats

    :param begining:
    :param quanti:
    :param psm:
    :param metacell_comp:
    :param meta:
    :param trans_quanti:
    :return:
    """
    # Merge all tab in one
    full_table = pd.merge(begining, quanti.round(2))
    full_table = pd.merge(full_table, psm)
    full_table = pd.merge(full_table, metacell_comp)
    if go_set is not None:
        full_table = pd.merge(full_table, go_set, how='left')
    table_stats = pd.merge(meta, trans_quanti, on='Sample')
    table_stats.set_index("Sample", inplace=True)
    print("Compilation OK")

    return full_table, table_stats


def comparison2x2_and_exp_bea(metacell_comp, looking_adjpv, threshold_FC, threshold_pvalue):

    # Simplified dataframe comparison + table for the beatrice export

    # 1) preparation of the 2 tables
    looking_2x2compare = metacell_comp[["accession", "description"]]
    export_beatrice = metacell_comp[["accession", "description"]]

    # 2) preparation of the lists in order to insert the adjusted pval next to their corresponding logFC
    looking_FC = metacell_comp.filter(regex=f"logFC$", axis=1)
    looking_FC = looking_FC.map(lambda x: float(x))
    adjpval_nolog10 = looking_adjpv  # keep the adjusted p-value without log10 transformation (for Beatrice export)
    nolog10_list = adjpval_nolog10.columns.tolist()
    looking_adjpv = looking_adjpv.map(lambda x: -np.log10(x) if x != 0 else -np.log10(0.001))  # transform to -log10
    list_adjpv = looking_adjpv.columns.tolist()

    # 3) Ordering the columns
    ordering_col = ["accession", "description"]
    ordering_col_nolog10 = ["accession", "description"]
    for col in looking_FC:
        ordering_col.append(col)
        ordering_col.append(list_adjpv.pop(0))
        ordering_col_nolog10.append(col)
        ordering_col_nolog10.append(nolog10_list.pop(0))

    # 4) Merges tables and reindex with the right order thanks to the ordering columns
    looking_2x2compare = pd.merge(looking_2x2compare, looking_FC, left_index=True, right_index=True)
    looking_2x2compare = pd.merge(looking_2x2compare, looking_adjpv, left_index=True, right_index=True)
    looking_2x2compare = looking_2x2compare.reindex(columns=ordering_col)
    export_beatrice = pd.merge(export_beatrice, looking_FC, left_index=True, right_index=True)
    export_beatrice = pd.merge(export_beatrice, adjpval_nolog10, left_index=True, right_index=True)
    export_beatrice = export_beatrice.reindex(columns=ordering_col_nolog10)
    export_beatrice.rename(columns=lambda x: x.replace("(", "").replace(")", "").replace("_", " "), inplace=True)


    # Add column for signficatifs values
    df_signif = pd.DataFrame()
    for col in looking_FC:
        col.replace("(", "").replace(")", "").replace("_", " ")
        corresponding_pval = looking_2x2compare.columns.get_loc(col) + 1
        corresponding_pval = looking_2x2compare.columns[corresponding_pval]
        df_signif = pd.DataFrame(looking_2x2compare[[col, corresponding_pval]])
        df_signif[f"{col}_Significant"] = df_signif.apply(
            lambda row: "Yes" if abs(row[col]) > threshold_FC and row[corresponding_pval] > threshold_pvalue
            else "No", axis=1)
        col_signif = df_signif[f"{col}_Significant"]
        looking_2x2compare.insert(looking_2x2compare.columns.get_loc(col) + 2, f"{col}_Significant", col_signif)

    return looking_FC, looking_2x2compare, export_beatrice