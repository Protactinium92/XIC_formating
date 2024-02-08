# -*- coding: utf-8 -*-
"""
    Part dedicated for the statistics and counting

"""

from import_package import *


def counting(looking_metacell,  looking_2x2compare, psm, spc):
    """
        Functions to count
            - the percentage  of imputation
            - the percentage of significatif
            - the number of protein and spectre in psm
            - the number of protein and spectre in spc if it exists
    :param looking_metacell:
    :return: dataframe
    """

    # Count the number of each imputation and give a percentage
    percentage_imputation = looking_metacell.apply(lambda col: col.value_counts() / len(col) * 100).round(0)
    percentage_imputation.reset_index(inplace=True)
    percentage_imputation["index"].replace({"D": "Direct Identification",
                                   "M": "Missing Entire Condition",
                                   "P": "Partially Observed Values",
                                   "R": "Recovery Cross Assignment"}, inplace=True)

    # Count the number of Significant 2 by 2
    looking_signif = looking_2x2compare.filter(regex=f"Significant$", axis=1)
    percentage_signif2x2 = looking_signif.apply(lambda col: col.value_counts() / len(col) * 100).round(0)
    percentage_signif2x2.drop("No", inplace=True)
    percentage_signif2x2.reset_index(inplace=True)

    # Count the number of protein in psm (spectral count after the normalization)
    looking_psm = psm.filter(regex=f"psm_count", axis=1)
    nbr_prot = looking_psm.apply(lambda col: np.count_nonzero(col))  # number of cell different to 0
    nbr_spectre = looking_psm.sum().astype(int)
    # Summarize this information in one table
    summary_psm = pd.DataFrame({'Number of protein': nbr_prot, 'Number of spectre': nbr_spectre})
    summary_psm = summary_psm.T
    summary_psm.reset_index(inplace=True)  # Put the name of the lign in the first column
    summary_psm.rename(columns={'index': ""}, inplace=True)

    if spc is not None:  # If spc exist
        looking_spc = spc.filter(regex=f"psm_count", axis=1)
        nbr_prot = looking_spc.apply(lambda col: np.count_nonzero(col))  # number of cell different to 0
        nbr_spectre = looking_spc.sum().astype(int)
        # Summarize this information in one table
        summary_spc = pd.DataFrame({'Number of protein': nbr_prot, 'Number of spectre': nbr_spectre})
        summary_spc = summary_spc.T
        summary_spc.reset_index(inplace=True)  # Put the name of the lign in the first column
        summary_spc.rename(columns={'index': ""}, inplace=True)
    else:
        summary_spc = None

    return percentage_imputation, percentage_signif2x2, summary_psm, summary_spc


def find_M(looking_metacell, quanti):

    # Find the value of M and P

    # 1) Search the position (n° lign) of the first occurrence of M & P for each column in looking_metacell
    search = ['M', 'P']
    position = looking_metacell.apply(lambda col: [col.eq(val).idxmax() for val in search])
    position.index = search

    # 2 ) Creation of a new df and copy the value in quanti following the position found
    # Names between position and quanti are not exactly the sames
    imputation_value = pd.DataFrame(index=search, columns=quanti.columns)
    for col in position.columns.tolist():
        indice_col = position.columns.get_loc(col)  # Name col position -> n° column
        for i in range(len(position)):
            value_lign = position.iat[i, indice_col]  # n° lign found
            imputation_value.iat[i, indice_col] = quanti.iat[value_lign, indice_col]
    imputation_value.drop("P", inplace=True)
    imputation_value.reset_index(inplace=True)
    imputation_value["index"].replace({"M": "Missing Entire Condition"}, inplace=True)

    return imputation_value


def my_pca(table_stats):
    """
        Make the PCA analysis on the XIC abundance
    :param table_stats: table of XIC abundance + conditions
    :return: table of th 5th first PCA and the scree values
    """

    x = table_stats.drop('Condition', axis=1)
    x_standardized = StandardScaler().fit_transform(x)

    pca = PCA()
    x_pca = pca.fit_transform(x_standardized)

    table_pca = table_stats['Condition']
    for i in range(5):
        table_pca[f'PCA{i + 1}'] = x_pca[:, i]

    scree = pd.DataFrame({"Dimension": ["Dim" + str(i + 1) for i in range(len(pca.explained_variance_ratio_))],
                          "% Explaine Variance": (pca.explained_variance_ratio_ * 100).round(1),
                          "% Cumulative": (np.cumsum(pca.explained_variance_ratio_ * 100)).round(1)
                          })

    return table_pca, scree, pca

