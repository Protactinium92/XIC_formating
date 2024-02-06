# -*- coding: utf-8 -*-

"""
    Part for the data vizualization
"""

from import_package import *


def sheet_exists(workbook, sheet_name):
    """ Test if the shit exist """
    try:
        var = workbook[sheet_name]
        rep = True
    except KeyError:
        rep = False
    return rep


def add_table(writer, sheets, dataframe, table_name, lig=3, col=1, stcol=0):
    """
        Create a table in an Excel sheet with a header

    :param sheets: name of desired sheet
    :param dataframe: desired table
    :param table_name: name of the desired table
    :param lig: lign number
    :param col: column number
    :param stcol: start column of the table
    """
    # Creation of the workbook
    wb = writer.book

    # Creation of the sheet if it not exists
    if not sheet_exists(wb, sheets):
        wb.create_sheet(sheets)

    # Add table
    dataframe.to_excel(writer, sheet_name=sheets, startrow=lig, startcol=stcol, header=True, index=False)
    # Add a header above the table
    feuille = wb[sheets]
    feuille.cell(row=lig, column=col).value = table_name
    print(f'Add {table_name} in {sheets}')


def output(new_excel_file, looking_metacell, full_table, percentage_imputation, imputation_value, summary_psm,
           summary_spc, percentage_signif2x2, threshold_FC, threshold_pvalue):
    if os.path.exists(new_excel_file):
        os.remove(new_excel_file)

    # Determination of indic where the metacell start to add the percentage table at the right column
    start_meta = looking_metacell.columns[0]
    # start_meta = full_table.columns.get_loc(start_meta)
    # Writing to the Excel file
    with pd.ExcelWriter(new_excel_file, engine='openpyxl') as writer:
        add_table(writer, "Global", full_table, "Synthesis")
        add_table(writer, "Statistics", percentage_imputation, table_name="Percentage of imputation")
        add_table(writer, "Statistics", imputation_value,
                  lig=len(percentage_imputation) + 6,
                  table_name="Imputed value for the MEC and the POV")
        add_table(writer, "Statistics", summary_psm,
                  lig=len(percentage_imputation) + len(imputation_value) + 9,
                  table_name="PROSTAR Spectral Count (FDR < 1%)")
        if summary_spc is not None:
            add_table(writer, "Statistics", summary_spc,
                      lig=len(percentage_imputation) + len(imputation_value) + len(summary_psm) + 12,
                      table_name="PROLINE Spectral Count (FDR < 1%)")
            add_table(writer, "Statistics", percentage_signif2x2,
                      lig=len(percentage_imputation) + len(imputation_value) + len(summary_psm) + len(summary_spc) + 15,
                      table_name=(f"Percentage of Significant, comparison 2 by 2:  log10 adjusted p-value threshold: {threshold_pvalue}, "
                                  f"log2 FC threshold: {threshold_FC}"))
        else:
            add_table(writer, "Statistics", percentage_signif2x2,
                      lig=len(percentage_imputation) + len(imputation_value) + len(summary_psm) + 12,
                      table_name=(f"Percentage of Significant, comparison 2 by 2: log10 adjusted p-value threshold: {threshold_pvalue}, "
                                  f"log2 FC threshold: {threshold_FC}"))


def styling(new_excel_file, looking_metacell, full_table):
    # Load file
    workbook = openpyxl.load_workbook(new_excel_file)
    sheet = workbook["Global"]

    # Reduce size of the metacell columns
    print("Formatting metacell columns size")
    i = 1
    for col in looking_metacell:
        index_col = full_table.columns.get_loc(col) + 1
        sheet.column_dimensions[openpyxl.utils.get_column_letter(index_col)].width = 3
        print(f"{i}/{looking_metacell.shape[1]}")
        i += 1

    # increase the width of the A1 cell on Statistics sheet
    sheet_stats = workbook["Statistics"]
    sheet_stats.column_dimensions["A"].width = 26

    # Saving
    print("Saving")
    workbook.save(new_excel_file)