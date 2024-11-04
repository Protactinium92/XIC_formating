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


def add_table(file, sheets, dataframe, table_name, lig=3, col=1, stcol=0):
    """
        Create a table in an Excel sheet with a header

    :param file: Path of the excel file
    :param sheets: name of desired sheet
    :param dataframe: desired table
    :param table_name: name of the desired table
    :param lig: lign number
    :param col: column number
    :param stcol: start column of the table
    """

    print(f'Addition of {table_name} in {sheets} on going')
    try:
        wb = openpyxl.load_workbook(file)
    except FileNotFoundError:
        wb = openpyxl.Workbook()
        wb.save(file)
        wb = openpyxl.load_workbook(file)

    # Check if the sheet exist, create it if not
    if not sheet_exists(wb, sheets):
        wb.create_sheet(sheets)

    # Add table inf the Excel File
    with pd.ExcelWriter(file, engine='openpyxl', mode='a', if_sheet_exists='overlay') as writer:
        dataframe.to_excel(writer, sheet_name=sheets, startrow=lig, startcol=stcol, header=True, index=False)

        # Add a header above the table
        worksheet = writer.sheets[sheets]
        worksheet.cell(row=lig, column=col).value = table_name
    print(f'Addition of {table_name} in {sheets} complete')


def output(new_excel_file, full_table, percentage_imputation, imputation_value, summary_psm,
           summary_spc, percentage_signif2x2, threshold_FC, threshold_pvalue):

    # Writing to the Excel file
    add_table(new_excel_file, "Quanti XIC", full_table, "Synthesis", lig=15)
    add_table(new_excel_file, "Statistics", percentage_imputation, table_name="Percentage of imputation")
    add_table(new_excel_file, "Statistics", imputation_value,
              lig=len(percentage_imputation) + 6,
              table_name="Imputed value for the MEC")
    add_table(new_excel_file, "Statistics", summary_psm,
              lig=len(percentage_imputation) + len(imputation_value) + 9,
              table_name="PROSTAR Spectral Count (FDR < 1%)")
    if summary_spc is not None:
        add_table(new_excel_file, "Statistics", summary_spc,
                  lig=len(percentage_imputation) + len(imputation_value) + len(summary_psm) + 12,
                  table_name="PROLINE Spectral Count (FDR < 1%)")
        add_table(new_excel_file, "Statistics", percentage_signif2x2,
                  lig=len(percentage_imputation) + len(imputation_value) + len(summary_psm) + len(summary_spc) + 15,
                  table_name=(
                      f"Percentage of Significant, comparison 2 by 2:  log10 adjusted p-value threshold: {threshold_pvalue}, "
                      f"log2 FC threshold: {threshold_FC}"))
    else:
        add_table(new_excel_file, "Statistics", percentage_signif2x2,
                  lig=len(percentage_imputation) + len(imputation_value) + len(summary_psm) + 12,
                  table_name=(
                      f"Percentage of significant protein, comparison 2 by 2: log10 adjusted p-value threshold: {threshold_pvalue}, "
                      f"log2 FC threshold: {threshold_FC}"))


def styling(new_excel_file, looking_metacell, full_table):
    # Load file
    workbook = openpyxl.load_workbook(new_excel_file)
    sheet = workbook["Quanti XIC"]

    # Reduce size of the metacell columns
    print("Formatting metacell columns size")
    count = 0
    progress_bar(count, looking_metacell.shape[1])
    for col in looking_metacell:
        index_col = full_table.columns.get_loc(col) + 1
        sheet.column_dimensions[openpyxl.utils.get_column_letter(index_col)].width = 3
        count += 1
        progress_bar(count, looking_metacell.shape[1], "Formatting metacell columns size")

    # increase the width of the A1 cell on Statistics sheet
    sheet_stats = workbook["Statistics"]
    sheet_stats.column_dimensions["A"].width = 26

    # Saving
    print("Saving")
    workbook.save(new_excel_file)