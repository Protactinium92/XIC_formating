# -*- coding: utf-8 -*-
"""
    Part dedicated for the importation of excel file

"""

from import_package import *

def import_excel(link_dir, file, sheet_number):
    """
        Import one sheet of the Excel file as a pandas table

    :param link_dir: link of the input directory
    :param file: name of the Excel file with the extension
    :param sheet_number: Desired sheet
    :return: pandas dataframe
    """
    print(f'Importation of the Excel {file} file, sheet {sheet_number} on going')
    directory = link_dir
    excel_file = file
    sheet_list = pd.ExcelFile(directory + excel_file).sheet_names
    print(f'Import of the Excel {file} file, sheet {sheet_list[sheet_number]} complete')
    return pd.read_excel(directory + excel_file, sheet_name=sheet_list[sheet_number])

def transform_into_df (link_dir_input, excel_XIC, excel_GO, excel_SPC):
    """
        Function to transform the excel into dataframe
            - quanti => XIC aboundance
            - feature => XIC all data
            - meta => XIC metadata (conditions)
            - protein set => SPC all data in protein set sheet
            - go_set => GO annotation

    :param link_dir_input: link of the input directory
    :param excel_XIC: excel file of XIC
    :param excel_GO: excel file of the annotation GO (can be just a .xlsx)
    :param excel_SPC: excel file of SPC (can be just a .xlsx)
    :return: return corresponding dataframe
    """


    # Choose the sheet in the excel file
    sheet_quanti = 0  # 1st sheet put in the python list
    sheet_meta = 1
    sheet_feature = 2
    sheet_protein = 3
    sheet_go = 0

    # Dataframe import
    quanti = import_excel(link_dir_input, excel_XIC, sheet_quanti)
    if quanti.columns[0] != 'ID':
        quanti.rename(columns={quanti.columns[0]: 'ID'}, inplace=True)
    meta = import_excel(link_dir_input, excel_XIC, sheet_meta)
    feature = import_excel(link_dir_input, excel_XIC, sheet_feature)

    # For GO and SPC excel file check if they exist (are more than .xlsx)
    # if not return None
    if excel_GO != ".xlsx":
        go_set = import_excel(link_dir_input, excel_GO, sheet_go)
        go_set.columns = [col.lower() for col in go_set.columns]  # lower case the name of each column
        if go_set.columns[0] != 'ID':
            go_set.rename(columns={go_set.columns[0]: 'ID'}, inplace=True)
    else:
        go_set = None
    if excel_SPC != ".xlsx":
        prot_set = import_excel(link_dir_input, excel_SPC, sheet_protein)
    else:
        prot_set = None

    return quanti, meta, feature, go_set, prot_set

