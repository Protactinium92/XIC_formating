# -*- coding: utf-8 -*-

"""
    Script to automatize the reporting of the XIC analysis

    Input : Excel File with multiple sheet
    Output : Excel file

    To perform the script, the following packages needs to be installed before using the command line
    pip install
        - pandas
        - numpy
        - scikit-learn as sklearn
        - openpyxl
        - plotly
        - matplotlib
        - seaborn


    This part is the main script

    Creation Date : 2024-01-23
    Author : Pierre-Alexandre
    Version : 1
    Edit with python 3.12.0
"""

# -------------------------------------------------------IMPORTS -------------------------------------------------------

try:
    print("Initialization")

    from import_package import *
    from importation_excel import transform_into_df
    import dataframe_manipulation as manip
    import my_stats as st
    import visualization as viz
    from output_excel import output, styling

    # ------------------------------------------------------------------------------------------------------------------
    #                                                       FUNCTIONS
    # ------------------------------------------------------------------------------------------------------------------


    def input_file():
        """
            Ask the user thanks to the input in the console
            1) The directory is the same for all files
            2) Need the proline file with the xlsx extension
            3) Need the go annotation file in xlsx

            To improve the code we should add some condition to check if the user write well and also if the files exist

        :return: variables in str characters
        """
        link = input('Write directory link (windows link): ')
        link = link.replace('\\', '/') + "/"
        excel_file_XIC = input('Write the name of the PROSTAR excel file (without .xlsx): ')
        excel_file_XIC = excel_file_XIC + ".xlsx"
        excel_file_SpC = input('Write the name of the PROLINE excel file (without .xlsx): ')
        excel_file_SpC = excel_file_SpC + ".xlsx"
        excel_file_go = input('Write the name of the annotation GO excel file (without .xlsx): ')
        excel_file_go = excel_file_go + ".xlsx"
        return link, excel_file_XIC, excel_file_SpC, excel_file_go


    # ------------------------------------------------------------------------------------------------------------------
    #                                                       MAIN
    # ------------------------------------------------------------------------------------------------------------------

    # ---------------------------------------------------- INPUT -------------------------------------------------------

    # Importation of the Prostar file SPC and the annotation file, they must be in the same directory
    link_dir_input, excel_XIC, excel_SPC, excel_GO = input_file()

    # in comment some lign for fast testing code
    # link_dir_input = r'C:\Users\hpier\Documents\code\Python\StageM1\XIC_DDA'
    # link_dir_input = link_dir_input.replace('\\', '/') + "/"
    # excel_XIC = '2023_S43_VMontoya_E22_Hyp_Test.xlsx'
    # excel_SPC = 'XIC_Proline_2023-S43_VMontoya_TIMS_export_jc_16012024_PPSE.xlsx'
    # excel_GO = 'TAIR_NEW_All-Annotations_For-Data-Analyzer-Proline_12022018.xlsx'


    # Signification following the adjusted p-val and the Fold-Change
    threshold_pvalue = float(input("Define the log10 adjusted p-value threshold: "))
    threshold_FC = float(input("Define the log2 Fold Change threshold in absolute value: "))

    # -------------------------------------------------- OUTPUT --------------------------------------------------------

    # name of the new file for export :
    new_excel_file = input('Write the name of the results Excel file (without .xlsx): ')
    # new_excel_file = 'Template Quanti XIC'

    # Complete name of the file with the path
    new_excel_file = os.path.join(link_dir_input, new_excel_file + ".xlsx")

    # Saving directory
    save_dir = link_dir_input + "export"

    # Check and creation of the directory if it not exists
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)


    # ------------------------------------------------ EXCEL FILE IMPORT -----------------------------------------------

    quanti, meta, feature, go_set, prot_set = transform_into_df(link_dir_input, excel_XIC, excel_GO, excel_SPC)


    # ------------------------------------------- DATAFRAME MANIPULATION -----------------------------------------------

    metacell_comp, psm, begining, spc, looking_adjpv = manip.manipulation(feature, prot_set)
    trans_quanti, meta = manip.quanti_meta_transform(quanti, meta)
    full_table, table_stats = manip.merging(begining, quanti, psm, go_set, metacell_comp, meta, trans_quanti)
    looking_FC, looking_2x2compare = manip.comparison2x2(metacell_comp, looking_adjpv,
                                                                                      threshold_FC, threshold_pvalue)

    # -------------------------------------------------- STATISTICS ----------------------------------------------------

    quanti.set_index("ID", inplace=True)
    looking_metacell = metacell_comp.filter(regex=f"^metacell", axis=1)  # metacell (beginning of the sentence)
    percentage_imputation, percentage_signif2x2, summary_psm, summary_spc = st.counting(looking_metacell,
                                                                                        looking_2x2compare, psm, spc)
    imputation_value = st.find_M(looking_metacell, quanti)
    table_pca, scree, pca = st.my_pca(table_stats)


    # ----------------------------------------------- DATA VISUALIZATION -----------------------------------------------

    plt.ioff()  # doesn't show the graphics interface
    viz.fullheatmap(quanti, save_dir)
    viz.volcanos(looking_FC, looking_2x2compare, threshold_pvalue, threshold_FC, save_dir)
    viz.scree_plot(scree, save_dir)
    viz.plot_PCA(table_pca, scree, save_dir, pca, table_stats)

    # ------------------------------- EXPORTATION IN NEW EXCEL FILE WITH DIFFERENT SHEET -------------------------------

    output(new_excel_file, full_table, percentage_imputation, imputation_value, summary_psm,
           summary_spc, percentage_signif2x2, threshold_FC, threshold_pvalue)
    styling(new_excel_file, looking_metacell, full_table)

    # ---------------------------------------------- END OF THE SCRIPT -------------------------------------------------

    print("The new excel file: ", new_excel_file.replace('/', '\\'), "has been updated successfully")

    # Open the Windows explorer and pass the system console to pause
    subprocess.Popen(['explorer', save_dir.replace('/', '\\')], shell=True)
    os.system("pause")

except Exception as e:
    print('\033[31m' + "ERROR:" + str(e) + '\033[0m')
    input("Press enter to close the window")