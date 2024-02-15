# -*- coding: utf-8 -*-

def progress_bar(progress, total, task='task'):
    """
        Generate a progress bar in the console
        yellow during the progression
        green when it's completed
    :param progress: current number
    :param total: representing the total number of tasks to be completed
    :param task: Name of the task
    :return: the progress bar in the console
    """
    percent = 100 * (progress/float(total))
    bar = '█'* int(percent) + '-' * (100-int(percent))
    color = '\033[33m'  # ANSI escape code for changing text color in the terminal
    print(color + f"\r|{progress}/{total}|{bar}| {percent:.0f}%", end="\r")
    if progress == total:
        color = '\033[92m'
        print(color + f"\r|{progress}/{total}|{bar}| {percent:.0f}%", end="\n")
        print(f"{task} Complete"+'\033[0m')

# Progress for each step
print("import package")

progress_bar(0, 12)
import os
progress_bar(1, 12)
import pandas as pd
progress_bar(2, 12)
import numpy as np
progress_bar(3, 12)
import openpyxl
progress_bar(4, 12)
import subprocess
progress_bar(5, 12)
from statsmodels.stats.multitest import multipletests
progress_bar(6, 12)
import seaborn as sns
progress_bar(7, 12)
from sklearn.decomposition import PCA
progress_bar(8, 12)
from sklearn.preprocessing import StandardScaler
progress_bar(9, 12)
import matplotlib.pyplot as plt
progress_bar(10, 12)
import plotly.express as px
progress_bar(11, 12)
import plotly.graph_objs as go
progress_bar( 12,12,'Initialisation')

