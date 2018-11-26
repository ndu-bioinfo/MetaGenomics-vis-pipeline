Functions for creating diversity plots using plotly

import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.spatial.distance import pdist,squareform
from skbio.stats.ordination import pcoa
from scipy.interpolate import griddata
from skbio import DistanceMatrix
from skbio import diversity
import seaborn as sns
from Bio import SeqIO,Phylo
from IPython.display import IFrame


import plotly.offline as py
import plotly.graph_objs as go
import ecopy as ep
from plotly import tools

import warnings
warnings.filterwarnings("ignore")

from Diversity_plots_0_3 import Diversity_bar_plot



1. Beta_Diversity_Contour_plot

2. Diversity_bar_plot - Generate AxB bar plots with selected alpha or beta diversity results and group the data using selected factor A and B; the significantly different groups are determined using  marked with different alphabets   

Method_choice =  'heip_e'
image_size = [800,500]

fig,layout = Diversity_bar_plot(df_alpha_meta = df_bar_alpha_16S,AxB_table = AxB_table,div_a_col=Method_choice, diversity_type='alpha', cols=1,Color_pattern = 'Spectral');
fig['layout'].update(height=image_size[0], width=image_size[1],showlegend=False);

py.plot(fig,filename='16S_'+Method_choice+'.html', auto_open=False,image='png',image_filename='16S_'+Method_choice,image_height=image_size[0],image_width=image_size[1])
IFrame(src='16S_'+Method_choice+'.html',height=image_size[0]+50, width=image_size[1]+50)

