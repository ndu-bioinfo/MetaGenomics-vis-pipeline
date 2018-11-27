## Functions for creating diversity plots using plotly

#### alpha or beta diversity bar plot
```python
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

from function_diversity_analysis import Diversity_bar_plot

# import alpha diversity dataframe: df_bar_alpha
```
| index	| sample_code	| heip_e	| season	| year	| region|
| --- | --- | --- | --- | --- | --- |
| 201501_086.7_110.0_86	| ND0184	| 0.113271	| Winter	| 2015	| Off|
| 201402_090.0_090.0_10	| ND0012	| 0.103633	| Winter	| 2014	| Off|
| 201607_093.3_110.0_17	| ND0398	| 0.065743	| Summer	| 2016	| Off|
| 201604_093.3_026.7_20	| ND0355	| 0.179908	| Spring	| 2016	| SCB|
| 201407_086.7_045.0_53	| ND0094	| 0.133548	| Summer	| 2014	| Up|
| ...|

```python
# import beta diversity matrix: dict_bar_beta['unifrac'] - dict of matrics calculated using different methods
dict_bar_beta['unifrac']
```
|     | 201501_086.7_110.0_86	| 201402_090.0_090.0_10 |	201607_093.3_110.0_17 |	201604_093.3_026.7_20 |	201407_086.7_045.0_53|
| --- | --- | --- | --- | --- | --- |
| 201501_086.7_110.0_86	| 0.000000	| 0.518750	| 0.624583	| 0.781583	| 0.651833| 
| 201402_090.0_090.0_10	| 0.518750	| 0.000000	| 0.480417	| 0.641250	| 0.536583| 
| 201607_093.3_110.0_17	| 0.624583	| 0.480417	| 0.000000	| 0.796750	| 0.761083| 
| 201604_093.3_026.7_20	| 0.781583	| 0.641250	| 0.796750	| 0.000000	| 0.566000| 
| 201407_086.7_045.0_53	| 0.651833	| 0.536583	| 0.761083	| 0.566000	| 0.000000| 
```python
# define A x B factors 
AxB_table = [['season','year'],['year','season'],['year','region'],['region','season']]

# 1. Diversity_bar_plot - Generate AxB bar plots with selected alpha or beta diversity results and group the data using selected factor A and B; the significantly different groups are determined using  marked with different alphabets   

# Example for alpha diversity
Method_choice =  'heip_e'
image_size = [800,500]

fig,layout = Diversity_bar_plot(df_alpha_meta = df_bar_alpha,AxB_table = AxB_table,div_a_col=Method_choice, diversity_type='alpha', cols=1,Color_pattern = 'Spectral');
fig['layout'].update(height=image_size[0], width=image_size[1],showlegend=False);

py.plot(fig,filename='16S_'+Method_choice+'.html', auto_open=False,image='png',image_filename='16S_'+Method_choice,image_height=image_size[0],image_width=image_size[1])
IFrame(src='16S_'+Method_choice+'.html',height=image_size[0]+50, width=image_size[1]+50)
export_values(fig).to_csv('16S_'+Method_choice+'.csv')

# Example for beta diversity

'''
