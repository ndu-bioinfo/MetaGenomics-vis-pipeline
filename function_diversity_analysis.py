########## Diversity plots - Functions for making diversity plots using plotly ##########

##########  Niu Du dniu at jcvi.org 08/31/2018 ##########

import pandas as pd
import numpy as np
import plotly.graph_objs as go
from scipy.interpolate import griddata
from plotly import tools
import colorlover as cl
from scipy.stats import f_oneway

def Hawk_smash(List):
    '''Flattern lists within a list '''
    return [item for sublist in List for item in sublist]

def plot_contour(metadata,factor,resolution = 50,contour_method='linear'):
    resolution = str(resolution)+'j'
    x = metadata['x'].values
    y = metadata['y'].values
    z = metadata[factor].values
    X,Y = np.mgrid[min(x):max(x):complex(resolution), min(y):max(y):complex(resolution)]
    points = [[a,b] for a,b in zip(x,y)]
    Z = griddata(points, z, (X, Y), method=contour_method)
    return X,Y,Z


def group_data_ANOVA(Subset,plot_list):
    '''Group data by alphabetic names using oneway ANOVA'''
    ANOVA_m = np.zeros((len(Subset),len(Subset)))
    for u,U in enumerate(plot_list):
        for v,V in enumerate(plot_list):
            ANOVA_m[u,v] = f_oneway(Subset[U],Subset[V])[1]
    
    df_ANOVA_m = pd.DataFrame(ANOVA_m,index = plot_list,columns=plot_list).fillna(1)
    df_ANOVA_m = df_ANOVA_m[df_ANOVA_m<0.05]
    df_ANOVA_m = df_ANOVA_m.replace(to_replace=0,value=1)
    df_ANOVA_m = df_ANOVA_m.fillna(1)
    df_ANOVA_test = df_ANOVA_m.copy() ### Test here
    i = 0
    
    for col in plot_list:
        if len(df_ANOVA_m.loc[df_ANOVA_m[df_ANOVA_m[col]==1].index])>0:
            df_ANOVA_m.loc[df_ANOVA_m[df_ANOVA_m[col]==1].index] = chr(65+i)
            i+=1
    return dict(df_ANOVA_m[df_ANOVA_m.columns[0]])

def export_values(fig):
    sample_ID = []
    y_values = []
    y_error = []
    label = []
    for g in fig['data']:
        sample_ID.append([A.split(' (')[0][3:]+'-'+B for A,B in zip(g['text'],g['x'])])
        y_values.append(g['y'])
        y_error.append(g['error_y']['array'])
        label.append([A.split('(')[-1][:1] for A in g['text']])
    sample_ID = Hawk_smash(sample_ID)
    y_values = Hawk_smash(y_values)
    y_error = Hawk_smash(y_error)
    label = Hawk_smash(label)
    return pd.DataFrame({'y':y_values,'y_error':y_error,'label':label},index = sample_ID)

def Beta_Diversity_Contour_plot(df_alpha_meta,df_PCoA_Matrix,a_div_col,group_factors,contour_factors,node_scale=0.15,contour_resolution=50,opacity_node=0.9,opacity_contour=0.6,contour_cmap='heatmap',contour_method = 'linear'):
    
    
    
    Site_size = dict(zip(df_alpha_meta.index,df_alpha_meta[a_div_col]))
    Scatter_data = []
    Scatter_data_length = []
    for Group_factor in group_factors:
        
        if Group_factor == 'season':
            List_Group_factor = ['Winter','Spring','Summer','Fall']
        elif Group_factor == 'region':
            List_Group_factor = ['Off','CC','SCB','Up']
        else:
            List_Group_factor = list(np.sort(list(set(df_alpha_meta[Group_factor]))))
        
        Index = []
        trace = []
        for trace_id,trace_value in enumerate(List_Group_factor):
                Index.append(list(df_alpha_meta[df_alpha_meta[Group_factor]==trace_value].index))
                trace.append(go.Scatter(
                    x=df_PCoA_Matrix['PC1'].loc[Index[trace_id]],
                    y=df_PCoA_Matrix['PC2'].loc[Index[trace_id]],
                    name=trace_value,
                    text=list(df_alpha_meta[df_alpha_meta[Group_factor]==trace_value].index),
                    hoverinfo='text',
                    mode='markers',
                    marker=dict(
                        size=[Site_size[m]*node_scale for m in df_alpha_meta[df_alpha_meta[Group_factor]==trace_value].index],
                        opacity=opacity_node
                        )
                    )
                )
        Scatter_data= Scatter_data+trace
        Scatter_data_length.append(len(trace))
    
    Vis_table = []
    for i in np.arange(len(Scatter_data_length)):
        Vis_list = []
        for j,elements in enumerate(Scatter_data_length):
            if j == i:
                Vis_list.append([True]*elements)
            else:
                Vis_list.append([False]*elements)
        Vis_table.append([item for sublist in Vis_list for item in sublist])
    dict_Vis = dict(zip(group_factors,Vis_table))

    M = {}
    for Factor in contour_factors:
        M[Factor] = plot_contour(df_alpha_meta,Factor,contour_resolution,contour_method)

     # data list for plotting
 # data list for plotting
    data = Scatter_data + [go.Contour(
            z= M[Factor][2],
            x= M[Factor][0][:,0],
            y= M[Factor][1][0,:],
            colorscale='Jet',
            showlegend = False,
            connectgaps= True,
            line=dict(smoothing=0.85),
            showscale = True,
            opacity = opacity_contour,
            hoverinfo = 'skip',
            contours=dict(
                coloring=contour_cmap,
            )
        )]
    
    updatemenus = list([
    dict(active=0,
         buttons=list([dict(label = x,
                 method = 'restyle',
                 args = ['z', [M[x][2]]]) for x in contour_factors]),
        direction = 'down',
        pad = {'r': 10, 't': 10},
        showactive = True,
        x = -0.23,
        xanchor = 'left',
        y = 0.9,
        yanchor = 'top'
        ),
    dict(active=0,
         buttons=list([   
            dict(label = Group_factor,
                 method = 'update',
                 args = [ {'visible': dict_Vis[Group_factor]+[True]},{'title': str(Group_factor)}]) for Group_factor in group_factors
                             ])
             )
    ])

    layout = go.Layout(legend=dict(x=-0.15, y=0),
            xaxis=dict(
            title='PCoA1',
            titlefont=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
                )
            ),

            yaxis=dict(
            title='PCoA2',
            titlefont=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
                )
            ),
            updatemenus=updatemenus
            )
    fig = go.Figure(data=data,layout = layout)
    return fig,layout

def Diversity_bar_plot(df_alpha_meta,AxB_table,div_a_col = 'faith_pd',df_matrix = None,diversity_type = 'alpha',cols = 1,Color_pattern = 'Spectral'):
    
    dict_Bar_data = {}
    
    for factor_A, factor_B in AxB_table:
        Subsets = {}
        plot_list = []
        Index = []
        trace = []

        
        if factor_A == 'season':
            List_A = ['Winter','Spring','Summer','Fall']
        elif factor_A == 'region':
            List_A = ['Off','CC','SCB','Up']
        else:
            List_A = list(np.sort(list(set(df_alpha_meta[factor_A]))))
        
        if factor_B == 'season':
            List_B = ['Winter','Spring','Summer','Fall']
        elif factor_B == 'region':
            List_B = ['Off','CC','SCB','Up']
        else:
            List_B = list(np.sort(list(set(df_alpha_meta[factor_B]))))
        
        
        Colors = cl.scales[str(len(List_B))]['div'][Color_pattern]
        
        ########## Function here to group data by ANOVA results ############            
        for A in List_A:
            for B in List_B:
                plot_list.append(A+'-'+B)    
        
        for B_id,B in enumerate(List_B):
            
                Index.append(list(df_alpha_meta[df_alpha_meta[factor_B]==B].index))
                meta_subset = df_alpha_meta.loc[Index[B_id]]
                
                Dict_AxB = dict(zip(List_A,[meta_subset[meta_subset[factor_A] == y].index for y in List_A]))
        
                if diversity_type == 'beta':
                    for y in List_A:
                        Subsets[y+'-'+B] = np.mean(df_matrix[Dict_AxB[y]].loc[Dict_AxB[y]]).values
                        
                if diversity_type == 'alpha':
                    for y in List_A:
                        Subsets[y+'-'+B] = df_alpha_meta[div_a_col].loc[Dict_AxB[y]].values
        
        
        group_label = group_data_ANOVA(Subsets,plot_list)

        
        ######### Plot data ##########
        for B_id,B in enumerate(List_B):
                Index.append(list(df_alpha_meta[df_alpha_meta[factor_B]==B].index))
                meta_subset = df_alpha_meta.loc[Index[B_id]]
                
                Dict_AxB = dict(zip(List_A,[meta_subset[meta_subset[factor_A] == y].index for y in List_A]))
                Color = Colors[B_id]
                
                if diversity_type == 'beta':
                    trace.append(go.Bar( 
                        x=List_A,
                        # beta diversity
                        y = [np.mean(np.mean(df_matrix[Dict_AxB[y]].loc[Dict_AxB[y]])) for y in List_A],
                        error_y=dict(
                        type='data',
                        array=[np.std(np.mean(df_matrix[Dict_AxB[y]].loc[Dict_AxB[y]])) for y in List_A],
                        visible=True,opacity=0.6
                        ),
                        name = B,
                        text =['   '+B+' ('+group_label[y+'-'+B]+')' for y in List_A],
                        textposition = 'auto',
                        hoverinfo = 'none',
                        opacity=1,
                        marker=dict(line=dict(
                            width=1),
                            color = Color,
                            opacity=0.6
                            )
                        )
                    )
                        
                if diversity_type == 'alpha':
                    trace.append(go.Bar( 
                        x=List_A,
                        # beta diversity
                        y = [np.mean(df_alpha_meta[div_a_col].loc[Dict_AxB[y]]) for y in List_A],
                        error_y=dict(
                        type='data',
                        array=[np.std(df_alpha_meta[div_a_col].loc[Dict_AxB[y]]) for y in List_A],
                        visible=True,opacity=0.6
                        ),
                        name = B,
                        text =['   '+B+' ('+group_label[y+'-'+B]+')' for y in List_A],
                        textposition = 'auto',
                        hoverinfo = 'none',
                        opacity=1,
                        marker=dict(line=dict(
                            width=1),
                            color = Color,
                            opacity=0.6
                            ),
                        )
                    )
                    
        
        dict_Bar_data[factor_A+' by '+factor_B]=trace
        
        
    layout = go.Layout(
        yaxis=dict(
        titlefont=dict(
            family='Courier New, monospace',
            size=18,
            color='#7f7f7f'
            )
            ),
        barmode='group',
        showlegend=False)

    fig = tools.make_subplots(rows=round(len(AxB_table)/cols), cols=cols,subplot_titles=(list(dict_Bar_data.keys())))
    
    dict_index = 0
    for col in np.arange(cols):
        for row in np.arange(round(len(AxB_table)/cols)):
            [fig.append_trace(i,int(row+1),int(col+1)) for i in dict_Bar_data[list(dict_Bar_data.keys())[dict_index]]]
            dict_index+=1
            
    return fig,layout
