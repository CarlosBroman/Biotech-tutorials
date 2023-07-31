# Import the required libraries
import pandas as pd # for working with data frames
import numpy as np # for numerical operations
import panel as pn # for creating interactive plots and dashboards
import hvplot.pandas # for creating visualizations with HoloViews
import holoviews as hv # for creating interactive visualizations

import random # for generating random colors

pn.extension('tabulator') # enable the tabulator extension for Panel
hv.extension('bokeh') # enable the Bokeh extension for HoloViews

# Define our functions
def get_unique_values(df, column_name): # this function takes a DataFrame and a column name as input and returns a list of unique values in the specified column.
    return df[column_name].unique().tolist()

def get_color_dict(df, column_name, n_colors): # this functions creates a dictionary to map a certain number of colors to the unique values of the column of a data frame
    therapies = df[column_name].unique()
    colors = []
    used_colors = set()
    while len(colors) < len(therapies):
        color = "#" + "".join([random.choice("0123456789ABCDEF") for j in range(6)])
        if color not in used_colors:
            colors.append(color)
            used_colors.add(color)
    color_dict = {}
    for i in range(len(therapies)):
        color_dict[therapies[i]] = colors[i]
    return color_dict

# Read the csv we created previously and store it as a data frame
df = pd.read_csv("IB_therapies.csv")
df # check if the data frame is correct

idf = df.interactive() # convert it into an interactive data frame

# Create the RadioBoxGroup for y selection
y_columns = list(df.columns[5:9]) # select the columns for the y-axis using list slicing
yaxis_CT = pn.widgets.RadioBoxGroup(
    name='Y axis',
    options=y_columns,
    value=y_columns[0]
)

# Create the RadioBoxGroup for x selection
x_columns = ['Target', 'Focus Area'] # select the colums for the x-axis
xaxis_CT = pn.widgets.RadioBoxGroup(
    name='X axis',
    options=x_columns,
    value=x_columns[0]
)

target = get_unique_values(df, 'Target') # create a list of unique therapy targets
FA = get_unique_values(df, 'Focus Area') # create a list of unique therapy Focus Area


# Create the colors for out plot
num_colors = len(df['Therapy'].unique()) # how many colors do we need (equals the number of therapies)
color_dict = get_color_dict(df, 'Therapy', num_colors) # create the dictionary with the colors assigned to the therapies
print(color_dict)


def plot(x, y, color_dict): 
    if x == 'Target':
        data = df[df.Target.isin(target)]
        x_label = 'Target'
    elif x == 'Focus Area':
        data = df[df['Focus Area'].isin(FA)]
        x_label = 'Focus Area'
    else:
        raise ValueError(f'Invalid value for x axis: {x}')
        
    CT_bar_pipeline = (
        data.groupby([x, 'Therapy'])[y].sum()
        .reset_index()
    )
    
    # Extract the unique therapies in the plot
    therapies = CT_bar_pipeline['Therapy'].unique().tolist()
    
    # Filter the color dictionary to include only the therapies in the plot
    plot_color_dict = {k: v for k, v in color_dict.items() if k in therapies}
    
    # Sort the Therapy column alphabetically
    CT_bar_pipeline = CT_bar_pipeline.sort_values(by=['Therapy'])
    
    return CT_bar_pipeline.hvplot(kind='bar', x=x, y=y, by='Therapy', 
                                  stacked=True, width=1200, height=600,
                                  hover_cols=['Therapy', y],
                                  xlabel=x_label, ylabel='Clinical Trials',
                                  color='Therapy', cmap=plot_color_dict,
                                  legend='right', fontsize={'title': 14, 'labels': 12})
# Create and show the plot
CT_plot = pn.panel(pn.bind(plot, x=xaxis_CT, y=yaxis_CT, color_dict=color_dict))

pn.Row(CT_plot, pn.Column(xaxis_CT, yaxis_CT)).show()