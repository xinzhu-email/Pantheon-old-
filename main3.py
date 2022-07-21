from cProfile import label
from bokeh.io import curdoc
from bokeh.models import FileInput, Button
from src.source import CreateTool, connect_figure, class2json
from src.transform import data_trans
from bokeh.layouts import row, column
import anndata
import scanpy as sc
import json

upload_button = FileInput()
curdoc().add_root(upload_button)
def eb(attr,old,new):
    try:
        adata = anndata.read_csv(upload_button.filename)
        print('csv')
    except:
        adata = anndata.read(upload_button.filename)
        print('h5ad')
    print(upload_button.filename)
    mainplot, panel1 = CreateTool(adata).base_tool()
    new_button = Button(label='Change Color')
    new_button.on_click(lambda: change_color(mainplot))
    panel1 = column(panel1, new_button)
    #new_button = connect_figure(mainplot).add_layout(callback_function=change_color, button=True, label='Change Color', parameter=(mainplot,panel1))
    #panel1 = column(panel1, new_button, row(*filter_cells(mainplot)))
    hl_figure, panel2 = CreateTool(adata).highlight_gene(mainplot)
    tab = CreateTool(adata).multi_panel([mainplot,hl_figure],[panel1,panel2], ['Main View', 'Highlight Gene'], update_view=True)
    curdoc().add_root(tab)
upload_button.on_change('filename',eb)



### Examples to add buttons ###
def change_color(Figure): 
    Figure.show_checked()  
    trans = class2json(Figure)
    to_json = trans.transfer()
    #print(to_json)
    data_dict = json.loads(to_json)
    color = data_dict['selected_color']
    selected_class = data_dict['checked_class']
    group = data_dict['selected_group']
    data = Figure.adata
    data.uns['category_dict'][group]['color'][[i for i in selected_class]] = color   
    indices = data_dict['selected_indices']
    data_color = data_dict['data_color']
    for i in indices:
        data_color[i] = color
    # Save change of data into the Figure
    data_dict['data_color'] = data_color
    trans.renew_data(Figure,data_dict)
    #add.set_change(data_color=data_color)
    Figure.text_color()

def filter_cells_callback(Figure, input):
    adata = Figure.adata.copy()
    sc.pp.filter_cells(adata, min_genes=int(input.value))
    print(input.value)
    print(len(adata.obs))
    cells = adata.obs.index
    change = connect_figure(Figure)
    change.set_change(filter_cells=cells)


def filter_cells(Figure):
    add = connect_figure(Figure)
    input = add.add_layout(input=True, label='Min genes expressed in cell:')
    comfirm_button = add.add_layout(button=True, label='Filter the Cells', callback_function=filter_cells_callback, parameter=(Figure,input))
    return input, comfirm_button

