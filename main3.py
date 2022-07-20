from bokeh.io import curdoc
from bokeh.models import FileInput
from src.source import CreateTool, connect_figure
from bokeh.layouts import row, column
import anndata
import scanpy as sc

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
    new_button = connect_figure(mainplot).add_layout(callback_function=change_color, button=True, label='Change Color', parameter=(mainplot,panel1))
    panel1 = column(panel1, new_button, row(*filter_cells(mainplot)))
    hl_figure, panel2 = CreateTool(adata).highlight_gene(mainplot)
    tab = CreateTool(adata).multi_panel([mainplot,hl_figure],[panel1,panel2], ['Main View', 'Highlight Gene'], update_view=True)
    curdoc().add_root(tab)
upload_button.on_change('filename',eb)



### Examples to add buttons ###
def change_color(Figure,panel):   
    # The needed data from Figure
    add = connect_figure(Figure)
    color = add.get(selected_color=True)
    selected_class = add.get(checked_class=True)
    group = add.get(selected_group=True)
    data = add.get(data=True)
    data.uns['category_dict'][group]['color'][[i for i in selected_class]] = color
    Figure.show_checked()
    indices = add.get(selected_indices=True)
    data_color = add.get(data_color=True)
    for i in indices:
        data_color[i] = color
    # Save change of data into the Figure
    add.set_change(data_color=data_color)
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

