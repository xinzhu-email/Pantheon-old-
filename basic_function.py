from bokeh.io import curdoc
from bokeh.models import FileInput
from src.core import CreateTool
import anndata

upload_button = FileInput()
curdoc().add_root(upload_button)
def eb(attr,old,new):
    print(upload_button.filename)
    try:
        adata = anndata.read_csv(upload_button.filename)
        print('csv')
    except:
        adata = anndata.read(upload_button.filename)
        print('h5ad')
    
    mainplot, panel1 = CreateTool(adata).base_tool()
    #print(mainplot)
    hl_figure, panel2 = CreateTool(adata).highlight_gene(mainplot)
    tab = CreateTool(adata).multi_panel([mainplot,hl_figure],[panel1,panel2], ['Main View', 'Highlight Gene'], update_view=True)
    curdoc().add_root(tab)
upload_button.on_change('filename',eb)