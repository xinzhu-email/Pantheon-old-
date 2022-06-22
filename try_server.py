from cProfile import label
from bokeh.io import show
from bokeh.models import Slider, ColumnDataSource, CDSView, IndexFilter, CustomJS, Circle, Div, Panel, Tabs, CheckboxGroup
from bokeh.models.widgets import Select, Button, ColorPicker,TextInput, DataTable, MultiSelect
from bokeh.events import ButtonClick
from bokeh.transform import linear_cmap
from bokeh.palettes import d3
from bokeh.layouts import row, column
from bokeh.io import curdoc
from bokeh.layouts import row
from bokeh.plotting import figure
from numpy.random import random
import pandas
import numpy as np
import anndata
from sklearn import metrics
import pandas as pd
import scipy.sparse as ss

# Loading data
#data_path = "E:/项目/图形化界面/ADT_gater-master/"
data_path = 'CD4_memory_Naive.h5ad'
adata = anndata.read(data_path)
#data_array = []
#data_array.append(pandas.read_csv(data_path+'ADT.csv'))
#data_df = pandas.concat(data_array,axis=1,join='inner')
#print(data)
data_df = adata.to_df()
generic_columns = data_df.columns.values.tolist()[1:]
#print(generic_columns)
# Initialize color attribute
data_df['color'] = pandas.Series(np.full(data_df.shape[0], 0).astype(str), index=data_df.index)
# Initialize highly variable gene
data_df['hl_gene'] = pandas.Series(np.full(data_df.shape[0], 0), index=data_df.index)

color_list = d3['Category20c'][20]
# Function
def tag_func(selector, effector, attr, plot):
    axis = getattr(plot, attr + "axis")
    axis.axis_label = selector.value
    setattr(effector, attr, selector.value)

def select_color_func(attr, old, new):
    global select_color, cur_color, Figure, source

    cur_color = int(select_color.value)
    Figure.refresh()

    # Just to trigger plot update 
    source.data["color"] = source.data["color"]


def color_func():
    global source, cur_color, data_df, Figure

    c_list = source.data["color"]

    for i in source.selected.indices:
        c_list[i] =  cur_color
    print('chosen color:',c_list)
    source.data["color"] = c_list
    data_df["color"] = c_list



def correct_func():
    global source
    source.selected.indices = []

def showall_func():
    global view
    view.filters = []

def selection_func():
    global view
    view.filters = [IndexFilter(source.selected.indices)]


def remove_func():
    global view, source
    remove_set = set(source.selected.indices)
    if len(view.filters) == 0:
        view.filters = [IndexFilter(range(source.data.shape[0]))]
    remain_indices = [x for x in view.filters[0].indices if x not in source.selected.indices]
    view.filters = [IndexFilter(remain_indices)]




TOOLTIPS = [
        ("index", "@index"),
        ("(x,y)", "($x, $y)"),
        ("color", "@color"),
]

class FlowPlot:
    def __init__(self, opts, source, view, columns, color_map, title = "", x_init_idx = 0, y_init_idx = 0, allow_select = True, select_color_change = True, legend = None):
        self.opts = opts
        self.source = source
        self.view = view
        self.columns = columns
        self.color_map = color_map
        self.p = figure(width=500, height=500, tools="pan,lasso_select,box_select,tap,wheel_zoom,save,hover",title="Surface Marker Gating Panel", tooltips=TOOLTIPS)
        #self.p.output_backend = "svg"
        print("backend is ", self.p.output_backend)
        self.p.xaxis.axis_label = self.columns[x_init_idx]
        self.p.yaxis.axis_label = self.columns[y_init_idx]
        self.r = self.p.circle(self.columns[x_init_idx], self.columns[y_init_idx],  source=self.source, view=self.view, fill_alpha=1,fill_color=color_map, line_color=None )
        self.p.legend.click_policy="hide"
        self.s_x = Select(title="x:", value=self.columns[x_init_idx], options=self.columns)
        self.s_y = Select(title="y:", value=self.columns[y_init_idx], options=self.columns)
        # Attach reaction
        self.s_x.on_change("value", lambda attr, old, new: tag_func(self.s_x, self.r.glyph, 'x', self.p) )
        self.s_y.on_change("value", lambda attr, old, new: tag_func(self.s_y, self.r.glyph, 'y', self.p) )
        # Set default fill color
        if select_color_change:
            self.r.selection_glyph = Circle(fill_alpha=1, fill_color=color_list[0], line_color=None)
        self.allow_select = allow_select

    def refresh(self):
        global cur_color
        print('color list of data: ',self.r.selection_glyph.fill_color)
        self.r.selection_glyph.fill_color = color_list[cur_color]
    
    def show_choosed(self,cate_name,class_name):
        
        #self.r.data_source.selected.indices = adata.uns[cate_name][adata.uns[cate_name]['class_name']==class_name].tolist()
        self.r.selection_glyph.line_color = 'black'


# Initialize data source
opts = dict(plot_width=500, plot_height=500, min_border=0, tools="pan,lasso_select,box_select,wheel_zoom,save")
source = ColumnDataSource(data=data_df)
view = CDSView(source=source, filters=[IndexFilter([i for i in range(data_df.shape[0])])])
color_map = linear_cmap('color',color_list, low=0, high=20)
Figure = FlowPlot(opts, source, view, generic_columns, color_map, "Surface Marker Gating Panel", legend='color')

# Change the color of selected parts
cur_color = color_list[0]
select_color = Select(title="Select color:", options=list(str(i) for i in range(20)), value=str(0))
select_color.on_change("value", select_color_func)

"""def select_clr_func(attr, old, new):
    global select_color, cur_color, Figure, source
    cur_color = select_color.color
    Figure.refresh()#in this function: change "self.r.selection_glyph.fill_color = color_list[cur_color]" to "self.r.selection_glyph.fill_color=cur_color"
    source.data["color"] = source.data["color"]
select_color = ColorPicker(title="Select color:", color=color_list[1], css_classes=color_list)
select_color.on_change("color", select_clr_func)"""

color_button = Button(label="Color")
color_button.on_click(color_func)

correct_button = Button(label="Correct Plot")
correct_button.on_click(correct_func)


gate_button = Button(label="Gate")
gate_button.on_click(selection_func)

remove_button = Button(label="Remove")
remove_button.on_click(remove_func)

showall_button = Button(label="Show All")
showall_button.on_click(showall_func)


### Category Functions ###
# New Category
category_options = {}
def new_category():
    global adata, category_options, cat_opt, category_options
    adata.uns[name.value] =  pandas.DataFrame(index=range(data_df.shape[0]), columns=['class_name','color'],dtype=object)
    category_options[name.value] = []
    cat_opt.options = list(category_options.keys())
    cat_opt.value = cat_opt.options[0]
    print(adata.uns[name.value])
    curdoc().remove_root(class_func)
    curdoc().remove_root(delclass_panels)
    curdoc().add_root(class_func)
    curdoc().remove_root(newclass_panels)
    curdoc().add_root(newclass_panels)

# Edit category
def edit_category():
    
    global cat_opt, adata, category_options
    old_name = cat_opt.value
    new_name = name.value
    adata.uns[new_name] = adata.uns.pop(old_name)
    category_options[new_name] = category_options.pop(old_name)
    cat_opt.options = list(category_options.keys())
    cat_opt.value = cat_opt.options[0]
    print(adata.uns[new_name])

# Delete category
def del_category():
    global cat_opt, category_options, adata
    del category_options[cat_opt.value]
    del adata.uns[cat_opt.value]   
    cat_opt.options = list(category_options.keys())
    if len(cat_opt.options) == 0:
        cat_opt.value = "No Category"
    else:
        cat_opt.value = cat_opt.options[0]
    
    

def catfunc_cb(attr,old,new):
    global name, save_button
    if cat_func.value == 'Create New Category':
        
        save_button.label = 'Confirm to Create new Category!'
    elif cat_func.value == 'Change Nme of Category':       
        edit_category()
    elif cat_func.value == 'Delete Selected Category':
        name.value = ''
        save_button.label = 'Confirm to Delete! (Please check the selected CATEGORY bellow!!)'

def save_cb():
    global save_button, cat_opt
    if cat_func.value == 'Create New Category':
        new_category()
    elif cat_func.value == 'Change Name of Category':
        edit_category()
    elif cat_func.value == 'Delete Selected Category':
        del_category()

name = TextInput(title='Input Category Name', value='')
name.js_on_change("value", CustomJS(code="""
    console.log('text_input: value=' + this.value, this.toString())
"""))
cat_func = Select(title="Select Function of Category:",options=['Create New Category','Change Name of Category','Delete Selected Category'],value='Create New Category')
cat_func.on_change("value",catfunc_cb)
save_button = Button(label="Save New Category")
save_button.on_click(save_cb)

# Choose Category
def choose_cat(attr,old,new):
    cat_choice = cat_opt.value
    class_select.title = 'Choose class:'
    if cat_choice == 'No Category':
        class_select.options = ['No Category']
        class_select.value = 'No Category'
    else:
        class_select.options = category_options[cat_choice]
    print('------------',class_select.value)

cat_opt = Select(title='Select Category', options=['No Category'], value='No Category')
cat_opt.on_change("value", choose_cat)

new_catboard = column(cat_func,row(cat_opt,name),save_button)

def class_change():
    choosed_list = class_checkbox.active
    


cls_label = []
class_checkbox = CheckboxGroup(labels=cls_label,active=[0,1])
class_checkbox.js_on_click(CustomJS(code="""
    console.log('checkbox_group: active=' + this.active, this.toString())
"""))
# New Class
def add_entry():
    global cls_label
    category_options[cat_opt.value] = list(category_options.get(cat_opt.value,[]) + [text_input.value])
    #class_options.append(text_input.value)
    class_select.options = category_options[cat_opt.value]
    class_select.value = text_input.value
    
    save_class(cat_opt.value, text_input.value)
    text_input.value = ''
    print(f'\n{class_select.options}\n')

def save_class(category_name, class_name):
    global adata, class_checkbox, cls_label
    class_list = adata.uns[category_name]['class_name']
    color_l = adata.uns[category_name]['color']
    for i in source.selected.indices:
        class_list[i] =  class_name
        color_l[i] = cur_color
    adata.uns[category_name]['class_name'] = class_list
    adata.uns[category_name]['color'] = color_l
    cls_label = list(category_options[category_name])
    for i in range(len(cls_label)):
        ind = adata.uns[category_name][adata.uns[category_name]['class_name'] == category_options[category_name][i]]
        cls_label[i] = cls_label[i] + ': color=' + str(ind['color'].values[0]) + ', cell_nums='+ str(len(ind))
    class_checkbox.labels = cls_label
    curdoc().remove_root(class_checkbox)
    curdoc().add_root(class_checkbox)
    print(adata.uns[category_name])

# Delete Class
def del_class():
    global adata, class_select, category_options
    cate_name = cat_opt.value
    class_name = class_select.value
    print('class_name',class_name)
    adata.uns[cate_name][adata.uns[cate_name]['class_name']==class_name] = np.NAN

    del_list = list(category_options[cate_name])
    del_list.remove(class_name)
    del category_options[cate_name]
    category_options[cate_name] = del_list

    class_select.options = category_options[cate_name]
    print(adata.uns[cate_name])
    

    
# Class Function
def clsfunc_cb(attr,old,new):
    global class_button
    if class_func.value == 'Create New Class':
        curdoc().remove_root(delclass_panels)
        curdoc().remove_root(newclass_panels)
        curdoc().add_root(newclass_panels)
        class_button.label = 'Confirm to Create New Class!'
    elif class_func.value == 'Edit Class':
        curdoc().remove_root(delclass_panels)
        curdoc().remove_root(newclass_panels)
        curdoc().add_root(delclass_panels)
        class_button.label = 'Save Selected Dots into the Class'
    elif class_func.value == 'Delete Selected Class':
        curdoc().remove_root(newclass_panels)
        curdoc().remove_root(delclass_panels)
        curdoc().add_root(delclass_panels)
        class_button.label =  'Confirm to Delete! (Please check the selected CLASS bellow!!)'
    elif class_func.value == 'Change Category Name':
        class_button.label = 'Comfirm to Change Category Name!'
def save_cls_cb(event):
    global save_button, cat_opt, Figure
    if class_func.value == 'Create New Class':
        add_entry()
    elif class_func.value == 'Edit Class':
        Figure.show_choosed(cat_opt.value,class_select.value)
        save_class(cat_opt.value,class_select.value)
    elif class_func.value == 'Delete Selected Class':
        del_class()
    elif class_func.value == 'Change Category Name':
        edit_category()
class_func = Select(title='Select Class Function',options=['Create New Class','Edit Class','Delete Selected Class','Change Category Name'],value='Create New Class')
class_func.on_change("value",clsfunc_cb)
# Site selector
class_select = Select(title = 'Choose Class: ', options = ['Choose Class'], value='Choose Class')
text_input = TextInput(value = '', title = 'Enter Class Name:')
class_button = Button(label = 'Add New Class')
# Listen for events
class_button.on_event(ButtonClick,save_cls_cb)

newclass_panels = column(text_input,class_button)
delclass_panels = column(class_select,class_button)



# Show class
def show_category(attr, old, new):
    global class_select
    #choice = class_select.value
    d.text = 'Now Category: '+ cat_opt.value + '. Now Class: ' + class_select.value

d = Div(text='Now class: not specified')
class_select.on_change("value",show_category)



def function_callback(attr, old, new):
    global layout, adata
    func_list = function_select.value

    for i in range(len(func_list)):
        if func_list[i] == 'Create Category':           
            curdoc().add_root(new_catboard)
            
        


function_options = ['Create Category','Clustering','Cell Type Annotation']
function_select = MultiSelect(title='Function Select',options=function_options)
function_select.on_change("value",function_callback)

def addpanel(event):
    global tab_list,figure_panel
    layout2 = figure_panel
    paneltab2 = Panel(child=layout2, title='Panel2')
    print('==========tab_list',tab_list)
    tab_list.append(paneltab2)
    print('------tab_list',tab_list)
    tabs.update(tabs=tab_list)

add_panel = Button(label='Add Panel')
add_panel.on_event(ButtonClick,addpanel)

control_panel = column(Figure.s_x, Figure.s_y,select_color, row(color_button, correct_button), gate_button, remove_button, showall_button,new_catboard,add_panel)
figure_panel = column(Figure.p)
layout = row(figure_panel,control_panel)

panle1 = Panel(child=layout,title='Original Panel')
tab_list = [panle1]

tabs = Tabs(tabs=tab_list)




curdoc().add_root(tabs)
