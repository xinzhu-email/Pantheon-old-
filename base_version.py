from cProfile import label
from bokeh.io import show
from bokeh.models import Slider, ColumnDataSource, CDSView, IndexFilter, CustomJS, Circle, Div, Panel, Tabs, CheckboxGroup
from bokeh.models.widgets import Select, Button, ColorPicker,TextInput, DataTable, MultiSelect
from bokeh.events import ButtonClick
from bokeh.palettes import d3
from bokeh.layouts import row, column
from bokeh.io import curdoc
from bokeh.layouts import row
from bokeh.plotting import figure
import pandas
import numpy as np
import anndata
import pandas as pd
import scipy.sparse as ss

# Loading data
#data_path = "E:/项目/图形化界面/ADT_gater-master/"
data_path = 'CD4_memory_Naive.h5ad'
adata = anndata.read(data_path)
#data_array = []
#data_array.append(pandas.read_csv(data_path+'ADT.csv'))
#data_df = pandas.concat(data_array,axis=1,join='inner')
data_df = adata.to_df()
generic_columns = data_df.columns.values.tolist()
#print(generic_columns)
# Initialize color attribute
data_df['color'] = pandas.Series(d3['Category20c'][20][0], index=data_df.index)
# Initialize highly variable gene
data_df['hl_gene'] = pandas.Series(np.full(data_df.shape[0], 0), index=data_df.index)
color_list = d3['Category20c'][20]

##############################
##### Function Definition ####
##############################
def tag_func(selector, effector, attr, plot):
    axis = getattr(plot, attr + "axis")
    axis.axis_label = selector.value
    setattr(effector, attr, selector.value)

# Callback of colorpicker(selection), update the selected dots with picked color
def select_color_func(attr, old, new):
    global select_color, cur_color, Figure, source
    cur_color = select_color.color
    Figure.refresh()
    # Just to trigger plot update 
    source.data["color"] = source.data["color"]

# Save the color of the dots
def color_func():
    global source, cur_color, data_df, Figure
    c_list = source.data["color"]
    for i in source.selected.indices:
        c_list[i] =  cur_color
    print('chosen color:',c_list)
    source.data["color"] = c_list
    data_df["color"] = c_list

# Show the saved color of dots
def correct_func():
    global source
    source.selected.indices = []

# Show all, gate, and remove function
def showall_func():
    global view
    view.filters = []

def selection_func():
    global view
    view.filters = [IndexFilter(source.selected.indices)]
    print(source.selected.indices)

def remove_func():
    global view, source
    remove_set = set(source.selected.indices)
    if len(view.filters) == 0:
        view.filters = [IndexFilter(range(source.data.shape[0]))]
    remain_indices = [x for x in view.filters[0].indices if x not in source.selected.indices]
    view.filters = [IndexFilter(remain_indices)]


### Category Functions ###

# New Category
def new_category():
    global adata, category_options, cat_opt, category_options
    adata.uns[name.value] =  pandas.DataFrame(index=range(data_df.shape[0]), columns=['class_name','color'],dtype=object)
    adata.uns[name.value]['color'] = color_list[0]
    category_options[name.value] = []
    cat_opt.options = list(category_options.keys())
    cat_opt.value = cat_opt.options[0]
    print(adata.uns[name.value])



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

# Choose Category
def choose_cat(attr,old,new):
    global source
    if class_checkbox.labels != []:
        cls_label = []
        cate = cat_opt.value
        for cls in list(category_options[cate]):
            l = adata.uns[cate][adata.uns[cate]['class_name'] == cls]
            s = cls + [cls + ': color=' + str(l['color'].values[0]) + ', cell_nums='+ str(l.shape[0])]
            cls_label = np.append(cls_label,s)
        class_checkbox.labels = cls_label
        source.data['color'] = adata.uns[cate]['color']

# Callback of function selection about category
def catfunc_cb(attr,old,new):
    global name, save_button
    curdoc().remove_root(cat_opt)
    curdoc().remove_root(name)
    curdoc().remove_root(save_button)
    if cat_func.value == 'Create New Category':       
        save_button.label = 'Confirm to Create new Category!'
        curdoc().add_root(name)
        curdoc().remove_root(class_func)
        curdoc().remove_root(input_t)
        curdoc().remove_root(class_button)        
    elif cat_func.value == 'Change Nme of Category':  
        save_button.label = 'Confirm to Change Category Name!'     
        curdoc().add_root(cat_opt)
        curdoc().add_root(name)
        curdoc().add_root(save_button)
    elif cat_func.value == 'Delete Selected Category':
        save_button.label = 'Confirm to Delete! (Please check the selected CATEGORY !!)'
        curdoc().add_root(cat_opt)
        curdoc().add_root(save_button)

# Callback of save button
def save_cb():
    global save_button, cat_opt
    curdoc().remove_root(cat_opt)
    curdoc().remove_root(name)
    curdoc().remove_root(save_button)
    curdoc().remove_root(show_color_button)
    curdoc().add_root(show_color_button)
    if cat_func.value == 'Create New Category':
        new_category()
        curdoc().add_root(cat_opt)
        curdoc().remove_root(class_func)
        curdoc().add_root(class_func)
        curdoc().add_root(input_t)
        curdoc().add_root(class_button)
    elif cat_func.value == 'Change Name of Category':
        edit_category()
        curdoc().add_root(cat_opt)
    elif cat_func.value == 'Delete Selected Category':
        del_category()
        curdoc().add_root(cat_opt)




# New Class
def add_entry():
    global cls_label
    category_options[cat_opt.value] = list(category_options.get(cat_opt.value,[]) + [input_t.value])
    save_class(cat_opt.value, input_t.value)
    input_t.value = ''

# Merge checked classes
def merge_class(toclass,color):
    global category_options,class_checkbox,class_select
    checked_class = []
    for i in class_checkbox.active:
        checked_class = np.append(checked_class,category_options[cat_opt.value][i])
    l = adata.uns[cat_opt.value]['class_name'].isin(checked_class)
    adata.uns[cat_opt.value]['color'][l] = color
    adata.uns[cat_opt.value]['class_name'][l] = toclass
    del_list = list(category_options[cat_opt.value])
    del_list2 = class_checkbox.labels
    print('category_option:',category_options[cat_opt.value])
    print('checked class:',checked_class)
    for i in range(len(checked_class)):
        del_list.remove(checked_class[i])
        del del_list2[class_checkbox.active[i]-i]
    del category_options[cat_opt.value]    
    category_options[cat_opt.value] = del_list + [toclass]
    #class_select.options = category_options[cat_opt.value]
    del_list2 = del_list2 + [toclass+ ': color=' + str(color) + ', cell_nums='+ str(len(adata.uns[cat_opt.value][l]))]
    class_checkbox.labels = del_list2

# Save change of classes
def save_class(category_name, class_name):
    global adata, class_checkbox, cls_label
    class_list = adata.uns[category_name]['class_name']
    color_l = adata.uns[category_name]['color']
    num = 0
    for i in source.selected.indices:
        class_list[i] =  class_name
        color_l[i] = cur_color
        num = num + 1
    adata.uns[category_name]['class_name'] = class_list
    adata.uns[category_name]['color'] = color_l
    cls_label = class_checkbox.labels
    cls_label = cls_label + [class_name + ': color=' + str(cur_color) + ', cell_nums='+ str(num)]
    class_checkbox.labels = cls_label
    curdoc().remove_root(class_checkbox)
    curdoc().add_root(class_checkbox)
    print(adata.uns[category_name])

# Delete Class
def del_class():
    global adata, class_checkbox, category_options
    checked_class = []
    for i in class_checkbox.active:
        checked_class = np.append(checked_class,category_options[cat_opt.value][i])
    adata.uns[cat_opt.value][adata.uns[cat_opt.value]['class_name'].isin(list(checked_class))] = np.NAN
    del_list = list(category_options[cat_opt.value])
    del_list2 = class_checkbox.labels
    for i in range(len(checked_class)):
        del_list.remove(checked_class[i])
        del del_list2[class_checkbox.active[i]-i]
    category_options[cat_opt.value] = del_list
    class_checkbox.labels = del_list2    
    print(adata.uns[cat_opt.value])

# Gate checked classes
def gate_class():
    global view
    checked_class = []
    print('===',class_checkbox.active)
    for i in class_checkbox.active:
        checked_class = np.append(checked_class,category_options[cat_opt.value][i])
    print('_____',adata.uns[cat_opt.value][adata.uns[cat_opt.value]['class_name'].isin(list(checked_class))].index)
    view.filters = [IndexFilter(adata.uns[cat_opt.value][adata.uns[cat_opt.value]['class_name'].isin(list(checked_class))].index)]
    
# Callback of edit function change of class
def cls_func_change(attr,old,new):
    global input_t
    curdoc().remove_root(input_t)
    curdoc().remove_root(class_button)
    if edit_classes.value == 'Merge checked classes':
        input_t.value = 'Merged Class Name'
        curdoc().add_root(input_t)
        curdoc().add_root(class_button)
    elif edit_classes.value == 'Gate checked classes':
        class_button.label = 'Gate checked classes'
        curdoc().add_root(class_button)

# Callback of Class Function Change 
def clsfunc_cb(attr,old,new):
    global class_button
    curdoc().remove_root(input_t)
    curdoc().remove_root(class_button)
    curdoc().remove_root(edit_classes)
    if class_func.value == 'Create New Class':
        curdoc().add_root(input_t)
        curdoc().add_root(class_button)
        class_button.label = 'Confirm to Create New Class!'
    elif class_func.value == 'Edit Class':
        curdoc().add_root(edit_classes)
        curdoc().add_root(input_t)
        curdoc().add_root(class_button)
        class_button.label = 'Save changes to the Class'
    elif class_func.value == 'Delete Selected Class':
        curdoc().add_root(class_button)
        class_button.label =  'Confirm to Delete! (Please check the selected CLASS!!)'

# Callback of class Save Button
def save_cls_button(event):
    global save_button, cat_opt, Figure
    if class_func.value == 'Create New Class':
        add_entry()
    elif class_func.value == 'Edit Class':
        if edit_classes.value == 'Merge checked classes':
            toclass = input_t.value
            color = cur_color
            merge_class(toclass,color)
        elif edit_classes.value == 'Gate checked classes':
            gate_class()
    elif class_func.value == 'Delete Selected Class':
        del_class()

# Show color of category
def show_color():
    source.data['color'] = adata.uns[cat_opt.value]['color']

# Add function
def function_callback(attr, old, new):
    global layout, adata
    func_list = function_select.value

    for i in range(len(func_list)):
        if func_list[i] == 'Create Category':           
            curdoc().add_root(cat_func)
            
# Add New panel
def addpanel(event):
    global tab_list,figure_panel
    layout2 = figure_panel
    paneltab2 = Panel(child=layout2, title='Panel2')
    print('==========tab_list',tab_list)
    tab_list.append(paneltab2)
    print('------tab_list',tab_list)
    tabs.update(tabs=tab_list)


TOOLTIPS = [
        ("index", "@index"),
        ("(x,y)", "($x, $y)"),
        ("color", "@color"),
]

class FlowPlot:
    def __init__(self, opts, source, view, columns, title = "", x_init_idx = 0, y_init_idx = 0, allow_select = True, select_color_change = True, legend = None):
        self.opts = opts
        self.source = source
        self.view = view
        self.columns = columns
        self.p = figure(width=500, height=500, tools="pan,lasso_select,box_select,tap,wheel_zoom,save,hover",title="Surface Marker Gating Panel", tooltips=TOOLTIPS)
        #self.p.output_backend = "svg"
        print("backend is ", self.p.output_backend)
        self.p.xaxis.axis_label = self.columns[x_init_idx]
        self.p.yaxis.axis_label = self.columns[y_init_idx]
        self.r = self.p.circle(self.columns[x_init_idx], self.columns[y_init_idx],  source=self.source, view=self.view, fill_alpha=1,fill_color='color', line_color=None )
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
        self.r.selection_glyph.fill_color = cur_color
    
    def show_choosed(self,cate_name,class_name):
        
        #self.r.data_source.selected.indices = adata.uns[cate_name][adata.uns[cate_name]['class_name']==class_name].tolist()
        self.r.selection_glyph.line_color = 'black'

##############################
### Initialize data source ###
##############################
opts = dict(plot_width=500, plot_height=500, min_border=0, tools="pan,lasso_select,box_select,wheel_zoom,save")
source = ColumnDataSource(data=data_df)
view = CDSView(source=source, filters=[IndexFilter([i for i in range(data_df.shape[0])])])
Figure = FlowPlot(opts, source, view, generic_columns, "Surface Marker Gating Panel", legend='color')


##############################
##### Buttons Definition #####
##############################

# Change the color of selected parts
cur_color = color_list[0]
select_color = ColorPicker(title="Select color:", color=color_list[1], css_classes=color_list)
select_color.on_change("color", select_color_func)

color_button = Button(label="Color")
color_button.on_click(color_func)

correct_button = Button(label="Correct Plot")
correct_button.on_click(correct_func)

# Gate, remove, and show all button
gate_button = Button(label="Gate")
gate_button.on_click(selection_func)

remove_button = Button(label="Remove")
remove_button.on_click(remove_func)

showall_button = Button(label="Show All")
showall_button.on_click(showall_func)

# Show color of category
show_color_button = Button(label='Show Color of Category')
show_color_button.on_click(show_color)


### Category Functions Buttons ###
category_options = {} #A DICT to save name of Categories and the name of classes belong to each Category

# Select functions about category
cat_func = Select(title="Select Function of Category:",options=['Create New Category','Change Name of Category','Delete Selected Category'],value='Create New Category')
cat_func.on_change("value",catfunc_cb)

# Input name of new category
name = TextInput(title='Input Category Name: ', value='New Category')
name.js_on_change("value", CustomJS(code="""
    console.log('text_input: value=' + this.value, this.toString())
"""))

# Confirm to create new category
save_button = Button(label="Create New Category")
save_button.on_click(save_cb)

# Select existed categories
cat_opt = Select(title='Select Category', options=['No Category'], value='No Category')
cat_opt.on_change("value", choose_cat)


### Class Functions Button ###
# Select Function of Class
class_func = Select(title='Select Class Function',options=['Create New Class','Edit Class','Delete Selected Class','Change Category Name'],value='Create New Class')
class_func.on_change("value",clsfunc_cb)

# Input of class name
input_t = TextInput(title='Enter class name: ',value='New Class')

# Select of class (use checkbox)
cls_label = [] # Checkbox label of class
class_checkbox = CheckboxGroup(labels=cls_label,active=[0])
class_checkbox.js_on_click(CustomJS(code="""
    console.log('checkbox_group: active=' + this.active, this.toString())
"""))
#class_select = Select(title = 'Choose Class: ', options = ['Choose Class'], value='Choose Class')

# Selection of editing classes
edit_classes = Select(title='Classes function',options=['Merge checked classes','Gate checked classes'], value='Merge checked classes')
edit_classes.on_change('value',cls_func_change)

# Save button of class
class_button = Button(label = 'Add New Class')
class_button.on_event(ButtonClick, save_cls_button)


# New function button (not used now)
function_options = ['Create Category','Clustering','Cell Type Annotation']
function_select = MultiSelect(title='Function Select',options=function_options)
function_select.on_change("value",function_callback)

# Add new panel (not used now)
add_panel = Button(label='Add Panel')
add_panel.on_event(ButtonClick,addpanel)


### Layout ###
figure_panel = column(Figure.p)
control_panel = column(Figure.s_x, Figure.s_y, select_color, gate_button, remove_button, showall_button, cat_func, name, save_button)
layout = row(figure_panel,control_panel)

# Panel
panle1 = Panel(child=layout,title='Original Panel')
tab_list = [panle1]
tabs = Tabs(tabs=tab_list)


curdoc().add_root(tabs)
