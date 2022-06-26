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
import scipy.sparse as ss
import scanpy as sc

# Loading data
data_path = "E:/项目/图形化界面/ADT_gater-master/"
#data_path = 'CD4_memory_Naive.h5ad'
#adata = anndata.read(data_path)
adata = anndata.read_csv(data_path + 'ADT.csv')
#data_array = []
#data_array.append(pandas.read_csv(data_path+'ADT.csv'))
#data_df = pandas.concat(data_array,axis=1,join='inner')
data_df = adata.to_df()
generic_columns = data_df.columns.values.tolist()
data_log = np.log1p(data_df)
#print(generic_columns)
# Initialize color attribute
#data_df['ind'] = pandas.Series(range(data_df.shape[0]), index=data_df.index)
data_df['color'] = pandas.Series(d3['Category20c'][20][0], index=data_df.index)
data_log['color'] = pandas.Series(d3['Category20c'][20][0], index=data_df.index)
# Initialize highly variable gene
data_df['hl_gene'] = pandas.Series(np.full(data_df.shape[0], 0), index=data_df.index)
color_list = d3['Category20c'][20]
adata.obs['ind'] = pandas.Series(range(data_df.shape[0]), index=data_df.index)


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

# Log
def log_cb(attr, old, new):
    global Figure
    if log_select.value == 'Raw data':
        data_df['color'] = data_log['color']
        Figure.source.data = data_df
    else:
        data_log['color'] = data_df['color']
        Figure.source.data = data_log

    

# Show all, gate, and remove function
def showall_func():
    global view
    view.filters = list([])

def selection_func():
    global view
    view.filters = [IndexFilter(source.selected.indices)]
    print(source.selected.indices)

def remove_func():
    global view, source
    remove_set = set(source.selected.indices)
    if len(view.filters) == 0:
        view.filters = [IndexFilter(np.object_(range(data_df.shape[0])))]
    remain_indices = [x for x in view.filters[0].indices if x not in source.selected.indices]
    view.filters = [IndexFilter(remain_indices)]

def select_class():
    global source,class_checkbox
    ind = source.selected.indices[0]
    class_id = adata.obs.iloc[ind][cat_opt.value]
    class_checkbox.active = [int(class_id)]
    #temp = pandas.DataFrame(list(adata.obs[cat_opt.value]),index=range(adata.obs.shape[0]),columns=[cat_opt.value])
    #print(class_id)   
    source.selected.indices = list(adata.obs[adata.obs[cat_opt.value]==class_id]['ind'])
    print(source.selected.indices)

    
def show_checked(new):
    global source
    #source.selected.indices = temp[temp[cat_opt.value]==str(class_checkbox.active[0])].index
    source.selected.indices = list(adata.obs[adata.obs[cat_opt.value].isin(list(str(i) for i in class_checkbox.active))]['ind'])
   

def save_profile():
    adata.write('./RESULT.h5ad')
    adata.obs[cat_opt.value].to_csv('RESULT.csv')
    """label = adata.obs[cat_opt.value]
    for i in range(data_df.shape[0]):
        ind  = int(adata.obs[cat_opt.value][i])
        label[i] = adata.uns['category_dict'][cat_opt.value]['class_name'][ind]
    label.to_csv('./RESULT__%s.csv'%cat_opt.value)"""

##########################
### Category Functions ###
##########################
# New Category
def new_category():
    global adata, cat_opt
    adata.uns['category_dict'][name.value] = pandas.DataFrame(columns=['class_name','color'])
    adata.obs[name.value] = pandas.Series(index=data_df.index,dtype=object)
    #adata.uns[name.value] =  pandas.DataFrame(index=range(data_df.shape[0]), columns=['class_name','color'],dtype=object)
    #adata.uns[name.value]['color'] = color_list[19]
    #category_options[name.value] = []
    cat_opt.options = list(adata.uns['category_dict'].keys())
    cat_opt.value = cat_opt.options[0]
    print(adata.uns['category_dict'])

# Edit category
def edit_category():   
    global cat_opt, adata
    old_name = cat_opt.value
    new_name = name.value
    adata.obs[new_name] = adata.obs.pop(old_name)
    adata.uns['category_dict'][new_name] = adata.uns['category_dict'].pop(old_name)
    #category_options[new_name] = category_options.pop(old_name)
    cat_opt.options = list(adata.uns['category_name'].keys())
    cat_opt.value = cat_opt.options[0]
    print(adata.uns['category_dict'])

# Delete category
def del_category():
    global cat_opt, adata
    del adata.uns['category_dict'][cat_opt.value]   
    del adata.obs[cat_opt.value]
    cat_opt.options = list(adata.uns['category_dict'].keys())
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
        for i in range(len(adata.uns['category_dict'][cate])):
            l = adata.obs[cate][adata.obs[cate] == i]
            s = adata.uns['category_dict'][cate]['class_name'][i] + ': color=' + str(adata.uns['category_dict'][cate]['color'][i]) + ', cell_nums='+ str(l.shape[0])
            cls_label = np.append(cls_label,s)
        class_checkbox.labels = list(cls_label)
        class_checkbox.active = []
        show_color()

# Callback of function selection about category
def catfunc_cb(attr,old,new):
    global name, save_button
    curdoc().remove_root(cat_opt)
    curdoc().remove_root(name)
    curdoc().remove_root(save_button)
    if cat_func.value == 'Create New Category':       
        save_button.label = 'Confirm to Create New Category!'
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
        #curdoc().add_root(class_func)
        #curdoc().add_root(input_t)
        #curdoc().add_root(class_button)
    elif cat_func.value == 'Change Name of Category':
        edit_category()
        curdoc().add_root(cat_opt)
    elif cat_func.value == 'Delete Selected Category':
        del_category()
        curdoc().add_root(cat_opt)

########################
#### Class Function ####
########################
# New Class
def add_entry():
    global cls_label, adata
    adata.uns['category_dict'][cat_opt.value].loc[len(adata.uns['category_dict'][cat_opt.value])] = {'class_name':input_t.value,'color':cur_color}
    save_class(cat_opt.value, input_t.value, cur_color, 0)
    input_t.value = ''

# Save change of classes
def save_class(category_name, class_name, color, n):
    global adata, class_checkbox, cls_label
    num = n
    ind = adata.uns['category_dict'][category_name][adata.uns['category_dict'][category_name]['class_name'] == class_name].index[0]   
    class_label = list(adata.obs[category_name])
    #print('ind=====',class_label)
    for i in source.selected.indices:
        class_label[i] = str(ind)
        num = num + 1
    adata.obs[category_name] = class_label
    if n == 0:
        cls_label = class_checkbox.labels
        cls_label = cls_label + [str(class_name) + ': color=' + str(color) + ', cell_nums=' + str(num)]
        class_checkbox.labels = cls_label
    else:
        class_checkbox.labels[class_checkbox.active[0]] = str(class_name) + ': color=' + str(color) + ', cell_nums=' + str(num)
    #print(adata.uns[category_name])

# Merge checked classes
def merge_class(toclass,color):
    global class_checkbox, adata
    adata.uns['category_dict'][cat_opt.value].drop(index=class_checkbox.active,inplace=True)
    temp = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value]) 
    temp['new'] = range(len(temp))
    count = 0
    for i in range(data_df.shape[0]):
        ind  = adata.obs[cat_opt.value][i]
        old = True
        for j in range(len(class_checkbox.active)):
            if adata.obs[cat_opt.value][i] == class_checkbox.active[j]:
                old = False
                adata.obs[cat_opt.value][i] = len(temp)
                count = count + 1
                break
        if old:
            adata.obs[cat_opt.value][i] = temp[temp.index == ind]['new']
    print('----',adata.uns['category_dict'][cat_opt.value])
    adata.uns['category_dict'][cat_opt.value] = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value],index=temp['new'])
    adata.uns['category_dict'][cat_opt.value].loc[i] = {'class_name':toclass,'color':color}
    del_list2 = class_checkbox.labels
    for i in range(len(class_checkbox.active)):
        del del_list2[class_checkbox.active[i]-i]
    del_list2 = del_list2 + [toclass+ ': color=' + str(color) + ', cell_nums='+ str(count)]
    class_checkbox.labels = del_list2
    class_checkbox.active = []

# Delete Class
def del_class():
    global class_checkbox, adata
    adata.uns['category_dict'][cat_opt.value].drop(index=class_checkbox.active,inplace=True)
    temp = adata.uns['category_dict'][cat_opt.value]
    temp['new'] = np.object_(range(temp.shape[0]))
    count = 0
    for i in range(data_df.shape[0]):
        ind  = adata.obs[cat_opt.value][i]
        old = True
        for j in range(len(class_checkbox.active)):
            if adata.obs[cat_opt.value][i] == class_checkbox.active[j]:
                old = False
                count = count + 1
                break
        if old:
            adata.obs[cat_opt.value][i] = temp[temp.index == ind]['new']
    adata.uns['category_dict'][cat_opt.value] = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value],index=temp['new'])

# Gate checked classes
def gate_class():
    global view
    checked_class = []
    print('===',class_checkbox.active)
    for i in class_checkbox.active:
        checked_class = np.append(checked_class,adata.uns['category_dict'][cat_opt.value][i])
    print('_____',adata.uns[cat_opt.value][adata.uns[cat_opt.value]['class_name'].isin(list(checked_class))].index)
    temp = pandas.DataFrame(list(adata.obs[cat_opt]),index=None)
    view.filters = [IndexFilter(temp[temp.isin(list(checked_class))].index)]
    
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
    elif edit_classes.value == 'Add dots to the class':
        class_button.label = 'Add selected dots to the checked class'
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
        show_color()
    elif class_func.value == 'Edit Class':
        if edit_classes.value == 'Merge checked classes':
            toclass = input_t.value
            color = cur_color
            merge_class(toclass,color)
        elif edit_classes.value == 'Gate checked classes':
            gate_class()
        elif edit_classes.value == 'Add dots to the class':
            class_name = adata.uns['category_dict'][cat_opt.value]['class_name'][class_checkbox.active[0]]
            color = adata.uns['category_dict'][cat_opt.value]['color'][class_checkbox.active[0]]
            print(adata.obs[cat_opt.value],str(class_checkbox.active[0]))
            l = adata[adata.obs[cat_opt.value]==str(class_checkbox.active[0])]
            print('===',l.shape[0],'==',np.object_(class_checkbox.active[0]),'==')
            save_class(cat_opt.value, class_name,color,len(l))
    elif class_func.value == 'Delete Selected Class':
        del_class()

# Show color of category
def show_color():
    col_list = source.data['color']
    for i in range(data_df.shape[0]):
        ind = adata.obs.iloc[i][cat_opt.value]
        try:
            ind = int(ind)
            col_list[i] = adata.uns['category_dict'][cat_opt.value]['color'][ind]
        except:
            col_list[i] = color_list[1]
    source.data['color'] = col_list

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

# Save labels
save_label = Button(label='Save profiles')
save_label.on_click(save_profile)

# Log
log_select = Select(title='Log or not',options=['Raw data','Log data'],value='Raw data')
log_select.on_change('value',log_cb)

# Select the class of the dot
dot_class = Button(label='Show the class')
dot_class.on_click(select_class)

# Show color of category
show_color_button = Button(label='Show Color of Category')
show_color_button.on_click(show_color)

##################################
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
cat_opt = Select(title='Select Category', options=['New Category'], value='New Category')
cat_opt.on_change("value", choose_cat)

##############################
### Class Functions Button ###
# Select Function of Class
class_func = Select(title='Select Class Function',options=['Create New Class','Edit Class','Delete Selected Class'],value='Create New Class')
class_func.on_change("value",clsfunc_cb)

# Input of class name
input_t = TextInput(title='Enter class name: ',value='New Class')


    

# Select of class (use checkbox)
cls_label = [] # Checkbox label of class
class_checkbox = CheckboxGroup(labels=cls_label,active=[0])
class_checkbox.on_click(show_checked)
#class_checkbox.js_on_click(CustomJS(code="""
#    console.log('checkbox_group: active=' + this.active, this.toString())
#"""))
#class_select = Select(title = 'Choose Class: ', options = ['Choose Class'], value='Choose Class')

# Selection of editing classes
edit_classes = Select(title='Edit Function:',options=['Merge checked classes','Gate checked classes','Add dots to the class'], value='Add dots to the class')
edit_classes.on_change('value',cls_func_change)

# Save button of class
class_button = Button(label = 'Add New Class')
class_button.on_event(ButtonClick, save_cls_button)

add_dots = Button(label='Add dots')
add_dots.on_event(ButtonClick,save_cls_button)


# New function button (not used now)
function_options = ['Create Category','Clustering','Cell Type Annotation']
function_select = MultiSelect(title='Function Select',options=function_options)
function_select.on_change("value",function_callback)

# Add new panel (not used now)
add_panel = Button(label='Add Panel')
add_panel.on_event(ButtonClick,addpanel)

def axis_cb(attr,old,new):
    global Figure, data_df
    if choose_panel.value != 'generic_columns':
        #print(adata.obsm[choose_panel.value])
        data_df = pandas.DataFrame(data=adata.obsm[choose_panel.value],columns=['0','1'])
        Figure.source = ColumnDataSource(data=data_df)     
        Figure.s_x.options = list(str(i) for i in range(data_df.shape[1]))
        Figure.s_y.options = list(str(i) for i in range(data_df.shape[1]))
        Figure.s_x.value = '0'
        Figure.s_y.value = '1'
        Figure.source.data['color'] = pandas.Series(d3['Category20c'][20][0],index=data_df.index)
    else:
        Figure.s_x.options = generic_columns
        Figure.s_y.options = generic_columns
        Figure.s_x.value = generic_columns[0]
        Figure.s_y.value = generic_columns[1]

# Read category
try:
    category_name = adata.uns['category_dict'].keys()
    cat_opt.options = list(category_name)
    cat_opt.value = list(category_name)[0]
    cate_panel = column(cat_opt)
    cls_label = []
    cate = cat_opt.value
    choose_panel = Select(title='Choose map:',value='generic_columns',options=list(adata.obsm.keys())+['generic_columns'])
    choose_panel.on_change('value',axis_cb)
    curdoc().add_root(choose_panel)
    for i in range(len(adata.uns['category_dict'][cate])):
        l = adata.obs[cate][adata.obs[cate] == str(i)]
        s = str(adata.uns['category_dict'][cate]['class_name'][i]) + ': color=' + str(adata.uns['category_dict'][cate]['color'][i]) + ', cell_nums='+ str(l.shape[0])
        cls_label = np.append(cls_label,s)
    show_color()
    class_checkbox.labels = list(cls_label)
    class_checkbox.active = []
except:
    adata.uns['category_dict'] = dict({'New Category':pandas.DataFrame(columns=['class_name','color'])})
    adata.obs['New Category'] = pandas.Series(index=data_df.index,dtype=object)
    d = Div(text='No Existed Category! So Create New Category!')
    curdoc().add_root(d)
    cate_panel = column(cat_opt)


### Layout ###
figure_panel = column(Figure.p)
control_panel = column(Figure.s_x, Figure.s_y, select_color, log_select, gate_button, remove_button, showall_button,save_label,dot_class)
class_panel = column(cate_panel, class_func,input_t,class_button,class_checkbox)
layout = row(figure_panel,column(control_panel),class_panel)

# Panel
panle1 = Panel(child=layout,title='Original Panel')
tab_list = [panle1]
tabs = Tabs(tabs=tab_list)


curdoc().add_root(tabs)
