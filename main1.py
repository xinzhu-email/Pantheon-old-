from cProfile import label
from bokeh.io import show
from bokeh.models import Slider, ColumnDataSource, CDSView, IndexFilter, CustomJS, Circle, Div, Panel, Tabs, CheckboxGroup, FileInput
from bokeh.models.widgets import Select, Button, ColorPicker,TextInput, DataTable, MultiSelect
from bokeh.events import ButtonClick
from bokeh.palettes import d3
from bokeh.layouts import row, column
from bokeh.io import curdoc
from bokeh.layouts import row
from bokeh.plotting import figure
from click import option
import pandas
import numpy as np
import anndata
import scipy.sparse as ss
import scanpy as sc
from torch import cat

# Loading data
#data_path = "E:/项目/图形化界面/ADT_gater-master/"
#data_path = 'CD4_memory_Naive.h5ad'
#adata = anndata.read(data_path)
adata = anndata.read_csv('ADT.csv')
#data_array = []
#data_array.append(pandas.read_csv(data_path+'ADT.csv'))
#data_df = pandas.concat(data_array,axis=1,join='inner')
data_df = adata.to_df()
generic_columns = data_df.columns.values.tolist()
data_log = np.log1p(data_df)
#print(generic_columns)
# Initialize color attribute
def color_define():
    global data_df, data_log
    #data_df['ind'] = pandas.Series(range(data_df.shape[0]), index=data_df.index)
    data_df['color'] = pandas.Series(d3['Category20c'][20][0], index=data_df.index)
    data_log['color'] = pandas.Series(d3['Category20c'][20][0], index=data_df.index)
    #data_df['label'] = pandas.Series(index=data_df.index)
    # Initialize highly variable gene
    data_df['hl_gene'] = pandas.Series(np.full(data_df.shape[0], 0), index=data_df.index)   
    adata.obs['ind'] = pandas.Series(range(data_df.shape[0]), index=data_df.index)
color_define()
color_list = d3['Category20c'][20]

upload_button = FileInput()
def eb(attr,old,new):
    print(upload_button.filename)

try:
    upload_button.filename = "ADT.csv"
except:
    d = Div(text="Please select file!")
    curdoc().add_root(d)
upload_button.on_change('filename',eb)


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
def log_cb(new):
    global source
    if log_check.active == []:
        data_df['color'] = source.data['color']
        source.data = data_df
    else:
        data_log['color'] = source.data['color']
        source.data = data_log

    

# Show all, gate, and remove function
def showall_func():
    global view
    view.filters = list([])


def gate_func():
    global view
    view.filters = [IndexFilter(source.selected.indices)]
    if len(class_checkbox.labels)==0:
        new_category()
    
    #print(source.selected.indices)

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
    show_color()
   

def save_profile():
    adata.write('./RESULT.h5ad')
    for cate in list(adata.uns['category_dict'].keys()):
        adata.obs[cate].to_csv('%s.csv'%cate)
    """label = adata.obs[cat_opt.value]
    for i in range(data_df.shape[0]):
        ind  = int(adata.obs[cat_opt.value][i])
        label[i] = adata.uns['category_dict'][cat_opt.value]['class_name'][ind]
    label.to_csv('./RESULT__%s.csv'%cat_opt.value)"""

def show_legend():
    global Figure
    Figure.source.data['label'] = list(adata.obs[cat_opt.value])
    print(source.data['label'])

##########################
### Category Functions ###
##########################
# New Category
def new_category():
    global adata, cat_opt
    marker = str(Figure.p.xaxis.axis_label) + '+' + str(Figure.p.yaxis.axis_label)
    adata.uns['category_dict'][marker] = pandas.DataFrame(columns=['class_name','color','cell_num'])
    adata.obs[marker] = pandas.Series(index=data_df.index,dtype=object)
    #adata.uns[name.value] =  pandas.DataFrame(index=range(data_df.shape[0]), columns=['class_name','color'],dtype=object)
    #adata.uns[name.value]['color'] = color_list[19]
    #category_options[name.value] = []
    cat_opt.options = [' '] + list(adata.uns['category_dict'].keys())
    cat_opt.value = marker

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
    if cat_opt.value == ' ':
        new_category()
    elif class_checkbox.labels != []:
        cls_label = []
        cate = cat_opt.value
        for i in range(len(adata.uns['category_dict'][cate])):
            class_name = adata.uns['category_dict'][cate]['class_name'][i]
            color = adata.uns['category_dict'][cate]['color'][i]
            cell_num = adata.uns['category_dict'][cate]['cell_num'][i]
            s = str(class_name) + ': color=' + str(color) + ', cell_nums=' + str(cell_num)
            cls_label = np.append(cls_label,s)
        class_checkbox.labels = list(cls_label)
        class_checkbox.active = []
        show_color()
    else:
        class_checkbox.labels = []


########################
#### Class Function ####
########################
# New Class
def add_entry():
    global cls_label, adata
    xaxis = str(Figure.p.xaxis.axis_label)
    yaxis = str(Figure.p.yaxis.axis_label)
    if cat_opt.value == ' ' or (str(cat_opt.value) != xaxis+'+'+yaxis and str(cat_opt.value) != yaxis+'+'+xaxis):
        print(str(cat_opt.value),xaxis+'+'+yaxis)
        class_checkbox.labels = []
        new_category()       
    cell_num = len(source.selected.indices)
    print(cell_num)
    adata.uns['category_dict'][cat_opt.value].loc[len(adata.uns['category_dict'][cat_opt.value])] = {'class_name':input_t.value,'color':cur_color,'cell_num':cell_num}
    save_class(cat_opt.value, input_t.value, cur_color, 0)
    input_t.value = ''

# Save change of classes
def save_class(cate, class_name, color, n):
    global adata, class_checkbox, cls_label, Figure
    
    if n == 0:
        ind = len(class_checkbox.labels)
    else:
        ind = class_checkbox.active[0]  
    class_label = list(adata.obs[cate])
    print('ind=====',ind)
    for i in source.selected.indices:
        class_label[i] = str(ind)
        #print('i:',i)
    adata.obs[cate] = class_label
    #print(class_label)
    num = n + len(source.selected.indices)
    cls_label = []
    cate = cat_opt.value
    for i in range(adata.uns['category_dict'][cate].shape[0]):
        class_name = adata.uns['category_dict'][cate]['class_name'][i]
        color = adata.uns['category_dict'][cate]['color'][i]
        cell_num = adata.uns['category_dict'][cate]['cell_num'][i]
        s = str(class_name) + ': color=' + str(color) + ', cell_nums=' + str(cell_num)
        cls_label = np.append(cls_label,s)
    class_checkbox.labels = list(cls_label)
    class_checkbox.active = [ind]
    show_color()
    #show_legend()
    #print(adata.uns[category_name])

# Merge checked classes
def merge_class():
    global class_checkbox, adata
    if input_t.value == '':
        toclass = adata.uns['category_dict'][cat_opt.value]['class_name'][class_checkbox.active[0]]
        color = adata.uns['category_dict'][cat_opt.value]['color'][class_checkbox.active[0]]
    else:
        toclass = input_t.value
        color = cur_color

    adata.uns['category_dict'][cat_opt.value].drop(index=class_checkbox.active,inplace=True)
    temp = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value],columns=['class_name','color'],index=adata.uns['category_dict'][cat_opt.value].index)
    temp['new'] = pandas.Series(range(temp.shape[0]),index=temp.index)
    count = 0
    clr = source.data['color']
    for i in range(data_df.shape[0]):
        ind  = str(adata.obs[cat_opt.value][i])
        old = True
        for j in range(len(class_checkbox.active)):
            if adata.obs[cat_opt.value][i] == str(class_checkbox.active[j]):
                old = False
                clr[i] = color_list[18]
                adata.obs[cat_opt.value][i] = str(len(temp))
                count = count + 1
                break
        if old:
            try:
                adata.obs[cat_opt.value][i] = str(temp[temp.index == int(ind)].loc['new'])
            except:
                count = count
    print('----',adata.uns['category_dict'][cat_opt.value])
    adata.uns['category_dict'][cat_opt.value] = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value],index=temp['new'])
    adata.uns['category_dict'][cat_opt.value].loc[len(temp)+1] = {'class_name':toclass,'color':color,'cell_num':count}
    del_list2 = class_checkbox.labels
    for i in range(len(class_checkbox.active)):
        del del_list2[class_checkbox.active[i]-i]
    del_list2 = del_list2 + [toclass+ ': color=' + str(color) + ', cell_nums='+ str(count)]
    class_checkbox.labels = del_list2
    class_checkbox.active = []

# Delete Class
def del_class():
    global class_checkbox, adata, source
    adata.uns['category_dict'][cat_opt.value].drop(index=class_checkbox.active,inplace=True)
    temp = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value],columns=['class_name','color','cell_num'],index=adata.uns['category_dict'][cat_opt.value].index)
    temp['new'] = pandas.Series(range(temp.shape[0]),index=temp.index)
    count = 0
    clr = source.data['color']
    print(clr)
    for i in range(data_df.shape[0]):
        ind  = str(adata.obs[cat_opt.value][i])
        old = True
        for j in range(len(class_checkbox.active)):
            if adata.obs[cat_opt.value][i] == str(class_checkbox.active[j]):
                old = False
                clr[i] = color_list[18]
                adata.obs[cat_opt.value][i] = np.nan
                break
        if old:
            try:
                adata.obs[cat_opt.value][i] = str(temp[temp.index == int(ind)].loc['new'])
            except:
                count = count
    del_list2 = class_checkbox.labels
    for i in range(len(class_checkbox.active)):
        del del_list2[class_checkbox.active[i]-i]
    class_checkbox.labels = del_list2
    class_checkbox.active = []
    adata.uns['category_dict'][cat_opt.value] = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value],index=temp['new'])
    source.data['color'] = clr

# Rename class
def rename():
    ind = class_checkbox.active[0]
    color = adata.uns['category_dict'][cat_opt.value]['color'][ind]
    cell_num = adata.uns['category_dict'][cat_opt.value]['cell_num'][ind]
    labels = class_checkbox.labels
    labels[ind] = str(input_t.value) + ': color=' + str(color) + ', cell_nums=' + str(cell_num)
    print(labels)
    class_checkbox.labels = labels


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


def save_cls_button(event):
    class_name = adata.uns['category_dict'][cat_opt.value]['class_name'][class_checkbox.active[0]]
    color = adata.uns['category_dict'][cat_opt.value]['color'][class_checkbox.active[0]]
    cell_num = len(source.selected.indices)
    save_class(cat_opt.value, class_name,color,cell_num)

# Show color of category
def show_color():
    global source
    col_list = source.data['color']
    print(adata.obs[cat_opt.value])
    for i in range(data_df.shape[0]):
        ind = adata.obs[cat_opt.value][i]
        #print(ind)
        try:
            ind = int(ind)
            col_list[i] = adata.uns['category_dict'][cat_opt.value]['color'][ind]
        except:
            col_list[i] = color_list[18]
    source.data['color'] = col_list

# change color of class
def change_color():
    color_l = source.data['color']
    for i in source.selected.indices:
        color_l[i] = cur_color
    ind = class_checkbox.active[0]
    adata.uns['category_dict'][cat_opt.value]['color'][ind] = cur_color
    source.data['color'] = color_l



TOOLTIPS = [
        ("(x,y)", "($x, $y)"),
        ("color", "@color"),
]

class FlowPlot:
    def __init__(self, opts, source, view, columns, title = "", x_init_idx = 0, y_init_idx = 0, allow_select = True, select_color_change = False, legend = None):
        self.opts = opts
        self.source = source
        self.view = view
        self.columns = columns
        self.p = figure(width=500, height=500, tools="pan,lasso_select,box_select,tap,wheel_zoom,save,hover",title="Surface Marker Gating Panel", tooltips=TOOLTIPS)
        #self.p.output_backend = "svg"
        print("backend is ", self.p.output_backend)
        self.p.xaxis.axis_label = self.columns[x_init_idx]
        self.p.yaxis.axis_label = self.columns[y_init_idx]
        self.r = self.p.circle(self.columns[x_init_idx], self.columns[y_init_idx],  source=self.source, view=self.view, fill_alpha=1,fill_color='color',line_color=None )
        self.p.legend.click_policy="hide"
        self.s_x = Select(title="x:", value=self.columns[x_init_idx], options=self.columns)
        self.s_y = Select(title="y:", value=self.columns[y_init_idx], options=self.columns)
        # Attach reaction
        self.s_x.on_change("value", lambda attr, old, new: tag_func(self.s_x, self.r.glyph, 'x', self.p) )
        self.s_y.on_change("value", lambda attr, old, new: tag_func(self.s_y, self.r.glyph, 'y', self.p) )
        # Set default fill color
        if select_color_change:
            self.r.selection_glyph = Circle(fill_alpha=1,fill_color=None, line_color='black')
        self.allow_select = allow_select

    def refresh(self):
        global cur_color
        #print('color list of data: ',self.r.selection_glyph.fill_color)
        self.r.selection_glyph = Circle(fill_alpha=1,fill_color=None)
        self.r.selection_glyph.fill_color = cur_color
    
    def show_legend(self):
        global source
        label_list = adata.uns['category_dict'][cat_opt.value]['class_name']
        source.data['label'] = list(adata.obs[cat_opt.value])
        print(source.data['label'])
        #self.p.add_glyph.legend_group = 'label'



##############################
### Initialize data source ###
##############################
def define(data_df):
    opts = dict(plot_width=500, plot_height=500, min_border=0, tools="pan,lasso_select,box_select,wheel_zoom,save")
    source = ColumnDataSource(data=data_df)
    view = CDSView(source=source, filters=[IndexFilter([i for i in range(data_df.shape[0])])])
    return opts, source, view

opts, source, view = define(data_df)
Figure = FlowPlot(opts, source, view, generic_columns, "Surface Marker Gating Panel", legend='color')


##############################
##### Buttons Definition #####
##############################

# Log
log_check = CheckboxGroup(labels=['Log axis'],active=[])
log_check.on_click(log_cb)

# Change the color of selected parts
cur_color = color_list[0]
select_color = ColorPicker(title="Select color:", color=color_list[0], css_classes=color_list)
select_color.on_change("color", select_color_func)

# Gate, remove, and show all button
gate_button = Button(label="Gate")
gate_button.on_click(gate_func)

remove_button = Button(label="Trim")
remove_button.on_click(remove_func)

showall_button = Button(label="Show All")
showall_button.on_click(showall_func)


# Select the class of the dot
dot_class = Button(label='Show the class')
dot_class.on_click(select_class)

# Show color of category
show_color_button = Button(label='Show Color of Classes')
show_color_button.on_click(show_color)



##################################
### Category Functions Buttons ###

# Marker
cat_opt = Select(title='Marker of Groups:',options=[' '],value=' ')
cat_opt.on_change('value',choose_cat)

# Input name of new category
name = TextInput(title='Input Name: ', value=' ')
name.js_on_change("value", CustomJS(code="""
    console.log('text_input: value=' + this.value, this.toString())
"""))


# Select existed categories
#cat_opt = Select(title='Select Category', options=['New Category'], value='New Category')
#cat_opt.on_change("value", choose_cat)

##############################
### Class Functions Button ###


# Input of class name
input_t = TextInput(title='Input name: ',value='')


    

# Select of class (use checkbox)
cls_label = [] # Checkbox label of class
class_checkbox = CheckboxGroup(labels=cls_label,active=[0])
class_checkbox.on_click(show_checked)
#class_checkbox.js_on_click(CustomJS(code="""
#    console.log('checkbox_group: active=' + this.active, this.toString())
#"""))
#class_select = Select(title = 'Choose Class: ', options = ['Choose Class'], value='Choose Class')

# Create new group
create_button = Button(label='Create Cluster')
create_button.on_click(add_entry)

# Delete Group
del_button = Button(label='Delete Cluster')
del_button.on_click(del_class)

# Rename Group
rename_button = Button(label='Rename Cluster')
rename_button.on_click(rename)

# Add dots button
add_dots = Button(label='Add dots')
add_dots.on_event(ButtonClick,save_cls_button)

# Merge Button
merge_button = Button(label='Merge Cluster')
merge_button.on_click(merge_class)

# Change Color
change_clr_button = Button(label='Change Color')
change_clr_button.on_click(change_color)

# Save button of class
class_button = Button(label = 'Add New Class')
class_button.on_event(ButtonClick, save_cls_button)


# Export Button
export_button = Button(label='Export Results')
export_button.on_click(save_profile)



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
file_panel = row(upload_button, export_button)
figure_panel = column(Figure.p)
control_panel = column(Figure.s_x, Figure.s_y, log_check, select_color, gate_button, remove_button, showall_button)
class_panel = column(cat_opt, input_t, create_button, rename_button, change_clr_button, merge_button, del_button, add_dots, class_checkbox)
layout = column(file_panel,row(figure_panel,column(control_panel),class_panel))

# Panel
panle1 = Panel(child=layout,title='Original Panel')
tab_list = [panle1]
tabs = Tabs(tabs=tab_list)


curdoc().add_root(layout)
