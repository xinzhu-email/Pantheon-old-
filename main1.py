from bokeh.models import ColumnDataSource, CDSView, IndexFilter, CustomJS, Circle, Div, Panel, Tabs, CheckboxGroup, FileInput,FixedTicker, ColorBar, LogColorMapper
from bokeh.models.widgets import Select, Button, ColorPicker,TextInput, DataTable, MultiSelect, AutocompleteInput
from bokeh.events import ButtonClick
from bokeh.transform import linear_cmap, log_cmap
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
import colorcet as cc
from torch import cat

# Loading data
#adata = anndata.read('CD4_memory_Naive.h5ad')
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
    data_df['hl_gene'] = pandas.Series(np.full(data_df.shape[0], 3), index=data_df.index)   
    adata.obs['ind'] = pandas.Series(np.array(range(data_df.shape[0])).astype(int).tolist(), index=data_df.index)
    print(adata.obs['ind'][0:20])
color_define()
color_list = d3['Category20c'][20]
cur_color = color_list[0]

upload_button = FileInput()
def eb(attr,old,new):
    print(upload_button.filename)

try:
    upload_button.filename = "ADT.csv"
except:
    d = Div(text="Please select file!")
    #curdoc().add_root(d)
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
        data_df['hl_gene'] = source.data['hl_gene']
        source.data = data_df
    else:
        data_log['color'] = source.data['color']
        data_log['hl_gene'] = source.data['hl_gene']
        source.data = data_log

    

# Show all, gate, and remove function
def showall_func():
    global view
    view.filters = list([])


def gate_func():
    global view
    view.filters = [IndexFilter(source.selected.indices)]
    
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

    
def show_checked():
    global source
    #source.selected.indices = temp[temp[cat_opt.value]==str(class_checkbox.active[0])].index
    source.selected.indices = list(adata.obs[adata.obs[cat_opt.value].isin(list(str(i) for i in class_checkbox.active))]['ind'])
    show_color()

   

def save_profile():
    adata.write('./RESULT.h5ad')
    for cate in list(adata.uns['category_dict'].keys()):
        adata.obs[cate].to_csv('%s.csv'%cate)
    #adata.uns['category_dict']('cluster name.csv') 

def show_legend():
    global Figure
    Figure.source.data['label'] = list(adata.obs[cat_opt.value])
    print(source.data['label'])

##########################
### Category Functions ###
##########################
# New Category
def new_category():
    global adata, cat_opt, name
    if name.value == '':
        marker = str(Figure.p.xaxis.axis_label) + '+' + str(Figure.p.yaxis.axis_label)
    else:
        marker = name.value
    adata.uns['category_dict'][marker] = pandas.DataFrame(columns=['class_name','color','cell_num'])
    adata.obs[marker] = pandas.Series(index=data_df.index,dtype=object)
    #adata.uns[name.value] =  pandas.DataFrame(index=range(data_df.shape[0]), columns=['class_name','color'],dtype=object)
    #adata.uns[name.value]['color'] = color_list[19]
    #category_options[name.value] = []
    cat_opt.options = list(adata.uns['category_dict'].keys())
    cat_opt.value = marker
    name.value = ''

    

# Edit category
def edit_category():   
    global cat_opt, adata
    old_name = cat_opt.value
    new_name = name.value
    adata.obs[new_name] = adata.obs.pop(old_name)
    adata.uns['category_dict'][new_name] = adata.uns['category_dict'].pop(old_name)
    cat_opt.options = list(adata.uns['category_dict'].keys())
    cat_opt.value = cat_opt.options[0]
    print(adata.uns['category_dict'])

# Delete category
def del_category():
    global cat_opt, adata
    del adata.uns['category_dict'][cat_opt.value]   
    del adata.obs[cat_opt.value]
    cat_opt.options = list(adata.uns['category_dict'].keys())
    if len(cat_opt.options) == 0:
        cat_opt.value = ''
    else:
        cat_opt.value = cat_opt.options[0]

# Choose Category
def choose_cat(attr,old,new):
    global source, class_checkbox
    try: 
        adata.uns['category_dict'][cat_opt.value]['class_name'][0]
        update_checkbox()
        class_checkbox.active = [0]
        show_color()
    except:
        class_checkbox.labels = ['Unassigned: color=grey, cell_nums=' + str(data_df.shape[0])]
        class_checkbox.active = []
    text_color()

########################
#### Class Function ####
########################
# New Class
def add_entry():
    global cls_label, adata
    xaxis = str(Figure.p.xaxis.axis_label)
    yaxis = str(Figure.p.yaxis.axis_label)
    #if str(cat_opt.value) != xaxis+'+'+yaxis and str(cat_opt.value) != yaxis+'+'+xaxis:
    if cat_opt.value == ' ':
        print(str(cat_opt.value),xaxis+'+'+yaxis)
        #class_checkbox.labels = ['no cluster: color=' + color_list[18] + ', cell_nums=' + str(data_df.shape[0])]
        new_category()       
    cell_num = len(source.selected.indices)
    #print(cell_num)
    adata.uns['category_dict'][cat_opt.value].loc[len(adata.uns['category_dict'][cat_opt.value])] = {'class_name':input_t.value,'color':cur_color,'cell_num':cell_num}
    save_class(cat_opt.value, input_t.value, cur_color, 0)
    input_t.value = ''
    
    #print(now_color,hide.value)

# Save change of classes
def save_class(cate, class_name, color, n):
    global adata, class_checkbox, cls_label, Figure, cls_label
    
    if n == 0:
        ind = len(class_checkbox.labels)-1
    else:
        ind = class_checkbox.active[0]  
    class_label = list(adata.obs[cate])
    #print('ind=====',ind)
    for i in source.selected.indices:
        class_label[i] = str(ind)
        #print('i:',i)
    adata.obs[cate] = class_label
    #print(class_label)
    cate = cat_opt.value
    num = 0
    cls_label = [] 
    update_checkbox()
    class_checkbox.active = [ind]
    show_color()
    correct_func()
    text_color()
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
    temp = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value],columns=['class_name','color','cell_num'],index=adata.uns['category_dict'][cat_opt.value].index)
    temp['new'] = pandas.Series(range(temp.shape[0]),index=temp.index)
    count = 0
    clr = source.data['color']
    for i in range(data_df.shape[0]):
        ind  = str(adata.obs[cat_opt.value][i])
        old = True
        for j in range(len(class_checkbox.active)):
            if adata.obs[cat_opt.value][i] == str(class_checkbox.active[j]):
                old = False
                clr[i] = color
                adata.obs[cat_opt.value][i] = str(len(temp))
                count = count + 1
                break
        if old:
            try:
                adata.obs[cat_opt.value][i] = str(temp[temp.index == int(ind)].loc['new'])
            except:
                count = count
    
    adata.uns['category_dict'][cat_opt.value] = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value],index=temp['new'])
    print('----',adata.uns['category_dict'][cat_opt.value])
    adata.uns['category_dict'][cat_opt.value].loc[len(temp)] = {'class_name':toclass,'color':color,'cell_num':count}
    del_list2 = class_checkbox.labels
    for i in range(len(class_checkbox.active)):
        del del_list2[class_checkbox.active[i]-i]
    tt = del_list2[-1]
    del_list2[-1] = str(toclass+ ': cell_nums='+ str(count))
    del_list2 = del_list2 + [tt]
    class_checkbox.labels = del_list2
    class_checkbox.active = []
    source.data['color'] = clr
    text_color()

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
                count  = count + 1
                break
        if old: 
            try:
                adata.obs[cat_opt.value][i] = str(temp[temp.index == int(ind)].loc['new'])
            except:
                count = count + 1
    del_list2 = class_checkbox.labels
    for i in range(len(class_checkbox.active)):
        del del_list2[class_checkbox.active[i]-i]
    del_list2[-1] = str('unassigned: color=grey, cell_nums=' + str(count))
    class_checkbox.labels = del_list2
    class_checkbox.active = []
    adata.uns['category_dict'][cat_opt.value] = pandas.DataFrame(adata.uns['category_dict'][cat_opt.value],index=temp['new'])
    source.data['color'] = clr
    text_color()

# Update class
def update_clus():
    global adata
    ind = class_checkbox.active[0]
    #adata.obs[adata.obs[cat_opt.value]==str(ind)].loc[cat_opt.value] = pandas.Series(index=[0],dtype=object)[0]
    cl_label = adata.obs[cat_opt.value]
    cl_label[adata.obs[cat_opt.value]==str(ind)] = np.NAN
    print(cl_label)
    for i in source.selected.indices:
        cl_label[i] = str(ind)
    adata.obs[cat_opt.value] = cl_label
    update_checkbox()
    show_color()
    text_color()


# Rename class
def rename():
    global class_checkbox, adata
    ind = class_checkbox.active[0]
    cell_num = adata.uns['category_dict'][cat_opt.value]['cell_num'][ind]
    labels = class_checkbox.labels
    labels[ind] = str(input_t.value) + ': cell_nums=' + str(cell_num)
    adata.uns['category_dict'][cat_opt.value]['class_name'] = input_t.value
    input_t.value = ''
    print(labels)
    class_checkbox.labels = labels
    text_color()



# Change color of checkbox
hide = Div(text='0',visible=False,css_classes=['hide'])
hide1 = Div(text='1',visible=False)
color_js = [color_list[18]]
length = 1

hide1.js_on_change('text',CustomJS(code="""
    const collection = document.getElementsByClassName("class_checkbox_label");
    var str = document.getElementsByClassName('hide')[0].children[0].innerHTML;
    console.log(document.getElementsByClassName('hide')[0].children[0].innerHTML);
    const color = str.split(' ');
    var k = color.length;
    console.log(k,color);
    for (var i=0;i<k;i++)
    {
        collection[0].children[0].children[i].style.color = color[i];
    }
    console.log('collection:' + collection[0].children[0].innerHTML);
    
"""))

def text_color():
    global hide
    color_js = ''
    try:
        length = len(adata.uns['category_dict'][cat_opt.value]['color']) 
        for i in range(length):
            color_js = color_js + adata.uns['category_dict'][cat_opt.value]['color'][i] + ' '
        hide.text = color_js + color_list[18]
    except:
        length = 1
        #color_js = [color_list[18]]
        hide.text = str(color_list[18])
    hide1.text = hide1.text + '1'
    print('CALL',color_js)
    
# Add dots to cluster
def save_cls_button(event):
    class_name = adata.uns['category_dict'][cat_opt.value]['class_name'][class_checkbox.active[0]]
    color = adata.uns['category_dict'][cat_opt.value]['color'][class_checkbox.active[0]]
    cell_num = len(source.selected.indices)
    save_class(cat_opt.value, class_name,color,cell_num)

# Remove dots from cluster
def remove_dot():
    global adata, class_checkbox
    cl_label = adata.obs[cat_opt.value]
    for i in source.selected.indices:
        for j in class_checkbox.active:
            if cl_label[i] == str(j):
                cl_label[i] = np.NAN
                break
    adata.obs[cat_opt.value] = cl_label
    update_checkbox()
    show_color()

# Show color of category
def show_color():
    global source
    col_list = source.data['color']
    #print(adata.obs[cat_opt.value])
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
    global adata, class_checkbox, source
    color_l = source.data['color']
    show_checked()
    for i in source.selected.indices:
        color_l[i] = cur_color
    source.data['color'] = color_l
    adata.uns['category_dict'][cat_opt.value]['color'][[i for i in class_checkbox.active]] = cur_color
    text_color()
    #print(hide.value,now_color)


def update_checkbox():
    global class_checkbox    
    cate = cat_opt.value
    cls_label = []
    num = 0
    for i in range(adata.uns['category_dict'][cate].shape[0]):
        class_name = adata.uns['category_dict'][cate]['class_name'][i]
        cell_num = len(data_df[adata.obs[cate]==str(i)])       
        s = str(class_name) +  ': cell_nums=' + str(cell_num)
        cls_label = np.append(cls_label,s)
        num = num + cell_num
    cls_label = np.append(cls_label,str('Unassigned: color=grey, cell_nums=' + str(data_df.shape[0]-num))) 
    class_checkbox.labels = list(cls_label)



TOOLTIPS = [
        ("(x,y)", "($x, $y)"),
        ("color", "@color"),
]

class FlowPlot:
    def __init__(self, opts, source, view, columns, color_map, title = "", x_init_idx = 0, y_init_idx = 0, allow_select = True, select_color_change = True, legend = None):
        self.opts = opts
        self.source = source
        self.view = view
        self.columns = columns
        #self.color_map = color_map
        self.p = figure(width=500, height=500, tools="pan,lasso_select,box_select,tap,wheel_zoom,save,hover",title="Surface Marker Gating Panel", tooltips=TOOLTIPS)
        #self.p.output_backend = "svg"
        print("backend is ", self.p.output_backend)
        self.p.xaxis.axis_label = self.columns[x_init_idx]
        self.p.yaxis.axis_label = self.columns[y_init_idx]
        self.r = self.p.circle(self.columns[x_init_idx], self.columns[y_init_idx],  source=self.source, view=self.view, fill_alpha=1,fill_color=color_map,line_color=None )
        self.p.legend.click_policy="hide"
        self.s_x = AutocompleteInput(title="x:", value=self.columns[x_init_idx], completions=self.columns,min_characters=1)
        self.s_y = AutocompleteInput(title="y:", value=self.columns[y_init_idx], completions=self.columns,min_characters=1)
        # Attach reaction
        self.s_x.on_change("value", lambda attr, old, new: tag_func(self.s_x, self.r.glyph, 'x', self.p) )
        self.s_y.on_change("value", lambda attr, old, new: tag_func(self.s_y, self.r.glyph, 'y', self.p) )
        # Set default fill color
        if select_color_change:
            self.r.selection_glyph = Circle(fill_alpha=1,fill_color=cur_color, line_color='black')
        self.allow_select = allow_select

    def refresh(self):
        global cur_color
        #print('color list of data: ',self.r.selection_glyph.fill_color)
        #self.r.selection_glyph = Circle(fill_alpha=1,fill_color=None)
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
Figure = FlowPlot(opts, source, view, generic_columns, 'color',"Surface Marker Gating Panel", legend='color')


##############################
##### Buttons Definition #####
##############################

# Show gene list
d_gene = Div(text='Gene/Marker List: '+str(generic_columns))

# Log
log_check = CheckboxGroup(labels=['Log-scaled axis'],active=[])
log_check.on_click(log_cb)

# Change the color of selected parts

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
cat_opt = Select(title='Select Cluster Group:',options=[' '],value=' ')
cat_opt.on_change('value',choose_cat)

# Input name of new category
name = TextInput(title='Input Group Name: ', value='')
name.js_on_change("value", CustomJS(code="""
    console.log('text_input: value=' + this.value, this.toString())
"""))

# Create New Category
new_view = Button(label='Create Group')
new_view.on_click(new_category)

# Rename Category
rename_view = Button(label='Rename Group')
rename_view.on_click(edit_category)

# Delete Category
del_view = Button(label='Delete Group')
del_view.on_click(del_category)

# Select existed categories
#cat_opt = Select(title='Select Category', options=['New Category'], value='New Category')
#cat_opt.on_change("value", choose_cat)

##############################
### Class Functions Button ###


# Input of class name
input_t = TextInput(title='Input Cluster Name: ', value='')


    

# Select of class (use checkbox)
cls_label = ['Unassigned: color=grey, cell_nums=' + str(data_df.shape[0])] # Checkbox label of class
class_checkbox = CheckboxGroup(labels=cls_label, active=[], css_classes=["class_checkbox_label"])
select_cluster = Button(label='Select Cluster')
select_cluster.on_click(show_checked)
show_text = Button(label='Show Color on Checkbox')
show_text.on_click(text_color)


#class_checkbox.on_change('labels',tex_color)
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
add_dots = Button(label='Add to')
add_dots.on_event(ButtonClick,save_cls_button)

# Remove dots from cluster
remove_dots = Button(label='Remove from')
remove_dots.on_click(remove_dot)

# Update cluster
update_cluster = Button(label='Update Cluster')
update_cluster.on_click(update_clus)

# Merge Button
merge_button = Button(label='Merge Cluster')
merge_button.on_click(merge_class)

# Change Color
change_clr_button = Button(label='Change Color')
change_clr_button.on_click(change_color)



# Export Button
export_button = Button(label='Export Results')
export_button.on_click(save_profile)



def axis_cb(attr,old,new):
    global Figure, data_df, opts, source, view 
    Figure.source = source
    Figure.view = view
    Figure.s_x.completions, Figure.s_y.completions = Figure.columns,Figure.columns
    Figure.r.view = view
    Figure.r.data_source = source
    show_color()

# Read category
try:
    category_name = adata.uns['category_dict'].keys()
    cat_opt.options = list(category_name)
    cat_opt.value = list(category_name)[0]
    cate_panel = column(cat_opt)
    ccount = 0
    cls_label = []
    cate = cat_opt.value
    views = list(adata.obsm.keys())
    for view_name in views:
        for i in range(adata.obsm[view_name].shape[1]):
            data_df[view_name+str(i)] = pandas.Series(adata.obsm[view_name][:,i],index=data_df.index)
            data_log[view_name+str(i)] = data_df[view_name+str(i)]
    opts, source, view = define(data_df)
    Figure.columns = list(data_df.columns)
    choose_panel = Select(title='Choose map:',value='generic_columns',options=list(adata.obsm.keys())+['generic_columns'])
    choose_panel.on_change('value',axis_cb)
    curdoc().add_root(choose_panel)
    for i in range(len(adata.uns['category_dict'][cate])):
        l = adata.obs[cate][adata.obs[cate] == str(i)]
        s = str(adata.uns['category_dict'][cate]['class_name'][i]) + ': color=' + str(adata.uns['category_dict'][cate]['color'][i]) + ', cell_nums='+ str(l.shape[0])
        cls_label = np.append(cls_label,s)
        ccount = ccount + l.shape[0]
    show_color()
    cls_label = cls_label + ['no cluster: color=grey, cell_nums=' + str(data_df.shape[0]-ccount)]
    class_checkbox.labels = list(cls_label)
    class_checkbox.active = []
except:
    d = Div(text='No Existed Clusters! ')
    adata.uns['category_dict'] = dict()
    curdoc().add_root(d)



### Layout ###
file_panel = row(upload_button, export_button)
figure_panel = column(Figure.p,d_gene,hide,hide1)
control_panel = column(Figure.s_x, Figure.s_y, log_check, select_color, gate_button, remove_button, showall_button)
class_panel = column(cat_opt,name,new_view,rename_view,del_view, input_t, create_button, show_text,class_checkbox)
edit_panel = column(select_cluster,update_cluster, add_dots,remove_dots,rename_button, change_clr_button, merge_button, del_button)
layout = column(file_panel,row(figure_panel,column(control_panel),class_panel,edit_panel))

# Panel of highlight gene
def change_view(attr,old,new):
    global hl_gene_plot, view
    hl_gene_plot.r.glyph.x = Figure.r.glyph.x
    hl_gene_plot.r.glyph.y = Figure.r.glyph.y
    hl_gene_plot.p.xaxis.axis_label = Figure.p.xaxis.axis_label
    hl_gene_plot.p.yaxis.axis_label = Figure.p.yaxis.axis_label
    view = Figure.r.view
    

def show_colorbar():
    global source1, Figure,hl_color_bar
    updated_color = source.data[hl_input.value]
    #updated_color = (updated_color-hl_bar_map.low)/(hl_bar_map.high - hl_bar_map.low)
    source1.data["hl_gene"] = updated_color
    print(source.data['hl_gene'])
    

def hl_filter():
    global source1, Figure
    if hl_filt.value == 'Gene Expression >':
        source1.selected.indices = list(adata.obs[source.data[hl_input.value] > float(hl_filt_num.value)]['ind'])
    elif hl_filt.value == 'Gene Expression <':
        source1.selected.indices = list(adata.obs[source.data[hl_input.value] < float(hl_filt_num.value)]['ind'])
    else:
        source1.selected.indices = list(adata.obs[source.data[hl_input.value] == float(hl_filt_num.value)]['ind'])
    print(len(source1.selected.indices))
    #new_r = show_colorbar()
    hl_gene_plot.r.selection_glyph = Circle(fill_alpha=1,fill_color='Black')
    #print(source.selected.indices)

def change_select():
    global source
    source.selected.indices = source1.selected.indices


hl_gene_map = log_cmap('hl_gene', cc.kbc[::-1], low=1, high=20)
source1 = ColumnDataSource(data=source.data)
hl_gene_plot = FlowPlot(opts, source1, view, Figure.columns, hl_gene_map, "Highlight Gene Viewing Window", select_color_change = False)
#hl_gene_plot.color_map = log_cmap('hl_gene', cc.kbc[::-1], low=0, high=40)
hl_bar_map = LogColorMapper(palette=cc.kbc[::-1], low=1, high=20)
hl_color_bar = ColorBar(color_mapper=hl_bar_map, label_standoff=8, border_line_color=None)
hl_gene_plot.p.add_layout(hl_color_bar,'right')
hl_input = AutocompleteInput(completions=generic_columns, title="Select Highlight Gene: ", min_characters=1)
hl_button = Button(label="Show Highlight Gene")
hl_button.on_click(show_colorbar)

hl_filt = Select(options=['Gene Expression >','Gene Expression =','Gene Expression <'],value='Gene Expression >')
hl_filt_num = TextInput()
hl_filt_button = Button(label='Filter')
hl_filt_button.on_click(hl_filter)
hl_select = Button(label='Change Select')
hl_select.on_click(change_select)
hl_comfirm = Button(label='Change Selected')
hl_comfirm.on_click(change_select)

control_panel2 = column(hl_input,hl_button,row(hl_filt, hl_filt_num),hl_filt_button,hl_comfirm)
layout2 = row(hl_gene_plot.p, control_panel2)
panel2 = Panel(child=layout2,title='Highlight Gene')


# Panel
panel1 = Panel(child=layout,title='Main View')
#
tab_list = [panel1,panel2]
tabs = Tabs(tabs=tab_list)
tabs.on_change('active',change_view)


curdoc().add_root(tabs)
