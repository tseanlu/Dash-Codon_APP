import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
from dash.dependencies import Input, Output
import pandas as pd
import numpy as np
from decimal import Decimal
import math
import urllib
import sqlite3
import os
from datetime import datetime

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
GENE_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}
GENE_CODE_REVERSE = {}
for x in GENE_CODE:
    symbol = GENE_CODE.get(x)
    group = sum(x.count(char) for char in ['G', 'C', 'g', 'c', 'S', 's'])
    if (not GENE_CODE_REVERSE.get(symbol)):
        new_dic = {}
        new_dic[group] = [x]
        GENE_CODE_REVERSE[symbol] = new_dic
    else:
        if (not GENE_CODE_REVERSE[symbol].get(group)):
            GENE_CODE_REVERSE[symbol][group] = [x]
        else:
            GENE_CODE_REVERSE[symbol][group].append(x)

def add_operation(src_l, des_l, adder):
#    print('add_operation...')
#    print('adder size:',len(adder))
#    print('src_l size:',len(src_l))
#    print('src_l:', src_l)
#    print('des_l size:',len(des_l))
#    print('des_l:',des_l)
    for k in range(len(src_l)):
        for s in range(len(adder)):
            des_l.append(src_l[k] + adder[s])
    return

def GC(seq):
    gc = sum(seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
    return (gc*100.0/len(seq))
def function1(aa_list, CGratio):
    ans=[]
    list_len = str(len(aa_list))
    cg_ratio_max=float(CGratio[1])
    cg_ratio_min=float(CGratio[0])
    max_cg = math.floor(3*Decimal(list_len)*Decimal(str(cg_ratio_max)))
    min_cg = math.ceil(3*Decimal(list_len)*Decimal(str(cg_ratio_min)))
##    print(max_cg)
##    print(min_cg)
    ## start doing DP
    tem = []     ## 0~max_cg for 2 line
    tem1 = []
    DP_list = []                     ## initial value
    for k in range(max_cg+1):
        init = [None]
        tem.append(init)
        init_1 = [None]
        tem1.append(init_1)
    tem1[0][0] = ''

    DP_list.append(tem)
    DP_list.append(tem1)
    aa_len = len(aa_list)
    for i in range(aa_len):       ## total need to do n round
#        print('i:',i)
        now = i%2
        past = now -1
        symbol = aa_list[i]             ## now dealing with symbol
#        print('symbol:',symbol)
        for h in range(len(DP_list[0])):
            for group in GENE_CODE_REVERSE[symbol]:
#                print('group',group)
                if (h-group > -1):
                    if (DP_list[past][h-group][0] != None):
                        if (DP_list[now][h][0] == None):
                            DP_list[now][h].pop()
                        DP_list[now][h].clear()
                        add_operation(DP_list[past][h-group], DP_list[now][h],
                                    GENE_CODE_REVERSE[symbol][group])
    count = 0
    for k in range(min_cg,max_cg+1):
        if( DP_list[(len(aa_list)-1)%2][k][0] == None):
            continue
        elif( DP_list[(len(aa_list)-1)%2][k][0] != None):
            if (len(DP_list[(len(aa_list)-1)%2][k][0])< len(aa_list)*3):
                DP_list[(len(aa_list)-1)%2][k].clear()
                DP_list[(len(aa_list)-1)%2][k].append(None)
                continue
        count += len(DP_list[(len(aa_list)-1)%2][k])
        ans.append(DP_list[(len(aa_list)-1)%2][k])
#        print('answer',DP_list[(len(aa_list)-1)%2][k])
        for string in DP_list[(len(aa_list)-1)%2][k]:
            if ( string == None):
                continue
            if ( not check(string,aa_list)):
                print("something wrong with:", string)
    ans1=[]
    for l in range(len(ans)):
        ans1.append(ans[l][0])
    GC_option=[]
    for t in range(len(ans1)):
        GC_option.append(GC(ans1[t]))
    return (ans1,GC_option)

def check(DNA,symbol):
    for k in range(len(symbol)):
        if (GENE_CODE.get(DNA[3*k:3*k+3]) != symbol[k]) :
#            print(GENE_CODE.get(DNA[3*k:3*k+3]),symbol[k])
            break
            return False
    return True
def download_result(data,codon):
    conn = sqlite3.connect('test1.sqlite.db')
    sql=pd.read_sql('select * from tab', conn)
    if codon in sql['CODON'].values:
        return ("It already exists")
    else :
        d=dict()
        d['AA']=[data]
        d['CODON']=[codon]
        d=pd.DataFrame(d)
        d=sql.append(d)
        d.to_sql('tab', conn, if_exists='replace', index=False)
        return("successfully register!!")
### App
#### Function
def generate_control_card():
    """
    :return: A Div containing controls for graphs.
    """
    return html.Div(
        id="control-card",
        children=[
            html.Br(),
            html.Br(),
            html.H3("Codon-to-Nucleotide Tool"),
            html.H5("Amino Acid Translation"),
            html.H6("By Andy Lu"),
            html.Br(),
            html.Br(),
            html.P("Input Amino Acid Sequence"),
            dcc.Input(id='AA', value='AGADG', type='text',placeholder='Enter amino acid sequence...',),
            html.Br(),
            html.Br(),
        ],
    )
######
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.Div(
            id="banner",
            className="banner",
            children=[html.Img(src=app.get_asset_url("dna_logo2.png"))],
        ),
    # Left column
    html.Div(
            id="left-column",
            className="four columns",
            children=[generate_control_card()],
        ),
    # Right column
    html.Div(
            id="right-column",
            className="eight columns",
            children=[
                # Selective Choices
                html.Div(
                    id="patient_volume_card",
                    children=[
                        html.B("Choose CG Ratio Range"),
                        html.Hr(),
                        html.Br(),
                        dcc.RangeSlider(
                            id='CG_range_slider',
                            min=0,
                            max=1,
                            step=0.01,
                            value=[0.5,0.7],
                            marks={
                                0:{'label':'0%'},
                                0.1:{'label':'10%'},
                                0.2:{'label':'20%'},
                                0.3:{'label':'30%'},
                                0.4:{'label':'40%','style':{'color':'#f50'}},
                                0.5:{'label':'50%'},
                                0.6:{'label':'60%','style':{'color':'#f50'}},
                                0.7:{'label':'70%'},
                                0.8:{'label':'80%'},
                                0.9:{'label':'90%'},
                                1:{'label':'100%'},
                            }),
                        html.Div(id='output_CG_range', style={'color':'blue','margin-top': 30,'fontSize': 15}
                        ),
                        html.Button('Update Choices',id='bttn0'),
                        html.Br(),
                        html.Hr(),
                        html.Br(),
                        html.B("Possible CG Ratio Choices"),
                        html.Br(),
                        dcc.Dropdown(
                            id="answer_dropdown"
                        ),
                        html.Br(),
                        html.Div(
                            id="NTchoice",style={'color':'red','margin-top': 30,'fontSize': 18}
                        ),
                        html.Br(),
                    ],
                ),
                # Results
                html.Hr(),
                html.Br(),
                html.Div(className="four columns",),
                html.Button('Check and Register',id='bttn1'),
                html.Div(id='register_alert'),
            ],
        ),

])
## Call back function
### Update CG value display
@app.callback(
    dash.dependencies.Output('output_CG_range', 'children'),
    [dash.dependencies.Input('CG_range_slider', 'value')])
def update_output(value):
    CGrange='Cureent selected CG range is {}'.format(value)
    return CGrange
### Update Answer Options
@app.callback(
    dash.dependencies.Output('answer_dropdown', 'options'),
    [dash.dependencies.Input(component_id='bttn0',component_property='n_clicks')],
    [dash.dependencies.State(component_id='AA', component_property='value'),
     dash.dependencies.State(component_id='CG_range_slider', component_property='value')])
def update_alert(bttn0,AA,CG_range_slider):
    result,cg_option=function1(AA,CG_range_slider)
    return ([{'label': i, 'value': j} for i,j in zip(cg_option,result)])
### Update Answer Value
@app.callback(
    Output('NTchoice', 'children'),
    [Input('answer_dropdown', 'value')])
def set_cities_value(available_options):
    Ans='Possible Nucleic Acid Sequence {}'.format(available_options)
    return Ans
## Alert Function
@app.callback(
    [dash.dependencies.Output(component_id='register_alert',component_property='children')],
    [dash.dependencies.Input(component_id='bttn1',component_property='n_clicks')],
    [dash.dependencies.State(component_id='AA', component_property='value'),
     dash.dependencies.State(component_id='answer_dropdown', component_property='value')]
)
def update_alert(bttn1,AA,answer_dropdown):
    result=answer_dropdown
    alert = download_result(result,AA)
    time=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    return ([alert+time])



if __name__ == '__main__':
    app.run_server(debug=True)
