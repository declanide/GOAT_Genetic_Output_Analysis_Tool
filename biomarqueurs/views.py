import random                                   # pour importer bokeh server
from bokeh.models import *
from bokeh.models import HoverTool
from bokeh.plotting import *
from bokeh.resources import CDN
from bokeh.embed import components, Resources
from django.shortcuts import render,redirect,HttpResponse,render_to_response
from models import *
from django.contrib.auth.models import User
from django.contrib.auth import authenticate, login, logout
from django.db import connection
from django.template.loader import render_to_string



#--- Panda ---#
import pandas as pd
import pandas.io.sql as psql
import numpy as np


#--- MatPlotLib ---#
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

#--- JSON for Ajax response ---#
import json
from django.http import HttpResponse


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

global threshold                                         #default threshold
global seuil                            #default threshold after transformation
colorManhattan= [( 'blue' ,' gray' ),( 'red' , 'blue' ),( 'black' , 'grey' )]            #for future development possible colors for manhattan

global sorted_data                                           # store dataframe in memory for future query
global phenotypes_found                                     # store significant phenotypes found
global phenotype_selected                                   # phenotype been displayed
global chr_max_pos                                          # STORE MAX position of each chromosome
global chr_min_pos                                          # STORE MIN position of each chromosome
global String_selection
global selection
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                            #Authentification page
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
def index ( request ) :

    if request.user.is_authenticated():
        return render(request, 'home.html', {'content':'modules/home', 'username': request.user})
    else:
        if request.method == 'POST':                  #Si methodes est post
            username = request.POST['username']       #demander user name et mot de passe
            pwd = request.POST['password']
            user = authenticate(username=username, password=pwd)  #authentifie l'utilisateur
            if user is not None:
                if user.is_active:
                    login( request, user)                     #retourner le nom de l'utilisateur
                    request.session.set_expiry(300)
                    return render(request, 'home.html', {'content':'modules/home', 'username': request.user})
                else:
                    return render(request, "modules/login.html")
            else:
                print ("erreur de login")
                return render(request, "modules/login.html",{'error':'Erreur de login'})
        else:
            return render(request, "modules/login.html")


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                            #Significant phenotypes by querying database and verifying
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

def connect( q, t, que ) :
    cursor = connection.cursor()
    phenotypes_found = []
    if len(q) == 1 and que == "rs":
            cursor.execute("select p.nom From marqueurs m join assoc a on a.rs_id_assoc = m.nom JOIN experiment xp on a.experiment = xp.idexperiment join phenotypes p on xp.phenotype=p.idphenotypes where pvalue_assoc <"+t+" and m.nom in ('"+q[0]+"') order by a.pvalue_assoc ASC ")
    if len(q) != 1 and que == "rs":
        cursor.execute("select p.nom From marqueurs m join assoc a on a.rs_id_assoc = m.nom JOIN experiment xp on a.experiment = xp.idexperiment join phenotypes p on xp.phenotype=p.idphenotypes where pvalue_assoc <"+t+" and m.nom in "+str(q)+" order by a.pvalue_assoc ASC ")
    if len(q) == 1 and que == "gene":
            cursor.execute("select DISTINCT p.nom From marqueurs m join assoc a on a.rs_id_assoc = m.nom JOIN experiment xp on a.experiment = xp.idexperiment join phenotypes p on xp.phenotype=p.idphenotypes where pvalue_assoc <"+t+" and (m.gene in ('"+q[0]+"') or m.gene_before in ('"+q[0]+"') or m.gene_after in ('"+q[0]+"') ) order by a.pvalue_assoc ASC ")
    if len(q) != 1 and que == "gene":
            cursor.execute("select DISTINCT p.nom From marqueurs m join assoc a on a.rs_id_assoc = m.nom JOIN experiment xp on a.experiment = xp.idexperiment join phenotypes p on xp.phenotype=p.idphenotypes where pvalue_assoc <"+t+" and (m.gene in "+str(q)+" or m.gene_before in "+str(q)+" or m.gene_after in "+str(q)+"  ) order by a.pvalue_assoc ASC ")

    solution = cursor.fetchall()                        #retourne une tuple avec les noms des phenotypes
    if len(solution) == 0:
        return "THERE ARE NO SIGNIFICANT SNP'S FOR THE THRESHOLD SELECTED. \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!"   # if no significant phenotype is found
    else:
        for n in range(len(solution)):                      #store phenotypes
            phenotypes_found.append(solution[n][0])
    return phenotypes_found
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                            #Create Manhattan in matplotlib according to phenotype selected
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

def manhattan( s, t ):
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    global chr_max_pos
    global chr_min_pos                          #GLOBAL variables call
    max_v=int(s.log10.max())
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    df = s.groupby('chromosome')                                    # CUTS GWAS DATA TO COMPUTE GWAS'S POSITIONS
    c = 0
    dictio = {}
    chr_max_pos=[]
    chr_min_pos=[]
    manhattan_max_pos = [0]
    df_concat=[]
    for name, group in df:
        dictio[str(name)] = group
        chr_min_pos.append(int(dictio[str(name)]['pos'][0:1]))                                #position min chromosome
        manhattan_max_pos.append(int(dictio[str(name)]['pos'][len(dictio[str(name)]['pos']) - 1:len(dictio[str(name)]['pos'])]) + manhattan_max_pos[c])
        chr_max_pos.append(int(dictio[str(name)]['pos'][len(dictio[str(name)]['pos']) - 1:len(dictio[str(name)]['pos'])]))              #position max chromosome
        dictio[str(name)]['position']=group['pos'] + manhattan_max_pos[c]
        c += 1
        df_concat.append(dictio[str(name)])
    c = 0
    f_data=pd.concat(df_concat)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    millieu = []                                                                  # HALF OF CHROMOSOME POSITIONS TO STORE

    for n in range(1,len(manhattan_max_pos)):
        millieu.append(manhattan_max_pos[n-1] + ((manhattan_max_pos[n]-manhattan_max_pos[n-1])//2))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                # CREATE MANHATTAN PLOT
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    source = ColumnDataSource(f_data)                    # SOURCE DATA FOR BOKEH PLOT
    TOOLS = [HoverTool(tooltips=[("SNP", "@rs_id_assoc"),("P-value","@pvalue_assoc")]), TapTool()]
    plot = figure(tools=["box_zoom", "box_select", "wheel_zoom", TOOLS[0], "crosshair", "reset", "save", TOOLS[1]] , x_axis_label='Chromosomes', y_axis_label='-log10(p)', plot_width=900, plot_height=500,y_range=(2.0,max_v+3))
    plot.circle('position', 'odd', source=source, size=3)
    plot.circle('position', 'even', source=source, size=3, color="black")
    plot.ray(x=[0],y=[t],length=0,angle=0, color='red')
    plot.ray(x=[0],y=[t],length=0, angle=np.pi,color='red')
    plot.axis.major_label_text_font_size='0pt'
    plot.text(millieu, [2.75]*(len(millieu)), text=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"], text_color='black', text_align='center', text_font_size='8pt', text_font_style='bold')
    plot.text(millieu[(len(millieu)/2)-4],[2.25],text=["Chromosomes"], text_color='black', text_align='center', text_font_size='10pt', text_font_style='bold')
    plot.xaxis.major_tick_line_color = None
    plot.xaxis.minor_tick_line_color = None
    plot.xaxis.visible=None
    plot.grid.grid_line_color = None            # TAKE GRID LINES OFF THE GRAPH
    graph, div1 = components(plot, CDN)
    return graph, div1

    '''
    plt.figure(figsize=(14,max_v+2))
    for n in dictio:
        if int(n) % 2 == 0:
            plt.plot(dictio[n]['position'],dictio[n]['log10_x'],".",mec="black", mfc="black",ms=4)
        else:
            plt.plot(dictio[n]['position'],dictio[n]['log10_x'],".",mec="gray", mfc="gray",ms=4)
    for n in dictio:
            plt.plot(dictio[n]['position'],dictio[n]['log10_y'],".",mec="red", mfc="red",ms=4)
    plt.axhline(y=(-np.log10(t)), xmin=0, xmax=1, color="black", ls="--")
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(p)")
    plt.ylim(0,max_v+2)
    man_min=manhattan_max_pos[0]-60000000
    man_max=manhattan_max_pos[len(manhattan_max_pos)-1]+60000000
    plt.xlim(man_min,man_max)
    plt.grid.grid_line_color = None
    plt.xticks(millieu, ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"], fontsize=8)
    '''
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    # SAVE MANHATTAN'S PICTURE
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    #plt.savefig("/Your/Path/Static/", dpi=300, facecolor='w', edgecolor='w',
    #    orientation='portrait', papertype=None, format=None,
    #    transparent=False, bbox_inches="tight", pad_inches=0.1,
    #    frameon=None)

    #return "manhattan.png"  # RETURN SAVED PICTURE
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                    #Create area Selection graph with table
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

def areaSelection(c,min,max,p,s):
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    global threshold                                # GLOBAL VARIABLE ACESS
    global seuil
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    sql="select distinct a.rs_id_assoc , a.chromosome,a.pos,a.info_assoc,a.pvalue_assoc,a.allele_A,a.allele_B,a.cohort_AA,a.cohort_BB,a.beta_assoc,a.maf, a.all_OR,xp.covariates,p.Risk_on_rise,dat.name from assoc a join experiment xp on a.experiment=xp.idexperiment join dataset dat on xp.dataset=dat.iddataset join phenotypes  p  on xp.phenotype=p.idphenotypes where p.nom='"+phenotype_selected+"' and a.chromosome="+c
    chrarea=pd.read_sql(sql,connection)
    #chrarea = s[s['chromosome'] == c]                      #SELECT CHROMOSOME
    area = chrarea[(chrarea['pos']>min) & (chrarea['pos']<max)]       #SELECT INTERVAL OF POSITIONS
    area["risk_allele"]=np.where((area['beta_assoc'] > 0 & area['Risk_on_rise']==1)|(area['beta_assoc'] < 0 & area['Risk_on_rise']==0), area['allele_B'], area['allele_A'])
    area["risk_af"]=np.where((area['beta_assoc'] > 0 & area['Risk_on_rise']==1)|(area['beta_assoc'] < 0 & area['Risk_on_rise']==0), ((2*area["cohort_BB"])+area["cohort_AB"])/((2*area["cohort_AA"])+(2*area["cohort_AB"])+(2*area["cohort_BB"])), ((2*area["cohort_AA"])+area["cohort_AB"])/((2*area["cohort_AA"])+(2*area["cohort_AB"])+(2*area["cohort_BB"])))
    area["risk_allele_beta"]=np.where((area['beta_assoc'] > 0 & area['Risk_on_rise']==1)|(area['beta_assoc'] < 0 & area['Risk_on_rise']==0), area['beta_assoc'], area['beta_assoc']*-1)
    long = len(area['rs_id_assoc'])                                   # SIZE OF DATASET
    sql = "select m.nom,m.end_gen_after,m.end_gen,m.start_gen,m.end_gen_before,m.func,m.position,m.start_gen_after,m.start_gen_before, m.observed  from marqueurs m where m.chromosome="+str(c)+" and position between "+str(min)+" and "+str(max)
    df = pd.read_sql(sql, connection)                               # SQL QUERY SAVED IN PANDAS DATAFRAME
    d1 = df.rename(columns = {'nom':'rs_id_assoc'})                   # RENAME COLUMN
    #sql2 = "select ex.idexperiment, dat.name from experiment ex join dataset dat on ex.dataset=dat.iddataset where ex.idexperiment="+str(int(area['experiment'][:1]))
    #df2 = pd.read_sql(sql2, connection)                               # SQL QUERY SAVED IN PANDAS DATAFRAME
    #d = df2.rename(columns = {'idexperiment':'experiment'})           # RENAME COLUMNS
    string_name = "CHROMOSOME :"+str(c)+" \nBETWEEN POSITIONS "+str(min)+" AND "+str(max)                # TITLE FOR GRAPH
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    pd.set_option('display.max_colwidth',-1)                        #important to make links appear in pandas dataframe
    data = pd.merge(area,d1,on='rs_id_assoc',how='outer')             # MERGE ALL DATAFRAME
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    c_data = data[0:long]                                        # SELECT ONLY DATA FROM MANHATTAN
    del c_data["pos"]
    sort_df = c_data.sort(['position'])              #sort BY position
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                        # CREATING TABLES COLUMNS IN ADDING INFORMATION
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#    pd.set_option('display.max_colwidth', -1)

    #sort_df.drop_duplicates(subset='rs_id_assoc', inplace=True,take_last=True)  # TAKE DUPLICATES OFF
    sort_df["log10"] = -np.log10(sort_df.pvalue_assoc)                     #transformation
    sort_df = sort_df.sort(['pvalue_assoc'])                                  # SORT BY P-VALUE
    max_pvalue = int(sort_df.log10[0:1])                                      # GET MINIMUM P-VALUE
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                    #NEW COLUMNS AND RENAMING
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    sort_df['imp'] = np.where(sort_df['info_assoc']==1, sort_df['log10'], 'NaN')
    sort_df['Imputed'] = np.where(sort_df['info_assoc']==1, True, False)
    c_data['Imputed'] = np.where(c_data['info_assoc']==1, True, False)
    sort_df['interest'] = np.where(sort_df['log10']>=(-np.log10(0.00000005)), sort_df['log10'], 'NaN')

    c_data['dbSNP'] = '<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs='+c_data.rs_id_assoc+'"target="_blank">dbSNP</a>'
    c_data['GWAS_Catalog'] = '<p><a href="http://www.genome.gov/gwastudies/index.cfm?snp='+c_data.rs_id_assoc+'"target="_blank">By rsID</a></p><p><a href="http://www.genome.gov/gwastudies/index.cfm?gene='+c_data['gene_before']+'"target="_blank">By gene before</a><p><a href="http://www.genome.gov/gwastudies/index.cfm?gene='+c_data['gene']+'"target="_blank">By gene</a></p><p><a href="http://www.genome.gov/gwastudies/index.cfm?gene='+c_data['gene_after']+'"target="_blank">By gene after</a>'
    c_data['Genecards'] = '<p><a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene='+c_data['gene_before']+'"target="_blank">By gene before</a></p><p><a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene='+c_data['gene']+'"target="_blank">By gene</a></p><p><a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene='+c_data['gene_after']+'"target="_blank">By gene after</a>'

    c_data = c_data.rename(columns = {'rs_id_assoc': 'rs_ID', 'position': 'Position', 'chromosome': 'Chr','func':'Region', 'gene_before': 'GeneBefore', 'pvalue_assoc': 'P-value', 'allele_A': 'Allele A', 'allele_B': 'Allele B', 'gene': 'Gene', 'name': 'Cohort', 'gene_after': 'GeneAfter'})
    c_data.sort_index(inplace=True)
    c_data["Allele\nmineure"] = np.where(c_data['cohort_AA'] > c_data['cohort_BB'], c_data['Allele B'], c_data['Allele A'])

    DF = c_data[['rs_ID', 'Chr', 'Position', 'GeneBefore', 'Gene', 'GeneAfter', 'Region','P-value', 'Imputed', 'Allele A', 'Allele B', 'Allele\nmineure','risk_allele','risk_af','risk_allele_beta', 'Cohort', 'dbSNP', 'GWAS_Catalog', 'Genecards']]
    #DF = c_data[['rs_ID', 'Chr', 'Position', 'Gene before', 'Gene', 'Gene after', 'Region','P-value   ', 'Imputed', 'Allele A', 'Allele B', 'Allele\nmineure', 'Cohort', 'dbSNP', 'GWAS_Catalog', 'Genecards']]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    #file = DF.to_csv("/Users/beatrizkanzki/PycharmProjects/GOAT_Genetic_Ouput_Analysis_Tool/static/GeneSelection.csv", cols=['rs_id_assoc','Chr','Position','Gene before','Gene','Gene after','Region','P-value','Allele A','Allele B','Allele\nmineure','Cohort']) #SELECT COLUMNS FOR CSV FILE
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    source = ColumnDataSource(sort_df)                    # SOURCE DATA FOR BOKEH PLOT
    TOOLS = [HoverTool(tooltips=[("SNP", "@rs_id_assoc"), ("Gene", "@gene"),("P-value","@pvalue_assoc"),("Region","@func")]), TapTool()]
    stringLegend = "pvalue < "+str(threshold)
    plot = figure(tools=["box_zoom", "box_select", "wheel_zoom", TOOLS[0], "crosshair", "reset", "save", TOOLS[1]] , x_axis_label='Position', y_axis_label='-log10(p)', plot_width=900, plot_height=500, y_range=(-3.2,max_pvalue+3))
    plot.circle('position', 'log10', source=source, size=3, legend='Genotyped')
    plot.square('position', 'imp', source=source, size=3, color="olive", legend='Imputed')
    plot.circle('position', 'interest', source=source, size=3, color="red", legend=stringLegend)


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    #html_table = DF.to_html(index=False,escape=False)         # PUT TABLE IN HTML
    sort_df = sort_df.sort(['position'])                      # SORT POSITIONS
    sort_df.drop_duplicates(subset=('gene'), inplace=True,take_last=True) # TAKE GENE NAME DUPLICATES OFF


    jsontable = DF.to_json(orient='records')#json table
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                            # GENE POSITION PLOTTING
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    sort_df['ligne1'] = sort_df.start_gen[0:len(sort_df):10]
    sort_df['Fligne1'] = sort_df.end_gen[0:len(sort_df):10]

    sort_df['ligne2'] = sort_df.start_gen[1:len(sort_df):10]
    sort_df['Fligne2'] = sort_df.end_gen[1:len(sort_df):10]

    sort_df['ligne3'] = sort_df.start_gen[2:len(sort_df):10]
    sort_df['Fligne3'] = sort_df.end_gen[2:len(sort_df):10]

    sort_df['ligne4'] = sort_df.start_gen[3:len(sort_df):10]
    sort_df['Fligne4'] = sort_df.end_gen[3:len(sort_df):10]

    sort_df['ligne5'] = sort_df.start_gen[4:len(sort_df):10]
    sort_df['Fligne5'] = sort_df.end_gen[4:len(sort_df):10]

    sort_df['ligne6'] = sort_df.start_gen[5:len(sort_df):10]
    sort_df['Fligne6'] = sort_df.end_gen[5:len(sort_df):10]

    sort_df['ligne7'] = sort_df.start_gen[6:len(sort_df):10]
    sort_df['Fligne7'] = sort_df.end_gen[6:len(sort_df):10]

    sort_df['ligne8'] = sort_df.start_gen[7:len(sort_df):10]
    sort_df['Fligne8'] = sort_df.end_gen[7:len(sort_df):10]

    sort_df['ligne9'] = sort_df.start_gen[8:len(sort_df):10]
    sort_df['Fligne9'] = sort_df.end_gen[8:len(sort_df):10]

    sort_df['ligne10'] = sort_df.start_gen[9:len(sort_df):10]
    sort_df['Fligne10'] = sort_df.end_gen[9:len(sort_df):10]

    plot.segment(sort_df.ligne1, [-0.30]*(len(sort_df)), sort_df.Fligne1,[-0.30]*(len(sort_df)), line_width=6, line_color="#8b4513",)          #ligne 1
    plot.text(sort_df.ligne1+((sort_df.Fligne1-sort_df.ligne1)/2), [-0.25]*(len(sort_df)), text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt', text_font_style='bold')

    plot.segment(sort_df.ligne2, [-0.55]*(len(sort_df)), sort_df.Fligne2,[-0.55]*(len(sort_df)), line_width=6, line_color="#8b4513",)           #ligne 2
    plot.text(sort_df.ligne2+((sort_df.Fligne2-sort_df.ligne2)/2), [-0.50]*(len(sort_df)), text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt', text_font_style='bold')

    plot.segment(sort_df.ligne3, [-0.85]*(len(sort_df)), sort_df.Fligne3,[-0.85]*(len(sort_df)), line_width=6, line_color="#8b4513",)          #ligne 3
    plot.text(sort_df.ligne3+((sort_df.Fligne3-sort_df.ligne3)/2), [-0.80]*(len(sort_df)), text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt', text_font_style='bold')

    plot.segment(sort_df.ligne4, [-1.15]*(len(sort_df)), sort_df.Fligne4,[-1.15]*(len(sort_df)), line_width=6, line_color="#8b4513",)          #ligne 4
    plot.text(sort_df.ligne4+((sort_df.Fligne4-sort_df.ligne4)/2), [-1.10]*(len(sort_df)), text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt', text_font_style='bold')

    plot.segment(sort_df.ligne10,[-1.45]*(len(sort_df)), sort_df.Fligne10,[-1.45]*(len(sort_df)), line_width=6, line_color="#8b4513",)          #ligne 5
    plot.text(sort_df.ligne10+((sort_df.Fligne10-sort_df.ligne10)/2), [-1.40]*(len(sort_df)), text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt', text_font_style='bold')

    plot.segment(sort_df.ligne9, [-1.75]*(len(sort_df)), sort_df.Fligne9,[-1.75]*(len(sort_df)), line_width=6, line_color="#8b4513",)          #ligne 6
    plot.text(sort_df.ligne9+((sort_df.Fligne9-sort_df.ligne9)/2), [-1.70]*(len(sort_df)), text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt', text_font_style='bold')

    plot.segment(sort_df.ligne8, [-2.05]*(len(sort_df)), sort_df.Fligne8, [-2.05]*(len(sort_df)), line_width=6, line_color="#8b4513",)          #ligne 7
    plot.text(sort_df.ligne8+((sort_df.Fligne8-sort_df.ligne8)/2), [-2.00]*(len(sort_df)), text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt', text_font_style='bold')

    plot.segment(sort_df.ligne7,[-2.35]*(len(sort_df)), sort_df.Fligne7,[-2.35]*(len(sort_df)), line_width=6, line_color="#8b4513",)          #ligne 8
    plot.text(sort_df.ligne7+((sort_df.Fligne7-sort_df.ligne7)/2), [-2.30]*(len(sort_df)), text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt', text_font_style='bold')

    plot.segment(sort_df.ligne6,[-2.65]*(len(sort_df)), sort_df.Fligne6,[-2.65]*(len(sort_df)), line_width=6, line_color="#8b4513",)          #ligne 9
    plot.text(sort_df.ligne6+((sort_df.Fligne6-sort_df.ligne6)/2), [-2.60]*(len(sort_df)),text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt',text_font_style='bold')

    plot.segment(sort_df.ligne5,[-2.95]*(len(sort_df)), sort_df.Fligne5,[-2.95]*(len(sort_df)), line_width=6, line_color="#8b4513",)          #ligne 10
    plot.text(sort_df.ligne5+((sort_df.Fligne5-sort_df.ligne5)/2), [-2.90]*(len(sort_df)), text=sort_df.gene, text_color='black', text_align='center', text_font_size='5pt', text_font_style='bold')

    plot.grid.grid_line_color = None            # TAKE GRID LINES OFF THE GRAPH
    graph, div1 = components(plot, CDN)
    return graph, div1, jsontable,string_name
    #return render(request,'modules/areaSelection.html',{'username':request.user,'phenotype':phenotype_selected,'pheno': phenotypes_found, 'json':sorted_sorted_data_json})
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

                                                        #MAIN FUNCTION THAT VERIFY WHICH BUTTON HAS BEEN CLICKED
                                                                 #GENERATES GRAPH AND TABLES

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


def SGene(request):
    request.session.set_expiry(0)

    
    global phenotypes_found
    global phenotype_selected
    global sorted_data                              #GLOBAL VARIABLES TO BE USED
    global threshold
    global String_selection
    global selection
    global seuil

    if request.method == 'POST':

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                #Gene Query button clicked
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


        if 'rs_ID' in request.POST:




            rs = request.POST['rs_ID']
            g = request.POST['gene']
            t = request.POST['threshold']

            selection=request.POST.getlist('col')
            String_selection=""
            if(len(selection))!=0:
                for n in selection:
                    String_selection += ", a."+str(n)


            if len(t) == 0:
                    t = 0.001
                    threshold = t




#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                            #USER INPUT VERIFICATION
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


            if float(t) > 1:
                #return HttpResponse("YOUR THRESHOLD IS INVALID, PLEASE SELECT A THRESHOLD UNDER 1", status=400)
                return render(request, 'modules/error.html', {'error': "YOUR THRESHOLD IS INVALID, \n PLEASE SELECT A THRESHOLD UNDER 1!"})

            if (len(rs) != 0 and len(g) != 0)or(len(rs) == 0 and len(g) == 0):                             #Si recherche avec les deux champs pas permettre
                #return HttpResponse("YOU CAN ONLY SEARCH BY RS_ID OR GENE NOT BOTH,\n PLEASE TRY AGAIN!", status=400)
                return render(request, 'modules/error.html', {'error': "YOU CAN ONLY SEARCH BY RS_ID OR GENE NOT BOTH,\n PLEASE TRY AGAIN!"})


            if len(rs) != 0 and rs[0:2] == 'rs' and ',' not in rs:
                query = (str(rs))

            if rs != 0 and g == 0:
                for n in range(len(query)):
                    if query[n][0:2] != 'rs' and r'rs\d.+' not in rs:
                        #return HttpResponse("ONE OR MORE OF YOUR RS_ID'S ARE INVALID\n PLEASE ENTER A VALID RS_ID!", status=400)
                        return render(request, 'modules/error.html', {'error':"ONE OR MORE OF YOUR RS_ID'S ARE INVALID\n PLEASE ENTER A VALID RS_ID!"})
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                                            #VERIFICATION OF PONCTUATION IN QUERY
                                                                            # VERIFY SIGNIFICANT PHENOTYPES
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
            if ',' not in rs and ' ' not in rs and ';' not in rs and len(rs) != 0:
                q = "rs"
                query = (str(rs),)
                ph = connect(query, str(t), q)
                if "SIGNIFICANT" in ph:
                    #return HttpResponse("THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!", status=400)
                    return render(request, 'modules/error.html', {'error': "THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!"})
                else:
                    phenotypes_found=ph
                    phenotype_selected=phenotypes_found[0]
            if ',' not in g and ' ' not in g and ';' not in g and len(g) !=0:
                q = "gene"
                query = (str(g),)
                ph = connect(query, str(t), q)
                if "SIGNIFICANT" in ph:
                    #return HttpResponse("THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!", status=400)
                    return render(request, 'modules/error.html', {'error': "THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!"})
                else:
                    phenotypes_found = ph
                    phenotype_selected = phenotypes_found[0]


            else:
                if ',' in rs:
                    query = str(rs).split(',')                                 #Separe la string d'apres le separateur
                    query = tuple(query)
                    q = "rs"
                    ph = connect(query,str(t),q)
                    if "SIGNIFICANT" in ph:
                        #return HttpResponse("THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!", status=400)
                        return render(request, 'modules/error.html', {'error':"THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!"})
                    else:
                        phenotypes_found = ph
                        phenotype_selected = phenotypes_found[0]
                if ' ' in rs and ','not in rs:
                    query = str(rs).split()
                    query = tuple(query)
                    q = "rs"
                    ph = connect(query, str(t), q)
                    if "SIGNIFICANT" in ph:
                        #return HttpResponse("THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!", status=400)
                        return render(request, 'modules/error.html', {'error': "THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!"})
                    else:
                        phenotypes_found = ph
                        phenotype_selected = phenotypes_found[0]
                if ';' in rs:
                    query = str(rs).split(';')
                    query = tuple(query)
                    q = "rs"
                    ph = connect(query,str(t),q)
                    if "SIGNIFICANT" in ph:
                        #return HttpResponse("THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!", status=400)
                        return render(request, 'modules/error.html', {'error': "THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!"})
                    else:
                        phenotypes_found = ph
                        phenotype_selected = phenotypes_found[0]
                if ',' in g:
                    query = str(g).split(',')
                    query = tuple(query)                                 #Separe la string d'apres le separateur
                    q = "gene"
                    ph = connect(query,str(t),q)
                    if "SIGNIFICANT" in ph:
                        #return HttpResponse("THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!", status=400)
                        return render(request, 'modules/error.html', {'error': "THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!"})
                    else:
                        phenotypes_found = ph
                        phenotype_selected = phenotypes_found[0]
                if ' ' in g and ','not in g:
                    query = str(g).split()
                    query = tuple(query)                                      #USER INPUT VERIFICATION
#---------------------
                    q = "gene"
                    ph = connect(query, str(t), q)
                    if "SIGNIFICANT" in ph:
                        #return HttpResponse("THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!", status=400)
                        return render(request, 'modules/error.html', {'error': "THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!"})
                    else:
                        phenotypes_found = ph
                        phenotype_selected = phenotypes_found[0]
                if ';' in g:
                    query = str(g).split(';')
                    query = tuple(query)
                    q = "gene"
                    ph = connect(query, str(t), q)
                    if "SIGNIFICANT" in ph:
                        #return HttpResponse("THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!", status=400)
                        return render(request, 'modules/error.html', {'error': "THERE ARE NO SIGNIFICANT SNP'S WITH THE THRESHOLD SELECTED FOR "+str(query)+". \nPLEASE CHANGE THE THRESHOLD FOR BETTER RESULTS!"})
                    else:
                        phenotypes_found = ph
                        phenotype_selected = phenotypes_found[0]

                if ',' not in rs and ' ' not in rs and ';' not in rs and len(rs)!=0 and r"rs[0-9].+" in rs:
                    #return HttpResponse("PLEASE SEPARATE YOUR QUERY BY \",\" OR \";\" OR ONE SPACE\"", status=400)
                    return render(request, 'modules/error.html', {'error': "PLEASE SEPARATE YOUR QUERY BY \",\" OR \";\" OR ONE SPACE\""})
                if ',' not in g and ' ' not in g  and ';' not in g and r"[\P{P}-]" in g and len(g)!=0:
                    #return HttpResponse("PLEASE SEPARATE YOUR QUERY BY \",\" OR \";\" OR ONE SPACE\"", status=400)
                    return render(request, 'modules/error.html', {'error': "PLEASE SEPARATE YOUR QUERY BY \",\" OR \";\" OR ONE SPACE\""})

            threshold = float(t)




#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                                    #View significant phenotype is clicked
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

        if 'view' in request.POST:
            phe = request.POST['view']
            phenotype_selected = phe


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                #SI LE BOUTON AREA SELECTION EST CLIQUEE
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
        if 'chromosome' in request.POST:
            chromosome = int(request.POST['chromosome'])
            position = int(request.POST['position'])
	    global chr_min_pos
	    global chr_max_pos
            global phenotypes_found
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                        # VERIFY USER INPUT
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
            if len(phenotypes_found) == 0:
                return render(request,'modules/error.html', {'error': 'YOU HAVE TO SEARCH FOR A GENE OR RS_ID \n BEFORE ACCESSING \"AREA SELECTION\"'})
            if chromosome >22 and chromosome <0 :
                return render(request, 'modules/error.html', {'error': "CHROMOSOME NUMBER MUST BE BETWEEN 1 AND 22 ONLY NOT "+str(chromosome)+""})
            if position < 0:
                return render(request, 'modules/error.html', {'error': "POSITION NUMBER MUST BE A POSITIVE NUMBER, SO "+str(chromosome)+" IS NOT VALID"})
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                    #VERIFY IF POSITION SELECTED IS WITHIN CHROMOSOME SELECTED
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
            position_max = position+ 1500000
            position_min = position-1500000
            if position_min < chr_min_pos[chromosome-1]:
                position_min = chr_min_pos[chromosome-1]
                position_max = position_min+3000000
            if position_max > chr_max_pos[chromosome-1]:
                position_max = chr_max_pos[chromosome-1]
                position_min = position_max-3000000
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
            graph, div1, jsonData,name = areaSelection(chromosome, position_min, position_max, str(phenotype_selected), sorted_data)
            return render(request, 'modules/areaSelection.html', {'phenotype': phenotype_selected,'name':name, 'graph': graph, 'pheno': phenotypes_found, 'div1':div1, 'link': "<a href=\"/static/GeneSelection.csv\">Download CSV File</a>", 'jsonData': jsonData, 'user': request.user})




#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                                    #GATHER INFORMATION FOR MANHATTAN AND MORE
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

        #extra_column = ", a.cohort_AB, a.cases_AA, a.cases_AB, a.cases_BB, a.allele_A_freq, a.allele_B_freq, a.cases_maf, a.controls_hwe, a.miss_dat_prop, a.cohort_hwe, a.cases_hwe, a.het_OR, a.het_OR_lower, a.het_OR_upper, a.all_OR, a.all_OR_lower, a.all_OR_upper"
        extra_column = ""
        #sql = "select distinct a.rs_id_assoc ,a.experiment, a.chromosome,a.pos,a.info_assoc,a.pvalue_assoc,a.allele_A,a.allele_B,a.cohort_AA,a.cohort_BB,a.beta_assoc,a.maf "+extra_column+", xp.covariates from assoc a join experiment xp on a.experiment=xp.idexperiment join phenotypes  p  on xp.phenotype=p.idphenotypes where p.nom='"+phenotype_selected+"' limit 100000"
        #sql = "select distinct a.rs_id_assoc ,a.experiment, a.chromosome,a.pos,a.info_assoc,a.pvalue_assoc,a.allele_A,a.allele_B,a.cohort_AA,a.cohort_BB,a.beta_assoc,a.maf "+extra_column+", xp.covariates from assoc a join experiment xp on a.experiment=xp.idexperiment join phenotypes  p  on xp.phenotype=p.idphenotypes where p.nom='"+phenotype_selected+"'"
        #sql="select distinct a.rs_id_assoc ,a.experiment, a.chromosome,a.pos,a.info_assoc,a.pvalue_assoc,a.allele_A,a.allele_B,a.cohort_AA,a.cohort_BB,a.beta_assoc,a.maf, a.all_OR "+extra_column+",xp.covariates,m.gene,m.gene_before,m.gene_after, case when a.beta_assoc>0 and p.Risk_on_rise is true then a.allele_B when a.beta_assoc<0 and p.Risk_on_rise is true then a.allele_A when a.beta_assoc>0 and p.Risk_on_rise is false then a.allele_A when a.beta_assoc<0 and p.Risk_on_rise is false then a.allele_B end as risk_allele, case when a.beta_assoc>0 and p.Risk_on_rise is true then ( (2*a.cohort_BB) + a.cohort_AB ) / ((2*a.cohort_AA) + (2*a.cohort_AB) +(2*a.cohort_BB) ) when a.beta_assoc<0 and p.Risk_on_rise is true then ( (2*a.cohort_AA) + a.cohort_AB ) / ((2*a.cohort_AA) + (2*a.cohort_AB) +(2*a.cohort_BB) ) when a.beta_assoc>0 and p.Risk_on_rise is false then ( (2*a.cohort_AA) + a.cohort_AB ) / ((2*a.cohort_AA) + (2*a.cohort_AB) +(2*a.cohort_BB) ) when a.beta_assoc<0 and p.Risk_on_rise is false then ( (2*a.cohort_BB) + a.cohort_AB ) / ((2*a.cohort_AA) + (2*a.cohort_AB) +(2*a.cohort_BB) ) end as risk_af, case when a.beta_assoc>0 and p.Risk_on_rise is true then a.beta_assoc when a.beta_assoc<0 and p.Risk_on_rise is true then a.beta_assoc*(-1) when a.beta_assoc>0 and p.Risk_on_rise is false then a.beta_assoc*(-1) when a.beta_assoc<0 and p.Risk_on_rise is false then a.beta_assoc end as risk_allele_beta from assoc a join experiment xp on a.experiment=xp.idexperiment join phenotypes  p  on xp.phenotype=p.idphenotypes join marqueurs m on a.rs_id_assoc=m.nom where p.nom='"+phenotype_selected+"'"
        sql="select distinct a.rs_id_assoc, a.chromosome,a.pos,a.info_assoc,a.pvalue_assoc,a.allele_A,a.allele_B,a.cohort_AA,a.cohort_BB,a.beta_assoc,a.maf, a.all_OR,xp.covariates,m.gene,m.gene_before,m.gene_after from assoc a join experiment xp on a.experiment=xp.idexperiment join phenotypes  p  on xp.phenotype=p.idphenotypes join marqueurs m on a.rs_id_assoc=m.nom where p.nom='"+phenotype_selected+"' and a.pvalue_assoc<=0.001"
        sorted_data = pd.read_sql(sql, connection)

        sorted_data.drop_duplicates(subset='rs_id_assoc', inplace=True,take_last=True)
        sorted_data["log10"] = -np.log10(sorted_data.pvalue_assoc)               #ADD COLUMN LOG10
        sorted_data = sorted_data.sort(['chromosome', 'pos'])
        #sorted_sorted_data = sorted_data[sorted_data.log10 >=(-np.log10(float(threshold)))]
        #sorted_sorted_data.sort(['chromosome', 'pvalue_assoc'],ascending=[True,True])
        sorted_data['even']=np.where(sorted_data['chromosome'] %2==0,sorted_data['log10'] , 'NaN')
        sorted_data["odd"]=np.where(sorted_data['chromosome'] %2!=0,sorted_data['log10'] , 'NaN')
        
        #h_data = pd.merge(sorted_data,sorted_data, on='rs_id_assoc', how='outer')
        col=['rs_id_assoc', 'chromosome', 'pos', 'pvalue_assoc', 'allele_A', 'allele_B', 'covariates', 'cohort_BB', 'cohort_AA', 'beta_assoc', 'maf']+selection

        DF = sorted_data[col]         # select COLUMN OF DATAFRAME TO SAVE IN CSV



        #DF.to_csv("/Your/Path/Static/", cols=col)
            #"a.av_max_post_call,a.info_imp,a.cohort_AB,a.cohort_NULL,a.cases_AA,a.cases_AB,a.cases_BB,a.controls_AA,a.controls_AB,a.controls_BB,a.allele_A_freq,a.allele_B_freq,a.cases_maf,a.controls_maf,a.miss_dat_prop,a.cohort_hwe,a.controls_hwe, a.het_OR,a.het_OR_lower,a.het_OR_upper,a.hom_OR,a.hom_OR_lower,a.hom_OR_upper,all_OR_lower,a.all_OR_upper,a.se_assoc,"

        #seuil = -np.log10(threshold)
        #df_HTML = DF.to_html(index=False,escape=False)                #Dataframe to html
        graph,div1=manhattan(sorted_data, -np.log10(0.00000005))                                       # CREATE MANHATTAN

        #return render(request, 'index.html',{'phenotype':phenotype_selected,'div1':"<div><a href=\"/static/manhattan.png\">Download GRAPH</a></div>",'graph':"<div><img src=\"/static/"+m+"\" alt=\"graph\" type='image/png' style=\"background-size:100%;width:80%;height:400px;\"\></div>",'link':"<div><a href=\"/static/Manhattan.csv\">Download CSV File</a></div>",'pheno':phenotypes_found,'table':df_HTML,'user':request.user})


        #JSON
        #print(sorted_sorted_data)
        sorted_data_json = sorted_data.to_json(orient='records')
        return render(request,'modules/sGene.html',{'username':request.user,'phenotype':phenotype_selected,'graph': graph, 'div1':div1,'pheno': phenotypes_found, 'jsonData':sorted_data_json,'h_data':sorted_data,'gene':g})



#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                                                    #Logout
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

def logoutUser(request):
    logout(request)
    return redirect('/')







#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                                                                                    #END
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

