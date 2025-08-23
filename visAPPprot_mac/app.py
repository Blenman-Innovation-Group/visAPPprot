from rpy2 import robjects
import rpy2.robjects as robjects
import sys
from flask import Flask, redirect, url_for, render_template, request, flash, jsonify
import numpy as np
import pandas as pd
import pickle
import os
import json
import rpy2
from rpy2.rinterface import BoolSexpVector
from rpy2.robjects import pandas2ri
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
from io import StringIO
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import copy
import matplotlib.patches as patches
import requests
from xml.etree import ElementTree
from os import listdir, system
from os.path import isfile, join
import math
import re
import imagesize
from skimage.io import imread
from exemplar_process.method1 import get_coordinates
from exemplar_process.object_extract import get_objects
from exemplar_process.text_extract import get_text
import inflect
import torch
from googleapiclient.discovery import build
from rpy2.rinterface_lib.embedded import RRuntimeError



sys.path.append('exemplar_process')




from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import DataFrame, FloatVector, IntVector, StrVector, ListVector, Matrix
from rpy2.robjects import conversion, default_converter



app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = ((2 ** 10) ** 2) * 900

r = robjects.r
r.source('deseq_functions.R')
r.source('functions.R')
array_vst_nocluster = ''
array_vst = ''
cur_dataset = ''
lfcshrink_type = ''
pathways_data = 'Pathways.RData'
exp_mat = ''
levels = []
cluster_rows = None
cluster_cols = None
raw = False
proteins = None
cols = None
submat = None
submat_proteins = None
submat_patients = None
heatmap_fig_count = 0
wgcna_fig_count = 0
pathway_fig_count = -1
volcano_fig_count = -1
resShrink = None
rowcol = {}
times_def = {"volcano":65, "pathways": 100, "raw": 20, "diff_exp": 130, "wgcna_all": 720}
times_cur = {"volcano":65, "pathways": 100, "raw": 20, "diff_exp": 130, "wgcna_all": 720}
name_to_uniprot = {}
uniprot_to_name = {}
gpr_header = []
column_name = None
gpr_proteins = 0
gpr_patients = 0
data_dir = "./processed_data/"
gpr_dir = "./microarray_data/"
filename = "./signatures/GAL_Fig2e-5PBMCs_scRNAseq_matrix.csv"
filename_terms = "./uniprot_keywords/function_corrV9.csv"
df = pd.read_csv(filename, sep=',', lineterminator='\n')
df_terms = pd.read_csv(filename_terms, sep=',', lineterminator='\n')
gene_cells = {}
gene_processes = {}
entity_type = {}
states = {}
states_index = 0
download_img_prefix = ""
template_images = []
words_list = {}
current = 0
api_key = ""
cx = ""



def prep_uniprot():
    global name_to_uniprot, uniprot_to_name

    name_to_uniprot_csv = open("ensembl_uniprot_only.csv", 'r')
    name_to_uniprot_csv.readline()

    # column 1: gene name. column 5: uniprot entry ID
    for row in name_to_uniprot_csv:
        row = row.strip().split(",")
        name_to_uniprot[row[1]] = row[5]
        uniprot_to_name[row[5]] = row[1]


def get_gpr_size():
    global gpr_dir

    onlyfiles = [f for f in listdir(gpr_dir) if isfile(join(gpr_dir, f))]
    gprfile = None
    if ".gpr" in onlyfiles[0]:
        gprfile = open(gpr_dir + onlyfiles[0], "r", encoding='latin-1')
    else:
        gprfile = open(gpr_dir + onlyfiles[1], "r", encoding='latin-1')
    
    # count lines in file
    num_lines = sum(1 for _ in gprfile)

    return num_lines, len(onlyfiles)



def check_gpr_dir():
    global gpr_dir

    if os.path.isdir(gpr_dir):
        if len(os.listdir(gpr_dir)) > 0:
            return True

        return False

    return False


def get_image_search_keys():
    global api_key, cx

    keys_file = pd.read_csv("./image_search.csv")

    cx = keys_file['CX'].values[0]
    api_key = keys_file['KEY'].values[0]


@app.route("/", methods=['GET', 'POST'])
def index():
    global cur_dataset, cluster_rows, cluster_cols, raw, exp_mat, levels, lfcshrink_type, times_cur, rowcol, pathways_data, uniprot_to_name, gpr_header, column_name, data_dir, gpr_dir, gpr_proteins, gpr_patients, download_img_prefix, template_images, words_list, base
    
    try:
        with conversion.localconverter(default_converter):
            base = importr('base')

        cluster_rows = None
        cluster_cols = None
        raw = False
        
        if not uniprot_to_name:
            prep_uniprot()

            # get keys for google image search 
            get_image_search_keys()

            if check_gpr_dir():
                with conversion.localconverter(default_converter):
                    # first check if there are any files in microarray_data
                    microarray_dir = [f for f in os.listdir("./microarray_data") if not f.startswith('.')]

                    if len(microarray_dir) > 0:
                        get_gpr_header = robjects.globalenv['get_gpr_header']
                        gpr_header_return = get_gpr_header(gpr_dir)

                        # get rid of [1] from R 
                        gpr_header = []
                        for header in gpr_header_return:
                            gpr_header.append(header)

                        # don't use the snr values from the gpr file, compute them on the fly
                        if "SNR" in gpr_header[-4]:
                            gpr_header = gpr_header[:-4]
                            if any("f635 mean" in s.lower() for s in gpr_header) and any("b635 mean" in s.lower() for s in gpr_header):
                                gpr_header = gpr_header + ["SNR mean 635"]
                            if any("f635 median" in s.lower() for s in gpr_header) and any("b635 median" in s.lower() for s in gpr_header):
                                gpr_header = gpr_header + ["SNR median 635"]
                            if any("532 mean" in s.lower() for s in gpr_header):
                                gpr_header = gpr_header + ["SNR mean 532"]
                            if any("f532 median" in s.lower() for s in gpr_header) and any("b532 median" in s.lower() for s in gpr_header):
                                gpr_header = gpr_header + ["SNR median 532"]
                        else:
                            if any("f635 mean" in s.lower() for s in gpr_header) and any("b635 mean" in s.lower() for s in gpr_header):
                                gpr_header = gpr_header + ["SNR mean 635"]
                            if any("f635 median" in s.lower() for s in gpr_header) and any("b635 median" in s.lower() for s in gpr_header):
                                gpr_header = gpr_header + ["SNR median 635"]
                            if any("532 mean" in s.lower() for s in gpr_header):
                                gpr_header = gpr_header + ["SNR mean 532"]
                            if any("f532 median" in s.lower() for s in gpr_header) and any("b532 median" in s.lower() for s in gpr_header):
                                gpr_header = gpr_header + ["SNR median 532"]
                    
                        # get gpr size to compute progress bar estimate
                        gpr_proteins, gpr_patients = get_gpr_size()
                        

        if request.method == 'POST':

            # before we have expression matrix file 
            if os.listdir("./microarray_data/") and (request.form.get('option') == 'GPR'):
                column_name = request.form['gpr_header']

                return redirect(url_for('compute_expmat'))

            if (column_name == None):
                column_name = "SNR"

            # after we have expression matrix file we can do these computations
            cur_dataset = request.form.get('rdata')
            exp_mat = request.form.get('csv')
            levels = [request.form.get('level1'), request.form.get('level2')]
            lfcshrink_type = request.form.get('lfcshrink_type')
            # set download image directory 
            download_img_prefix = "./static/download_imgs_" + cur_dataset.split(".")[0] + "/"

            # set up template image and words_list
            if len(template_images) == 0:
                template_dir = './static/exemplar_references/' + cur_dataset.split('.')[0] + '/'

                pathway_dir = template_dir + 'pathway/'

                for filename in os.listdir(pathway_dir):
                    if (".csv" not in filename) and (".txt" not in filename) and ("DS_Store" not in filename): # skip over words_list file
                        img_path = os.path.join(pathway_dir, filename)
                        template_images.append(img_path)

                # get list of words from each template image
                words_list_f = open(pathway_dir + "words_list.txt", "r")
                words_list = {} # image_name : {word: (x, y), word: (x, y), ...}

                for line in words_list_f:
                    line = line.strip().split("||")
                    filename = line[0]
                    words = line[1:]
                    words_pos = {}
                    for i in range(0, len(words), 3):
                        word = words[i]
                        word_x = float(words[i+1])
                        word_y = float(words[i+2])
                        words_pos[word] = (word_x, word_y)

                    words_list[filename] = words_pos


                words_list_f.close()


            if request.form.get('option') == 'Volcano Plot':
                return redirect(url_for('volcano'))
            elif request.form.get('option') == 'Pathway Map':
                pathways_data = request.form.getlist('pathways_data')
                return redirect(url_for('pathways'))
            elif request.form.get('option') == 'Heatmap':
                raw_norm = request.form['raw_norm']

                if "raw" in raw_norm.lower():
                    raw = True

                clustering = request.form.get('clustering')
                cluster_rows = False
                cluster_cols = False
                if clustering == 'rows':
                    cluster_rows = True
                    cluster_cols = False
                elif clustering == 'columns':
                    cluster_rows = False
                    cluster_cols = True
                else:
                    cluster_rows = True
                    cluster_cols = True

                return redirect(url_for('heatmap', r_func=cur_dataset)) # do something else
            
        elif request.method == 'GET':

            # detect RData and csv files in local dir for user to select as input
            rdata = []
            pathways_list = []
            csv = []
            if len(exp_mat) > 0:
                csv = [exp_mat]

            onlyfiles = [f for f in listdir(data_dir) if isfile(join(data_dir, f))]

            for file in onlyfiles:
                if file.endswith('.txt'):
                    pathways_list.append(file)
                if file.endswith('.csv') and not ("expmat" in file.lower()):
                    rdata.append(file)
                if file.endswith('.csv') and ("expmat" in file.lower()) and (file not in csv):
                    csv.append(file)


            # prep all the proteins/patient lengths to use for progress bar compute in JS
            rowcol = {}

            for csv_file in csv:
                df = pd.read_csv(open(data_dir + csv_file, "r"))
                rowcol[csv_file] = df.shape    

            # fill in levels from previous
            level1 = ""
            level2 = ""
            if len(levels) > 0:
                level1 = levels[0]
                level2 = levels[1]


            # make download_img directories for each of the datasets they use in advance
            for rdataset in rdata:
                dataset_name = rdataset.split(".")[0]
                download_dir_name = "./static/download_imgs_" + dataset_name.split(".")[0] + "/"
                cmd = "mkdir " + download_dir_name
                system(cmd)

                # create directories for diff exp table files
                download_dir_name = "./differential_expression_tables/" + dataset_name.split(".")[0] + "/"
                cmd = "mkdir " + download_dir_name
                system(cmd)

                # also create directories for the exemplar images
                download_dir_name = "./static/exemplar_references/" + dataset_name.split(".")[0] + "/"
                cmd = "mkdir " + download_dir_name
                system(cmd)
                download_dir_name = "./static/exemplar_references/" + dataset_name.split(".")[0] + "/pathway/"
                cmd = "mkdir " + download_dir_name
                system(cmd)
                

                # if words_list not in pathway directory create it
                if os.path.exists(download_dir_name + "words_list.txt"):
                    continue
                else:
                    cmd = "touch " + download_dir_name + "words_list.txt"
                    system(cmd)

            lfcshrink_types = {"GLM": "Generalized Linear Model (GLM)", "normal": "GLM lfcShrink normal", "apeglm": "GLM lfcShrink apeglm", "ashr": "GLM lfcShrink ashr"}


            return render_template('index_horizontal.html', rdata=rdata, csv=csv, pathways_data=pathways_list, rowcol=rowcol, gpr_header=gpr_header, current_column=column_name, gpr_proteins=gpr_proteins, gpr_patients=gpr_patients, current_expmat=exp_mat, current_dataset=cur_dataset, level1=level1, level2=level2, lfcshrink_types=lfcshrink_types)
        
        return render_template("index_horizontal.html")

    except RRuntimeError as e:
        error_str = "R function error. " + str(e)
        return render_template('error.html', message=error_str)

    except Exception as e:
        error_str = "Error initialzing visAPPprot system. You likely have an issue with incorrect data formatting. " + str(e)
        return render_template('error.html', message=error_str)


def get_uniprot_info(gene_name):
    global name_to_uniprot

    try:
        # if "frag" within gene name this is not going to be within uniprot so get rid of "frag"
        gene_name_update = gene_name

        if "frag" in gene_name_update:
            ind = gene_name_update.find("frag")
            if (ind == len(gene_name_update) - 4):
                gene_name_update = gene_name_update[:ind]

        # only look for exact gene names within human proteins
        # https://bioinformatics.stackexchange.com/questions/22020/fetching-protein-sequences-through-uniprot-gives-connectionerror
        url = f'https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene_name_update}+AND+taxonomy_id:9606&format=json'

        response = requests.get(url)

        function_text = []

        if response.status_code == 200:
            response = response.json()
            results = response["results"]

            for entry in results:
                uniprot_entry_name = entry['primaryAccession']
                if "comments" in entry:
                    comments = entry["comments"]
                    for comment in comments:
                        if comment['commentType'] == "FUNCTION": 
                            texts = comment['texts']
                            for text in texts:
                                if text not in function_text:
                                    function_text.append(text['value'])

        if len(function_text) > 0:
            return "\n".join(function_text)

        # if no function text about this protein return
        return "No function information available"

    except Exception as e:
        error_str = "Error getting UniProt function information. " + str(e)
        return render_template('error.html', message=error_str)



@app.route('/compute_expmat/', methods=['GET', 'POST']) 
def compute_expmat():
    global column_name, exp_mat, data_dir, gpr_dir

    try:
        with conversion.localconverter(default_converter):
            if "SNR" in column_name and (("mean" in column_name) or ("median" in column_name)):
                print("SNR")
                compute_expmat_snr = robjects.globalenv['compute_expmat_snr']
                res = compute_expmat_snr(column_name, gpr_dir, data_dir)
                print(res)
                exp_mat = res[0]
            else:
                print("else")
                compute_expmat_column = robjects.globalenv['compute_expmat_column']
                res = compute_expmat_column(column_name, gpr_dir, data_dir)
                print(res)
                exp_mat = res[0]

            return render_template("complete_expmat.html")

    except RRuntimeError as e:
        error_str = "R function error for computing expression matrix. " + str(e)
        return render_template('error.html', message=error_str)

    except Exception as e:
        error_str = "Error computing expression matrix. " + str(e)
        return render_template('error.html', message=error_str)



def find_cell_corr(diffexp_genes):
    global df, gene_cells, df_terms, entity_type

    try:
        columns_all = df.columns
        columns = list(df.columns)[3:] # only cell names

        columns_all_terms = df_terms.columns
        columns_terms = list(df_terms.columns)[2:] # only cell names

        gene_cells = {} # gene : [cell]
        entity_type = {} # entity : colname


        for gene in diffexp_genes:
            gene_row_index = df.index[df['GAL'] == gene] # get index of row matching gene name
            gene_row = df.loc[gene_row_index, columns] # get row of values

            if len(gene_row) > 0:
                gene_row_corr = gene_row.gt(0) # make all values > 0 True, otherwise False
                corr_cells = gene_row.columns[gene_row_corr.all()].tolist() # column names where cell-gene correlation exists

                # get uniprot and original gene name
                uniprot = df.loc[gene_row_index, "uniprot"].values[0]
                orig_gene = df.loc[gene_row_index, "gene"].values[0]

                #ignore nans
                if not isinstance(uniprot, float):
                    
                    gene_cells[gene] = corr_cells

                    for e in corr_cells:
                        entity_type[e.lower()] = "Cell Subclass"


            # repeat with terms df
            gene_row_index = df_terms.index[df_terms['GAL'] == gene] # get index of row matching gene name
            gene_row = df_terms.loc[gene_row_index, columns_terms] # get row of values


            if len(gene_row) > 0:
                gene_row_corr = gene_row.dropna(axis=1,how='all')
                corr_cols = list(gene_row_corr.columns) # column names where cell-gene correlation exists

                # get uniprot and original gene name
                uniprot = df_terms.loc[gene_row_index, " Uniprot"].values[0]

                #ignore nans
                if not isinstance(uniprot, float):
                    
                    for col in corr_cols:
                        entities = df_terms.loc[gene_row_index, col].values[0] # string
                        entities = entities.split("*")
                        if gene not in gene_cells:
                            gene_cells[gene] = []
                            
                        gene_cells[gene] += entities

                        for e in entities:
                            entity_type[e.lower()] = col.strip()


    except Exception as e:
        error_str = "Error finding cell correlation. " + str(e)
        return render_template('error.html', message=error_str)



def find_cell_corr_condensed_terms(diffexp_genes):
    global df, gene_cells, df_terms, entity_type, gene_processes

    try:
        engine = inflect.engine()

        columns_all = df.columns
        columns = list(df.columns)[3:] # only cell names

        columns_all_terms = df_terms.columns
        columns_terms = list(df_terms.columns)[2:] # only cell names

        gene_cells = {} # gene : [cell]
        entity_type = {} # entity : colname
        gene_processes = {} # gene : [process]

        # read in condensed cellxgene terms file for entity type: entity match
        entity_type_dict = {}

        entity_type_file = open("./uniprot_keywords/cellxgene_terms_v14.csv", mode="r", encoding='utf-8-sig')

        entity_type_header = entity_type_file.readline() # get rid of header
        entity_type_header = entity_type_header.strip().split(",")

        for line in entity_type_file:
            line = line.strip().split(",")
            
            for i in range(len(entity_type_header)):
                typ = entity_type_header[i]
                e = line[i].lower().replace('-', ' ')
                
                if typ not in entity_type_dict:
                    entity_type_dict[typ] = []
                entity_type_dict[typ].append(e)


        for gene in diffexp_genes:
            gene_row_index = df.index[df['GAL'] == gene] # get index of row matching gene name
            gene_row = df.loc[gene_row_index, columns] # get row of values

            if len(gene_row) > 0:
                gene_row_corr = gene_row.gt(0) # make all values > 0 True, otherwise False
                corr_cells = gene_row.columns[gene_row_corr.all()].tolist() # column names where cell-gene correlation exists

                # get uniprot and original gene name
                uniprot = df.loc[gene_row_index, "uniprot"].values[0]
                orig_gene = df.loc[gene_row_index, "gene"].values[0]

                #ignore nans
                if not isinstance(uniprot, float):
                    gene_cells[gene] = corr_cells

                    for e in corr_cells:
                        entity_type[e.lower().replace('-', ' ')] = "Cell Type"


            # repeat with terms df
            gene_row_index = df_terms.index[df_terms['GAL'] == gene] # get index of row matching gene name
            gene_row = df_terms.loc[gene_row_index, "Annotation"] # get row of values


            if len(gene_row) > 0:
                gene_row_corr = gene_row.dropna(axis=0,how='all')

                # get uniprot and original gene name
                uniprot = df_terms.loc[gene_row_index, "Uniprot"].values[0]

                #ignore nans
                if not isinstance(uniprot, float):

                    entities = gene_row.values[0] # string
                    # ignore nans AKA no correlations found
                    if isinstance(entities, float):
                        continue
                        
                    entities_raw = entities.strip().lower().replace('-', ' ').split("*")
                    if gene not in gene_cells:
                        gene_cells[gene] = []

                    if gene not in gene_processes:
                        gene_processes[gene] = []
                    
                    # if entity is not part of a cellxgene terms don't include it in vis
                    entities = []
                    for e in entities_raw:    
                        e = e.strip() # get rid of extra spaces if they used to be within quotes

                        e_singular = e
                        e_plural = engine.plural(e)
                        if (engine.singular_noun(e)):
                            e_singular = engine.singular_noun(e)
                        e_saved = ""

                        exists = False
                        for typ in entity_type_dict:
                            if (e in entity_type_dict[typ]):
                                exists = True
                                e_saved = e
                            elif(e_singular in entity_type_dict[typ]):
                                exists = True
                                e_saved = e_singular
                            elif(e_plural in entity_type_dict[typ]):
                                exists = True
                                e_saved = e_plural
                        if exists:
                            entities.append(e_saved)
                       
                    # check whether entity belongs to entity or process
                    for e in entities:
                        if e in entity_type_dict["Processes"]:
                            gene_processes[gene].append(e)
                        else:
                            gene_cells[gene].append(e)

                    for e in entities:
                        # if it's a process it's not going to be part of the cellxgene terms
                        for typ in entity_type_dict:
                            if e in entity_type_dict[typ]:
                                entity_type[e] = typ

    except Exception as e:
        error_str = "Error processing cell correlation. " + str(e)
        return render_template('error.html', message=error_str)



def img_name_from_url(url):
    img_name = url.split("/")[-1]

    if (".png" in img_name):
        ind = img_name.find('.png')
        img_name = img_name[0:ind] + '.png'
    elif (".jpg" in img_name):
        ind = img_name.find('.jpg')
        img_name = img_name[0:ind] + '.jpg'
    elif (".jpeg" in img_name):
        ind = img_name.find('.jpeg')
        img_name = img_name[0:ind] + '.jpeg'
    elif (".gif" in img_name):
        ind = img_name.find('.gif')
        img_name = img_name[0:ind] + '.gif'
    elif (".svg" in img_name):
        ind = img_name.find('.svg')
        img_name = img_name[0:ind] + '.svg'
    else:
        img_name += ".png"


    return img_name


def is_image_accessible(link):

    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/111.0.0.0 Safari/537.36'}

    try:
        # https://stackoverflow.com/questions/43048099/image-download-mime-type-validation-python-requests
        rsp = requests.get(link)
        content_type_received = rsp.headers['Content-Type'] # mime type
        if content_type_received.startswith('image'): # image/jpeg, image/png etc
            return True
        return False

    except requests.RequestException:
        return False



@app.route('/google_image_search', methods=['GET', 'POST']) 
def google_image_search():
    global current, api_key, cx

    try:
        service = build("customsearch", "v1", developerKey=api_key)

        img_data = []

        query = ""

        if request.method == 'POST':
            query = request.form['query']
            current = int(request.form['current'])
            nex = int(request.form['next'])
            prev = int(request.form['prev'])
            start = current

            if not query:
                query = "None"

            if ((nex > -1)):
                start = current + nex
                current = start
                if (start < 0):
                    start = 0
                    current = 0
            if ((prev > -1)):
                start = current - prev
                current = start
                if (start < 0):
                    start = 0
                    current = 0
            

            result = service.cse().list(
                q=query,
                cx=cx,
                searchType='image',
                start=start,
                num=10
                ).execute()


            for item in result['items']:
                link = item['link']
                # skip .svg
                if (".svg" in link):
                    continue

                # skip over cloudflare protected links
                if is_image_accessible(link):
                    # https://stackoverflow.com/questions/48255936/get-call-for-getting-list-of-google-image-picture-urls#:~:text=q=tbn:ANd9GcR6DGPUW2U6pDAUmuVz_DR7cxU%2DLPUdMrOzbN2E0GHE8_FegyyHhvdOSX4V,.gstatic.com/images?
                    thumbnailLink = item['image']['thumbnailLink']
                    contextLink = item['image']['contextLink']
                    height = item['image']['height']
                    width = item['image']['width']

                    img_data.append({'link': link, 'thumbnailLink': thumbnailLink, 'contextLink': contextLink, 'height':height, 'width': width})

        
        return render_template("google_search.html", img_data=img_data, query_string=query, current=current)


    except Exception as e:
        error_str = "Error running Google image search. " + str(e)
        return render_template('error.html', message=error_str)


@app.route('/template_references', methods=['GET', 'POST']) 
def template_references():
    global cur_dataset, exp_mat, levels, lfcshrink_type, volcano_fig_count, rowcol, resShrink, data_dir, column_name

    try:
        query_string = cur_dataset.split('-')[0] + ' pathway'

        # sort through images to see which ones are viable templates
        download_dir = 'static/template_references/' + query_string.split(" ")[0] + "_" + query_string.split(" ")[1] + '/'
        cmd = 'mv ' + '\'static/template_references/' + query_string + '\' ' + download_dir
        os.system(cmd)

        pathway_dir = download_dir + 'pathway/'

        if (not os.path.exists(pathway_dir)):

            cmd = 'mkdir ' + download_dir + 'pathway/'
            os.system(cmd)

            count = 0
            num_objects_total = []
            num_words_total = []
            img_names = []

            for filename in os.listdir(download_dir):
                img_path = os.path.join(download_dir, filename)
                # checking if it is a file

                if os.path.isfile(img_path) and ("DS_Store" not in img_path):                
                    shadow_img_path = img_path
                    image = imread(img_path)
                    shadow_img = imread(shadow_img_path)
                    word_list = get_coordinates(img_path, "grayscale", 1, 1, 5)
                    svg_dir = './'
                    num_objects = len(get_objects(svg_dir, image, word_list))
                    num_words = get_text(word_list, shadow_img, img_path)
                    num_objects_total.append(num_objects)
                    num_words_total.append(num_words)
                    img_names.append("#" + filename.split(".")[0].split("_")[1])
                    count += 1
                    
                    if ((num_words == 3) or (num_words >= (0.5*num_objects))):
                        continue
                    else:
                        # probably more shapes than text
                        cmd = "mv " + img_path + ' ' + pathway_dir + filename
                        os.system(cmd)


        # collect names of template images
        imgs = []
        for filename in os.listdir(pathway_dir):
            img_path = os.path.join(pathway_dir, filename)
            imgs.append(img_path)

        # add image names to list and send to html
        return render_template("template_references.html", template_images=imgs)

    except Exception as e:
        error_str = "Error processing exemplar images. " + str(e)
        return render_template('error.html', message=error_str)



@app.route('/volcano', methods=['GET', 'POST']) 
def volcano():
    global cur_dataset, exp_mat, levels, lfcshrink_type, volcano_fig_count, rowcol, resShrink, data_dir, column_name

    try:
        if request.method == 'GET':
            with conversion.localconverter(default_converter):
                volcano_compute = robjects.globalenv['volcano']

            with conversion.localconverter(default_converter):
                res = volcano_compute(data_dir, cur_dataset, exp_mat, levels, lfcshrink_type, cur_dataset.split(".")[0], volcano_fig_count)


            with conversion.localconverter(default_converter):
                resShrink = res[0] # resShrink

                resShrink = pd.DataFrame(resShrink)
                resShrink = resShrink.T
                # clear nans from output
                # https://stackoverflow.com/questions/75223099/na-character-not-identidied-as-nan-after-importing-it-into-python-with-rpy2
                resShrink[0] = resShrink[0].apply(lambda val: np.nan if isinstance(
                    val, rpy2.rinterface_lib.sexp.NACharacterType) 
                    else val
                )
                resShrink[1] = resShrink[1].apply(lambda val: np.nan if isinstance(
                    val, rpy2.rinterface_lib.sexp.NACharacterType) 
                    else val
                )
                resShrink[2] = resShrink[2].apply(lambda val: np.nan if isinstance(
                    val, rpy2.rinterface_lib.sexp.NACharacterType) 
                    else val
                )
                resShrink[3] = resShrink[3].apply(lambda val: np.nan if isinstance(
                    val, rpy2.rinterface_lib.sexp.NACharacterType) 
                    else val
                )
                resShrink[4] = resShrink[4].apply(lambda val: np.nan if isinstance(
                    val, rpy2.rinterface_lib.sexp.NACharacterType) 
                    else val
                )
                resShrink[5] = resShrink[5].apply(lambda val: np.nan if isinstance(
                    val, rpy2.rinterface_lib.sexp.NACharacterType) 
                    else val
                )

                resShrink = resShrink.T.to_dict('list')

                volcano_fig_count +=1 

                diffexp_genes = []

                for key in resShrink:
                    gene = resShrink[key]
                    label_interest = gene[3]
                    gene_name = gene[5]
                    
                    if label_interest:
                        # ignore nans
                        if not isinstance(label_interest, float):
                            diffexp_genes.append(gene_name)


                find_cell_corr_condensed_terms(diffexp_genes)

                # add uniprot column to the res2 file

                dir_prefix = "./differential_expression_tables/" + cur_dataset.split(".")[0] + "/"
                res_filename = dir_prefix + "res2Shrink_" + cur_dataset.split(".")[0] + "_volcano_" + lfcshrink_type + "_" + str(volcano_fig_count) + ".csv"
                res_file = open(res_filename, "r")

                res_filename_out = dir_prefix + "res2Shrink_" + cur_dataset.split(".")[0] + "_volcano_" + str(volcano_fig_count) + "_temp.csv"
                res_file_out = open(res_filename_out, "w")

                # header
                header = res_file.readline()
                header = header.strip() + ",uniprot\n"
                res_file_out.write(header)

                for line in res_file:
                    line = line.strip().split(",")
                    # get label_interest
                    label_interest = line[4].strip("\"")
                    if len(label_interest) > 0:
                        # get uniprot information
                        function_text = get_uniprot_info(line[4]).strip()
                        line = ",".join(line) + "," + "\"" + function_text + "\"" + "\n"

                        res_file_out.write(line)
                    else:
                        line = ",".join(line) + "\n"
                        res_file_out.write(line)

                res_file.close()
                res_file_out.close()
                cmd = "mv " + res_filename_out + " " + res_filename
                os.system(cmd)


                return render_template("volcano.html", resShrink=resShrink, levels=levels, cur_dataset=cur_dataset, volcano_fig_count = volcano_fig_count, rowcol = rowcol, exp_mat=exp_mat, column_name=column_name)

        if request.method == 'POST':
            form_args = request.form.to_dict()

            function_text = get_uniprot_info(form_args['protein_name'])

            protein_info = {'function_text': function_text}
            return protein_info

    except RRuntimeError as e:
        error_str = "R function error for volcano plot. " + str(e)
        return render_template('error.html', message=error_str)


def get_words_from_exemplar(template_image_names):
    global cur_dataset, exp_mat, levels, volcano_fig_count, rowcol, resShrink, data_dir, column_name, words_list

    try:
        query_string = cur_dataset.split('-')[0] + ' pathway'

        # sort through images to see which ones are viable templates
        pathway_dir = './static/exemplar_references/' + cur_dataset.split('.')[0] + '/pathway/'

        # reopen outfile for appending this time
        words_list_f = open(pathway_dir + "words_list.txt", "a")

        for filename in template_image_names:

            # skip files that already exist 
            if filename in words_list:
                continue
            if ("DS_Store" in filename) or (".txt" in filename) or (".csv" in filename):
                continue

            img_path = os.path.join(pathway_dir, filename)
            print("starting to get words")
            print(img_path)

            # word detection
            shadow_img_path = img_path
            print("imread")
            image = imread(img_path)
            print("imread")
            shadow_img = imread(shadow_img_path)
            print("get_coordinates")
            word_list = get_coordinates(img_path, "grayscale", 1, 1, 5)
            svg_dir = './'
            print("get_objects")
            num_objects = len(get_objects(svg_dir, image, word_list))
            print("get_text")
            num_letters, words_pos = get_text(word_list, shadow_img, img_path)


            new_line = filename
            words_list[filename] = {}

            for i in range(0, len(words_pos), 3):
                word = words_pos[i]
                word_x = words_pos[i+1]
                word_y = words_pos[i+2]

                new_line = new_line + "||" + word + "||" + word_x + "||" + word_y
                words_list[filename][word] = (word_x, word_y)

            new_line += "\n"
            words_list_f.write(new_line)

        return

    except Exception as e:
        error_str = "Error extracting words from exemplar images. " + str(e)
        return render_template('error.html', message=error_str)


@app.route('/context_map', methods=['GET', 'POST']) 
def cell_correlation():
    global gene_cells, gene_processes, entity_type, resShrink, levels, states, states_index, cur_dataset, download_img_prefix, template_images, words_list, s3

    try:
        if request.method == 'POST':

            if request.is_json:
                # handle progress save
                form = request.get_json()

                if 'states' in form:
                    states = form["states"]
                    states_index = form["states_index"]
                    shapes_counter = form["shapes_counter"]
                    text_counter = form["text_counter"]


                    dataset_prefix = cur_dataset.split('.')[0]
                    with open('states' + dataset_prefix + '.pickle', 'wb') as f:
                        pickle.dump(states, f)
                    with open('states_index' + dataset_prefix + '.pickle', 'wb') as f:
                        pickle.dump(states_index, f)
                    with open('shapes_counter' + dataset_prefix + '.pickle', 'wb') as f:
                        pickle.dump(shapes_counter, f)
                    with open('text_counter' + dataset_prefix + '.pickle', 'wb') as f:
                        pickle.dump(text_counter, f)

            # download images from google, detect words, add to existing templates in front end
            else:
                if 'template_urls' in request.form:
                    template_urls = json.loads(request.form["template_urls"])

                    # get image names to extract words from
                    template_image_names = []

                    # download each image from the list of urls
                    for url in template_urls:
                        print(url)
                        img_name = img_name_from_url(url)
                        download_path = './static/exemplar_references/' + cur_dataset.split('.')[0] + '/pathway/' + img_name
                        
                        # skip images already downloaded
                        if download_path in template_images:
                            continue


                        template_image_names.append(img_name)

                        try:
                            f = open(download_path,'wb')

                            headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/111.0.0.0 Safari/537.36'}

                            r = requests.get(url, stream=True)
                            
                            # https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests
                            for chunk in r.iter_content(chunk_size=8192): 
                                f.write(chunk)

                            f.close()
                            template_images.append(download_path)

                            # update words_list to include new image 
                            get_words_from_exemplar(template_image_names)
                        except Exception as e:
                            print ("-- exception --")
                            print(e)


        # load state files if available
        states = {}
        states_index = 0
        shapes_counter = 0
        text_counter = 0
        # template_words_match = []

        dataset_prefix = cur_dataset.split('.')[0]

        if os.path.exists("./states" + dataset_prefix + ".pickle"):
            with open('states' + dataset_prefix + '.pickle', 'rb') as f:
                states = pickle.load(f)
            with open('states_index' + dataset_prefix + '.pickle', 'rb') as f:
                states_index = pickle.load(f)
            with open('shapes_counter' + dataset_prefix + '.pickle', 'rb') as f:
                shapes_counter = pickle.load(f)
            with open('text_counter' + dataset_prefix + '.pickle', 'rb') as f:
                text_counter = pickle.load(f)


        # get icon names from static folder
        onlyfiles = [f for f in listdir("./static/icons") if isfile(join("./static/icons", f))]
        icons = []

        for file in onlyfiles:
            if file.endswith('.png'):
                icons.append(file.split(".")[0])

        onlyfiles = [f for f in listdir("./static/cell_icons") if isfile(join("./static/cell_icons", f))]
        # sort names alphabetically
        onlyfiles.sort()

        rc_icons = ["Add entity"]
        rc_icons_size = {}

        for file in onlyfiles:
            if file.endswith('.png'):
                rc_icons.append(file.split(".")[0])
                width, height = imagesize.get("./static/cell_icons/" + file)
                rc_icons_size[file.split(".")[0]] = (width, height)

        # read in list of organ-gene correlations
        organ_gene_file = open("./static/tissue.tsv", "r")

        columns = organ_gene_file.readline().split("\t")

        organ_genes = {} # organ: [genes]
        gene_organs = {} # gene: [organs]

        for line in organ_gene_file:
            line = line.split("\t")
            organ = line[0]
            genes = line[1].split("*")
            organ_genes[organ] = genes
            for gene in genes:
                if gene not in gene_organs:
                    gene_organs[gene] = []
                gene_organs[gene].append(organ)

        # sort through images to see which ones are viable templates
        progress_thumbnail_dir = download_img_prefix

        progress_thumbnails = []
        for filename in os.listdir(progress_thumbnail_dir):
            if (".png" in filename and "progress" in filename):
                img_path = os.path.join(progress_thumbnail_dir, filename)
                progress_thumbnails.append(img_path)

        progress_thumbnails.sort(key=lambda x: os.path.getmtime(x))
        progress_thumbnails.reverse() # sort thumbnails based on most recent


        return render_template("layout.html", cur_dataset=cur_dataset, gene_cells=gene_cells, gene_processes=gene_processes, resShrink=resShrink, levels=levels, entity_type=entity_type, icons=icons, rc_icons=rc_icons, rc_icons_size=rc_icons_size, template_images=template_images, words_list=words_list, progress_thumbnails=progress_thumbnails, organ_genes=organ_genes, gene_organs=gene_organs, states=states, states_index=states_index, shapes_counter=shapes_counter, text_counter=text_counter)

    except Exception as e:
        error_str = "Error with context map functionality. " + str(e)
        return render_template('error.html', message=error_str)


@app.route('/pathways/')
def pathways():
    global cur_dataset, exp_mat, levels, lfcshrink_type, pathway_fig_count, rowcol, data_dir, column_name

    try:
        with conversion.localconverter(default_converter):
            pathways_compute = robjects.globalenv['pathways']

            res = pathways_compute(data_dir, cur_dataset, exp_mat, levels, pathways_data, lfcshrink_type, cur_dataset.split(".")[0], pathway_fig_count)

            fig = res[0] # vertical pathway heatmap
            index = res[1] # stars

            fig = pd.DataFrame(fig)
            fig = fig.to_dict('list')

            index = list(index)

            pathway_fig_count += 1

            # add uniprot column to the res2 file
            dir_prefix = "./differential_expression_tables/" + cur_dataset.split(".")[0] + "/"
            res_filename = dir_prefix + "res2Shrink_" + cur_dataset.split(".")[0] + "_pathway_" + lfcshrink_type + "_" + str(pathway_fig_count) + ".csv"
            res_file = open(res_filename, "r")

            res_filename_out = dir_prefix + "res2Shrink_" + cur_dataset.split(".")[0] + "_pathway_" + str(pathway_fig_count) + "_temp.csv"
            res_file_out = open(res_filename_out, "w")

            # header
            header = res_file.readline()
            header = header.strip() + ",uniprot\n"
            res_file_out.write(header)

            for line in res_file:
                line = line.strip().split(",")
                # get label_interest
                label_interest = line[4].strip("\"")
                if len(label_interest) > 0:
                    # get uniprot information
                    function_text = get_uniprot_info(line[4]).strip()
                    line = ",".join(line) + "," + "\"" + function_text + "\"" + "\n"

                    res_file_out.write(line)
                else:
                    line = ",".join(line) + "\n"
                    res_file_out.write(line)

            res_file.close()
            res_file_out.close()
            cmd = "mv " + res_filename_out + " " + res_filename
            os.system(cmd)

            
            return render_template("pathways.html", fig=fig, index=index, levels=levels, cur_dataset=cur_dataset, pathway_fig_count = pathway_fig_count, rowcol = rowcol, exp_mat=exp_mat, column_name=column_name)

    except RRuntimeError as e:
        error_str = "R function error for pathway map. " + str(e)
        return render_template('error.html', message=error_str)

    except Exception as e:
        error_str = "Error generating pathway map. " + str(e)
        return render_template('error.html', message=error_str)


@app.route('/render_heatmap_hm')
def render_heatmap_hm():
    global array_vst, raw, proteins, cols, cur_dataset, heatmap_fig_count, cluster_rows, cluster_cols, column_name, download_img_prefix

    try:
        # convert array_vst to np array, then re-orient array to work for matplotlib
        # nope need to keep same orientation to match pheatmap orientation
        np_expmat = copy.deepcopy(array_vst).T
        np_expmat = np_expmat.to_numpy()
        rows = copy.deepcopy(proteins)

        nodes = [0, 0.5, 1.0]
        colors = ["blue", "white", "red"]
        v_max = np.percentile(np_expmat, 99)
        v_min = 0

        if (not raw):
            v_max = np.max(np_expmat)
            one_cutoff_lower = 0.5 - 0.5/v_max
            one_cutoff_upper = 0.5 + 0.5/v_max

            nodes = [0, one_cutoff_lower, one_cutoff_upper, 1.0]
            colors = ["blue", "white", "white", "red"]
            v_min = -1*v_max

        
        cm = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
        cm.set_under("black")

        fig_w = min(round(len(cols)/3), 430)
        fig_h = min(round(len(rows)/5), 430)

        fig, ax = plt.subplots(figsize=(fig_w,fig_h))

        im = ax.imshow(np_expmat, cmap=cm, vmin=v_min, vmax=v_max, aspect='auto')

        start_y = -0.5

        min_lim, max_lim = ax.get_xlim()
        xmin = min(min_lim, max_lim)  
        offset = abs(max_lim - min_lim) * 0.4

        start_bar = xmin - offset

        patch_row = None
        patch_col = None

        if cluster_rows:
            patch_row = ax.add_patch(patches.Rectangle((start_bar,-0.5), 0.02, len(rows),color="green",
                                      clip_on=False))
        if cluster_cols: 
            patch_col = ax.add_patch(patches.Rectangle((-0.5,-3), len(cols), 0.02, color="green",
                                      clip_on=False))


        # https://stackoverflow.com/questions/61498858/using-setp-to-hide-axes-spines
        plt.setp(ax.spines.values(), visible=False) 

        # https://www.pythoncharts.com/matplotlib/customizing-grid-matplotlib/
        ax.set_xticks(np.arange(np_expmat.shape[1]+1)-.5, minor=True)
        ax.set_yticks(np.arange(np_expmat.shape[0]+1)-.5, minor=True)
        ax.grid(which="minor", color="black", linestyle='-', linewidth=0.5)
        ax.tick_params(which="minor", bottom=False, left=False)


        fs_cols = fig_w*0.8
        fs_rows = fig_h*0.02

        if fig_w == 430:
            fs_cols = 2
        if fig_h == 430:
            fs_rows = 2


        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(0, len(cols)))
        ax.set_xticklabels(cols, fontsize=fs_cols)
        ax.set_yticks(np.arange(0, len(rows)))
        ax.set_yticklabels(rows, fontsize=fs_rows)

        plt.setp(ax.get_xticklabels(), rotation=90, ha="center",
                 rotation_mode="default")

        # https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
        # https://matplotlib.org/stable/gallery/axes_grid1/demo_colorbar_with_axes_divider.html
        divider = make_axes_locatable(ax)

        cax = divider.append_axes("top", size="0.02%", pad=0.3)

        cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
        cbar.ax.tick_params(labelsize=5)


        chart_title = column_name
        if (not raw):
            chart_title = 'VST ' + chart_title
        cbar.ax.set_title(chart_title)

        column_name_file = re.sub(r'[\\/*?:"<>|]',"_",column_name)

        fig_title = download_img_prefix + "heatmap_" + cur_dataset.split(".")[0] + "_unnormalized_"  + str(heatmap_fig_count) + ".svg"

        if (not raw):
            fig_title = download_img_prefix + "heatmap_" + cur_dataset.split(".")[0] + "_normalized_" + str(heatmap_fig_count) + ".svg"


        bbox_artist = []
        if patch_row:
            bbox_artist.append(patch_row)
        if patch_col:
            bbox_artist.append(patch_col)

        fig.subplots_adjust(left=0.5)  
        plt.savefig(fig_title, dpi=150, bbox_extra_artists=bbox_artist)
        plt.close()

        heatmap_fig_count += 1

        return redirect(url_for('heatmap', r_func=cur_dataset))


    except Exception as e:
        error_str = "Error saving heatmap. " + str(e)
        return render_template('error.html', message=error_str)



@app.route('/render_heatmap_wgcna', methods=["POST", "GET"])
def render_heatmap_wgcna():
    global raw, cur_dataset, wgcna_fig_count, column_name

    try:
        data = request.get_json()

        submat_hm = np.array(data.get("submat"))
        submat_proteins_hm = data.get("submat_proteins")
        submat_patients_hm = data.get("submat_patients")
        cluster_all = data.get("cluster_all")
        mods = data.get("mods")
        colorscheme = data.get("colorscheme")

        # convert array_vst to np array, then re-orient array to work for matplotlib
        np_expmat = copy.deepcopy(submat_hm)
        np_expmat = np.flip(np_expmat, axis=0)
        rows = copy.deepcopy(submat_proteins_hm)
        rows.reverse()
        cols = copy.deepcopy(submat_patients_hm)

        nodes = [0, 0.5, 1.0]
        colors = ["blue", "white", "red"]
        v_max = np.percentile(np_expmat, 99)
        v_min = 0

        if (not raw):
            v_max = np.max(np_expmat)

            one_cutoff_lower = 0.5 - 0.5/v_max
            one_cutoff_upper = 0.5 + 0.5/v_max

            nodes = [0, one_cutoff_lower, one_cutoff_upper, 1.0]
            colors = ["blue", "white", "white", "red"]
            v_min = -1*v_max

        
        cm = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
        cm.set_under("black")


        fig_w = min(round(len(cols)/3), 430)
        fig_h = min(round(len(rows)/5), 430)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        im = ax.imshow(np_expmat, cmap=cm, vmin=v_min, vmax=v_max, aspect='auto')

        start_y = -0.5

        min_lim, max_lim = ax.get_xlim()
        xmin = min(min_lim, max_lim)  
        offset = abs(max_lim - min_lim) * 0.2

        start_bar = xmin - offset

        for mod in mods:
            mod_ind = mod[0] 
            num_proteins = mod[1] 

            # sometimes R returns 1 index, sometimes 0
            if mods[len(mods)-1][0] == 1:
                mod_ind -=1

            rect = patches.Rectangle((start_bar, start_y),offset, num_proteins, color=colorscheme[mod_ind], clip_on=False) 
            ax.add_patch(rect)
            start_y += num_proteins


        plt.setp(ax.spines.values(), visible=False) 

        linew = 3
        fs_cols = fig_w*0.8
        fs_rows = fig_h*0.02

        # if max limit of height or width reached then the cells are really compressed so make font size small
        if fig_w == 430:
            fs_cols = 2
        if fig_h == 430:
            fs_rows = 2
            linew = 0.5


        ax.set_xticks(np.arange(np_expmat.shape[1]+1)-.5, minor=True)
        ax.set_yticks(np.arange(np_expmat.shape[0]+1)-.5, minor=True)
        ax.grid(which="minor", color="black", linestyle='-', linewidth=linew)
        ax.tick_params(which="minor", bottom=False, left=False)

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(0, len(cols)))
        
        ax.set_xticklabels(cols, fontsize=fs_cols)
        ax.set_yticks(np.arange(0, len(rows)))
        
        ax.set_yticklabels(rows, fontsize=fs_rows)

        plt.setp(ax.get_xticklabels(), rotation=90, ha="center",
                 rotation_mode="default")


        divider = make_axes_locatable(ax)

        cax = divider.append_axes("top", size="0.03%", pad=0.2)

        cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
        cbar.ax.tick_params(labelsize=7)


        chart_title = column_name
        if (not raw):
            chart_title = 'VST ' + column_name
        cbar.ax.set_title(chart_title)

        column_name_file = re.sub(r'[\\/*?:"<>|]',"_",column_name)

        fig_title = download_img_prefix + "heatmap_" + cur_dataset.split(".")[0] + "_wgcna_" + str(wgcna_fig_count) + ".svg"

        if (not raw):
            fig_title = download_img_prefix + "heatmap_" + cur_dataset.split(".")[0] + "_wgcna_" + str(wgcna_fig_count) + ".svg"

        fig.subplots_adjust(left=0.2)  
        plt.savefig(fig_title, dpi=150)
        plt.close()

        wgcna_fig_count += 1

        return jsonify({"redirect": url_for("heatmap", r_func=cur_dataset)})


    except Exception as e:
        error_str = "Error saving WGCNA heatmap. " + str(e)
        return render_template('error.html', message=error_str)


@app.route('/heatmap') 
def heatmap():
    global array_vst_nocluster, array_vst, cluster_rows, cluster_cols, raw, levels, lfcshrink_type, cur_dataset, exp_mat, proteins, cols, data_dir, column_name, heatmap_fig_count

    try:
        array_vst_js = {}
        proteins = []
            
        # run deseq for heatmap
        with conversion.localconverter(default_converter):
            deseq = robjects.globalenv["deseq"]

            if raw:
                deseq = robjects.globalenv["rawdata"]

            res = None
           
            # NOTE rpy2 will convert all boolean values to TRUE unless explcitly stated this way
            if (cluster_rows and cluster_cols):
                res = deseq(data_dir, cur_dataset, exp_mat, levels, lfcshrink_type, BoolSexpVector([True]), BoolSexpVector([True]), cur_dataset.split(".")[0], heatmap_fig_count)
            elif (cluster_rows):
                res = deseq(data_dir, cur_dataset, exp_mat, levels, lfcshrink_type, BoolSexpVector([True]), BoolSexpVector([False]), cur_dataset.split(".")[0], heatmap_fig_count)
            elif (cluster_cols):
                res = deseq(data_dir, cur_dataset, exp_mat, levels, lfcshrink_type, BoolSexpVector([False]), BoolSexpVector([True]), cur_dataset.split(".")[0], heatmap_fig_count)


            array_vst = robjects.conversion.rpy2py(res[0])
            proteins = list(res[1])
            array_vst_nocluster = robjects.conversion.rpy2py(res[2])

            # we will lose column and row names when converting to dictionary later so save these and put it back later
            cols = list(array_vst.names)
            array_vst_js = pd.DataFrame(array_vst)

            array_vst_js.columns = proteins
            array_vst_js.index = cols

            array_vst = array_vst_js

            array_vst_js = array_vst_js.T.to_dict('list')


            # adjust height and width of heatmap based on total number of proteins & patients included in heatmap
            cols = list(array_vst.index)
            h = len(proteins) * 12 # 25
            w = len(cols) * 25 # 50

           
            return render_template("heatmap.html", array_vst = array_vst_js, proteins = proteins, patients = cols, w=w, h=h, raw=raw, cluster_rows=cluster_rows, cluster_cols=cluster_cols, column_name=column_name)

    except RRuntimeError as e:
        error_str = "R function error for heatmap. " + str(e)
        return render_template('error.html', message=error_str)

    except Exception as e:
        error_str = "Error generating overview heatmap. " + str(e)
        return render_template('error.html', message=error_str)


@app.route('/heatmap_return')
def heatmap_return():
    global cur_dataset
    return redirect(url_for('heatmap', r_func=cur_dataset)) # do something else

@app.route("/wgcna",methods=["POST"]) 
def wgcna():
    global array_vst_nocluster, array_vst, raw, levels, lfcshrink_type, cur_dataset, exp_mat, submat, submat_proteins, submat_patients, wgcna_fig_count, data_dir, column_name

    try:
        proteins = request.form["wgcna_proteins"]
        patients = request.form["wgcna_patients"]
        y_vals = request.form["wgcna_y"]
        proteins = proteins.split(",")
        patients = patients.split(",")
        y_vals = y_vals.split(",")
        y_vals = [int(i) for i in y_vals]
       
        array_vst = array_vst.T

        # again get column and row names before we lose them in conversion to dictionary
        submat = array_vst.loc[proteins, patients]
        # submat = array_vst
        submat_patients = list(submat.columns)
        submat_proteins = list(submat.index)
        submat_js = submat.to_dict('list')

        # total height of heatmap, based on total number of selected proteins
        heatmap_h = len(submat_proteins) * 25

        with conversion.localconverter(default_converter):
            wgcna_compute = robjects.globalenv['wgcna_compute']

            res = wgcna_compute(proteins, patients, y_vals, array_vst_nocluster, levels, lfcshrink_type, cur_dataset, exp_mat, data_dir)

            diss1_clustered = res[0] # don't use
            gene_names = res[1] # don't use
            dynamicMods = list(res[5]) # colors of modules
            dynamicMods_dd = res[6] # colors of modules
            diss_l = res[7] # list of modules
            genenames_l = res[8] # list of top 30 proteins per module

            cols = list(diss1_clustered.names)
            diss1_clustered = pd.DataFrame(diss1_clustered)

            diss1_clustered.columns = gene_names
            diss1_clustered.index = cols

            diss1_clustered = diss1_clustered.T.to_dict('list')

            gene_names = list(gene_names)
            dynamicMods_dd = list(dynamicMods_dd)


            diss_l = list(diss_l)
            genenames_l = list(genenames_l)

            for i in range(len(diss_l)):
                diss_l[i] = robjects.conversion.rpy2py(pd.DataFrame(diss_l[i]))
                diss_l[i] = pd.DataFrame(diss_l[i])
                genenames_l[i] = list(genenames_l[i])
                diss_l[i].columns = genenames_l[i]
                diss_l[i].index = genenames_l[i]
                diss_l[i] = diss_l[i].T.to_dict('list')

            # order, merge, height hclust attributes to draw dendrogram
            order = list(res[2])
            merge = np.array(res[3]).tolist()
            height = list(res[4])


            return render_template("wgcna_sle.html", diss1_clustered = [], gene_names = [], dynamicMods_dd = dynamicMods_dd, dynamicMods = dynamicMods, diss_l = diss_l, genenames_l = genenames_l, order = order, merge = merge, height = height, submat = submat_js, patients = submat_patients, proteins = submat_proteins, heatmap_h = heatmap_h, raw=raw, total_proteins_n=len(array_vst.index), wgcna_fig_count=wgcna_fig_count, cluster_all=False, cur_dataset=cur_dataset, column_name=column_name)


    except RRuntimeError as e:
        error_str = "R function error in WGCNA clustering. " + str(e)
        return render_template('error.html', message=error_str)

    except Exception as e:
        error_str = "Error generating WGCNA figures. " + str(e)
        return render_template('error.html', message=error_str)




@app.route("/wgcna_all",methods=["POST"]) 
def wgcna_all():
    global array_vst_nocluster, array_vst, raw, levels, lfcshrink_type, cur_dataset, exp_mat, wgcna_fig_count, data_dir, column_name

    try:
        array_vst = array_vst.T

        # again get column and row names before we lose them in conversion to dictionary
        # submat = array_vst.loc[proteins]
        submat = array_vst
        submat_patients = list(submat.columns)
        submat_proteins = list(submat.index)
        submat_js = submat.to_dict('list')

        # total height of heatmap, based on total number of selected proteins
        heatmap_h = len(submat_proteins) * 25

        with conversion.localconverter(default_converter):
            wgcna_compute = robjects.globalenv['wgcna_compute_all']

            res = wgcna_compute(submat_proteins, submat_patients, array_vst_nocluster, levels, lfcshrink_type, cur_dataset, exp_mat, data_dir)

            diss1_clustered = res[0] # don't use
            gene_names = res[1] # don't use
            dynamicMods = list(res[5]) # colors of modules
            dynamicMods_dd = res[6] # colors of modules
            diss_l = res[7] # list of modules
            genenames_l = res[8] # list of top 30 proteins per module

            cols = list(diss1_clustered.names)
            diss1_clustered = pd.DataFrame(diss1_clustered)

            diss1_clustered.columns = gene_names
            diss1_clustered.index = cols

            diss1_clustered = diss1_clustered.T.to_dict('list')

            gene_names = list(gene_names)
            dynamicMods_dd = list(dynamicMods_dd)

            diss_l = list(diss_l)
            genenames_l = list(genenames_l)

            for i in range(len(diss_l)):
                diss_l[i] = robjects.conversion.rpy2py(pd.DataFrame(diss_l[i]))
                diss_l[i] = pd.DataFrame(diss_l[i])
                genenames_l[i] = list(genenames_l[i])
                diss_l[i].columns = genenames_l[i]
                diss_l[i].index = genenames_l[i]
                diss_l[i] = diss_l[i].T.to_dict('list')

            # order, merge, height hclust attributes to draw dendrogram
            order = list(res[2])
            merge = np.array(res[3]).tolist()
            height = list(res[4])

            return render_template("wgcna_sle.html", diss1_clustered = [], gene_names = [], dynamicMods_dd = dynamicMods_dd, dynamicMods = dynamicMods, diss_l = diss_l, genenames_l = genenames_l, order = order, merge = merge, height = height, submat = submat_js, patients = submat_patients, proteins = submat_proteins, heatmap_h = heatmap_h, raw=raw, total_proteins_n=len(array_vst.index), wgcna_fig_count=wgcna_fig_count, cluster_all=True, cur_dataset=cur_dataset, column_name=column_name)

    except RRuntimeError as e:
        error_str = "R function error in WGCNA clustering. " + str(e)
        return render_template('error.html', message=error_str)

    except Exception as e:
        error_str = "Error generating WGCNA figures. " + str(e)
        return render_template('error.html', message=error_str)


if __name__ == '__main__':
    app.run(host="localhost", port=8888, debug=True)

