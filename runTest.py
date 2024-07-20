import markcorr
import os
import shutil
import pandas as pd
import numpy as np
import altair as alt
from window import *
from pointPattern import *
import math
import pandas as pd
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from breakpts import *
from sm_density import *
from closepairs import *
from window import *
from unnormdensity import *
from multiprocessing import Pool
import random
import betacells

def same_cell_type():
    random.seed(12)
    xrange = [0,330]
    yrange = [0,330]
    W = window(xrange,yrange)
    x = []
    y = []
    for i in range(0,500):
        n = random.randint(80, 250)
        x.append(n)
        n = random.randint(0, 100)
        y.append(n)
    for i in range(0,300):
        n = random.randint(0, 120)
        x.append(n)
        n = random.randint(200, 300)
        y.append(n)
    mark1 = []
    for i in range(0,len(x)):
        n = random.randint(1,10)
        mark1.append(n)
    mark2 = []
    d = [3]*len(x)
    adipocyte_x = [35,98,30,150,200,230,280,80,50,200,300,130]
    adipocyte_y = [27,60,100,150,300,180,50,190,260,80,280,230]
    rm_idx = []
    for i in range(len(x)):
        for j in range(len(adipocyte_x)):
            d2_0 = (x[i]-adipocyte_x[j])**2 + (y[i]-adipocyte_y[j])**2
            if d2_0 < 961:
                rm_idx.append(i)
    x = np.delete(x,rm_idx).tolist()
    y = np.delete(y,rm_idx).tolist()
    mark1 = np.delete(mark1,rm_idx).tolist()
    d = np.delete(d,rm_idx).tolist()
    mark2.extend(random.sample(range(1, 10), 3))
    p_adipo = pointPattern(adipocyte_x,adipocyte_y,[60]*len(adipocyte_x),W,mark2)
    n=len(x)
    fig, ax = plt.subplots()
    plt.xticks(range(0, 360, 30))   
    plt.yticks(range(0, 360, 30))
    plt.axis('square')
    for i in range(n):
        x_p = x[i]
        y_p = y[i]
        d_p = d[i]
        c=plt.Circle((x_p, y_p), d_p/2)
        ax.add_artist(c)
    for i in range(len(adipocyte_x)):
        x_p = adipocyte_x[i]
        y_p = adipocyte_y[i]
        c=plt.Circle((x_p, y_p), 30)
        ax.add_artist(c)
    plt.title("Spatial Plot of Simulated Data - Same Cell Types", fontsize=15)
    plt.xlabel("X-corrdinate", fontsize=12)
    plt.ylabel("Y-corrdinate", fontsize=12)
    plt.savefig("./testResult/sameType.png")
    plt.show()
    mark = pd.DataFrame(np.array([mark1]).transpose(),columns=["small cell"])
    d = None
    p = pointPattern(x,y,d,W,mark)
    r, funs = markcorr.markcorr(p,pp=None, saveImage=True, savefolder="./testResult/sameType/", remove_zeros=False)
    
    iso, trans = {}, {}
    for i in funs:
        trans[i] = (funs[i][0])
        iso[i] = (funs[i][1])
    iso = pd.DataFrame(iso)
    trans = pd.DataFrame(trans)

    iso.to_csv("./testResult/sameType/iso.csv")
    trans.to_csv("./testResult/sameType/trans.csv")


def diff_cell_type():
    random.seed(12)
    xrange = [0,330]
    yrange = [0,330]
    W = window(xrange,yrange)
    x = []
    y = []
    for i in range(0,500):
        n = random.randint(80, 250)
        x.append(n)
        n = random.randint(0, 100)
        y.append(n)
    for i in range(0,300):
        n = random.randint(0, 120)
        x.append(n)
        n = random.randint(200, 300)
        y.append(n)
    mark1 = []
    mark2 = []
    for i in range(500):
        n = random.randint(1,10)
        mark1.append(n)
    for i in range(300):
        n = random.randint(1,10)
        mark2.append(n)
    mark3 = []
    d = [3]*len(x)
    adipocyte_x = [35,98,30,150,200,230,280,80,50,200,300,130]
    adipocyte_y = [27,60,100,150,300,180,50,190,260,80,280,230]
    rm_idx = []
    for i in range(len(x)):
        for j in range(len(adipocyte_x)):
            d2_0 = (x[i]-adipocyte_x[j])**2 + (y[i]-adipocyte_y[j])**2
            if d2_0 < 961:
                rm_idx.append(i)
    x = np.delete(x,rm_idx).tolist()
    y = np.delete(y,rm_idx).tolist()
    mark1 = np.delete(mark1,np.where(np.array(rm_idx) <= 500)).tolist()
    d = np.delete(d,rm_idx).tolist()
    mark2 = np.delete(mark2, np.where(np.array(rm_idx) > 500)).tolist()

    mark3.extend(random.sample(range(1, 10), 3))
    p_adipo = pointPattern(adipocyte_x,adipocyte_y,[60]*len(adipocyte_x),W,mark3)
    n=len(x)
    fig, ax = plt.subplots()
    plt.xticks(range(0, 360, 30))   
    plt.yticks(range(0, 360, 30))
    plt.axis('square')
    for i in range(n):
        x_p = x[i]
        y_p = y[i]
        d_p = d[i]
        c=plt.Circle((x_p, y_p), d_p/2)
        ax.add_artist(c)
    for i in range(len(adipocyte_x)):
        x_p = adipocyte_x[i]
        y_p = adipocyte_y[i]
        c=plt.Circle((x_p, y_p), 30)
        ax.add_artist(c)
    plt.title("Spatial Plot of Simulated Data - Different Cell Types", fontsize=15)
    plt.xlabel("X-corrdinate", fontsize=12)
    plt.ylabel("Y-corrdinate", fontsize=12)
    plt.savefig("./testResult/diffTypes.png")
    plt.show()
    mark_1 = mark1+[0]*len(mark2)
    mark_2 = [0]*len(mark1)+mark2
   
    mark = pd.DataFrame()
    mark["smallcell1"] = mark_1
    mark["smallcell2"] = mark_2
    d = None
    p = pointPattern(x,y,d,W,mark)
    r, funs = markcorr.markcorr(p,pp=None, saveImage=True, savefolder="./testResult/diffTypes/", remove_zeros= False)

    iso, trans = {}, {}
    for i in funs:
        trans[i] = (funs[i][0])
        iso[i] = (funs[i][1])
    iso = pd.DataFrame(iso)
    trans = pd.DataFrame(trans)

    iso.to_csv("./testResult/diffTypes/iso.csv")
    trans.to_csv("./testResult/diffTypes/trans.csv")


def general_case():
    random.seed(12)
    xrange = [0,330]
    yrange = [0,330]
    W = window(xrange,yrange)
    x = []
    y = []
    for i in range(0,2000):
        n = random.randint(0, 330)
        x.append(n)
        n = random.randint(0, 330)
        y.append(n)
    # for i in range(0,500):
    #     n = random.randint(80, 250)
    #     x.append(n)
    #     n = random.randint(0, 100)
    #     y.append(n)
    mark1 = []
    for i in range(0,len(x)):
        n = random.randint(1,10)
        mark1.append(n)
    mark2 = []
    d = [3]*len(x)
    adipocyte_x = [35,98,30,150,200,230,280,80,50,200,300,130]
    adipocyte_y = [27,60,100,150,300,180,50,190,260,80,280,230]
    rm_idx = []
    for i in range(len(x)):
        for j in range(len(adipocyte_x)):
            d2_0 = (x[i]-adipocyte_x[j])**2 + (y[i]-adipocyte_y[j])**2
            if d2_0 < 961:
                rm_idx.append(i)
    x = np.delete(x,rm_idx).tolist()
    y = np.delete(y,rm_idx).tolist()
    mark1 = np.delete(mark1,rm_idx).tolist()
    d = np.delete(d,rm_idx).tolist()

    mark2.extend(random.sample(range(1, 10), 3))
    # d.extend([60,60,60])
    p_adipo = pointPattern(adipocyte_x,adipocyte_y,[60]*len(adipocyte_x),W,mark2)
    n=len(x)
    fig, ax = plt.subplots()
    plt.xticks(range(0, 360, 30))   
    plt.yticks(range(0, 360, 30))
    plt.axis('square')
    for i in range(n):
        x_p = x[i]
        y_p = y[i]
        d_p = d[i]
        c=plt.Circle((x_p, y_p), d_p/2)
        ax.add_artist(c)
    for i in range(len(adipocyte_x)):
        x_p = adipocyte_x[i]
        y_p = adipocyte_y[i]
        c=plt.Circle((x_p, y_p), 30)
        ax.add_artist(c)
    plt.title("Spatial Plot of Simulated Data - General Case", fontsize=15)
    plt.xlabel("X-corrdinate", fontsize=12)
    plt.ylabel("Y-corrdinate", fontsize=12)
    plt.savefig("./testResult/general.png")
    plt.show()
    mark = pd.DataFrame(np.array([mark1]).transpose(),columns=["small cell"])
    d=None
    p = pointPattern(x,y,d,W,mark)
    r, funs = markcorr.markcorr(p,pp=None, saveImage=True, savefolder="./testResult/general/", remove_zeros=False)

    iso, trans = {}, {}
    for i in funs:
        trans[i] = (funs[i][0])
        iso[i] = (funs[i][1])
    iso = pd.DataFrame(iso)
    trans = pd.DataFrame(trans)

    iso.to_csv("./testResult/general/iso.csv")
    trans.to_csv("./testResult/general/trans.csv")

def run_betacells():
    xrange = [28.08, 778.08]
    yrange = [16.2, 1007.02]
    W = window(xrange,yrange) 
    x = betacells.x
    y = betacells.y
    d = [1] * len(x)
    mark = {}
    mark["type"] = betacells.type
    m = pd.DataFrame(mark)
    m["type"] = m["type"].astype("category")
    p = pointPattern(x,y,d,W,m)
    r, funs = markcorr.markcorr(p, savefolder="./testResult/betacells/", remove_zeros=False, saveImage=False)
    
    iso, trans = {}, {}
    for i in funs:
        trans[i] = (funs[i][0])
        iso[i] = (funs[i][1])
    iso = pd.DataFrame(iso)
    trans = pd.DataFrame(trans)

    iso.to_csv("./testResult/betacells/iso.csv")
    trans.to_csv("./testResult/betacells/trans.csv")


def run_testcases(i):
    if i == 0:
        general_case()
    elif i == 1: 
        same_cell_type()
    elif i == 2:
        diff_cell_type()

def runAML(i):
    aml = pd.read_csv("../data/output/AML2.csv")
    imageSize_x = 1056
    imageSize_y = 642
    W = window((0, imageSize_x), (0, imageSize_y))

    AML_image = aml.loc[aml["ImageNumber"] == i + 1]
    x = AML_image["x"].tolist()
    y = AML_image["y"].tolist()
    mark = AML_image.drop(["x", "y", "Area", "ImageNumber", "Unnamed: 0", "CellType"], axis = 1)
    pp = AML_image.loc[AML_image["CellType"] ==  "Adipocytes"]
    if len(pp) == 0:
        pp = None
    else:
        pp_x = pp["x"].tolist()
        pp_y = pp["y"].tolist()
        pp_d = (np.sqrt((pp["Area"] / 3.14)) * 2).tolist()
        pp = pointPattern(pp_x, pp_y, pp_d, W)
    d = (np.sqrt(AML_image["Area"] / 3.14) * 2).tolist()

    # d = None
    # pp = None

    points = pointPattern(x, y, d, W, mark)

    folder_name = "./result/AML/image_%d/"%i
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    else:
        shutil.rmtree(folder_name)
        os.mkdir(folder_name)
  
    r, funs = markcorr.markcorr(points, savefolder = folder_name, pp= pp, remove_zeros=False)
    iso, trans = {}, {}
    for i in funs:
        trans[i] = (funs[i][0])
        iso[i] = (funs[i][1])
    iso = pd.DataFrame(iso)
    trans = pd.DataFrame(trans)

    iso.to_csv(folder_name + "iso.csv")
    trans.to_csv(folder_name + "trans.csv")

def runNormal (i):
    normal = pd.read_csv("../data/output/Normal2.csv")
    imageSize_x = 1056
    imageSize_y = 642
    W = window((0, imageSize_x), (0, imageSize_y))

    Normal_image = normal.loc[normal["ImageNumber"] == i + 1]
    x = Normal_image["x"].tolist()
    y = Normal_image["y"].tolist()
    mark = Normal_image.drop(["x", "y", "Area", "ImageNumber", "Unnamed: 0", "CellType"], axis = 1)
    pp = Normal_image.loc[Normal_image["CellType"] ==  "Adipocytes"]
    if len(pp) == 0:
        pp = None
    else:
        pp_x = pp["x"].tolist()
        pp_y = pp["y"].tolist()
        pp_d = (np.sqrt((pp["Area"] / 3.14)) * 2).tolist()
        pp = pointPattern(pp_x, pp_y, pp_d, W)
    d = (np.sqrt(Normal_image["Area"] / 3.14) * 2).tolist()

    # pp = None
    # d = None

    points = pointPattern(x, y, d, W, mark)
    
    folder_name = "./result/Normal/image_%d/"%i
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    else:
        shutil.rmtree(folder_name)
        os.mkdir(folder_name)
    
    
    r, funs = markcorr.markcorr(points, savefolder = folder_name, pp = pp, remove_zeros=False, saveImage=False)
    iso, trans = {}, {}
    for i in funs:
        trans[i] = (funs[i][0])
        iso[i] = (funs[i][1])
    iso = pd.DataFrame(iso)
    trans = pd.DataFrame(trans)

    iso.to_csv(folder_name + "iso.csv")
    trans.to_csv(folder_name + "trans.csv")

if __name__ == '__main__':

    # a = [0, 1, 5, 7, 13, 18, 19, 22, 23, 24, 25, 35]
    # with Pool(15) as p:
    #     p.map(runNormal, range(15))
    # same_cell_type()
    # diff_cell_type()
    # general_case()
    run_betacells()
