import sys, argparse
import logging
import os
import scanpy as sc
from markcorr import *
import pandas as pd

def main(args_input = sys.argv[1:]):
    anndata_input = ""
    output_path = "./"
    cell_type_number = 2
    USAGE='''
        Main script for mark cross correlation calculation.
        Usage: python main.py -i <AnndataInput> -o <OutputPath> -n <CellTypeNumber>
        required argument:
            -i | --AnndataInput : path to Anndata input file 
            -o | --OutputPath : path to store output result
        optional argument:
            -n | --CellTypeNumber : the number of cell types to calculate mark cross correlation
        '''
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i','--AnndataInput',
        help="path to Anndata input file"
    )
    parser.add_argument(
        '-o','--OutputPath',
        help="path to store output result"
    )
    parser.add_argument(
        '-n','--CellTypeNumber',
        help="the number of cell types to calculate mark cross correlation"
    )

    args = parser.parse_args(args_input)
    if not args.AnndataInput:
        logging.info(USAGE)
        logging.error("Input Anndata file is required")
        sys.exit(2)
    # if not args.OutputPath:
    #     logging.error("Output path is required")
    #     sys.exit(2)
    if args.CellTypeNumber:
        cell_type_number=args.CellTypeNumber
    anndata_input = args.AnndataInput
    output_path = args.OutputPath

    ext = os.path.splitext(anndata_input)[-1].lower()
    if ext == ".csv":
        adata = sc.read_csv(anndata_input)
        adata_df = adata.to_df()
    elif ext == ".h5ad":
        adata = sc.read_h5ad(anndata_input)
        adata_df = adata.to_df()
    else:
        logging.error(anndata_input, "is an unknown file format.")

    spatial = adata.obsm['spatial']
    size = adata.obsm['size']
    if size == None:
        logging.error("Size information should be stored in the obsm of Anndata file")
    adata_df['x spatial coordinate'] = pd.Series(spatial[:,0]).values
    adata_df['y spatial coordinate'] = pd.Series(spatial[:,1]).values
    adata_df['marks'] = adata.obs
    adata_df['diameter'] =  pd.Series(size[:,0]).values

    x = adata_df['x spatial coordinate']
    y = adata_df['y spatial coordinate']
    m = adata_df['marks']
    d = adata_df['diameter']
    xrange = [0,max(x)]
    yrange = [0,max(y)]
    W = window(xrange,yrange)
    p = pointPattern(adata_df['x spatial coordinate'],adata_df['y spatial coordinate'],d,W,m,cell_type_number)
    r, funs = markcorr(p)

    # plot the figures
    output_pdf = os.path.join(output_path,"markcrosscorr.pdf")
    fig=plt.figure()
    i = 1
    for f in funs:
        fig.add_subplot(cell_type_number,cell_type_number,i)
        plt.imshow(r,f)
        i+=1
    fig.savefig(output_pdf, format="pdf", bbox_inches="tight")

if __name__ == '__main__':
    main()