# (C) Robin Fallegger
import scanpy as sc
import pandas as pd
import muon as mu
import numpy as np

import warnings

def preprocess_data(mdata, target_sum = 10000, exclude_highly_expressed = False,
                    remove_mitochondrial = False, view_separator = ':'):
    """
    Preprocesses data in a MuData object.

    Parameters:
    ----------
    mdata :class:`MuData`
        The MuData object.
    target_sum :class:`int`, optional:
        The target sum of total counts after normalization. Parameter passed on to scanpy.pp.normalize_total. Defaults to 10000.
    exclude_highly_expressed :class:`bool`, optional:
        Whether to exclude highly expressed genes during normalization. Parameter passed on to scanpy.pp.normalize_total. Defaults to False.
    remove_mitochondrial :class:`bool`, optional:
        Whether to remove mitochondrial genes. Defaults to False.
    view_separator :class:`str`, optional:
        Separator used to identify view-specific genes. Defaults to ':'.
    """

    for view in mdata.mod.keys():
        mdata.mod[view].layers['counts'] = mdata.mod[view].X.copy()
        print('INFO: Preprocessing {0}'.format(view))
        mdata.mod[view].X = np.nan_to_num(mdata.mod[view].X, copy=True) # remove nan values

        if remove_mitochondrial:
            #identify mitochondrial genes
            mdata.mod[view].var['mt'] = mdata.mod[view].var_names.str.startswith(view + view_separator + 'MT-')
            #remove mitochondrial genes from analysis
            mu.pp.filter_var(mdata.mod[view], 'mt', lambda x: ~x)

        # normalise data using median of ratios (deseq2)
        # deseq2_counts, size_factors = deseq2_norm(mdata.mod[view].X)
        # mdata.mod[view].X = deseq2_counts

        sc.pp.normalize_total(mdata.mod[view], target_sum=target_sum, exclude_highly_expressed=exclude_highly_expressed)
        sc.pp.log1p(mdata.mod[view])

        if mdata.mod[view].shape[0] > 1:
            sc.pp.highly_variable_genes(mdata.mod[view])
            print('Highly variable genes in {0}: {1}'.format(view, mdata.mod[view].var.highly_variable.sum()))
        else:
            print('Not enough cells in {0}'.format(view))
            mdata.mod[view].var['highly_variable'] = False
    # return mdata

def remove_hvg_marker_genes(mdata, markers, view_separator=':', hvg_column='highly_variable', hvg_filtered='filtered_highly_variable'):
    """In each view, sets highly variable genes to False if they are in the markers dict for another view, but not if they are in the markers for the same view.
    Used for removing potential cell type marker genes found in the background of other views and thought to be contamination.

    Parameters
    ----------
    mdata :class:`~mudata.MuData`
        MuData object. Highly variable genes should be computed already in .var for each modality.
    markers :class:`dict`
        Dictionary with markers for each view. Keys are the views and values are lists of markers. Can contain markers for views that are not in mdata.mod.keys().
    view_separator :class:`str`, optional
        Separator between view and gene names. Defaults to ':'.
    hvg_column :class:`str`, optional
        Column in mdata.mod['some_view'].var that contains the highly variable genes. Defaults to 'highly_variable'.
    hvg_filtered :class:`str`, optional
        Column in mdata.mod['some_view'].var where filtered highly variable genes will be stored. Defaults to 'filtered_highly_variable'.
    """
    # check if markers is a dict
    if not isinstance(markers, dict):
        raise TypeError('markers is not a dict')

    # check that all keys in markers are lists
    if not all(isinstance(markers[mod], list) for mod in markers.keys()):
        raise TypeError('not all values in markers are lists')

    # check that hvg_column is in var for all modalities
    if not all(hvg_column in mdata.mod[mod].var.columns for mod in mdata.mod.keys()):
        raise ValueError('{0} is not in the columns of .var for all modalities'.format(hvg_column))

    for current_mod in mdata.mod.keys():
        # markers in markers dict for each modality except for current_mod
        negative_markers = [marker for mod in markers.keys() if mod != current_mod for marker in markers[mod]]

        if current_mod not in list(markers.keys()):
            warnings.warn('no markers in dict for view: {0}'.format(current_mod), Warning)
        else:
            #keep negative_markers not in markers[current_mod] and add view_separator
            negative_markers = [current_mod + view_separator + marker for marker in negative_markers if marker not in markers[current_mod]]
        
        # duplicate hvg_column to hvg_filtered
        mdata.mod[current_mod].var[hvg_filtered] = mdata.mod[current_mod].var[hvg_column]
        # set negative_markers to False in current_mod
        mdata.mod[current_mod].var.loc[mdata.mod[current_mod].var_names.isin(negative_markers), hvg_filtered] = False


def remove_low_variability_views(mdata, highly_variable_col='highly_variable', min_variable_genes=100, to_exclude=None):
    """
    Remove views with low variability from a MuData object.

    Parameters
    ----------
    mdata :class:`~mudata.MuData`
        The MuData object.
    highly_variable_col :class:`str`, optional
        Column name containing the highly variable genes. (default: 'highly_variable')
    min_variable_genes :class:`int`, optional
        Minimum number of highly variable genes required. (default: 100)
    to_exclude :class:`list`, optional
        List of views to exclude. (default: None)
    """

    if to_exclude is None:
        to_exclude = []
    
    for view in to_exclude:
        if view in list(mdata.mod.keys()):
            print('Removing {0} due to exclusion list'.format(view))
            mdata.mod[view].var['to_remove'] =  mdata.mod[view].var_names.str.startswith(view)
            mu.pp.filter_var(mdata.mod[view], 'to_remove', lambda x: ~x)
            del mdata.mod[view]

    for view in list(mdata.mod.keys()):

        if mdata.mod[view].shape[0] > 1:
            if highly_variable_col not in mdata.mod[view].var.columns:
                print('Removing {0} because it does not contain any highly variable genes'.format(view))
                mdata.mod[view].var['to_remove'] =  mdata.mod[view].var_names.str.startswith(view)
                mu.pp.filter_var(mdata.mod[view], 'to_remove', lambda x: ~x)
                mdata.mod[view].var = mdata.mod[view].var.drop(columns='to_remove')
            elif mdata.mod[view].var[highly_variable_col].sum() < min_variable_genes:
                print('Removing {0} due to low number of highly variable genes ({1})'.format(view,  mdata.mod[view].var[highly_variable_col].sum()))
                mdata.mod[view].var['to_remove'] =  mdata.mod[view].var_names.str.startswith(view)
                mu.pp.filter_var(mdata.mod[view], 'to_remove', lambda x: ~x)
                mdata.mod[view].var = mdata.mod[view].var.drop(columns='to_remove')
            else:
                print('Highly variable genes in {0}: {1}'.format(view, mdata.mod[view].var[highly_variable_col].sum()))

    #update mdata after cleaning
    mdata.var.drop(columns=mdata.var.columns, inplace=True)
    # for all views concatenate the var columns
    mdata.var = pd.concat([mdata.mod[view].var for view in mdata.mod.keys()], axis=0)

    # return mdata

def join_mdata(mdata1, mdata2, genes='outer', views='outer', obsm_keys= [], uns_keys=[]):
    
    #check that views is either outer, inner, or left
    assert views in ['outer', 'inner', 'left'], "views must be either 'outer', 'inner', or 'left'"
    assert genes in ['outer', 'inner', 'left'], "genes must be either 'outer', 'inner', or 'left'"
    
    mdatas = [mdata1, mdata2]

    adatas = {}

    if views == 'outer':
        #get all unique views in mdatas
        views = sorted(list(set(view for mdata in mdatas for view in mdata.mod.keys())))
    elif views == 'inner':
        #get all views that are in both mdatas
        views = sorted(list(set(mdata1.mod.keys()).intersection(set(mdata2.mod.keys()))))
    elif views == 'left':
        views = sorted(list(set(mdata1.mod.keys())))
    
    if len(views) == 0:
        return None, None

    for view in views:
        if genes in ['outer', 'inner']:
            adatas[view] = sc.concat([mdata.mod[view] for mdata in mdatas if view in list(mdata.mod.keys()) and 0 not in mdata.mod[view].shape], join=genes)
        elif genes == 'left':
            temp_list = []
            if view in mdata1.mod.keys() and 0 not in mdata1.mod[view].shape:
                temp_list.append(mdata1.mod[view])
            else:
                continue
            if view in mdata2.mod.keys():
                mu.pp.filter_var(mdata2.mod[view], mdata1.mod[view].var_names)
                if 0 not in mdata2.mod[view].shape:
                    temp_list.append(mdata2.mod[view])
            
            if len(temp_list) == 1:
                adatas[view] = temp_list[0]
            else:
                adatas[view] = sc.concat(temp_list, join='outer')

    mdata = mu.MuData(adatas)
    metadata = pd.concat([mdata.obs for mdata in mdatas]).drop_duplicates()
    temp = mdata.obs.reset_index().merge(metadata.reset_index(), on='index', how='outer').set_index('index')
    temp.index.name = None
    mdata.obs = temp
    
    #concatenate obsm keys
    for obsm_key in obsm_keys:
        # check that obsm_key is in both mdatas
        if obsm_key in mdata1.obsm.keys() and obsm_key in mdata2.obsm.keys():

            if all([isinstance(data.obsm[obsm_key], pd.DataFrame) for data in mdatas]):
                obsm_dfs = [data.obsm[obsm_key] for data in mdatas]
                
            elif all([isinstance(mdata.obsm[obsm_key], np.ndarray) for mdata in mdatas]):
                obsm_dfs = [pd.DataFrame(data.obsm[obsm_key], index=data.obs.index) for data in mdatas]
                assert all([df.shape[1] == obsm_dfs[0].shape[1] for df in obsm_dfs]), 'both numpy arrays must have same number of columns'
            else:
                raise TypeError('obsm key: ' + obsm_key + ' must be either pandas dataframe or numpy array')
                
            df = pd.concat(obsm_dfs, axis = 0)
            #order df by index as in mdata.obs
            df = df.loc[mdata.obs.index]
            mdata.obsm[obsm_key] = df.values if isinstance(mdata1.obsm[obsm_key], np.ndarray) else df
    
    for uns_key in uns_keys:
        if all([uns_key in data.uns.keys() for data in mdatas]):
            if all([isinstance(data.uns[uns_key], pd.DataFrame) for data in mdatas]):
                mdata.uns[uns_key] = pd.concat([data.uns[uns_key] for data in mdatas], axis=0)
                continue

        for data, suffix in zip(mdatas, ['_x','_y']):
            if uns_key in data.uns.keys():
                mdata.uns[uns_key + suffix] = data.uns[uns_key]

    return mdata, mdata.obs

def get_psbulk_stats(mdata, stat, view_separator = ':', uns_key = 'psbulk_stats'):
    """Get cell counts for each patient and omic in a MuData object.

    Parameters
    ----------
    mdata : MuData
        MuData object containing the data.
    stat : str
        Name of the stat to get. Must be in mdata.uns[uns_key].columns. e.g. 'psbulk_n_cells' or 'psbulk_counts'
    view_separator : str, optional
        Separator used to separate the view name from column names in the MuData object (e.g. of var, or psbulk_stats), by default ':'
    uns_key : str, optional
        Key in mdata.uns where the cell counts will is stored, by default 'stats_psbulk'

    Returns
    -------
    pd.DataFrame
        Dataframe containing the cell counts for each patient and omic.
    """
    
    assert uns_key in mdata.uns.keys(), 'Key {0} not found in mdata.uns'.format(uns_key)
    stats = mdata.uns[uns_key]
    
    #get all columns with psbulk_n_cells in the name
    stats = stats[[col for col in stats.columns if stat in col]]
    
    #split column names by : and get the first element
    stats.columns = [col.split(view_separator)[0] for col in stats.columns]
    
    #keep only columns with name in mdata.mod.keys()
    stats = stats[[col for col in stats.columns if col in mdata.mod.keys()]]
    
    #order columns alphabetically
    stats = stats[sorted(stats.columns)]
    
    return stats