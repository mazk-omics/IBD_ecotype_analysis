#!/usr/bin/env python3
"""Integration fixtures covering source alignment, shared holdout and MTX scan."""
import argparse
import hashlib
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from export_three_reference_web_pilots import (
    ONTOLOGY, SCP_TYPES, BUNDLES, run_gse, run_scp,
    stream_mtx_selected, coalesce_selected,
)


def args(root, mode):
    return argparse.Namespace(mode=mode, input=str(root/'source.h5ad'), raw=str(root/'raw'),
        candidate=str(root/('gse.parquet' if mode.startswith('gse') else 'scp.parquet')),
        donor_split=str(root/'split.tsv'), outdir=str(root/mode), cells_per_type=10, seed=20261006)


def validate_web(root, mode, expected_types, expected_values=None):
    out=root/mode
    m=pd.read_csv(out/'selected_cell_manifest.tsv',sep='\t')
    assert len(m)==10*len(expected_types) and m.groupby('pilot_type').size().eq(10).all()
    assert set(m.pilot_type)==set(expected_types) and m.groupby(['pilot_type','individual_id']).size().max()<=5
    filename={'gse_joint':'GSE282122_Joint20_scRNA_reference.txt',
              'gse_ileum':'GSE282122_Ileum19_scRNA_reference.txt',
              'scp259':'SCP259_ColonRectum18_scRNA_reference.txt'}[mode]
    web=pd.read_csv(out/filename,sep='\t',index_col=0)
    assert web.columns.tolist()==m.web_column.tolist()
    assert 'ZERO' not in web.index
    if expected_values is not None:
        assert np.array_equal(web.loc['EPCAM'].to_numpy(),expected_values[m.source_cell_row,0])
        assert not m.individual_id.isin(['CD0','UC0','Healthy0']).any()
        assert 'RGS5' not in web.index
    else:
        assert m.total_counts_matches_obs.all() and m.n_genes_matches_obs.all()
        split=pd.read_csv(out/'donor_split.tsv',sep='\t')
        assert not m.individual_id.isin(split.loc[split.split.eq('holdout'),'individual_id']).any()
    return m


def main():
    with tempfile.TemporaryDirectory(prefix='three_reference_tests_') as tmp:
        root=Path(tmp);records=[]
        for disease,n in [('CD',12),('UC',12),('Healthy',3)]:
            for donor in range(n):
                for name in ONTOLOGY:
                    for rep in range(4):
                        i=len(records); site='Terminal_Ileum' if (name=='Ileal_absorptive' or (name!='Colonic_absorptive' and rep<2)) else 'Descending_Colon'
                        records.append({'cell_id':f'cell{i}','obs_name':f'cell{i}','Patient':f'{disease}{donor}',
                            'Disease':disease,'Site':site,'sample_id':f'{disease}{donor}_{site}',
                            'Treatment':'ND' if disease=='Healthy' else 'Pre','Inflammation':'Healthy' if disease=='Healthy' else 'Inflamed',
                            'reference_eligible':True,'coarse_identity':name,'identity_status':'CONFIDENT',
                            'primary_reference_action':'HOLD_CROSS_DATASET' if name in ['Glial','Tuft'] else 'KEEP',
                            'final_analysis':name+' subtype','predicted_doublets':rep==3,
                            'n_genes_by_counts':4,'total_counts':84+i%7})
        c=pd.DataFrame(records);n=len(c)
        c.sample(frac=1,random_state=9).to_parquet(root/'gse.parquet',index=False)
        d=c[['Patient','Disease']].drop_duplicates();d['split']=np.where(d.Patient.isin(['CD0','UC0','Healthy0']),'holdout','training')
        d.to_csv(root/'split.tsv',sep='\t',index=False)
        values=np.column_stack([70+np.arange(n)%7,np.ones(n)*12,np.ones(n),np.ones(n)]).astype(np.int32)
        strings=h5py.string_dtype('utf-8')
        with h5py.File(root/'source.h5ad','w') as h:
            obs=h.create_group('obs');obs.attrs['_index']='_index'
            obs.create_dataset('_index',data=c.obs_name.to_numpy(dtype=object),dtype=strings)
            for name in ['sample_id','Patient','Disease','Site','final_analysis']:
                obs.create_dataset(name,data=c[name].to_numpy(dtype=object),dtype=strings)
            var=h.create_group('var')
            var.create_dataset('gene_id',data=np.asarray(['g0','g1','g2','g3','g4'],dtype=object),dtype=strings)
            var.create_dataset('gene_symbol',data=np.asarray(['EPCAM','COL1A1','RGS5','RGS5','ZERO'],dtype=object),dtype=strings)
            x=h.create_group('X');x.attrs['encoding-type']='csr_matrix';x.attrs['shape']=[n,5]
            x.create_dataset('data',data=values.ravel());x.create_dataset('indices',data=np.tile(np.arange(4),n));x.create_dataset('indptr',data=np.arange(n+1)*4)
        before=hashlib.sha256((root/'source.h5ad').read_bytes()).hexdigest()
        run_gse(args(root,'gse_joint'));joint=validate_web(root,'gse_joint',ONTOLOGY,values)
        run_gse(args(root,'gse_ileum'));ileum=validate_web(root,'gse_ileum',[t for t in ONTOLOGY if t!='Colonic_absorptive'],values)
        assert ileum.Site.eq('Terminal_Ileum').all()
        assert (root/'gse_joint/donor_split.tsv').read_bytes()==(root/'gse_ileum/donor_split.tsv').read_bytes()
        assert hashlib.sha256((root/'source.h5ad').read_bytes()).hexdigest()==before
        duplicate=args(root,'gse_joint')
        try:run_gse(duplicate)
        except ValueError as e:assert 'overwrite' in str(e)
        else:raise AssertionError('Output overwrite accepted')
        # Build three released-format MTX bundles with overlapping gene symbols.
        epithelial={'Colonic_absorptive','Goblet','Stem_TA_progenitor','Enteroendocrine','Tuft'}
        stromal={'Fibroblast','Pericyte','Endothelial','Glial'}
        rec=[];by_bundle={b:[] for b in BUNDLES}
        for disease in ['Healthy','UC']:
            for donor in range(8):
                for name in SCP_TYPES:
                    for rep in range(3):
                        bundle='Epithelial' if name in epithelial else 'Stromal' if name in stromal else 'Immune'
                        cell=f'scp{len(rec)}';by_bundle[bundle].append(cell)
                        rec.append({'cell_id':cell,'NAME':cell,'Cluster':name,'Subject':f'{disease}{donor}',
                            'Health':'Healthy' if disease=='Healthy' else 'Inflamed','Location':'Epi' if bundle=='Epithelial' else 'LP',
                            'Sample':f'{disease}{donor}sample','expression_bundle':bundle,'donor_id':f'{disease}{donor}',
                            'disease':disease,'inflammation':'Healthy' if disease=='Healthy' else 'Inflamed',
                            'coarse_identity':name,'primary_reference_action':'HOLD_CROSS_DATASET' if name in ['Glial','Tuft'] else 'KEEP',
                            'nGene':2,'nUMI':82})
        scp=pd.DataFrame(rec);scp.sample(frac=1,random_state=3).to_parquet(root/'scp.parquet',index=False)
        meta=root/'raw/metadata';meta.mkdir(parents=True)
        metadata=scp[['NAME','Cluster','nGene','nUMI','Subject','Health','Location','Sample']].copy()
        pd.concat([pd.DataFrame([{k:'TYPE' if k=='NAME' else 'group' for k in metadata.columns}]),metadata]).to_csv(meta/'all.meta2.txt',sep='\t',index=False)
        paths=[]
        for b,(folder,prefix) in BUNDLES.items():
            p=root/'raw/expression'/folder;p.mkdir(parents=True)
            cells=by_bundle[b];(p/(prefix+'.barcodes2.tsv')).write_text('\n'.join(cells)+'\n')
            (p/(prefix+'.genes.tsv')).write_text('EPCAM\nCOL1A1\nZERO\n')
            path=p/('gene_sorted-'+prefix+'.matrix.mtx');paths.append(path)
            with path.open('w') as f:
                f.write('%%MatrixMarket matrix coordinate integer general\n% test comment\n'+f'3 {len(cells)} {2*len(cells)}\n')
                for j in range(len(cells)):f.write(f'1 {j+1} 70\n2 {j+1} 12\n')
        run_scp(args(root,'scp259'));validate_web(root,'scp259',SCP_TYPES)
        # Chunk boundaries and repeated coordinates must preserve exact totals.
        path=root/'duplicate.mtx';path.write_text('%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 2 3\n1 2 4\n2 2 5\n')
        selected=pd.DataFrame({'cell_id':['b']})
        ci,gi,vv,nnz=stream_mtx_selected(path,['g0','g1'],['a','b'],selected,{'g0':0,'g1':1},chunksize=1)
        result=coalesce_selected(ci,gi,vv,1,2)
        assert nnz==3 and result[1].tolist()==[7,5] and result[3].tolist()==[12]
        path.write_text('%%MatrixMarket matrix coordinate integer general\n2 2 1\n1 2 -3\n')
        try:stream_mtx_selected(path,['g0','g1'],['a','b'],selected,{'g0':0,'g1':1},chunksize=1)
        except ValueError as e:assert 'nonnegative' in str(e)
        else:raise AssertionError('Negative MTX count accepted')
    print('PASS: GSE joint/ileum keyed alignment, shared donor holdout, Glial/Tuft promotion, SCP bundle/barcode alignment, bounded streaming, duplicate-coordinate coalescing, gene exclusion, no overwrite and negative-count rejection')


if __name__=='__main__':main()
