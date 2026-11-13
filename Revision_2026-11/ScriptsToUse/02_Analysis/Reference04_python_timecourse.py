import pandas as pd,numpy as np
from scipy import stats
r=pd.read_csv('97df6d39-Expt1234_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_RAWvalues.csv').rename(columns={'TimePoints':'TimePoint'})
r['Treatment']=r.Sample_Condition.str.replace(r'_[^_]+$','',regex=True)
meta=['Unnamed: 0','File_Name','Sample_Condition','Reordered','File_Num','TimePoint','Cluster','Replicate','Treatment']
mk=[c for c in r.columns if c not in meta]
dg=pd.read_csv('9c56c7e6-Statistics_diagnostics.csv').set_index(['Marker','Treatment'])
def holm(p):
    p=np.asarray(p,float);o=np.argsort(p);m=len(p);out=np.empty(m);run=0
    for rank,i in enumerate(o):
        run=max(run,(m-rank)*p[i]);out[i]=min(1,run)
    return out
def bh(p):
    p=np.asarray(p,float);ok=~np.isnan(p);q=np.full(len(p),np.nan);pp=p[ok];m=len(pp);o=np.argsort(-pp)
    cm=np.minimum.accumulate((pp[o]*m/(m-np.arange(m))));res=np.empty(m);res[o]=np.minimum(1,cm);q[ok]=res;return q
def welch_anova(groups):
    k=len(groups);n=np.array([len(g) for g in groups]);m=np.array([g.mean() for g in groups]);v=np.array([g.var(ddof=1) for g in groups])
    w=n/v;mw=(w*m).sum()/w.sum();A=(w*(m-mw)**2).sum()/(k-1);B=1+2*(k-2)/(k**2-1)*((1-w/w.sum())**2/(n-1)).sum()
    F=A/B;df2=(k**2-1)/(3*((1-w/w.sum())**2/(n-1)).sum());return stats.f.sf(F,k-1,df2)
rows=[];fam=[]
for T in ['BDNF','CompMed','ContDepr']:
  d1=r[r.Treatment==T];d0=r[(r.Treatment=='NoTreat')&(r.TimePoint==0)]
  tps=sorted(d1.TimePoint.unique())
  for m in mk:
    ref=d0[m].values; gs=[d1[d1.TimePoint==t][m].values for t in tps]
    d=dg.loc[(m,T)]
    reasons=[]
    if d.Min_N<3: reasons.append('n')
    if d.Variance_Ratio>3: reasons.append('vr')
    if d.BrownForsythe_p<0.05: reasons.append('bf')
    mode='welch' if reasons else 'classic'
    allg=[ref]+gs
    om=stats.f_oneway(*allg).pvalue if mode=='classic' else welch_anova(allg)
    pu=[stats.ttest_ind(g,ref,equal_var=False).pvalue for g in gs]
    if mode=='classic': pp=stats.dunnett(*gs,control=ref,random_state=1).pvalue
    else: pp=holm(pu)
    fam.append(dict(Marker=m,Treatment=T,mode=mode,omni=om))
    for t,g,a,b in zip(tps,gs,pu,pp): rows.append(dict(Marker=m,Treatment=T,TimePoint=t,p_unadjusted=a,p_primary=b,mode=mode,omnibus=om))
L=pd.DataFrame(rows);F=pd.DataFrame(fam)
L['q_BH']=L.groupby('Treatment').p_unadjusted.transform(lambda s: pd.Series(bh(s.values),index=s.index))
F['q_om']=F.groupby('Treatment').omni.transform(lambda s: pd.Series(bh(s.values),index=s.index))
L.to_csv('/tmp/w/ref04_long.csv',index=False);F.to_csv('/tmp/w/ref04_fam.csv',index=False)
print(F['mode'].value_counts().to_dict())
print('unadj<.05',(L.p_unadjusted<.05).sum(),'primary<.05',(L.p_primary<.05).sum(),'BH<.05',(L.q_BH<.05).sum(),'of',len(L))
print('omnibus<.05',(F.omni<.05).sum(),'q<.05',(F.q_om<.05).sum(),'of',len(F))
print(L.groupby('Treatment').apply(lambda g: pd.Series({'unadj':(g.p_unadjusted<.05).sum(),'primary':(g.p_primary<.05).sum(),'BH':(g.q_BH<.05).sum(),'n':len(g)})))
print(L[L.Marker.isin(['pERK_171Yb','pTrkB_170Er'])].sort_values(['Marker','Treatment','TimePoint']).round(4).to_string())
